import wrf
import numpy as np
import pandas as pd
from datetime import datetime as dt
from netCDF4 import Dataset #type: ignore
from termcolor import cprint
from numpy import unravel_index
import collections.abc as c
import numpy.typing as npt


class Base_point:
    '''
    Base_point: parent class for validation. Sets location name.
    '''
    def __init__(self: object, loc: str, **kwargs) -> None:
        self.loc = loc

class UA_point(Base_point):
    def __init__(self, loc: str, **kwargs) -> None:
        super().__init__(loc, **kwargs)
    def __eq__(self, other: object) -> bool:
        try: #! come back to this...
            results: npt.NDArray[np.bool_] = np.empty(5, np.bool_)
            results[0] = np.allclose(self.p, other.p, atol=10.)
            results[1] = np.allclose(self.t, other.t, atol=1.)
            results[2] = np.allclose(self.td, other.td, atol=1.)
            results[3] = np.allclose(self.wdir, other.wdir, atol=5.)
            results[4] = np.allclose(self.wspd, other.wspd, atol=1.)
            result: bool = bool(results.all())
        except AttributeError:
            if isinstance(self, WRF_point):
                result: bool = other.__eq__(self)
            else:
                raise NotImplementedError('Compairison not implemented')
        finally:
            return result

class Sat_point(Base_point):
    def __init__(self, loc: str, **kwargs) -> None:
        super().__init__(loc, **kwargs)
    def sat_loc(self: object, ulat: float, ulon: float, lats: npt.ArrayLike, lons: npt.ArrayLike) -> tuple[int, int] | int:
        R: int = 6371000
        lat1: float = np.radians(ulat)
        lat2: npt.NDArray = np.radians(lats)
        delta_lat: npt.NDArray = np.radians(lats-ulat)
        delta_lon: npt.NDArray = np.radians(lons-ulon)
        a: npt.NDArray = (np.sin(delta_lat/2))*(np.sin(delta_lat/2))+(np.cos(lat1))*(np.cos(lat2))*(np.sin(delta_lon/2))*(np.sin(delta_lon/2))
        c: npt.NDArray = 2*np.arctan2(np.sqrt(a),np.sqrt(1-a))
        d: npt.NDArray = R*c
        if d.ndim == 1:
            return d.argmin()
        else:
            x: int
            y: int
            x, y = unravel_index(d.argmin(),d.shape) #type: tuple[int, int]
            return x, y
    def __eq__(self, other: object) -> bool:
        result = None
        try:
            if isinstance(self, Tropomi_point) or (isinstance(self, WRF_point) and isinstance(other, Tropomi_point)):
                results = np.empty(2, bool)
                results[0] = np.abs(self.xch4 - other.xch4) <= 1.e-1
                results[1] = np.abs(self.xco - other.xco) <= 1.e-1
                result = results.all()
            #! elif for OCO-2 here
        except AttributeError:
            if isinstance(self, WRF_point):
                result = other.__eq__(self)
            else:
                result = NotImplemented
        else:
            if result is None:
                if isinstance(self,WRF_point):
                    result = other.__eq__(self)
                else:
                    result = NotImplemented
            else:
                pass
        finally:
            return result
        

class Surface_point(Base_point):
    def __init__(self, loc, **kwargs):
        super().__init__(loc, **kwargs)
    def __eq__(self, other: object) -> bool:
        try:
            results = np.empty(5, bool)
            results[0] = abs(self.T2 - other.T2) <= 0.5
            results[1] = abs(self.td2 - other.td2) <= 0.5
            results[2] = abs(self.slp - other.p) <= 3.0
            results[3] = abs(self.wspd10 - other.wspd10) <= 0.2
            results[4] = abs(self.wdir10 - other.wdir10) <= 5
            result = results.all()
        except AttributeError:
            if isinstance(self, WRF_point):
                result = other.__eq__(self)
            else:
                result = NotImplemented
        finally:
            return result

class WRF_point(Surface_point,Sat_point,UA_point):
    '''
    WRF_point: sets up validation point using WRF data. Reads T2, TD2, SLP, and 10m Wind speed/direction. Inherits from Base_point.
    '''
    def __init__(self, wrffile: Dataset, lat: float, lon: float, loc: str, chem: bool = None, **kwargs):
        super().__init__(loc,**kwargs)
        self.lat = lat
        self.lon = lon
        self.x, self.y = wrf.ll_to_xy(wrffile,self.lat,self.lon)
        self.vars = ['T2', 'td2', 'slp','uvmet10_wspd_wdir','p','temp','td','uvmet_wspd_wdir']
        for var in self.vars:
            if var == 'T2':
                self.T2 = wrf.getvar(wrffile, var, meta=False)[self.y, self.x]
            elif var == 'td2':
                self.td2 = wrf.getvar(wrffile, var, meta=False, units='K')[self.y, self.x]
            elif var == 'slp':
                self.slp = wrf.getvar(wrffile, var, meta=False, units='hPa')[self.y, self.x]
            elif var == 'uvmet10_wspd_wdir':
                self.wspd10, self.wdir10 = wrf.getvar(wrffile, var, meta=False)[:, self.y, self.x]
            elif var == 'p':
                self.p = wrf.getvar(wrffile, var, meta=False, units='hpa')[:, self.y, self.x]
            elif var == 'temp':
                self.t = wrf.getvar(wrffile, var, meta=False)[:, self.y, self.x]
            elif var == 'td':
                self.td = wrf.getvar(wrffile, var, meta=False, units='K')[:, self.y, self.x]
            elif var == 'uvmet_wspd_wdir':
                self.wspd, self.wdir = wrf.getvar(wrffile, var, meta=False)[..., self.y, self.x]
        ## ! Next if/elif block needs to be edited depending on WRF-GHG output structure (include converting units) ! ## 
        self.sfc_pres = wrffile['PSFC'][0, self.y, self.x].tolist() #! may need to calculate or may be in WRF output
        if chem is True:
            self.xch4 = self._extract_ghg(wrffile, 'xch4')
            self.xco = self._extract_ghg(wrffile, 'xco')
            self.xco2 = self._extract_ghg(wrffile, 'xco2')
    def _extract_ghg(self, wrffile: Dataset, chem: str):
        if chem == 'xch4':
            _ant = wrffile['CH4_ANT'][0, :, self.y, self.x]
            _bck = wrffile['CH4_BCK'][0, :, self.y, self.x]
            _tst = wrffile['CH4_TST'][0, :, self.y, self.x]
        elif chem == 'xco2':
            _ant = wrffile['CO2_ANT'][0, :, self.y, self.x]
            _bck = wrffile['CO2_BCK'][0, :, self.y, self.x]
            _tst = wrffile['CO2_TST'][0, :, self.y, self.x]
        elif chem == 'xco':
            _ant = wrffile['CO_ANT'][0, :, self.y, self.x]
            _bck = wrffile['CO_BCK'][0, :, self.y, self.x]
            #_tst = wrffile['CO_BIO'][0, :, self.y, self.x] ##?
            _tst = np.zeros_like(_bck)
        _ghg = _tst + _ant -_bck
        if len(_ghg) == len(self.p):
            pres_bound = np.empty_like(self.p)
            for i, pres in enumerate(self.p):
                if i == 0:
                    pres_bound[i] = self.sfc_pres
                    pres_bound[i+1] = pres_bound[i] + (2*(pres-pres_bound[i]))
                else:
                    try:
                        pres_bound[i+1] = pres_bound[i] + (2*(pres-pres_bound[i]))
                    except IndexError:
                        pass
        p_layer_diff = np.array([pres_bound[i]-pres_bound[i-1] for i in range(1,len(pres_bound))]) #! This may need a value at beginning for xch4[0]
        p_diff = pres_bound[0] - pres_bound[-1]
        return np.sum(_ghg*p_layer_diff)/p_diff
    def __str__(self) -> str:
        return f'{self.loc} WRF Point has a temperature of {self.T2} K, a dewpoint of {self.td2} K, a slp of {self.slp} hPa, and the wind is {self.wspd10} m s^-1 at {self.wdir10} degrees.'
    
class Obs_point(Surface_point):
    '''
    Obs_point: sets up validation point for observation. Reads T2, TD2, slp, and wind speed/direction from ASOS data. Inherits from Base_point.
    '''
    def __init__(self, loc, obsfile, obstime, **kwargs):
       super.__init__(loc, **kwargs)
       data = pd.read_csv(obsfile,na_values='M',parse_dates=['valid'],date_format='%Y-%m-%d %H:%M')
       for i, time in enumerate(data.valid):
           if abs(obstime.timestamp() - time.timestamp()) == 420:
               idx = i
               break
       self.T2 = data.tmpc[idx] + 273.15
       self.td2 = data.dwpc[idx] + 273.15
       self.slp = data.mslp[idx]
       self.wdir10 = data.drct[idx]
       self.wspd10 = data.sped[idx] * 0.44704
       del data
    def __str__(self) -> str:
        return f'{self.loc} Observation Point has a temperature of {self.T2} K, a dewpoint of {self.td2} K, a slp of {self.slp} hPa, and the wind is {self.wspd10} m s^-1 at {self.wdir10} degrees.'
    

class UObs_point(UA_point):
    def __init__(self, loc, ua_file, wrffile, lat, lon, **kwargs):
        super().__init__(loc, **kwargs)
        x, y = wrf.ll_to_xy(wrffile, lat, lon)
        wrf_p = wrf.getvar(wrffile, 'p', meta=False, units='hPa')[:, y, x]
        data = pd.read_csv(ua_file, na_values='M', parse_dates=['validUTC'], date_format='%Y-%m-%d %H:%M')
        p = data.pressure_mb.to_numpy()
        idx = np.digitize(wrf_p, p)
        data = data.iloc[idx]
        self.p = data.pressure_mb.to_numpy() #* mb == hPa
        self.t = data.tmpc.to_numpy() + 273.15 #* deg C -> K
        self.td = data.dwpc.to_numpy() + 273.15 #* deg C -> K
        self.wdir = data.drct.to_numpy() #* deg
        self.wspd = data.speed_kts.to_numpy() * 0.514444 #* kt -> m s^-1
        del data, p, idx, wrf_p, x, y


class Tropomi_point(Sat_point): #! Find some way to do both xch4 and xco, will probably need two files
    def __init__(self, loc, xch4_f, xco_f, ulat, ulon, **kwargs):
        super().__init__(loc, **kwargs)
        self.xch4 = self._ghg(xch4_f, 'xch4', ulat, ulon)
        self.xco = self._ghg(xco_f, 'xco', ulat, ulon)
    def _ghg(self, tropomi_f, chem, ulat, ulon):
        ds = Dataset(tropomi_f, 'r')
        grp = 'PRODUCT'
        if chem in ['ch4','CH4','xch4','XCH4','methane','Methane','METHANE']: #!! find a way to do both
            self.xch4_sds = sds = 'methane_mixing_ratio'
            self.xch4_unit = 'ppb'
        elif chem in ['co','CO','xco','XCO','carbon monoxide','Carbon Monoxide', 'CARBON MONOXIDE']:
            self.xco_sds = sds = 'carbonmonoxide_total_column_corrected' # ! is this right?
            self.xco_unit = 'ppb' # * will have to convert for co
        lats = ds.groups[grp].variables['latitude'][0][:][:]
        lons = ds.groups[grp].variables['longitude'][0][:][:]
        qas = np.array(ds.groups[grp].variables['qa_value'][0][:][:])
        data = ds.groups[grp].variables[sds] #units: ppb (ch4) mol m^-2 --> ppb (co)
        fv = data._FillValue
        dA = np.array(data[0][:][:])
        dA[(dA==fv) & (qas<=0.5)] = np.nan #? Do I want to filter for qa here or later?
        if chem in ['co','CO','xco','XCO','carbon monoxide','Carbon Monoxide', 'CARBON MONOXIDE']:
            dA = dA * 28.01 * 0.0001 * 1000 # converts (mol m^-2) * (g mol^-1) * (m^-1) == g m^-3 == ppm --> ppb
        #TODO: check about averging kernel
        min_lat = np.min(lats)
        max_lat = np.max(lats)
        min_lon = np.min(lons)
        max_lon = np.max(lons)
        if not min_lat <= ulat <= max_lat:
            raise RuntimeError(f'User Latitude is not within TROPOMI file. {ulat}, {min_lat}, {max_lat}')
        if not min_lon <= ulon <= max_lon:
            raise RuntimeError(f'User Longitude is not within TROPOMI file. {ulon}, {min_lon}, {max_lon}')
        x, y = self.sat_loc(ulat, ulon, lats, lons)
        if np.isnan(dA[x,y]): #? qa check here?
            raise RuntimeError('No value at desired location.')
        if x < 1:
            x += 1
        if x > dA.shape[0]-2:
            x -= 2
        if y < 1:
            y += 1
        if y > dA.shape[1]-2:
            y -= 2
        t_b_t = dA[x-1:x+2,y-1:y+2].astype(float)
        # ? is this necessary? t_b_t[t_b_t==float(fv)] = np.nan
        nnan = np.count_nonzero(~np.isnan(t_b_t))
        if nnan == 0:
            raise RuntimeError('No valid pixels in 3x3 grid.')
        grid_avg = np.nanmean(t_b_t)
        grid_std = np.nanstd(t_b_t)
        if np.abs(grid_avg-dA[x,y]) <= grid_std:
             return grid_avg
        else:
             return dA[x,y]
    def __str__(self) -> str:
        return f'TROPOMI Satelite products at {self.loc} are {self.xch4:.2e} {self.xch4_unit} of CH4 and {self.xco:.2e} {self.xco_unit} of CO'


class OCO_Point(Sat_point):
    def __init__(self, loc, oco_f, ulat, ulon, **kwargs):
        super().__init__(loc, **kwargs)
        ds = Dataset(oco_f, 'r')
        lats = ds['latitude'][:]
        lons = ds['longitude'][:]
        idx = self.sat_loc(ulat, ulon, lats, lons)
        _xco2 = ds['']



def gprint(x): return cprint(x, 'white', 'on_green', attrs=['bold'],end=' ')
def rprint(x): return cprint(x, 'white', 'on_red', attrs=['blink','bold'], end=' ')

def main():
    wrf_file = 'insert path here'
    obs_files = [] #list of file paths
    tropomi_f = '' # TODO: add file path
    lats = [] #insert lat here#
    lons = [] #insert lon here#
    locs = [] #insert location names here
    chem = ''
    eval_time = dt(2024,7,3,13,0) # !edit as needed


    wrf_data = Dataset(wrf_file)
    sfc_results = np.empty(len(locs),bool)
    print('Beginning surface validation...')
    for i, lat, lon, loc, obs_data in enumerate(zip(lats, lons, locs, obs_files)):
        print(f'Evaluating {loc}...')
        wrf_p = WRF_point(wrf_data, lat, lon, loc, chem)
        obs_p = Obs_point(loc, obs_data, eval_time)
        print(wrf_p)
        print(obs_p)
        result = (wrf_p == obs_p)
        sfc_results[i] = result
        if result:
            gprint('PASS:'); print(f'Surface meteorology at {loc} is within tolerance.')
        else:
            rprint('FAIL:'); print(f'Surface meteorology at {loc} is',end=' '); rprint('NOT'); print('within tolerance.')
    print('Surface validation complete!')

    print('Beginning satelite validation...')
    sat_results = np.empty_like(sfc_results, bool)
    for i, lat, lon, loc in enumerate(zip(lats, lons, locs)):
        print(f'Evaluating {loc}...')
        wrf_p = WRF_point(wrf_data, lat, lon, loc, chem)
        tro_p = Tropomi_point(loc, tropomi_f, chem, lat, lon)
        print(wrf_p)
        print(tro_p)
        result = (wrf_p == tro_p)
        sat_results[i] = result
        if result:
            gprint('PASS:'); print(f'Satelite chemistry {chem} at {loc} is within tolerance.')
        else:
            rprint('FAIL:'); print(f'Satelite chemistry {chem} at {loc} is',end=' '); rprint('NOT'); print('within tolerance.')
    print('Satelite validation complete!')

    results = np.append(sfc_results, sat_results)
    if results.all():
        cprint('This run has PASSED validation.', 'white', 'on_green', attrs=['blink','bold'])
    else:
        cprint('This run has FAILED validation. Please check results and settings and rerun.', 'white', 'on_red', attrs=['blink','bold'])


if __name__ == '__main__':
    main()