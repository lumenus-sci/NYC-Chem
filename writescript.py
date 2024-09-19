#import argparse as ap
'''Writes the various submit scripts for submition to Stampede3 (TACC)'''
#parser = ap.ArgumentParser(description='short sample')
#parser.add_argument('--seq', action='store', dest='seq', default=0, type=int)
#parser.add_argument('--day', action='store', dest='day', default=None, type=int)
#parser.add_argument('--mon', action='store', dest='mon', default=None, type=int)
#parser.add_argument('--hr', action='store', dest='hr', default=None, type=int)
#args = parser.parse_args()
seqs = [x for x in range(1,49)]
days = [x for x in range(20,32) for _ in (0,1)]
days2 = [x for x in range(1,13) for _ in (0,1)]
days.extend(days2)
mons = [7] * 24
mons2 = [8] * 24
mons.extend(mons2)
hrs = [0, 12] * 24

for seq, mon, day, hr in zip(seqs, mons, days, hrs):
    with open(f'wrf_ghg_{seq:02}.sh','a') as fl:
        fl.write('#!/bin/bash -l\n')
        fl.write('#\n')
        fl.write(f'#SBATCH -J WRF_GHG_{seq:02}\n')
        fl.write(f'#SBATCH -e WRF_GHG_{seq:02}.e.%j\n')
        fl.write(f'#SBATCH -o WRF_GHG_{seq:02}.o.%j\n')
        fl.write('#SBATCH --ntasks-per-node 48\n')
        fl.write('#SBATCH -N 4       # skx-dev\n') 
        fl.write('#SBATCH -p skx # Queue name\n')
        fl.write('#SBATCH -t 48:00:00       # Run time (hh:mm:ss) - 1.5 hours\n')
        fl.write('#SBATCH --mail-user=steve@belumenus.com\n')
        fl.write('#SBATCH --mail-type=all\n')
        fl.write('\n')
        if seq == 1:
            fl.write('ln -s $HOME/work/WRF-4.5.2/WRF/WRFV4.5.2/run/* .\n')
            fl.write('ln -s ../../wps/2023-nudge/met_em* .\n')
            fl.write('ln -s $HOME/work/chem-files/wrf* .\n')
            fl.write('ln -s $HOME/work/chem-files/hist* .\n')
            fl.write('ln -sf namelist.input.real namelist.input\n')
            fl.write('ibrun ./real.exe >& real.log\n')
            fl.write('mv rsl.out.0000 rsl.out.real\n')
            fl.write('mv rsl.error.0000 rsl.error.real\n')
            fl.write('\n')
            fl.write('cp wrfinput_d01 wrfinput_d01_beforeCT\n')
            fl.write('cp wrfbdy_d01 wrfbdy_d01_beforeCT\n')
            fl.write('\n')
            fl.write('ln -s $HOME/mozbc/*2024*.inp .\n')
            fl.write('ln -s $HOME/mozbc/CO2CAMS*.inp .\n')
            fl.write('ln -s $HOME/mozbc/mozbc .\n')
            fl.write('ln -s $HOME/work/chem-files/CH4/CH4_CAMS/cams73_latest_ch4_conc_surface_satellite_inst_2023_00??.nc .\n')
            fl.write('./mozbc < CH4_2024yearly08.inp >& CH4_d01.log\n')
            fl.write('./mozbc < CH4_2024yearly08_d02.inp >& CH4_d02.log\n')
            fl.write('\n')
            fl.write('ln -s $HOME/work/chem-files/CO2_CAMS/cams73*2023* .\n')
            fl.write('\n')
            fl.write('./mozbc < CO2CAMS_2024yearly08.inp >& CO2_d01.log\n')
            fl.write('./mozbc < CO2CAMS_2024yearly08_d02.inp >& CO2_d02.log\n')
            fl.write('\n')
        fl.write(f'ln -sf namelist.input.{seq:02} namelist.input\n')
        if seq > 1:
            fl.write(f'ncatted -O -h -a MMINLU,global,m,c,"MODIFIED_IGBP_MODIS_NOAH" wrfrst_d02_2023-{mon:02}-{day:02}_{hr:02}:00:00 wrfrst_d02_2023-{mon:02}-{day:02}_{hr:02}:00:00\n')
        fl.write(f'\n')
        fl.write('ibrun ./wrf.exe >& wrf.log\n')
        fl.write(f'mv rsl.out.0000 rsl.out.wrf.{seq:02}\n')
        fl.write(f'mv rsl.error.0000 rsl.error.wrf.{seq:02}\n')
        fl.write('\n')
        



