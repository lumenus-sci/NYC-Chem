#! /usr/bin/env python3
#
# python script to download selected files from rda.ucar.edu
# after you save the file, don't forget to make it executable
#   i.e. - "chmod 755 <name_of_script>"
#
import requests
#
print('Beginning download...')
files = [
    "2023/20230812/gfs.0p25.2023081200.f000.grib2",
    "2023/20230812/gfs.0p25.2023081206.f000.grib2",
    "2023/20230812/gfs.0p25.2023081212.f000.grib2",
    "2023/20230812/gfs.0p25.2023081218.f000.grib2",
]

print(f'File list contains {len(files)} files')
#
# download the data file(s)
for file in files:
    print(f'Downloading {file}')
    idx = file.rfind("/")
    if (idx > 0):
        ofile = file[idx+1:]
    else:
        ofile = file

    response = requests.get("https://data.rda.ucar.edu/d084001/" + file)
    with open(ofile, "wb") as f:
        f.write(response.content)

print('Download complete!')
