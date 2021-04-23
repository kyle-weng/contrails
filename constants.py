import sys
import os

version = 1.11

local = sys.platform == 'win32'
h_res = 960
v_res = 1080

remote_source = [
	"/net/fusi/raid03/yzw/CESM/CESM2.0/FHISH/",
	"/net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST_Contrail/2019_2020/",
	"/net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST_Contrail/2020_COVID/",
    "/home/yzw/contrail/new_emission_run2/",
	"/net/fusi/raid03/yzw/CESM/CESM2.1.1/F2000_Contrail/",
	"/net/fusi/raid03/yzw/CESM/CESM2.1.1/F2000/",
    "/net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST_Contrail/2020_nonCOVID/"]

if local:
    arbitrary_filepath = "Y:\contrails\Ctrl.cam.h0.2005-01.nc"
	# see readme
    os.environ["PROJ_LIB"] = "C:\\Users\\Kyle\\Anaconda3\\Library\\share";
    src = "Y:\contrails"
    dest = ""
else:
    arbitrary_filepath = "/home/kaw/contrails/Ctrl.cam.h0.2005-01.nc"
    src = ""
    dest = "/home/kaw/contrails/average_all_months/"
    dpi = 96

#default_res = [1920, 1080]

# six plots to a page
default_res = [640, 540]
default_lev = 17

g = 9.8

# coordinates for plotting specific regions/windows
locations = {
    0: { 'desc': "Hawaii to LA",
        'lat_min': 15,
        'lat_max': 40,
        'lon_min': 200,
        'lon_max': 250 },
    1: { 'desc': "Eastern US",
        'lat_min': 23,
        'lat_max': 50,
        'lon_min': 275,
        'lon_max': 300 },
    2: { 'desc': "Western Europe",
        'lat_min': 33,
        'lat_max': 57,
        'lon_min': 344,
        'lon_max': 393 },
    3: { 'desc': "Eastern Asia",
        'lat_min': 15,
        'lat_max': 51,
        'lon_min': 445,
        'lon_max': 512 },
}