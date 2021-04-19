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
# remote_source = [
# 	"/net/fusi/raid03/yzw/CESM/CESM2.0/FHISH/",
# 	"/net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST/",
# 	"/net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST_Contrail/",
#     "/home/yzw/contrail/new_emission_run2/",
# 	"/net/fusi/raid03/yzw/CESM/CESM2.1.1/F2000_Contrail/",
# 	"/net/fusi/raid03/yzw/CESM/CESM2.1.1/F2000/"]
# remote_source = {
# 	"FHISH": "/net/fusi/raid03/yzw/CESM/CESM2.0/FHISH/",
# 	"FHIST": "/net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST/",
# 	"FHIST_Contrail": "/net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST_Contrail/",
# 	"F2000_Contrail": "/net/fusi/raid03/yzw/CESM/CESM2.1.1/F2000_Contrail/",
# 	"F2000": "/net/fusi/raid03/yzw/CESM/CESM2.1.1/F2000/"}

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

# for hawaii to LA: lat = lat[15:40], lon = [200:250]
# for Eastern US: lat = [23:50], lon = [275:300]
# for Western Europe: lat = [33:57], lon = [344:393]
# for eastern asia: lat = [15:51], lon = [445:512]