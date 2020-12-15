import sys
import os

local = sys.platform == 'win32'

remote_source = [
	"/net/fusi/raid03/yzw/CESM/CESM2.0/FHISH/",
	"/net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST/",
	"/net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST_Contrail/",
	"/net/fusi/raid03/yzw/CESM/CESM2.1.1/F2000_Contrail/",
	"/net/fusi/raid03/yzw/CESM/CESM2.1.1/F2000/"]
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
    dest = ""
    dpi = 96