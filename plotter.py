"""
netCDF file cruncher.
"""
# TODO: Transition code to Cartopy instead of using Basemap.
import os

local = 0

remote_source = ["/net/fusi/raid03/yzw/CESM/CESM2.0/FHISH/",
		 "/net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST/",
		 "/net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST_Contrail/",
		 "/net/fusi/raid03/yzw/CESM/CESM2.1.1/F2000_Contrail/",
		 "/net/fusi/raid03/yzw/CESM/CESM2.1.1/F2000/"]

# setting PROJ_LIB to espg as per https://stackoverflow.com/a/53751941
# and other settings related to local or remote execution
if local:
	os.environ["PROJ_LIB"] = "C:\\Users\\Kyle\\Anaconda3\\Library\\share";
	src = "Y:\contrails"
	dest = ""
else:
	src = ""
	dest = ""
	dpi = 96

import sys
from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import shiftgrid
from mpl_toolkits.basemap import Basemap
from os import listdir
from os.path import isfile, join
import argparse

# a lot of this is modeled off
# https://joehamman.com/2013/10/12/plotting-netCDF-data-with-Python/

'''
Initialize and return the parser.
'''
def setupParser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(description='NetCDF file cruncher for contrail simulation datasets.')

	# required arguments
	parser.add_argument('datasets', type=str, nargs='+', help='Up to two of: FHIST, FHIST_Contrail, diff, F2000_Contrail, F2000')
	parser.add_argument('-v', '--vars', required=True, type=str, nargs='+', help='At least one of: CLDICE, AREI, FREQI, ICIMR, IWC, QRL')
	
	# optional arguments
	parser.add_argument('-d', '--diff', action='store_const', const=True, help='Analyze two datasets.')
	parser.add_argument('-f', '--frac', action='store_const', const=True, help='Analyze fractional difference. -diff must be supplied.')
	parser.add_argument('-l', '--lev', type=int, help='level (pressure, hPa). Integer from [0, 31].')
	
	return parser

'''
Verify that potentially program-breaking arguments don't break the program.
'''
def validateArguments(n: argparse.Namespace):
	assert True if n.frac == (not n.diff) else False, "If -f is supplied, -d must be supplied."
	assert len(n.datasets) > 2 or (len(n.datasets) == 2 and not n.diff), "Too many datasets supplied."
	assert len(n.datasets) == 2, "Please specify both datasets."

'''
Handle arguments by calling the right functions.
'''
def handleArguments(n: argparse.Namespace):
	pass

# given a list of substrings and the length of the shortest substring, find the
# farthest/greatest index at which a slash occurs before the strings diverge
def searchCommonSubstring(shortest_len: int, l: list) -> int:
	# finally, the actual code
	slash = 0
	for pos in range(0, shortest_len):
		for i in range(0, len(l) - 1):
			# print("Checking", l[i][pos], "and", l[i + 1][pos], "at position", pos)
			if not (l[i][pos] == l[i + 1][pos]):
				# return pos to return the specific divergence point
				return slash
			elif l[i][pos] == '/':
				slash = pos
	return slash

# Given a list of strings, remove the common substring (beginning from the
# first element) from each element of the list and return that list. the
# substrings must all begin and end at the same indices.
def removeCommonSubstring(l: list) -> list:
	# length checking
	if len(l) == 1:
		print("The longest common substring is the entirety of the single\
			   element of the input list. Returning the input list.")
		return l
	
	# type checking within list
	# we can also use this opportunity to find the length of the shortest
	# string in the input list
	shortest_len = sys.maxsize
	for i in range(0, len(l)):
		if not isinstance(l[i], str):
			print("Each element in the list must be a string. Returning the\
				   input list.")
			return l
		elif len(l[i]) < shortest_len:
			shortest_len = len(l[i])
	
	ind = searchCommonSubstring(shortest_len, l)
	return [e[ind:] for e in l]
		

def basemapPlotMerc(m, lon, lat, var, var_units):
	xi, yi = m((lon), lat)
	
	# Plot Data
	cool = cm.cool
	cs = m.pcolormesh(xi,yi,np.squeeze(var),cmap=cool)
	
	# Add Grid Lines
	m.drawparallels(np.arange(-90., 91., 10.), labels=[1,0,0,0])
	m.drawmeridians(np.arange(0., 181. + 180., 10.), labels=[0,0,0,1])
	
	# Add Coastlines, States, and Country Boundaries
	m.drawcoastlines()
	m.drawstates()
	m.drawcountries()
	
	# Add Colorbar
	cbar = m.colorbar(cs, location='bottom', pad="10%")
	cbar.set_label(var_units)

# the key difference between this and the above method is that this doesn't
# have the offset the above has
# interesting mirror issue I ran into while plotting with m.lon_0 > 0:
# https://github.com/matplotlib/basemap/issues/463
def basemapPlotRob(m, lon, lat, var, var_units, frac_clim = False):
	xi, yi = m((lon), lat)
	
	# Plot Data
	cs = m.pcolormesh(xi,yi,np.squeeze(var),cmap=cm.RdBu)
	#cs = m.pcolormesh(xi,yi,np.squeeze(var))
	print(var)
	#print(cs)
	
	# Add Grid Lines
	m.drawparallels(np.arange(-90., 91., 10.), labels=[1,0,0,0])
#	m.drawmeridians(np.arange(0., 181. + 180., 20.), labels=[0,0,0,1],
#				 fmt=alternatingLongitudeLabels)
	m.drawmeridians(np.arange(0., 181. + 180., 20.), labels=[0,0,0,1])
	
	# Add Coastlines, States, and Country Boundaries
	m.drawcoastlines()
	m.drawstates()
	m.drawcountries()
	
	# Add Colorbar
	#cbar = m.colorbar(cs, location='bottom', pad="10%")
	#cbar.set_label(var_units)
	#if frac_clim:
	#	cbar.clim(-1, 1)

def plot():
	if local: # single file
		ds = Dataset("Y:\contrails\Ctrl.cam.h0.2005-01.nc")
		lats = ds.variables['lat'][:]
		lons = ds.variables['lon'][:]
		var = ds.variables['TGCLDIWP'][:] # total grid-box cloud ice water path
		var_units = ds.variables['TGCLDIWP'].units
		
		# set up stereographic projection
		lon_0 = lons.mean()
		lat_0 = lats.mean()
		
		lon, lat = np.meshgrid(lons, lats)
		
		m = Basemap(width=40000000,height=20000000,
		            resolution='l',projection='stere',
		            lat_ts=40,lat_0=lat_0,lon_0=lon_0)
		
		m1 = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
		             llcrnrlon=0,urcrnrlon=360,lat_ts=20,resolution='l')
		basemapPlotMerc(m1, lon, lat, var, var_units)
	else:
		m = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,    
				    llcrnrlon=0,urcrnrlon=360,lat_ts=20,resolution='l')
		for src in remote_source:
			files = [f for f in listdir(src) if isfile(join(src, f))]
			dest = "/home/kaw/contrails" + src
			if not os.path.exists(dest):
				os.makedirs(dest)
			plt.figure(figsize =(1920/dpi, 1080/dpi), dpi=dpi)
			for x in range(0, len(files)):
				ds = Dataset(src + files[x])
				lats = ds.variables['lat'][:]
				lons = ds.variables['lon'][:]
				var = ds.variables['TGCLDIWP'][:] # total grid-box cloud ice
												  # water path
				var_units = ds.variables['TGCLDIWP'].units
				
				# set up stereographic projection
				lon_0 = lons.mean()
				lat_0 = lats.mean()
				
				lon, lat = np.meshgrid(lons, lats)
				
				basemapPlotMerc(m, lon, lat, var, var_units)
				# add title (based on filepath)
				plt.title('total grid-box cloud ice water path - ' + files[x][:-3])
				plt.savefig(dest + files[x][:-3] + ".png", dpi=dpi)
				plt.clf()

def setupAverage(var_choice: str, base_file="/home/kaw/contrails/Ctrl.cam.h0.2005-01.nc"):
	# arbitrary netcdf file
	
	# todo: deprecate "local"
	if not local:
		ds = Dataset(base_file)
	else:
		ds = Dataset("Y:\contrails\Ctrl.cam.h0.2005-01.nc")
	
	lats = ds.variables['lat'][:]
	lons = ds.variables['lon'][:]
	
	# empty dataset - ice is a numpy masked array, not a netcdf4 dataset
	variable = var_choice
	all_avg = ds.variables[variable][:].copy()
	all_avg[:] = 0 # zero out values to provide starting point
	all_avg_units = ds.variables[variable].units
	
	# now, we have an empty-equivalent that retains the same shape. we can start
	# collecting all of the values we need.
	
	rmsrc_truncated = removeCommonSubstring(remote_source)
	
	print("Setup done.")
	return all_avg, all_avg_units, lats, lons, rmsrc_truncated

# give this an index corresponding to a filepath in remote_source
def plotAverage(src: str, var: str, lev: int = 16):
	var_choice = var
	diff = False
	if src == "diff":
		diff = True
	
	all_avg, all_avg_units, lats, lons, rmsrc_truncated = setupAverage(var_choice)
	
	# todo: do this more efficiently
	if diff:
		rmsrc_truncated_file = src
	else:
		if src == "FHIST":
			directory_index = 1
		elif src == "FHIST_Contrail":
			directory_index = 2
		# expand functionality here
		elif src == "F2000_Contrail":
			directory_index = 3
		elif src == "F2000":
			directory_index = 4
		
		rmsrc_truncated = rmsrc_truncated[directory_index]
		rmsrc_truncated_file = rmsrc_truncated.replace('/','')
	
	# else, the following code won't work
	assert (local == 0), "can't average values when only working with one file"
	
	time = 0
	no_files = 0
	
	all_avg = all_avg[time][lev][:][:]
	
	dest = "/home/kaw/contrails/average_all_months/"
	if not os.path.exists(dest):
		os.makedirs(dest)
		
	# take as input the remote source
	# [1:2] corresponds to ['/net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST/']
	# so [2:3] is the other value we're looking for

	if not diff:
		src = remote_source[directory_index]
		files = [f for f in listdir(src) if isfile(join(src, f))]
		files.sort() # This is a key step that allows for relatively lazy ordering without strict
					 # month/year checking. This matters for src = "diff".

	else:
		# TODO: MAKE THIS AN ACTUAL GENERALIZABLE APPROACH
		src, src2 = remote_source[1], remote_source[2]
		#src, src2 = remote_source[3], remote_source[4]
		# files = [f for f in listdir(src) if isfile(join(src, f))]
		# files2 = [f for f in listdir(src2) if isfile(join(src2, f))]
		files = [f for f in listdir(src) if isfile(join(src, f)) and 'h0' in f and '2014' in f]
		files2 = [f for f in listdir(src2) if isfile(join(src2, f)) and 'h0' in f]
		files.sort()
		files2.sort()
		#print(files)
		#print()
		#print(files2)
	
	print("Averaging begun.")
	# todo: this is technically not the best way to do it. after this, we'll have to iterate
	# through all_avg (lat, lon) again. shouldn't be a huge deal, though.
	manual_access = 0
	if manual_access:
		for x in range(0, len(files)):
			ds = Dataset(src + files[x])
			ds_var = ds.variables[var_choice]
			for lat in range(0, all_avg.shape[0]):
				for lon in range(0, all_avg.shape[1]):
					all_avg[lat][lon] += ds_var[time][lev][lat][lon]
			no_files += 1
	else:
		if not diff:
			for x in range(0, len(files)):
				ds = Dataset(src + files[x])
				ds_var = ds.variables[var_choice]
				all_avg[:][:] += ds_var[time][lev][:][:]
				no_files += 1
		else: # diff time
			for x in range(0, len(files)):
				ds = Dataset(src + files[x])
				ds_var = ds.variables[var_choice]
				
				ds2 = Dataset(src2 + files2[x])
				#print(src2 + files2[x])
				ds_var2 = ds2.variables[var_choice]
				
				all_avg[:][:] += \
					(ds_var2[time][lev][:][:] - ds_var[time][lev][:][:])
				no_files += 1

	for lat in range(0, all_avg.shape[0]):
		for lon in range(0, all_avg.shape[1]):
			all_avg[lat][lon] /= no_files
	print("Done averaging.")
	# plottime
	
	# this neat trick is courtesy of the example usage for the robinson
	# projection from the basemap api
	# for the robinson projection, when this > 0, mirror image results

	# 12/4/20 - removing the 180 degree shift from the longitudes

	lon_0 = 0.5 * (lons[0] + lons[-1]) - 180
	# 12/4/20 edits
	m = Basemap(projection = 'robin', lon_0=lon_0, resolution='l')
	#m = Basemap(projection = 'robin', lon_0=0, resolution='l')

	# for hawaii to LA
	#m = Basemap(projection='merc',llcrnrlat=13.66,urcrnrlat=41,llcrnrlon=-165,urcrnrlon=-106.25,resolution='h')
	# lats: 110:140, 13.66492147 to 40.9947644
	# lons: 12:60, -165 to -106.25 (before shiftgrid and -360)
	# var[0][12:60, 110:140]

	# quick fix to re-order longitude values
	# print("Attempting shifting.")
	all_avg_shifted, lons_shifted = shiftgrid(lon0=180, datain=all_avg, lonsin=lons, start=True)

	plt.figure(figsize =(1920/dpi, 1080/dpi), dpi=dpi)

	# for hawaii to LA
	#lon = lon[12:60]
	#lat = lat[110:140]
	#all_avg_shifted = all_avg_shifted[110:140, 12:60]

	#m = Basemap(projection='cyl',llcrnrlat= 15.,urcrnrlat= 40.,\
    #        	resolution='l',  llcrnrlon= 200.,urcrnrlon=250.)
	#map.contourf()

	#lon, lat = np.meshgrid(lons, lats)
	lon, lat = np.meshgrid(lons_shifted, lats)
	
	basemapPlotRob(m, lon, lat, all_avg_shifted, all_avg_units)
	# add title (based on filepath)
	if not diff:
		plt.title(var_choice + " averaged for " + rmsrc_truncated + ", level = " + str(lev))
		plt.savefig(dest + rmsrc_truncated_file + " " + var_choice + ".png", dpi=dpi)
	if diff:
		plt.title(var_choice + " averaged for diff" + ", level = " + str(lev))
		plt.savefig(dest + "diff " + var_choice + "_lev_" + str(lev) + ".png", dpi=dpi)
	plt.clf()
    
# # 1) You should do next is to calculate the fractional changes of each properties for level 18 only (~330 hPa). 
# The fractional change is defined as the absolute difference (you are plotting now) divided by the 10-year mean value (use no-contrail run) for each grid. 
# # 2) Please also calculate the standard deviation first for each grid point, and then divide the differences by the standard deviation.
# # 3) Please check if you can do the student-t test for those two group using Python. I cc Tianhao here who has lots of experience in performing statistical analysis using Python. 
# After you learn how to do it, perform t test at each grid point, and for the absolute difference plots, only plot the grids with p-value less than 0.05.

# okay, we're going to simplify plotAverage() into two smaller methods-- one for diff and one for single

# no_contrail_run = true -> use F2000 data
# else -> use F2000_Contrail data
# var <-> property
# At Yuan's request, this examines the 18th level by default.
# Also note: this is basically the same as plotAverage without diff. However, this also returns values instead of just plotting.
def mean(var: str, lev: int, src=remote_source[4], time=0):
	all_avg, all_avg_units, lats, lons, rmsrc_truncated = setupAverage(var)
	
	files = [f for f in listdir(src) if isfile(join(src, f))]
	files.sort()
	
	no_files = 0
	for x in range(0, len(files)):
		ds = Dataset(src + files[x])
		ds_var = ds.variables[var]
		all_avg[:][:] += ds_var[time][lev][:][:]
		no_files += 1

	for lat in range(0, all_avg.shape[0]):
		for lon in range(0, all_avg.shape[1]):
			all_avg[lat][lon] /= no_files
	
	return all_avg, all_avg_units, lats, lons, rmsrc_truncated

def difference(var: str, lev: int, src=(remote_source[2],remote_source[1]), time=0):
	all_avg, all_avg_units, lats, lons, rmsrc_truncated = setupAverage(var)
	
	files = ([f for f in listdir(src[0]) if isfile(join(src[0], f))],
		     [f for f in listdir(src[1]) if isfile(join(src[1], f))])
	files[0].sort()
	files[1].sort()
	
	no_files = 0
	for x in range(0, len(files)):
		ds = (Dataset(src[0] + files[0][x]), Dataset(src[1] + files[1][x]))
		ds_var = (ds[0].variables[var], ds[1].variables[var])
	
		all_avg[:][:] += (ds_var[1][time][lev][:][:] - ds_var[0][time][lev][:][:])
		no_files += 1

	return all_avg, all_avg_units, lats, lons, rmsrc_truncated

def fractionalChange(var: str, lev=17, contrail_filepath=remote_source[3], no_contrail_filepath=remote_source[4]):
	print("Finding fractional change for", var)
	mean_all_avg, all_avg_units, lats, lons, rmsrc_truncated = mean(var, lev, src=no_contrail_filepath)
	difference_all_avg, _, _, _, _ = difference(var, lev, src=(contrail_filepath, no_contrail_filepath))
	
	frac = mean_all_avg / difference_all_avg
	print("Average calculated.")
	
	# todo: move this to a distinct plotting method that takes data as input
	lon_0 = 0.5 * (lons[0] + lons[-1]) - 180
	m = Basemap(projection = 'robin', lon_0=lon_0, resolution='l')
	print("Attempting shifting.")
	all_avg_shifted, lons_shifted = shiftgrid(lon0=180, datain=frac[0][lev], lonsin=lons, start=True)

	plt.figure(figsize =(1920/dpi, 1080/dpi), dpi=dpi)
	lon, lat = np.meshgrid(lons_shifted, lats)
	basemapPlotRob(m, lon, lat, all_avg_shifted, all_avg_units, frac_clim=True)

	plt.title(var + " - fractional change" + ", level = " + str(lev))
	plt.savefig(dest + "frac " + var + "_lev_" + str(lev) + ".png", dpi=dpi)
	plt.clf()

# TODO: Actual proper input validation w/ argparse
	
if __name__ == "__main__":
	if local:
		print("Please don't run it like this.")
		sys.exit(0)
	
	if len(sys.argv) < 3:
		print("\
Usage: python3 plotter.py <operation: str> <variable: str> [-arg argument]\n\
- operation: One of: FHIST, FHIST_Contrail, diff, F2000_Contrail, F2000\n\
- variable : One of: CLDICE, AREI, FREQI, ICIMR, IWC, QRL, all\n\
args:\n\
 -lev : integer in [0, 31]\n\
 -frac: ")
		sys.exit(0)
	
	try:
		arg1, arg2 = str(sys.argv[1]), str(sys.argv[2])
	except:
		print("The two parameters must be valid strings.")
		sys.exit(0)
	
	# frac implemented as stopgap measure. please complete setupParser()
	ops = ["FHIST", "FHIST_Contrail", "diff", "F2000_Contrail", "F2000"]
	#ds_vars = ["CLDICE", "AREI", "FREQI", "ICIMR", "IWC", "QRL"]
	ds_vars = ["AREI", "CLOUD", "IWC", "T"]
	
	if arg1 == "frac":
		for arg in ds_vars:
			fractionalChange(arg)
			print()
		sys.exit(0)
	
	#if (arg1 not in ops) or (arg2 not in ds_vars + ["all"]):
	if False:
		print("Check usage statement for valid inputs.")
		sys.exit(0)
		
	arg4 = 16 # level
	
	if arg2 == "all":
		if sys.argv[3] == "-lev":
			try:
				arg4 = int(sys.argv[4])
			except:
				print("-lev must be an integer from [0, 31].")
				sys.exit(0)
		for arg in ds_vars:
			plotAverage(arg1, arg, lev=arg4)
	else:
		if sys.argv[3] == "-lev":
			try:
				arg4 = int(sys.argv[4])
			except:
				print("-lev must be an integer from [0, 31].")
				sys.exit(0)
		plotAverage(arg1, arg2, lev=arg4)