"""
netCDF file cruncher.
"""
import os
from constants import *
import sys
from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import shiftgrid
from mpl_toolkits.basemap import Basemap
from os import listdir
from os.path import isfile, join
from parse import *
from typing import List, Tuple

def searchCommonSubstring(shortest_len: int, l: list) -> int:
	"""
	Given a list of substrings and the length of the shortest substring, find the farthest/greatest index at which a 
	slash occurs before the strings diverge.
	"""
	slash = 0
	for pos in range(0, shortest_len):
		for i in range(0, len(l) - 1):
			if not (l[i][pos] == l[i + 1][pos]):
				# pos is the specific divergence point
				return slash
			elif l[i][pos] == '/':
				slash = pos
	return slash

def removeCommonSubstring(l: list) -> list:
	"""
	Given a list of strings, remove the common substring (beginning from the first element) from each element of the list 
	and return that list. The substrings must all begin and end at the same indices.
	"""
	# length checking
	if len(l) == 1:
		print("The longest common substring is the entirety of the single element of the input list. Returning the input list.")
		return l
	
	# type checking within list + finding the length of the shortest string
	shortest_len = sys.maxsize
	for i in range(0, len(l)):
		if not isinstance(l[i], str):
			print("Each element in the list must be a string. Returning the input list.")
			return l
		elif len(l[i]) < shortest_len:
			shortest_len = len(l[i])
	
	ind = searchCommonSubstring(shortest_len, l)
	return [e[ind:] for e in l]

def basemapPlot(m, lon, lat, var, var_units, min_val=None, max_val=None):
	xi, yi = m((lon), lat)
	
	# Plot Data
	if min_val != None and max_val != None:
		cs = m.pcolormesh(xi,yi,np.squeeze(var),cmap=cm.RdBu, vmin=min_val, vmax=max_val)
	else:
		cs = m.pcolormesh(xi,yi,np.squeeze(var),cmap=cm.RdBu)
	
	# Add Grid Lines
	m.drawparallels(np.arange(-90., 91., 10.), labels=[1,0,0,0])
	m.drawmeridians(np.arange(0., 181. + 180., 20.), labels=[0,0,0,1])
	
	# Add Coastlines, States, and Country Boundaries
	m.drawcoastlines()
	m.drawstates()
	m.drawcountries()
	
	# Add Colorbar
	cbar = m.colorbar(cs, location='bottom', pad="10%")
	cbar.set_label(var_units)

def setupAverage(var_choice: str):
	ds = arbitrary_file
	
	lats = ds.variables['lat'][:]
	lons = ds.variables['lon'][:]
	
	# empty/zeroed dataset - ice is a numpy masked array, not a netcdf4 dataset
	variable = var_choice
	all_avg = ds.variables[variable][:].copy()
	all_avg[:] = 0
	all_avg_units = ds.variables[variable].units
	
	# now, we have an empty-equivalent that retains the same shape. we can start
	# collecting all of the values we need.
	
	rmsrc_truncated = removeCommonSubstring(remote_source)
	
	print("Setup done for variable {0}.".format(var_choice))
	return all_avg, all_avg_units, lats, lons, rmsrc_truncated

# give this an index corresponding to a filepath in remote_source
def plotAverage(datasets, src: str, var: str, lev: int = 16):
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
	
	all_avg = all_avg[time, lev, :, :]
	
	dest = "/home/kaw/contrails/average_all_months/"
	if not os.path.exists(dest):
		os.makedirs(dest)
		
	# take as input the remote source
	# [1:2] corresponds to ['/net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST/']
	# so [2:3] is the other value we're looking for

	if not diff:
		src = remote_source[directory_index]
		files = [f for f in listdir(src) if isfile(join(src, f))]
		# This is a key step that allows for relatively lazy ordering without strict month/year checking.
		# This matters for src = "diff".
		files.sort() 

	else:
		# TODO: MAKE THIS AN ACTUAL GENERALIZABLE APPROACH
		#src, src2 = remote_source[datasets[0]], remote_source[datasets[1]]
		src, src2 = remote_source[1], remote_source[2]
		# files = [f for f in listdir(src) if isfile(join(src, f))]
		# files2 = [f for f in listdir(src2) if isfile(join(src2, f))]
		files = [f for f in listdir(src) if isfile(join(src, f)) and 'h0' in f and '2014' in f]
		files2 = [f for f in listdir(src2) if isfile(join(src2, f)) and 'h0' in f]
		files.sort()
		files2.sort()
	
	print("Averaging begun for variable {0} and level {1}.".format(var, str(lev)))
	# this is technically not the best way to do it. after this, we'll have to iterate
	# through all_avg (lat, lon) again. shouldn't be a huge deal, though.

	# todo: make indexing more efficient (are for loops necessary?)
	manual_access = 0
	if manual_access:
		for x in range(0, len(files)):
			ds = Dataset(src + files[x])
			ds_var = ds.variables[var_choice]
			for lat in range(0, all_avg.shape[0]):
				for lon in range(0, all_avg.shape[1]):
					all_avg[lat, lon] += ds_var[time, lev, lat, lon]
			no_files += 1
	else:
		if not diff:
			for x in range(0, len(files)):
				ds = Dataset(src + files[x])
				ds_var = ds.variables[var_choice]
				all_avg[:, :] += ds_var[time, lev, :, :]
				no_files += 1
		else:
			for x in range(0, len(files)):
				ds = Dataset(src + files[x])
				ds_var = ds.variables[var_choice]
				
				ds2 = Dataset(src2 + files2[x])
				ds_var2 = ds2.variables[var_choice]
				
				all_avg[:, :] += (ds_var2[time, lev, :, :] - ds_var[time, lev, :, :])
				no_files += 1

	for lat in range(0, all_avg.shape[0]):
		for lon in range(0, all_avg.shape[1]):
			all_avg[lat, lon] /= no_files
	print("Done averaging.")

	lon_0 = 0.5 * (lons[0] + lons[-1]) - 180
	# 12/4/20 edits
	#m = Basemap(projection = 'robin', lon_0=lon_0, resolution='l')
	#m = Basemap(projection = 'robin', lon_0=0, resolution='l')

	# quick fix to re-order longitude values - robinson projection
	print("Attempting shifting.")
	all_avg_shifted, lons_shifted = shiftgrid(lon0=180, datain=all_avg, lonsin=lons, start=True)

	plt.figure(figsize =(res[0]/dpi, res[1]/dpi), dpi=dpi)

	# for hawaii to LA
	#lon = lon[200:250]
	#lat = lat[15:40]

	# for Eastern US: lat = [23:50], lon = [275:300]
	# for Western Europe: lat = [33:57], lon = [344:393]
	# for eastern asia: lat = [15:51], lon = [445:512]

	lat_min = 15
	lat_max = 40
	lon_min = 200
	lon_max = 250

	# calculate minimum and maximum within plotted area-- for colorbar limits
	a = np.logical_and(lats >= lat_min, lats < lat_max)
	b = np.logical_and(lons >= lon_min - 180, lons <= lon_max - 180)
	min_val = all_avg_shifted[a, :][:, b].min()
	max_val = all_avg_shifted[a, :][:, b].max()

	m = Basemap(projection='cyl', llcrnrlat=lat_min, urcrnrlat=lat_max, \
            	resolution='h', llcrnrlon=lon_min, urcrnrlon=lon_max)
	# m = Basemap(projection='cyl',llcrnrlat= 15.,urcrnrlat= 40.,\
	# 	resolution='h',  llcrnrlon= 200.,urcrnrlon=250.)

	lon, lat = np.meshgrid(lons_shifted, lats)
	
	basemapPlot(m, lon, lat, all_avg_shifted, all_avg_units, min_val, max_val)

	# add title (based on filepath)
	if not diff:
		plt.title(var_choice + " averaged for " + rmsrc_truncated + ", level = " + str(lev))
		plt.savefig(dest + rmsrc_truncated_file + " " + var_choice + ".png", dpi=dpi)
	if diff:
		plt.title(var_choice + " averaged for diff" + ", level = " + str(lev))
		plt.savefig(dest + "diff " + var_choice + "_lev_" + str(lev) + ".png", dpi=dpi)
	plt.clf()
	
if __name__ == "__main__":
	if local:
		print("Please don't run it like this.")
		sys.exit(0)

	datasets, diff, frac, lev, var, res = handleArguments()

	for v in var:
		for l in range(lev[0], lev[1] + 1):
			plotAverage(datasets, "diff", v, l)