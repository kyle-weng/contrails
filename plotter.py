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
	""" Given a list of substrings and the length of the shortest substring, find the farthest/greatest index at which a 
	slash occurs before the strings diverge. """
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
	""" Given a list of strings, remove the common substring (beginning from the first element) from each element of the list 
	and return that list. The substrings must all begin and end at the same indices. """
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
		cs = m.pcolormesh(xi,yi,np.squeeze(var),cmap=cm.RdBu_r, vmin=min_val, vmax=max_val)
	else:
		cs = m.pcolormesh(xi,yi,np.squeeze(var),cmap=cm.RdBu_r)
	
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

def average(datasets: List[str], var: str, lev: List[int], region: int, month_year_constraint: str = None):
	difference = len(datasets) == 2

	src_a, src_b, lev_lower, lev_upper, all_avg, all_avg_units, lats, lons, levs, files, files2, time = \
		setup_temp(difference, datasets, var, lev, month_year_constraint)

	# summing
	if difference:
		for x in range(0, len(files)):
			print("File no. {0}, filepath: {1}".format(x, files[x]))
			for i in range(lev_lower, lev_upper):
				ds_a = Dataset(src_a + files[x]).variables[var]
				ds_b = Dataset(src_b + files2[x]).variables[var]
				term = (ds_b[time, i, :, :] - ds_a[time, i, :, :])
				all_avg[:, :] += term
	else:
		for x in range(0, len(files)):
			print("File no. {0}, filepath: {1}".format(x, files[x]))
			for i in range(lev_lower, lev_upper):
				ds_a = Dataset(src_a + files[x]).variables[var]
				term = (ds_a[time, i, :, :])
				all_avg[:, :] += term

	# averaging - check: is there a better way of doing this?
	for lat in range(0, all_avg.shape[0]):
		for lon in range(0, all_avg.shape[1]):
			all_avg[lat, lon] /= len(levs) * len(files)

	print("Done averaging.")
			

	lon_0 = 0.5 * (lons[0] + lons[-1]) - 180

	print("Attempting shifting.")
	all_avg_shifted, lons_shifted = shiftgrid(lon0=180, datain=all_avg, lonsin=lons, start=True)

	location = locations[region]
	lat_min = location['lat_min']
	lat_max = location['lat_max']
	lon_min = location['lon_min']
	lon_max = location['lon_max']

	# calculate minimum and maximum within plotted area-- for colorbar limits
	a = np.logical_and(lats >= lat_min, lats < lat_max)
	b = np.logical_and(lons >= lon_min - 180, lons <= lon_max - 180)
	min_val = all_avg_shifted[a, :][:, b].min()
	max_val = all_avg_shifted[a, :][:, b].max()

	m = Basemap(projection='cyl', llcrnrlat=lat_min, urcrnrlat=lat_max, \
            	resolution='h', llcrnrlon=lon_min, urcrnrlon=lon_max)

	lon, lat = np.meshgrid(lons_shifted, lats)

	basemapPlot(m, lon, lat, all_avg_shifted, all_avg_units, min_val=min_val, max_val=max_val)
	plt.title(var + " average" + (difference * " difference") + ", " + month_year_constraint)
	plt.savefig(dest + var + " average" + (difference * " difference") + ", " + month_year_constraint + ".png", dpi=dpi)
	plt.clf()

def setup_temp(difference: bool, datasets: List[str], var: str, lev: List[int], month_year_constraint: str=None):
	h0_constraint = 'h0'

	if difference:
		# python3 plotter.py /filepath2 /filepath1 -> filepath2 - filepath1
		src_b = datasets[0]
		src_a = datasets[1]
	else:
		src_a = datasets[0]

	# determining bounds on level iteration-- single value or range?
	# [a] -> level a, [a, b] -> levels a through b, inclusive
	lev_lower = lev[0]
	lev_upper = lev[1] + 1 if len(lev) == 2 else lev[0] + 1 # mmm a crunchy ternary operator
	
	ds = arbitrary_file
	
	lats = ds.variables['lat'][:]
	lons = ds.variables['lon'][:]
	levs = ds.variables['lev'][:]
	
	# empty/zeroed dataset - ice is a numpy masked array, not a netcdf4 dataset
	all_avg = ds.variables[var][:].copy()
	all_avg[:] = 0
	all_avg_units = ds.variables[var].units
	
	# now, we have an empty-equivalent that retains the same shape. we can start
	# collecting all of the values we need.
	
	rmsrc_truncated = removeCommonSubstring(remote_source)
	
	print("Setup done for variable {0}.".format(var))

	time = 0
	
	# note that 'lev' as an input variable here is functionally useless. this sums over all levs, but
	# we'll just stick everything in a single level
	all_avg = all_avg[time, 16, :, :]
	
	if not os.path.exists(dest):
		os.makedirs(dest)

	# this code is pretty spaghetti. there's probably a better way to rewrite this. you can remove a bit
	# of the redundancy (files = ...), but it'll make the code that much more unreadable.
	if difference:
		if not month_year_constraint:
			files = [f for f in listdir(src_a) if isfile(join(src_a, f)) and h0_constraint in f]
			files2 = [f for f in listdir(src_b) if isfile(join(src_b, f)) and h0_constraint in f]
		else:
			print("Month/year constraint is {0}".format(month_year_constraint))
			files = [f for f in listdir(src_a) if isfile(join(src_a, f)) and h0_constraint in f and month_year_constraint in f]
			files2 = [f for f in listdir(src_b) if isfile(join(src_b, f)) and h0_constraint in f and month_year_constraint in f]
		files.sort()
		files2.sort()
		assert len(files) == len(files2), "can't compute difference between uneven file sets"
	else:
		if not month_year_constraint:
			files = [f for f in listdir(src_a) if isfile(join(src_a, f)) and h0_constraint in f]
		else:
			files = [f for f in listdir(src_a) if isfile(join(src_a, f)) and h0_constraint in f and month_year_constraint in f]
		files.sort()
	
	# ugliest return statement I've ever written
	return src_a, src_b if difference else _, lev_lower, lev_upper, all_avg, all_avg_units, \
		lats, lons, levs, files, files2 if difference else _, time

# goal: rewrite integrate() but allow for single files as well as differences (double files)
def integrate(datasets: List[str], var: str, lev: List[int], region: int, month_year_constraint: str = None):
	difference = len(datasets) == 2

	src_a, src_b, lev_lower, lev_upper, all_avg, all_avg_units, lats, lons, levs, files, files2, time = \
		setup_temp(difference, datasets, var, lev, month_year_constraint)

	# IWC integration -> ice water path
	if var == "IWC":
		if difference:
			for x in range(0, len(files)):
				print("File no. {0}, filepath: {1}".format(x, files[x]))
				for i in range(lev_lower, lev_upper):
					print("\tAdding: level {0}\t{1} hPa".format(i, levs[i]))
					ds_a = Dataset(src_a + files[x]).variables[var]
					ds_b = Dataset(src_b + files2[x]).variables[var]
					
					delta_p = levs[i] - levs[i - 1]
					adjustment = -1 * delta_p / g
					term = (ds_b[time, i, :, :] - ds_a[time, i, :, :]) * adjustment
					all_avg[:, :] += term
		else:
			# single file time
			for x in range(0, len(files)):
				print("File no. {0}, filepath: {1}".format(x, files[x]))
				for i in range(lev_lower, lev_upper):
					print("\tAdding: level {0}\t{1} hPa".format(i, levs[i]))
					ds_a = Dataset(src_a + files[x]).variables[var]
					
					delta_p = levs[i] - levs[i - 1]
					adjustment = -1 * delta_p / g
					term = (ds_a[time, i, :, :]) * adjustment
					all_avg[:, :] += term
	else:
		print("integration undefined for non-IWC variables.")
		sys.exit(0) # throw an error instead later

	lon_0 = 0.5 * (lons[0] + lons[-1]) - 180

	print("Attempting shifting.")
	all_avg_shifted, lons_shifted = shiftgrid(lon0=180, datain=all_avg, lonsin=lons, start=True)

	# coordinate selection
	location = locations[region]
	lat_min = location['lat_min']
	lat_max = location['lat_max']
	lon_min = location['lon_min']
	lon_max = location['lon_max']

	# calculate minimum and maximum within plotted area-- for colorbar limits
	a = np.logical_and(lats >= lat_min, lats < lat_max)
	b = np.logical_and(lons >= lon_min - 180, lons <= lon_max - 180)
	min_val = all_avg_shifted[a, :][:, b].min()
	max_val = all_avg_shifted[a, :][:, b].max()

	m = Basemap(projection='cyl', llcrnrlat=lat_min, urcrnrlat=lat_max, \
            	resolution='h', llcrnrlon=lon_min, urcrnrlon=lon_max)

	lon, lat = np.meshgrid(lons_shifted, lats)
	if var == "IWC":
		all_avg_units = 'g/(m^2)'
	basemapPlot(m, lon, lat, all_avg_shifted, all_avg_units, min_val=-5e-6, max_val=10e-6)
	plt.title(var + " integration" + (difference * " difference") + ", " + month_year_constraint)
	plt.savefig(dest + var + " integration" + (difference * " difference") + ", " + month_year_constraint + ".png", dpi=dpi)

	# add title (based on filepath)
	plt.clf()

def fraction(datasets: List[str], var: str, lev: List[int], region: int, month_year_constraint: str = None):
	print("test")
	
if __name__ == "__main__":
	if local:
		print("Please don't run it like this.")
		sys.exit(0)

	datasets, action, lev, vars, res, region = handleArguments()

	plt.figure(figsize=(res[0]/dpi, res[1]/dpi), dpi=dpi)

	# sample script call:
	# python3 plotter.py /net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST_Contrail/2020_COVID/ /net/fusi/raid03/yzw/CESM/CESM2.1.1/FHIST_Contrail/2020_nonCOVID/ -v IWC -l 10 15 -i --region 0
	for v in vars:
		for i in range (1, 6 + 1): # month
			action_choices[action](datasets, v, lev=lev, region=region, month_year_constraint="2020-0{0}".format(i))
	plt.close()