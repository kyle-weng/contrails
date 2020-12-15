def plot():
	if local:
		# single file
		ds = Dataset(arbitrary_filepath)
		lats = ds.variables['lat'][:]
		lons = ds.variables['lon'][:]
		var = ds.variables['TGCLDIWP'][:]
		var_units = ds.variables['TGCLDIWP'].units
		
		# set up stereographic projection
		lon_0 = lons.mean()
		lat_0 = lats.mean()
		
		lon, lat = np.meshgrid(lons, lats)
		m = Basemap(width=40000000, height=20000000, resolution='l', projection='stere', \
			lat_ts=40, lat_0=lat_0, lon_0=lon_0)
		m1 = Basemap(projection='merc', llcrnrlat=-80, urcrnrlat=80, llcrnrlon=0, \
			urcrnrlon=360, lat_ts=20, resolution='l')
		basemapPlot(m1, lon, lat, var, var_units)
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
				var = ds.variables['TGCLDIWP'][:]
				var_units = ds.variables['TGCLDIWP'].units
				
				# set up stereographic projection
				lon_0 = lons.mean()
				lat_0 = lats.mean()
				lon, lat = np.meshgrid(lons, lats)
				basemapPlot(m, lon, lat, var, var_units)

				# add title (based on filepath)
				plt.title('total grid-box cloud ice water path - ' + files[x][:-3])
				plt.savefig(dest + files[x][:-3] + ".png", dpi=dpi)
				plt.clf()

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
	basemapPlot(m, lon, lat, all_avg_shifted, all_avg_units, frac_clim=True)

	plt.title(var + " - fractional change" + ", level = " + str(lev))
	plt.savefig(dest + "frac " + var + "_lev_" + str(lev) + ".png", dpi=dpi)
	plt.clf()