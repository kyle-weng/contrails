# contrails
Some stuff I worked on during summer and fall 2020. Storing information here in an attempt to declutter the actual code. Will update later.
Skeleton code taken from here: [Joe Hamman](https://joehamman.com/2013/10/12/plotting-netCDF-data-with-Python/)  
### Things worth noting
* When working locally, set `PROJ_LIB` to espg's location: [Stackoverflow](https://stackoverflow.com/a/53751941 "Hello!")
* Plots are mirrored when working with Robinson projections and positive `lon_0` values: [Github issue](https://github.com/matplotlib/basemap/issues/463 "Hi!").  
### Todo
* Transition code to Cartopy and abandon Basemap.
* Implement proper arguments (argparse).
* Move setup functions to other files as necessary.
