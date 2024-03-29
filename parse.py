import argparse
import sys
from typing import List, Tuple
from constants import arbitrary_filepath, default_res, default_lev, version, locations
from netCDF4 import Dataset

arbitrary_file = Dataset(arbitrary_filepath)

def setupParser() -> argparse.ArgumentParser:
	""" Initialize and return the parser.

	Returns:
	    The argument parser object.
	    
	"""
	parser = argparse.ArgumentParser(prog='NetCDF File Cruncher', description='NetCDF file cruncher for contrail simulation datasets.')

	# required arguments
	parser.add_argument('datasets', type=str, nargs='+', help='Up to two of: FHIST, FHIST_Contrail, F2000_Contrail, F2000')
	parser.add_argument('-v', '--vars', required=True, type=str, nargs='+', help='At least one of: CLDICE, AREI, FREQI, ICIMR, IWC, QRL')
	
	# optional arguments
	parser.add_argument('-f', '--frac', action='store_true', help='Analyze fractional difference. -diff must be supplied.')
	parser.add_argument('-l', '--lev', type=int, nargs='+', default=[0, 31], \
		help='level (pressure, hPa). Integer from [0, 31]. Supply two ints in ascending order for an inclusive range of plots by level.')
	parser.add_argument('-r', '--res', type=int, nargs=2, default=default_res, help='Resolution (horizontal, vertical). Defaults to 1920 x 1080.')
	parser.add_argument('-i', '--integration', action='store_true', help='Short-circuit normal decision logic and integrate IWC values.')
	parser.add_argument('--region', type=int, choices=list(locations.keys()), \
		help='Plot a specific area. See constants.py for areas.') # perhaps this should be a required argument

	# miscellaneous
	parser.add_argument('--version', action='version', version='%(prog)s {0}'.format(version))

	return parser

def validateArguments(n: argparse.Namespace):
	""" Verify that potentially program-breaking arguments don't break the program-- 
	check interactions that argparse can't catch and/or specific values. 

	Args:
	    n: The namespace of the parser after it has parsed arguments.
	
	"""
	def check(condition, message):
		if not condition:
			raise argparse.ArgumentTypeError(message)

	check(len(n.datasets) == 2 if n.frac else True, "Two datasets necessary for fractional difference analysis.")
	check(len(n.datasets) <= 2, "Too many datasets supplied.")
	check(len(n.lev) <= 2, "Please specify one or two levels.")
	check(len(n.lev) == 1 or n.lev[1] > n.lev[0], "The level bounds must be in ascending order.")
	check(n.res[0] > 0 and n.res[1] > 0, "Resolution dimensions must be positive.")
	check(n.lev[0] > 0 if n.integration else True, "Can't integrate level 0 because level -1 doesn't exist.")
	check(all(key in arbitrary_file.variables.keys() for key in n.vars), "Invalid variable names.")

def handleArguments() -> Tuple[List[str], bool, List[int], List[str], List[int], bool, int]:
	""" Overall argument handler.

	Returns:
		name.datasets: A list of dataset filepaths.
		name.frac: A boolean indicating whether fractional difference is to be calculated.
		name.lev: A list describing either a specific level to examine (if length 1) or a range thereof (if length 2).
		name.vars: A list of variables to analyze.
		name.res: A list indicating the resolution of the output graphs (width by height).
		name.integration: A boolean indicating whether a variabie is to be integrated.
		name.region: An int corresponding to a specific area to plot.

	"""
	parser = setupParser()
	name = parser.parse_args(sys.argv[1:])
	validateArguments(name)
	return name.datasets, name.frac, name.lev, name.vars, name.res, name.integration, name.region