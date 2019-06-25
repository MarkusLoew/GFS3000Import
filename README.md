[![Build Status](https://travis-ci.org/MarkusLoew/GFS3000Import.svg?branch=master)](https://travis-ci.org/MarkusLoew/GFS3000Import)

# GFS3000Import
R package to import files generated by the WALZ GFS3000 photosynthesis system.

For documentation see the online help:

	help(package = "GFS3000Import")

Installation:

	devtools::install_github("MarkusLoew/GFS3000Import")


## ImportGFS
GFS3000 imports zero points and measurement points as separate data frames and the units of the measured parameters as a vector. The "import" argument can be set to either "all", for everything, or "MP", for measurement points only. Default is "all".

It returns list containing three data frames: 
* "ZP" for Zeropoints, i.e. the calibration data, 
* "MP" for the measurement points, i.e. the hopefully meaningful data, and 
* "Units" with a vector of the units to accompany the data. 

Or, if the "import" argument is set to "MP", just a dataframe with the measurement points MP. Depends on the value for "import".

### Examples:

* ImportGFS("light_response/2019-01-16_Pyc_1.csv", load = "none")

* ImportGFS("light_response/2019-01-16_Pyc_1.csv", load = "all")

* ImportGFS("light_response/2019-01-16_Pyc_1.csv", load = "MP")


## recalcGFS
Recalculates gas exchange parameters based on a new leaf area. Equations from the GFS3000 manual are used. Rounding errors and deviations occur compared to recalculating data with the original GFS3000 software. Use with care. Equations from manual are included in the source code.

