# photomass_ls.py


## About the script

photomass_ls computes an estimate of the stellar mass for a given galaxy based on a photometry performed by GALFIT on g & r band data from DESI Legacy Imaging Surveys. The script automatically downloads images of the galaxy from the Legacy Surveys (DR10) database, creates image masks, generates GALFIT
input files with well-assessed initial values, performs the GALFIT photometry, and calculates the stellar mass estimate. More details are given in Ebrová, Bílek, & Eliášek (2025).


## Citation

If you use this code in your research or projects, please cite the following accompanied paper:

Ebrová, Bílek, & Eliášek: "Photometric stellar masses for galaxies in DESI Legacy Imaging
Surveys" (2025; accepted for publication in Astronomy & Astrophysics; arXiv:2510.02257); 
https://ui.adsabs.harvard.edu/abs/2025arXiv251002257E/abstract 


## Requirements

GALFIT : https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html

python libraries:
requests, urllib, numpy, scipy, SEP, astropy, [matplotlib (required for use of option `--plot`)]


You can install required packages into virtual environment:
```
python3 -m venv venv
. venv/bin/activate
pip3 install "numpy<2.0" scipy SEP astropy requests
```

> SEP 1.2.1 does not work with numpy 2.0 or newer


## Limitations
The script uses [NED](ned.ipac.caltech.edu) and [`get_icrs_coordinates` from `astropy`](cds.unistra.fr) to get galactic extinctions and coordinates for the given galaxy.
If the object is not present in these databases the corresponding exception will be raised.
You can use custom coordinates with `--coordinates` option to circumvent this limitation.
K correction is limited to redshift < 0.5 as well as other limitations of the [K-correction calculator](http://kcor.sai.msu.ru).


## Running the script
You can run the script to get list of options with descriptions like this:
```
$ ./photomass_ls.py -h
usage: photomass_ls.py [-h] [--refetch] [--galpath GALPATH] [--dist DIST] [--coordinates COORDINATES] [--additional-filters ADDITIONAL_FILTERS] [--plot] [--K] [--redshift REDSHIFT] OBJECT radius

computes a stellar mass estimate for a given galaxy based on GALFIT photometry with Legacy Surveys data

positional arguments:
  OBJECT                name of the galaxy for which the stellar mass estimate is done
  radius                approximate galaxy angular size in arcseconds (e.g. semi-major axis at SB=25.5mag/arcsec^2 in the Spitzer 3.6 micron filter)

options:
  -h, --help            show this help message and exit
  --refetch             re-download and re-analyse image even if a local copy exists
  --galpath GALPATH     path to galfit64 binary; if not set, $PATH is searched
  --dist DIST           galaxy distance in Mpc
  --coordinates COORDINATES
                        comma seperated coordinates: RA,DEC with unints (default units are deg); requires --dist
  --additional-filters ADDITIONAL_FILTERS
                        other filters without delimiter (e.g. 'zi')
  --plot                plot png image of source data and masked data
  --K, -K               apply K correction
  --redshift REDSHIFT   galaxy redshift - overrides the one dowloaded from NED
```

The numerical outpus are printed to stdout and to `OBJECT_photomass.out`. For example:
```
$./photomass_ls.py NGC474 106.3 --dist 30.88
Downloading https://www.legacysurvey.org/viewer/fits-cutout?ra=20.02786271688&dec=3.41551721475&height=1024&width=1024&layer=ls-dr10&pixscale=0.524&bands=gr
|===========================================================================================================================| 8.3M/8.3M (100.00%)         2s
Galaxy: NGC474
RA: 20.0279 Dec: 3.416
log10(M*[Sun]): E: 10.79 M: 10.84 Q: 10.88
Ext[mag]: g : 0.112 r : 0.075
Mag[mag]: g : -20.5 r : -21.29
Sersic index: g : 3.049 r : 3.872
R_e[px]: g : 44.73 r : 48.04
R_e[arcsec]: g : 23.44 r : 25.17
Axis ratio: g : 0.8185 r : 0.8188
Position angle[deg]: g : 21.18 r : 19.1
Distance[Mpc]: 30.88 (from input)
Plate scale[arcsec / px]: 0.524 0.524
Zero point[mag]: 20.99
Redshift: 0.007722
```
where:
First two lines are outputed to `stderr`.

`Galaxy`: `OBJECT` (converted to upper case if no coordinate are specified and stript of characters " []:")

`RA`, `Dec`: coordinates of the `OBJECT` in deg

`log(M*[Msun])`: the estimate of the logarithm of the galaxy stellar mass computed by Eq(1) of Ebrova et al. (2025), where Msun is the sollar masses, using calibrations - E: Eskew et al. (2012), M: Meidt et al. (2014), Q: Querejeta et al. (2015)

The next set of variables is calculated for different filters:
 - `Ext`: Galactic extinctions in DES filters from Schlafly & Finkbeiner (2011) using NASA/IPAC Extragalactic Database (NED)
 - `Mag`: extinction-correlated absolute magnitude calculated from the galaxy distance and the GLAFIT Integrated magnitude (includes K-correction with option --K)
 - `Sersic index`: Sersic index of the GALFIT model
 - `R_e`: effective radius of the GALFIT model - in pixels and arcsec
 - `Axis ratio`: axis ratio of the GALFIT model
 - `Position angle`: positioning angle of the GALFIT model (North=0, East=90)
 - `Mag_without_K-correction`: extinction-correlated absolute magnitude calculated from the galaxy distance and the GLAFIT Integrated magnitude before K-correction

`Distance`: value from input or from NED redshift (using WMAP9 model from astropy)

`Redshift`: redshift used for K correction. From the command line or NED.


Example with using coordinate:
```
photomass_ls.py ngc474 --coordinates "20.0233 deg, 3.4143 deg" 106.3  --dist 30.88 --galpath `pwd`
https://www.legacysurvey.org/viewer/fits-cutout?ra=20.0233&dec=3.4143&height=1024&width=1024&layer=ls-dr10&pixscale=0.524&bands=gr
Downloading https://www.legacysurvey.org/viewer/fits-cutout?ra=20.0233&dec=3.4143&height=1024&width=1024&layer=ls-dr10&pixscale=0.524&bands=gr
|===============================================================================================================================================================================================================================================================| 8.3M/8.3M (100.00%)         1s
Galaxy: ngc474
RA: 20.0233 Dec: 3.414
log10(M*[Sun]): E: 10.72 M: 10.77 Q: 10.81
Ext[mag]: g : 0.112 r : 0.075
Mag[mag]: g : -20.32 r : -21.12
Sersic index: g : 1.273 r : 1.929
R_e[px]: g : 68.91 r : 66.22
R_e[arcsec]: g : 36.11 r : 34.7
Axis ratio: g : 0.8257 r : 0.8142
Position angle[deg]: g : 14.79 r : 16.3
Distance[Mpc]: 30.88 (from input)
Plate scale[arcsec / px]: 0.524 0.524
Zero point[mag]: 20.99
Redshift: N/A
```
Here the `OBJECT` is not converted to uppercase, redshift is not received from the database and so it is not available. 
> Warning: the script checks if the file (filename) with data from legacy survey is present - there is no idication of coordinates in the filename, so you have to use `--refetch`, if you want to use same object name with different coordinates.

Example with K correction:
```
$photomass/photomass_ls.py ngc474 106.3  --dist 30.88 --galpath /path/to/galfit/ --K --redshift 0.00794
Using existing FITS from the local directory
Galaxy: NGC474
Mag_without_K-correction[Mag]: g: -20.5 r: -21.29
Applying K correction: g: 0.0085 r: 0.0078
RA: 20.0279 Dec: 3.416
log10(M*[Sun]): E: 10.79 M: 10.84 Q: 10.89
Ext[mag]: g : 0.112 r : 0.075
Mag[mag]: g : -20.51 r : -21.3
Sersic index: g : 3.049 r : 3.872
R_e[px]: g : 44.73 r : 48.04
R_e[arcsec]: g : 23.44 r : 25.17
Axis ratio: g : 0.8185 r : 0.8188
Position angle[deg]: g : 21.18 r : 19.1
Distance[Mpc]: 30.88 (from input)
Plate scale[arcsec / px]: 0.524 0.524
Zero point[mag]: 20.99
Redshift: 0.007722
```
this example is with data from LS already present in the diractory.

The script also saves these files:
 - downloaded FITS - with file name : `<object_name>_<downsize>_<size>px_ls-dr10_<filters>.fits`
 - jpeg with same filters downloaded from LS - with file name : `<object_name>_<downsize>_<size>px_ls-dr10_<filters>.jpg`
 - single-filter data: `<object_name>_<filter>.fit`
 - image mask for the respective filter: `<object_name>_<filter>_mask.fit`
 - input file for GALFIT: `<object_name>_<filter>_gal.inp`
 - GALFIT FITS image output: `<object_name>_<filter>_out.fits`
 - GALFIT ASCII text file output: `<object_name>_<filter>.galfit`

## Acknowledgements

Funded by the European Union under the Marie Skłodowska-Curie grant agreement No. 101067618 (GalaxyMergers). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or European Research Executive Agency (REA). Neither the European Union nor REA can be held responsible for them.

This research has made use of the *K-corrections calculator* service
http://kcor.sai.msu.ru and the NASA/IPAC Extragalactic Database, which is funded by the National Aeronautics and Space Administration and operated by the California Institute of Technology.
