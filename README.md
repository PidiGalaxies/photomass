# photmass_ls.py


## About the script

photomass_ls computes an estimate of the stellar mass for a given galaxy based on a photometry performed by GALFIT on g & r band data from DESI Legacy Imaging Surveys. The script automatically downloads images of the galaxy from the Legacy Surveys (DR10) database, creates image masks, generates GALFIT
input files with well-assessed initial values, performs the GALFIT photometry, and calculates the stellar mass estimate. More details are given in Ebrová, Bílek, & Eliášek (2025).


## Citation

If you use this code in your research or projects, please cite the following accompanied paper:

Ebrová, Bílek, & Eliášek: "Photometric stellar masses for galaxies in DESI Legacy Imaging
Surveys" (2025; submitted)


## Requirements

GALFIT : https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html

python libraries:
requests, numpy, scipy, SEP, astropy, [matplotlib (required for use of option `--plot`)]


You can install required packages into virtual environment:
```
python3 -m venv venv
. venv/bin/activate
pip3 install "numpy<2.0" scipy SEP astropy requests
```

> SEP 1.2.1 does not work with numpy 2.0 or newer


## Limitations
The script uses [NED](ned.ipac.caltech.edu) and [`get_icrs_coordinates` from `astropy`](cds.unistra.fr) to get galactic extinctions and coordinates for the given galaxy.
If the object is not present in these databases corresponding exception will be raised.


## Running the script
You can run the script to get list of options with descriptions like this:
```
$ ./photomass_ls.py -h
usage: photomass_ls.py [-h] [--refetch] [--galpath GALPATH] [--dist DIST] [--additional-filters ADDITIONAL_FILTERS] [--plot] OBJECT radius

computes a stellar mass estimate for a given galaxy based on GALFIT photometry with Legacy Surveys data

positional arguments:
  OBJECT                name of the galaxy for which the stellar mass estimate is done
  radius                approximate galaxy angular size in arcseconds (e.g. semi-major axis at SB=25.5mag/arcsec^2 in the Spitzer 3.6 micron filter)

options:
  -h, --help            show this help message and exit
  --refetch             re-download and re-analyse image even if a local copy exists
  --galpath GALPATH     path to galfit64 binary; if not set, $PATH is searched
  --dist DIST           galaxy distance in Mcp
  --additional-filters ADDITIONAL_FILTERS
                        other filters without delimiter (e.g. 'xi')
  --plot                plot png image of source data and masked data
```

Thea numerical outpus are printed to stdout. For example:
```
$./photomass_ls.py NGC474 106.3 --dist 30.88
Downloading https://www.legacysurvey.org/viewer/fits-cutout?ra=20.02786271688&dec=3.41551721475&height=1024&width=1024&layer=ls-dr10&pixscale=0.524&bands=gr
|===========================================================================================================================| 8.3M/8.3M (100.00%)         2s
Galaxy: NGC474
RA: 20.0233 deg Dec: 3.4143 deg
log10(M*[Sun]): 10.7895
Ext[mag]: g : 0.112 r : 0.075
Mag[mag]: g : -20.5046 r : -21.2934
Sersic index: g : 3.0489 r : 3.872
R_e[px]: g : 44.7267 r : 48.0426
R_e[arcsec]: g : 23.4368 r : 25.1743
Axis ratio: g : 0.8185 r : 0.8188
Position angle[deg]: g : 21.1833 r : 19.0952
Distance[Mpc]: 30.88 (from input)
Redshift: 0.007722
```
where:

`Galaxy`: `OBJECT`

`RA`, `Dec`: coordinates of the `OBJECT`

`log(M*[Msun]`: the estimate of the logarithm of the galaxy stellar mass computed by Eq(1) of Ebrova et al. (2025), where Msun is the sollar masses.

The next set of variables is calculated for different filters:
 - `Ext`: Galactic extinctions in DES filters from Schlafly & Finkbeiner (2011) using NASA/IPAC Extragalactic Database (NED)
 - `Mag`: extinction-correlated absolute magnitude calculated from the galaxy distance and the GLAFIT Integrated magnitude
 - `Sersic index`: Sersic index of the GALFIT model
 - `R_e`: effective radius of the GALFIT model - in pixels and arcsec
 - `Axis ratio`: axis ratio of the GALFIT model
 - `Position angle`: positioning angle of the GALFIT model (North=0, East=90)

`Distance`: value from input or from redshift (using WMAP9 model from astropy)

`Redshift`: redshift from NED


The script also saves these files:
 - downloaded FITS - with file name : `<object_name>_<downsize>_<size>px_ls-dr10_<filters>.fits`
 - single-filter data: `<object_name>_<filter>.fit`
 - image mask for the respective filter: `<object_name>_<filter>_mask.fit`
 - input file for GALFIT: `<object_name>_<filter>_gal.inp`
 - GALFIT FITS image output: `<object_name>_<filter>_out.fits`
 - GALFIT ASCII text file output: `<object_name>_<filter>.galfit`

## Acknowledgements

This project has received funding from the European Union's Horizon Europe Research and Innovation programme under the Marie Skłodowska-Curie grant agreement No. 101067618, GalaxyMergers.
