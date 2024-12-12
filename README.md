# photmass_ls.py


## About the script

photomass computes an estimate of the stellar mass for a given galaxy based on a photometry performed by GALFIT on g & r band data from DESI Legacy Imaging Surveys. The script automatically downloads images of the galaxy from the Legacy Surveys (DR10) database, create image masks, generate GALFIT
input files with well-assessed initial values, perform the GALFIT photometry, and calculate the stellar mass estimate. More details are given in Ebrová, Bílek, & Eliášek (2025).


## Citation 

If you use this code in your research or projects, please cite the following accompanied paper:

Ebrová, Bílek, & Eliášek: "Photometric stellar masses for galaxies in DESI Legacy Imaging
Surveys" (2025; submitted)


## Requirements 

GALFIT : https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html

python libraries:
requests, numpy, scipy, SEP, astropy, matplotlib (with option --plot)


You can install required packaged into virtual environment:
```
python3 -m venv venv
. venv/bin/activate
pip3 install "numpy<2.0" scipy SEP astropy requests
```

> SEP 1.2.1 does not work with numpy 2.0 or newer


## Limitations
The script uses [NED](ned.ipac.caltech.edu) and [`get_icrs_coordinates` from `astropy`](cds.unistra.fr) to get galactic extinctions and coordinates.
If the object is not present in these databases corresponding exception will be raised.


## Running the script
You can use `python3 photomass_ls.py` to get usage:
```
usage: photomass_ls.py [-h] [--overwrite] [--local] [--galpath GALPATH] [--dist DIST] [--additional-filters ADDITIONAL_FILTERS] [--plot] OBJECT radius

Computes a stellar mass estimate for a given galaxy based on GALFIT photometry with Legacy Surveys data.

positional arguments:
  OBJECT                Name of the galaxy for which the stellar mass estimate is done.
  radius                Approximate galaxy angular size in arcseconds (e.g. semi-major axis at SB=25.5mag/arcsec^2 in the Spitzer 3.6 micron filter).

options:
  -h, --help            Show this help message and exit.
  --overwrite           If set, the fits will be redownloaded and analysis redone - otherwise nothing will be done.
  --local               Use the local file if it exists.
  --galpath GALPATH     Path to galfit64 binary. If not set, $PATH is searched.
  --dist DIST           Galaxy distance in Mcp.
  --additional-filters ADDITIONAL_FILTERS
                        Other filters without delimiter.
  --plot                Plot png image of source data and masked data.
```

The outpus are to stdout. For example:
```
$ python3 ../photomass/photomass_ls.py  NGC4656 4 --galpath ../ --local
galaxy: NGC4656
RA,Dec =  
log(M*[Msun]) = 9.247399396604427
Ext[mag]: g : 0.043 r : 0.029
Mag[mag]: g : -18.493440451964197 r : -18.680040451964196
Sersic index: g : 1.1615 r : 1.1101
R_e[px]: g : 439.2828 r : 420.0659
R_e[arcsec]: g : 115.09209360000001 r : 110.05726580000001
Axis ratio: g : 0.1761 r : 0.1885
Distance[Mpc]: 9.315541188290236
redshift: 0.002155
```
where and log(M*[Msun]) is the estimate of the logarithm of the galaxy stellar mass computed by Eq(1) of Ebrova et al. (2025), and Msun is the sollar masses.
 
These variables are calculated for different filters:
 - `Ext`: Galactic extinctions in DES filters from Schlafly & Finkbeiner (2011) using NASA/IPAC Extragalactic Database (NED)
 - `Mag`: extinction-correlated absolute magnitude calculated from the galaxy distance and the GLAFIT Integrated magnitude
 - `Sersic index`: Sersic index of the GALFIT model
 - `R_e`: effective radius of the GALFIT model - in pixels and arcsec
 - `Axis ratio`: axis ratio of the GALFIT model

`Distance`: value from input or from redshift (using WMAP9 model from astropy)

`redshift`: redshift from NED 


The script also saves these files:
 - downloaded FISTS - with file name : `<object_name>_<downsize>_<size>px_ls-dr10_<filters>.fits`
 - fits with only one filter data: `<object_name>_<filter>.fit`
 - fits with mask for one filter data: `<object_name>_<filter>_mask.fit`
 - fits outputted from GALFIT: `<object_name>_<filter>_out.fits`
 - input for GALFIT: `<object_name>_<filter>_gal.inp`
 - output from GALFIT: `<object_name>_<filter>.galfit`

## Acknowledgements

This project has received funding from the European Union's Horizon Europe Research and Innovation programme under the Marie Skłodowska-Curie grant agreement No. 101067618, GalaxyMergers.
