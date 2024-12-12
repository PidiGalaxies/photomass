*** photmass_ls.py ***

** Requirements **
------------------
galfit : https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html

python libraries:
requests, numpy, scipy, SEP, astropy, matplotlib (with option --plot)


You can install required packaged into virtual environment:
```
python3 -m venv venv
. venv/bin/activate
pip3 install "numpy<2.0" scipy SEP astropy requests
```

> SEP 1.2.1 does not work with numpy 2.0 or newer

** Limitations **
The script uses [NED](ned.ipac.caltech.edu) and [`get_icrs_coordinates` from `astropy`](cds.unistra.fr) to get galactic extinctions and coordinates.
If the object is not present in these databases corresponding exception will be raised.


** Running script ** 
You can use `python3 photomass_ls.py` to get usage:
```
usage: photomass_ls.py [-h] [--overwrite] [--local] [--galpath GALPATH] [--dist DIST] [--additional-filters ADDITIONAL_FILTERS] [--plot] OBJECT radius

Calculates things from other things!

positional arguments:
  OBJECT                Name of the object for which to calculate things.
  radius                Semi-majir axis[arcsec] at level {mu}(3.6um)=25.5mag(AB)/arcsec^2^[ucd=phys.angSize.smajAxis].

options:
  -h, --help            show this help message and exit
  --overwrite           If set the fits will be redownloaded and analysis redone - otherwise nothing will be done.
  --local               use local file if it exists
  --galpath GALPATH     Path to galfit64 binary. It should be in $PATH if not set.
  --dist DIST           Distance in Mcp
  --additional-filters ADDITIONAL_FILTERS
                        Other filters without delimiter
  --plot                Plot png image of source data and masked data
```

The outpus are to stdout. For example:
```
$ python3 ../photomass/photomass_ls.py  NGC4656 4 --galpath ../ --local
M[Sun] = 9.247399396604427
Ext: g : 0.043 r : 0.029
Mag: g : -18.493440451964197 r : -18.680040451964196
Sersic: g : 1.1615 r : 1.1101
R_e[px]: g : 439.2828 r : 420.0659
R_e[arcsec]: g : 115.09209360000001 r : 110.05726580000001
Axis ratio: g : 0.1761 r : 0.1885
Distance: 9.315541188290236 Mpc
redshift: 0.002155
```
Where M[Sun] is mass estimate in Sun masses, 
These variables are calculated for different filters:
 - `Ext`: galactic extinctions from NED 
 - `Mag`: absolute magnitude (corrected using galactic extinction)
 - `Sersic`: Sersic radius
 - `R_e`: Effective radius - in pixels and arcsec
 - `Axis ratio`
`Distance`: value from input or from redshift (using WMAP9 model from astropy)
`redshift`: redshift from NED 


The script also saves these files:
 - downloaded FISTS - with file name : `<object_name>_<downsize>_<size>px_ls-dr10_<filters>.fits`
 - fits with only one filter data: `<object_name>_<filter>.fit`
 - fits with mask for one filter data: `<object_name>_<filter>_mask.fit`
 - fits outputted from galfit: `<object_name>_<filter>_out.fits`
 - input for galfit: `<object_name>_<filter>_gal.inp`
 - output from galfit: `<object_name>_<filter>.galfit`
