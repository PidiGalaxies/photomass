#import requests
#from numpy import *
import numpy as np
import astropy.io.fits as fits
from astropy import units as u
from astropy.coordinates import get_icrs_coordinates
from astropy.coordinates import SkyCoord
import sep
from scipy.ndimage import gaussian_filter
import os

def get_image_fits(ra, dec, wpx, downsize, LSversion='ls-dr9',bands="grz" ):
    '''
    Gets fits from legacy survay - details at https://www.legacysurvey.org/dr9/description/
    Args:
        ra, dec (float): center coordinates
        wpx (int) : size of fits in pixels
        downsize (int) : dowscale with respect to default plate
        LSversion (str) : Version of Legacy Survay, default='ls-dr9'
        band (str) : string with bands to be downloaded - available are g, r, z, i; default="grz"
    '''
    # https://www.legacysurvey.org/dr9/description/
    ps = 0.262 * downsize # Pixscale=0.262 will return (approximately) the native pixels

    url = "https://www.legacysurvey.org/viewer/fits-cutout?"
    url += "ra="        + str(ra)
    url += "&dec="      + str(dec)
    url += "&height="   + str(wpx)
    url += "&width="    + str(wpx)
    url += "&layer="    + LSversion
    url += "&pixscale=" + str(ps)
    url += "&bands="    + bands

    print(url)

    image_fits = fits.open(url, cache=False)
    return image_fits

def harvest_ned(_id):
    '''
    Get extinction for r, g, z filters and redshift.
    Args:
        _id (string): Id of object
    Returns:
        {'r': (float), 'z': (float), 'i': (float), 'redshift': (float)}
    '''
    import requests
    oid = _id.encode('ascii')   
    x = requests.get('https://ned.ipac.caltech.edu/cgi-bin/gmd?uplist=' + oid.decode('ascii') + '&delimiter=bar&nondb=user_objname&position=z&SchEXT=DES+r&SchEXT=DES+g&SchEXT=DES+z')
    data = [z for z in x.content.split(b'\n') if oid in z ]
    if len(data) == 0:
        raise(Exception("Failed to get data for "+ _id + " from NED"))
    data = data[0].split(b'|') # the line with values starts with oid and is deliminated by '|'
    return {'g':float(data[1]), 'r':float(data[2]),'z':float(data[3]),'redshift':float(data[4])} 

import sys
import os
from argparse import ArgumentParser

parser = ArgumentParser(description="Calculates things from other things!")
parser.add_argument("OBJECT", help="Name of the object for which to calculate things.")
parser.add_argument("radius", help="Semi-majir axis[arcsec] at level {mu}(3.6um)=25.5mag(AB)/arcsec^2^[ucd=phys.angSize.smajAxis].", type=float)
parser.add_argument("--overwrite", help="If set the fits will be redownloaded and analysis redone - otherwise nothing will be done.", action="store_true")
parser.add_argument("--local", help="use local file if it exists", action="store_true")
parser.add_argument("--galpath", help="Path to galfit64 binary. It should be in $PATH if not set.", action="store", default="")
parser.add_argument("--dist", help="Distance in Mcp", action="store", default=None)
parser.add_argument("--additional-filters", help="Other filters without delimiter", action="store", default='')
parser.add_argument("--plot", help="Plot png image of source data and masked data", action="store_true")

args = parser.parse_args()
filters = ['g', 'r']  + list(args.additional_filters)
obj_name =  args.OBJECT #'ngc4993'

ned = harvest_ned(obj_name)
if args.dist is None:
    from astropy.cosmology import WMAP9
    dist = WMAP9.comoving_distance(ned['redshift']).value #distance in Mpc
else:
    dist = float(args.dist)
obj_r =  args.radius # amaj(arcsec) Semi-major axis at level {mu}(3.6um)=25.5mag(AB)/arcsec^2^[ucd=phys.angSize.smajAxis]
r_fac = 2.5
ls_scale = 0.262 #LS default plate scale 0.262 arcsec/pixel
wpx = 1024 # The maximum size for cutouts (in number of pixels) is currently 512 # currently more!!    

downsize = int(2 * r_fac * obj_r / ls_scale / wpx + 0.5) # with respect to the default plate scale of 0.262 arcsec/pixel
downsize = downsize if downsize > 0 else 1

ZP = 22.5 - 2.5 * np.log10(downsize ** 2) # zero point [22.5 for the LS fiducal pixel scale of 0.262 arcsec/px]

center = get_icrs_coordinates(obj_name)

file_name = obj_name +'_' + str(downsize) + '_' + str(wpx) + 'px_' + 'ls-dr10_' + "".join(filters) +  '.fits'

if os.path.isfile(file_name) and args.local:
    pass
elif os.path.isfile(file_name) and not args.overwrite:
    print("File", file_name, "already exists!")
    print("If you want do download it and do the analysis again, please add option --overwrite or use already downloaded file with --local")
    sys.exit(0)
else:
    im_fits = get_image_fits(center.ra.value, center.dec.value, wpx, downsize=downsize, LSversion = 'ls-dr10', bands = "".join(filters))
    im_fits.writeto(file_name, overwrite=True)  # save / overwrite galaxy image
    im_fits.close()

magnitudes = {}
Ser_index  = {}
axis_ratio = {}
R_e        = {}
R_e_arcsec = {}

for i, band in zip(range(len(filters)), filters):
    im_fits = fits.open(file_name)
    
    out_image_name = obj_name + '_' + band + '.fit'
    out_mask_name  = obj_name + '_' + band + '_mask.fit'

    # create image with only one band

    im_fits[0].data = im_fits[0].data[i]
    im_fits[0].header['BANDS'] = band
    im_fits[0].header['BAND0'] = band
    for i in range(1, len(filters)):
            im_fits[0].header.remove('BAND' + str(i)) 
    im_fits.writeto(out_image_name, overwrite=True)  # save galaxy image
    
    
    im = im_fits[0].data
    if len(im) == 0 or np.all(im==0):
        print('do data for band:', band)
        filters.remove(band)
        continue
    
    
    ######
    ## SOURCE DETECTON, MASK CREATION
    
    # sep manual : https://sep.readthedocs.io/en/v1.1.x/reference.html
    
    sep.set_extract_pixstack(17e6)
    sep.set_sub_object_limit(1e4)
    
    # detect extended objects
    im = im.byteswap().newbyteorder()
    bkg = sep.Background(im, bw=wpx // 4, bh=wpx // 4)
    
    data_sub = im - bkg
    objects, segmap = sep.extract(data_sub, deblend_cont=1.0, thresh=1, err=bkg.globalrms, segmentation_map=True) #0.03
    
    mainobj_no = segmap[wpx // 2,wpx // 2] ## object number of the main galaxy (is in the image center)
    mask = (segmap > 0) * (segmap != mainobj_no) # objects to mask
    
    #detect small objects (stars)
    
    bkg_fine = sep.Background(im, bw = 2, bh = 2)
    
    data_sub_fine = im - bkg_fine
    objects_fine, segmap_fine = sep.extract(data_sub_fine, deblend_cont=1, thresh=2, err=bkg_fine.globalrms, segmentation_map=True) #0.03
    
    mainobj_no_fine = segmap_fine[wpx // 2,wpx // 2] ## object number of the main galaxy (is in the image center)
    
    maxratio = 2
    #elongated objects to exclude from mask
    l = [i + 1 for i in range(len(objects_fine)) if objects_fine['a'][i] / objects_fine['b'][i] > maxratio]
    l.append(mainobj_no_fine)
    segmap_fine[np.isin(segmap_fine, l)] = 0
    
    mask_fine = (segmap_fine > 0) # small objects to mask
    
    
    mask_tot = (mask_fine + mask) > 0 # masked by at least one mask
    mask_tot_smooth = gaussian_filter(mask_tot * 1., sigma=5) > 0.1 #expand mask
    mask_tot = (mask_tot + mask_tot_smooth) > 0 # masked pixels are True, unmasked False
    
    hdu = fits.PrimaryHDU((mask_tot).astype('int'))  # save mask
    hdu.writeto(out_mask_name, overwrite=True)
    
    
    
    
    #####
    # set initial estimates for Galfit
    # xmin, xmax, ymin, ymax, plate_scale, ZP
    # xc, xy, mag, Re, n, q, PA
    
    
    xc = wpx / 2
    yc = wpx / 2
    
    xmin = 0
    xmax = wpx - 1
    ymin = 0
    ymax = wpx - 1

    plate_scale = 0.262 * downsize # arcsec/pixel
     
    mag = -2.5 * np.log10(objects['flux'][mainobj_no - 1]) + ZP
    a     = objects['a']    [mainobj_no - 1]
    b     = objects['b']    [mainobj_no - 1]
    theta = objects['theta'][mainobj_no - 1]
    Re = np.sqrt(a * b)
    n = 4
    q = b / a
    PA = -theta / np.pi * 180 + 90
    sky_bkg = np.nanmedian(bkg)
    
    #print('xmin = ', xmin, '\nxmax = ', xmax, '\nymin = ', ymin, '\nymax = ', ymax, '\nplate_scale = ', plate_scale, '\nZeropoint = ', ZP,'\n\n') 
    #print('xc = ', xc, '\nyc = ', yc,'\nmag = ',mag,'\nRe = ', Re, '\nn = ', n, '\nq = ',q, '\nPA = ', PA) 
    
    
    
    galimp = f'''===============================================================================
    # IMAGE and GALFIT CONTROL PARAMETERS
    A) {out_image_name}            # Input data image (FITS file)
    B) {obj_name}_{band}_out.fits       # Output data image block
    C) none                # Sigma image name (made from data if blank or "none") 
    D) none   #        # Input PSF image and (optional) diffusion kernel
    E) 1                   # PSF fine sampling factor relative to data 
    F) {out_mask_name}                # Bad pixel mask (FITS image or ASCII coord list)
    G) none                # File with parameter constraints (ASCII file) 
    H) {xmin} {xmax} {ymin} {ymax}   # Image region to fit (xmin xmax ymin ymax)
    I) 100    100          # Size of the convolution box (x y)
    J) {ZP}              # Magnitude photometric zeropoint 
    K) {plate_scale} {plate_scale}        # Plate scale (dx dy)    [arcsec per pixel]
    O) regular             # Display type (regular, curses, both)
    P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps
    
    # INITIAL FITTING PARAMETERS
    #
    #   For object type, the allowed functions are: 
    #       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat, 
    #       ferrer, powsersic, sky, and isophote. 
    #  
    #   Hidden parameters will only appear when they're specified:
    #       C0 (diskyness/boxyness), 
    #       Fn (n=integer, Azimuthal Fourier Modes),
    #       R0-R10 (PA rotation, for creating spiral structures).
    # 
    # -----------------------------------------------------------------------------
    #   par)    par value(s)    fit toggle(s)    # parameter description 
    # -----------------------------------------------------------------------------
    
    # Object number: 1
     0) sersic                 #  object type
     1) {xc} {yc}  1 1  #  position x, y
     3) {mag}     1          #  Integrated magnitude
     4) {Re}        1          #  R_e (half-light radius)   [pix]
     5) {n}         1          #  Sersic index n (de Vaucouleurs n=4) 
     6) 0.0000      0          #     ----- 
     7) 0.0000      0          #     ----- 
     8) 0.0000      0          #     ----- 
     9) {q}      1          #  axis ratio (b/a)  
    10) {PA}    1          #  position angle (PA) [deg: Up=0, Left=90]
     Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 
    
    # Object number: 2
     0) sky                    #  object type
     1) {sky_bkg}      1          #  sky background at center of fitting region [ADUs]
     2) 0.0000      1          #  dsky/dx (sky gradient in x)
     3) 0.0000      1          #  dsky/dy (sky gradient in y)
     Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 
    
    # downsize = {downsize}
    
    ================================================================================ '''
    
    #Write galfit input
    gal_im = obj_name + '_' + band + '_gal.imp'
    with open(gal_im, 'w') as f:
        f.write(galimp)
    
    #Plot png image with masked data vs full data
    if args.plot:
        import matplotlib.pyplot as plt

        plt.subplot(121)
        plt.imshow(np.log10(im) * (~mask_tot))
        plt.subplot(122)
        plt.imshow(np.log10(im))
        plt.savefig(obj_name + '_' + band + '_image.png')
    
    
    #Tests if galfit ouputs exists - there should not be any galfit.* files!
    if os.path.isfile('galfit.*'):
        print('''File galfit.* exists!
Please delete all galfit ouputs! - (rm galfit.*)''')
        raise
    if os.system(args.galpath + 'galfit64 ' + gal_im + ' >/dev/null'):
        raise Exception('Problem executing galfit64!')
    
    #Parse galfit output
    with open('galfit.01') as f:
        for line in f:
            if 'Integrated magnitude' in line:
                magnitudes[band] = float(line.strip().split(')')[1].strip().split(' ')[0])
            if 'Sersic index' in line:
                Ser_index[band] = float(line.strip().split(')')[1].strip().split(' ')[0])
            if 'R_e' in line:
                R_e[band] = float(line.strip().split(')')[1].strip().split(' ')[0])
            if 'Axis ratio' in line:
                axis_ratio[band] = float(line.strip().split(')')[1].strip().split(' ')[0])
            if 'Plate scale' in line:
                _ = line.strip().split(')')[1].strip().split()
                px_to_ang_dx = float(_[0])
                px_to_ang_dy = float(_[1])
    if px_to_ang_dx == px_to_ang_dx and px_to_ang_dx == plate_scale:
        R_e_arcsec[band] = R_e[band] * px_to_ang_dx
    
    #Fix magnitudes - only for the filters with extinction from NED!!!
    if band in ned.keys():
        magnitudes[band] += - 25 - ned[band] - 5 * np.log10(dist)
    
    #Rename galfit output to unique name!
    os.rename('galfit.01', obj_name + '_' + band + '.galfit')


print("M[Sun] =", 0.673 * magnitudes['g'] - 1.108 * magnitudes['r'] + 0.996 )
print('Ext:',    " ".join([band + ' : '      + str(ned[band])        for band in filters]))
print('Mag:',    " ".join([band + ' : '      + str(magnitudes[band]) for band in filters]))
print('Sersic:', " ".join([band + ' : '      + str(Ser_index[band])  for band in filters]))
print('R_e[px]:', " ".join([band + ' : '     + str(R_e[band])        for band in filters]))
print('R_e[arcsec]:', " ".join([band + ' : ' + str(R_e_arcsec[band])        for band in filters]))
print('Axis ratio:', " ".join([band + ' : '  + str(axis_ratio[band]) for band in filters]))
print('Distance:', dist, 'Mpc')
print('redshift:', ned['redshift'])
