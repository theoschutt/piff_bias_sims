# get catalog of exposures #TODO: decide on exposures (maybe (part of) latest 1000 tile run)
# columns needed in the exp catalog
#    - expnum
#    - band
#    - boresight: TELRA, TELDEC
#    - 
# galaxy catalog: from Y6 Gold #TODO: make catalog
#    - RA/DEC
#    - X,Y
#    - size, e1, e2
#    - griz mags
#    - redshift
#    - match to template galaxy SED by color #TODO: get galaxy SEDs
#    - 
# star catalog: from Y6 Gold (or could use Piff output catalog?) #TODO: make catalog
#    - RA/DEC
#    - X,Y
#    - size, e1, e2
#    - griz mags
#    - match to template stellar SED by color
#    - 

import os
import numpy as np
import fitsio
import pixmappy
import galsim_extra
from toFocal import toFocalDegree

def getCCDCorners(ccdnum, telra, teldec):
    """
    Get lower left and upper right corners of sensor in RA/DEC degrees
    """

    # The pixel centers go from 1-2048, 1-4096,
    # so the edges are 0.5-2048.5, 0.5-4096.5
    u_ll, v_ll = toFocalDegree(ccdnum, 0.5, 0.5)
    u_ur, v_ur = toFocalDegree(ccdnum, 2048.5, 4096.5)

    min_ra = telra + u_ll
    min_dec = teldec + v_ll
    max_ra = telra + u_ur
    max_dec = teldec + v_ur

    return min_ra, min_dec, max_ra, max_dec

def makeCCDCatalog(catalog, cattype, ccdnum, expnum=None,
    min_ra=None, min_dec=None, max_ra=None, max_dec=None):
    """
    Get catalog of objects within the edges of a sensor.
    Can use RA/DEC coordinates to filter object (cattype=='radec'). 
    Can also use the Piff summary star catalog and just filter by expnum
    and ccdnum (cattype == 'piff').
    
    @return filtered catalog
    """

    if cattype == 'radec':
        ccd_obj_mask = ((catalog['RA'] > min_ra)
                        & (catalog['RA'] < max_ra)
                        & (catalog['DEC'] > min_dec)
                        & (catalog['DEC'] < max_dec))

    # if using piff summary star catalog, can just filter by exp&ccd
    elif cattype == 'piff':
        ccd_obj_mask = ((catalog['EXPNUM'] == expnum)
                        & (catalog['CCDNUM'] == ccdnum))
    else:
        print('Invalid cattype: ',cattype)
        exit

    return catalog[ccd_obj_mask]

def getCCDCatalogs(catalogs, cattypes, ccdnum, expnum,
    min_ra=None, min_dec=None, max_ra=None, max_dec=None):
    """
    Cuts arbitrary number of catalogs to coordinates of a sensor.
    @return: List of filtered catalogs 
    """

    filtered_cat_list = []
    for ccat, cattype in zip(catalogs, cattypes):

        ccd_cat = makeCCDCatalog(ccat, cattype, ccdnum, expnum,
                      **(min_ra, min_dec, max_ra, max_dec))
        filtered_cat_list.append(ccd_cat)

    return filtered_cat_list

def main(argv):

    #set version number
    ver = 1

    # setup directories
    if not os.path.isdir('output'):
        os.mkdir('output')

    # load catalogs that will be used repeatedly
    catpath = 'catalogs'
    exppath = os.path.join(catpath, 'expcat_v%i.fits'%ver)
    galpath = os.path.join(catpath, 'galcat_v%i.fits'%ver)
    starpath = os.path.join(catpath, 'starcat_v%i.fits'%ver)
    expcat = fitsio.read(exppath)
    galcat = fitsio.FITS(galpath)
    starcat = fitsio.FITS(starpath)
    bad_ccds = [2, 61] # We don't use CCDs 2 and 61

    # set config values that are true for all exposures
    config = {}

    config['image'] = {
        'pixel_scale' : 0.263, # DECam arcsec/pixel
        'x_size' : 2048,
        'y_size' : 4096,
        'wcs' : {
            'type' : 'Pixmappy',
            'dir' : 'pixmappy/data/',
            'file_name' : 'y6a1.guts.astro', #TODO: check
        },
        'noise' : {
            'type' : 'CCD',
            'sky_level_pixel' : 4000, #ADU/px
            'gain' : 1.0, #e-/ADU
            'read_noise' : 7 #ADU/px
        }
    }

    config['stamp'] = {
        'type' : 'MixedScene',
        'draw_method' : 'phot',
        'objects' : {'gal' : 0.8, 'star' : 0.2}, #TODO: necessary with gal/star catalogs?
    }    

    config['output']['dir'] = 'output'

    # loop over exposure catalog
    for exposure in expcat:
        # exposure-specific variables
        expnum = exposure['EXPNUM']
        band = exposure['BAND']
        telra = exposure['TELRA']
        teldec = exposure['TELDEC']

        # set exposure-specfic config values
        config['image'] = {
            'wcs': {'exp' : expnum},
            'bandpass' : {
                'file_name' : 'decam_bandpass_%s.dat'%band,
                'wave_type': 'angstrom'
            }
        }

        # loop over CCDs in exposure
        for ccdnum in range(1,63):
            # skip bad ccds
            if ccdnum not in bad_ccds:
                # get CCD RA/DEC limits
                min_ra, min_dec, max_ra, max_dec = getCCDCorners(ccdnum, telra, teldec)
                # get input catalogs for the CCD
                [ccdgalcat, ccdstarcat] = getCCDCatalogs([galcat, starcat],
                                              ['radec', 'radec'], ccdnum, expnum,
                                              **(min_ra, min_dec, max_ra, max_dec))
                # set filenames
                fn_galcat = os.path.join(catpath, 'sim_v%i_%i_c%i_%s_gal_catalog.fits'%(ver,
                    expnum,ccdnum,band))
                fn_starcat = os.path.join(catpath, 'sim_v%i_%i_c%i_%s_star_catalog.fits'%(ver,
                    expnum,ccdnum,band))

                # write CCD object catalogs
                fitsio.write(fn_galcat, ccdgalcat)
                fitsio.write(fn_starcat, ccdstarcat)

                # set CCD-specific config values
                config['input'] = {'catalog' : [{'file_name' : fn_galcat},
                                                {'file_name' : fn_starcat}]
                }
                config['image'] = {'wcs': {'ccdnum': ccdnum}},
                config['stamp'] = {
                    'image_pos' : {
                        'type' : 'XY',
                        'x' : {
                            'type' : 'Catalog', 
                            'col': 'X',
                            'num': '$current_obj_type_index'
                        },
                        'y' : {
                            'type' : 'Catalog', 
                            'col': 'Y',
                            'num': '$current_obj_type_index'
                        } #TODO: ask re: image_pos, world_pos, sky_pos and XY vs
                          # RADec type
                    }
                }

                config['psf'] = {'type' : 'InterpolatedImage'} #TODO: or custom GSObject type?
                config['gal'] = {
                    'type' : 'DeVaucouleurs', #TODO: ask about profile type
                    'half_light_radius' : {
                        'type' : 'Catalog',
                        'col' : 'FLUX_RADIUS_%s'%upper(band), #TODO: check if correct column to use
                        'num' : '$current_obj_type_index'
                    },
                    'sed' : {
                        'file_name' : {
                            'type' : 'Catalog',
                            'col' : 'SED_TEMPLATE_ID', #TODO: will need to match to gal SEDs
                            'num' : '$current_obj_type_index'
                        },
                        'wave_type' : 'Ang',
                        'flux_type' : 'flambda', #TODO: check
                        'norm_flux' : 1., #TODO: is this necessary?
                        'norm_bandpass' : '@image.bandpass'
                    }
                    'redshift' : 0.,#TODO: is this already baked into the SED matching?
                    'ellip' : {
                        'type' : 'G1G2', #Y6 Gold uses reduced shear
                        'g1' : {
                            'type' : 'Catalog', 
                            'col' : 'SOF_BDF_G_1',
                            'num' : '$current_obj_type_index'
                        },
                        'g2' : {
                            'type' : 'Catalog', 
                            'col' : 'SOF_BDF_G_2',
                            'num' : '$current_obj_type_index'
                        }
                    }
                    'rotate' : {'type' : 'Random'},
                    'flux' : #TODO: Do flux and norm_flux have to be specified?           
                }

                config['star'] = {
                    'type': 'DeltaFunction', #TODO: or also InterpolatedImage
                    'sed': {},#TODO: will need to match catalog star with SED template,
                    'flux': 1. #TODO: prob from catalog, or do norm_flux in SED?
                }

                config['output'] = {
                    'file_name' : 'sim_v%i_%i_c%i_%s_image.fits'%(ver,expnum,ccdnum,band),
                    'psf': {'file_name': 'sim_v%i_%i_c%i_%s_psf.fits'%(ver,expnum,ccdnum,band)}
                }

                # run galsim on this config
                galsim.config.Process(config)

if __name__ == "__main__":
    main(sys.argv)
