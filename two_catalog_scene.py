import galsim

class TwoCatalogSceneBuilder(galsim.config.StampBuilder):

    def buildProfile(self, config, base, psf, gsparams, logger):
        
        gal_cat = getInputObj('catalog', config, base, 'twoCatalogScene', num=0)
        star_cat = getInputObj('catalog', config, base, 'twoCatalogScene', num=1)
        
        ngal = gal_cat.nobjects
        nstar = star_cat.nobjects
        
        # setting sequence for Galsim's built-in indexing
        SetDefaultIndex(config, ngal+nstar)
        
        # grabs required config params (with type for parsing)
        req = {'index': int}
        kwargs, safe = GetAllParams(config, base, req=req)
        index = kwargs['index']
        
        # indexing for gal_cat and star_cat
        # set obj_type
        if index < ngal:
            obj_type = 'gal'
        else:
            obj_type = 'star'
            index -= ngal
        
        base[obj_type]['index'] = index # TODO: might need PropagateIndexKeyRNGNum
        
        # make correct object
        obj = galsim.config.BuildGSObject(base, obj_type)
        base['current_obj'] = obj
        
        return galsim.Convolve(obj,psf)

  
galsim.config.stamp.RegisterStampType('TwoCatalogScene', TwoCatalogSceneBuilder())

