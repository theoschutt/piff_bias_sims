select
    distinct(model.expnum),
    e.band,
    e.mjd_obs,
    e.telra,
    e.teldec,
    model.psf_fwhm
from prod.piff_hsm_model_qa model
    join proctag t on model.pfw_attempt_id=t.pfw_attempt_id
    join exposure e on model.expnum=e.expnum
    join QA_SUMMARY qa on model.expnum=qa.expnum
where
    t.tag = 'Y6A2_PIFF_V3_TEST'
    and model.flag = 0
    and model.nstar >= 30;
