import lightkurve as lk
import numpy as np
import sys

# this uses lightkurve to extract QLP/SPOC TESS data and format it for EXOFASTv2
# call it like this (python 3)
# python getdata.py 76989773
# where 76989773  is the TIC ID

# see https://docs.lightkurve.org/ 
# Please acknowledge
# This research made use of Lightkurve, a Python package for Kepler and TESS data analysis \citep{lightkurve:2018}.
#@MISC{lightkurve:2018,
#   author = {{Lightkurve Collaboration} and {Cardoso}, J.~V.~d.~M. and
#             {Hedges}, C. and {Gully-Santiago}, M. and {Saunders}, N. and
#             {Cody}, A.~M. and {Barclay}, T. and {Hall}, O. and
#             {Sagear}, S. and {Turtelboom}, E. and {Zhang}, J. and
#             {Tzanidakis}, A. and {Mighell}, K. and {Coughlin}, J. and
#             {Bell}, K. and {Berta-Thompson}, Z. and {Williams}, P. and
#             {Dotson}, J. and {Barentsen}, G.},
#    title = "{Lightkurve: Kepler and TESS time series analysis in Python}",
# keywords = {Software, NASA},
#howpublished = {Astrophysics Source Code Library},
#     year = 2018,
#    month = dec,
#archivePrefix = "ascl",
#   eprint = {1812.013},
#   adsurl = {http://adsabs.harvard.edu/abs/2018ascl.soft12013L},
#}


ticid = 'TIC ' + str(sys.argv[1])
search_result = lk.search_lightcurve(ticid,mission='TESS',author=('SPOC','QLP'))

lcs = search_result.download_all()
i=0
for lc in lcs:

    lc = lc.remove_nans()
    lc = lc.remove_outliers()
    lc = lc.flatten()
    lc = lc.normalize()
    
    time = np.array(lc.time.value) + 2457000.0
    flux = np.array(lc.flux.value)
    err = np.array(lc.flux_err.value)

    filename = 'n20000101.TESS.TESS.' + 'TIC' + str(search_result[i].target_name.data[0]) + '.S' + str(lc.sector).zfill(2) + '.' + str(int(search_result[i].exptime[0].value)).zfill(4) + '.' + lc.author + '.dat'

    np.savetxt(filename,np.column_stack([time,flux,err]))

    i = i+1
