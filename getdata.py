import lightkurve as lk
import numpy as np
import sys
import datetime
#import ipdb

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

t0 = datetime.datetime(2000,1,1)
jd0 = 2451544.5

#ticid = 'TIC ' + str(sys.argv[1])
id = str(sys.argv[1])
search_result = lk.search_lightcurve(id)#,mission='TESS',author=('SPOC','QLP'))

if len(search_result) == 0:
    print("No light curves found for " + id + ". Try a different name.")
    sys.exit()

lcs = search_result.download_all()
i=0
for lc in lcs:

    lc = lc.remove_nans()
#    lc = lc.remove_outliers() # this often clips the transits
#    lc = lc.flatten() # this can remove transit signals
    lc = lc.normalize()

    author = search_result[i].author[0] # SPOC, QLP, etc
    exptime = str(int(search_result[i].exptime[0].value)).zfill(4)

    if author == "Kepler":
        bjd_offset = 2454833.0
        sector = 'Q' + search_result[i].mission[0][-2:]
        filter = 'Kepler'
        telescope = 'Kepler'
    elif author == "TESS-SPOC" or author=='QLP' or author=='SPOC':
        bjd_offset = 2457000.0
        sector = 'S' + str(lc.sector).zfill(2)
        filter = 'TESS'
        telescope = 'TESS'
    else:
        print("WARNING: Skipping lightcurve with unrecognized author: " + author)
        continue
    
    time = np.array(lc.time.value) + bjd_offset
    flux = np.array(lc.flux.value)
    err = np.array(lc.flux_err.value)

    datestr = (t0 + datetime.timedelta(days=time[0]-jd0)).strftime('n%Y%m%d')

    # create the filename in EXOFASTv2 format
    filename = datestr + '.' + filter + '.' + telescope + '.' + id + '.' + sector + '.' + exptime + '.' + author + '.dat' 

    np.savetxt(filename,np.column_stack([time,flux,err]))

    i = i+1
