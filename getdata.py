'''
this uses lightkurve to extract TESS/Kepler data and format it for EXOFASTv2

The first argument is the SIMBAD-resolvable star name

The second (optional) argument is the transit depth, in
fraction. Values >= 1 - depth will never be clipped. 
Default = 0.03 = 3%.

The third (optional) argument is the sigma
clipping. Negative values will skip clipping.
Default=5. 

examples:

python getdata.py TIC76989773

# more aggressive clipping
python $EXOFAST_PATH/getdata.py WASP-43 0.0255 3

# disable sigma clipping
python $EXOFAST_PATH/getdata.py WASP-43 0 -1

see https://docs.lightkurve.org/ 

Please acknowledge

This research made use of Lightkurve, a Python package for Kepler and TESS data analysis \citep{lightkurve:2018}.

@MISC{lightkurve:2018,
   author = {{Lightkurve Collaboration} and {Cardoso}, J.~V.~d.~M. and
             {Hedges}, C. and {Gully-Santiago}, M. and {Saunders}, N. and
             {Cody}, A.~M. and {Barclay}, T. and {Hall}, O. and
             {Sagear}, S. and {Turtelboom}, E. and {Zhang}, J. and
             {Tzanidakis}, A. and {Mighell}, K. and {Coughlin}, J. and
             {Bell}, K. and {Berta-Thompson}, Z. and {Williams}, P. and
             {Dotson}, J. and {Barentsen}, G.},
    title = "{Lightkurve: Kepler and TESS time series analysis in Python}",
 keywords = {Software, NASA},
howpublished = {Astrophysics Source Code Library},
     year = 2018,
    month = dec,
archivePrefix = "ascl",
   eprint = {1812.013},
   adsurl = {http://adsabs.harvard.edu/abs/2018ascl.soft12013L},
}

'''

import lightkurve as lk
import numpy as np
import sys,os
import datetime
import glob
import argparse
#import ipdb

# parse command line arguments
parser = argparse.ArgumentParser(
    prog='getData',
    description='Downloads TESS/Kepler data and formats it for EXOFASTv2')
parser.add_argument('id',help='SIMBAD-resolvable star name')
parser.add_argument('-d','--depth',default=0.03,dest='depth',help='Fractional transit depth. Flux >= 1-depth will not be clipped.')
parser.add_argument('-n','--nsigma',default=5.0,dest='nsigma',help='N sigma clipping. Negative values will skip clipping.')
parser.add_argument('-p','--path',default='.',dest='path',help='path to output files')
parser.add_argument('-u','--undeblend',default=False,action='store_true',dest='undeblend',help='Undo deblending from lightcurves')
parser.add_argument('-v','--verbose',default=False,action='store_true',dest='verbose',help='Display clipping stats')
parser.add_argument('-a','--all',default=False,action='store_true',dest='download_all',help='Download all lightcurves')
parser.add_argument('-o','--overwrite',default=False,action='store_true',dest='overwrite',help='Overwrite previously downloaded files')
args = parser.parse_args()

# undo the deblending applied to TESS/Kepler lightcurves
if args.undeblend:
    # this feature probably won't be used often. don't make it a dependency if not used
    from astroquery.vizier import Vizier
    v = Vizier(columns=['TIC','Ncont','Rcont','_r'],
               catalog='IV/39/tic82').query_object(args.id)[0]
    ndx = np.argmin(v['_r'])
    v = v[ndx]

    if str(v['Rcont']) == '--':
        if args.verbose: print("No deblending to undo")
        contratio = 0.0
    else:
        if args.verbose: 
            print("Undoing deblending for " + str(v['Ncont']) + ' stars (' + str(v['Rcont']) + ')')
        contratio = float(v['Rcont'])
    file_ext = '.undeblended.dat'
else:
    file_ext='.dat'
    contratio = 0.0

t0 = datetime.datetime(2000,1,1)
jd0 = 2451544.5

search_results = lk.search_lightcurve(args.id,author=('SPOC','TESS-SPOC','QLP','TASOC','CDIPS','Kepler'))

if len(search_results) == 0:
    print("No light curves found for " + args.id + ". Name must be SIMBAD-resolveable")
    sys.exit()

# sometimes IDs match to multiple TIC IDs
# warn the user to select by TIC ID
unique_ids = list(set(search_results.target_name))
if len(unique_ids) > 1 and 'TIC' in args.id: 
    match = np.where(search_results.target_name == args.id[3:])
    if len(match) > 0:
        search_results = search_results[match]
        unique_ids = list(set(search_results.target_name))

if len(unique_ids) > 1:
    print("Multiple TIC IDs match " + args.id)
    print(unique_ids)
    print("Specify target by TIC ID")
    sys.exit()

unique_sectors = list(set(search_results.mission))

# only get one light curve per sector
if args.download_all: 
    to_download = list(range(len(search_results)))
else:
    to_download = []
    for sector in unique_sectors:
        match = np.where(search_results.mission == sector)[0]
        if len(match) == 1: to_download.append(match[0])
        if len(match) > 1:
            # prioritize by exptime: 120, 200, 20, 600, then 1800
            match2 = np.where(search_results[match].exptime.value == 120)[0]
            if len(match2) == 0: match2 = np.where(search_results[match].exptime.value == 200)[0]
            if len(match2) == 0: match2 = np.where(search_results[match].exptime.value == 20)[0]
            if len(match2) == 0: match2 = np.where(search_results[match].exptime.value == 600)[0]
            if len(match2) == 0: match2 = np.where(search_results[match].exptime.value == 1800)[0]
            if len(match2) == 1: to_download.append(match[match2[0]])
            if len(match2) > 1:
                # prioritize by author: TESS, TESS-SPOC, QLP, Kepler, CDIPS, then TASOC
                match3 = np.where(search_results[match[match2]].author == 'TESS')[0]
                if len(match3) == 0: match3 = np.where(search_results[match[match2]].author == 'TESS-SPOC')[0]
                if len(match3) == 0: match3 = np.where(search_results[match[match2]].author == 'QLP')[0]
                if len(match3) == 0: match3 = np.where(search_results[match[match2]].author == 'Kepler')[0]
                if len(match3) == 0: match3 = np.where(search_results[match[match2]].author == 'CDIPS')[0]
                if len(match3) == 0: match3 = np.where(search_results[match[match2]].author == 'TASOC')[0]
                if len(match3) == 1: to_download.append(match[match2[match3[0]]])

for search_result in search_results[to_download]:

    author = search_result.author[0] # SPOC, QLP, etc
    exptime = str(int(search_result.exptime[0].value)).zfill(4)
    ticid = 'TIC' + search_result.target_name[0]

    if author == "Kepler":
        bjd_offset = 2454833.0
        sector = 'Q' + search_result.mission[0][-2:]
        filter = 'Kepler'
        telescope = 'Kepler'
    elif author == "TESS-SPOC" or author=='QLP' or author=='SPOC':
        bjd_offset = 2457000.0
        sector = 'S' + str(int(str(search_result.mission[0]).split()[-1])).zfill(2)
        filter = 'TESS'
        telescope = 'TESS'
    elif author=='CDIPS' or author=='TASOC':
        # they don't have flux in the same place. Or errors. Are there any cases where these LCs are better?
        print("WARNING: CDIPS and TASOC LCs are not supported (yet?)")
        continue
        bjd_offset = 2457000.0
        sector = 'S' + str(int(str(search_result.mission[0]).split()[-1])).zfill(2)
        filter = 'TESS'
        telescope = 'TESS'
    else:
        print("WARNING: Skipping lightcurve with unrecognized author: " + author)
        continue
    file_suffix = '.' + filter + '.' + telescope + '.' + args.id + '.' + sector + '.' + exptime + '.' + author + file_ext

    # skip if I've already got it
    if not args.overwrite:
        files = glob.glob(os.path.join(args.path,'*' + file_suffix))
        if len(files) != 0: continue

    lc = search_result.download()
    lc = lc.remove_nans()
    lc *= (1.0 + contratio)
    lc = lc.normalize()

    time = np.array(lc.time.value) + bjd_offset
    flux = np.array(lc.flux.value)
    err = np.array(lc.flux_err.value)

    # replace nans in err with median absolute deviation
    nan = np.where(np.isnan(err) | np.isinf(err))[0]
    maderr = np.median(abs(flux[1:] - flux[0:-1])) * 1.48 /np.sqrt(2.0)
    err[nan] = maderr

    # lopsided 5 sigma clipping (keeping values that are low by transit depth)
    # Should not clip transits or stellar variability
    if args.nsigma < 0: 
        nbad = 0 
        ngood = len(flux)
    else: nbad = 1   

    while nbad != 0:
        rms = np.std(flux)
        median = np.median(flux)
        good = np.where((flux > (median-args.depth-args.nsigma*rms)) & (flux < (median + args.nsigma*rms)))[0]
        ngood = len(good)
        nbad = len(flux)-ngood

        time = time[good]
        flux = flux[good]/median
        err = err[good]/median
        if args.verbose: print((nbad,rms,median))

    # are they all bad?
    if ngood == 0:
        ipdb.set_trace()
        print("WARNING: no good points after sigma clipping")
        print("Skipping lightcurve " + sector + ' ' + exptime + ' ' + author)
        continue

    datestr = (t0 + datetime.timedelta(days=time[0]-jd0)).strftime('n%Y%m%d')

    # create the filename in EXOFASTv2 format
    filename = os.path.join(args.path,datestr + file_suffix)
  
    np.savetxt(filename,np.column_stack([time,flux,err]))
    print("Downloaded " + filename)

#import ipdb
#ipdb.set_trace()
