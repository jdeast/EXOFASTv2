;+
; NAME:
;   EXOFASTv2
;
; PURPOSE: 
;   This is the high-level code for EXOFASTv2 to perform global
;   models for exoplanetary systems. 
;
;   Please cite Eastman et al., 2019
;   (https://arxiv.org/abs/1907.09480) if you make use of this routine
;   in your research.
;
;   Report errors or bugs to jason.eastman@cfa.harvard.edu
;
; CALLING SEQUENCE:
;
;   exofastv2, PRIORFILE= [RVPATH=, TRANPATH=, FLUXFILE=, PREFIX=,
;             /CIRCULAR, NTHIN=, MAXSTEPS=, /DEBUG, RANDOMFUNC=,
;             SEED=, /SPECPRIORS, /BESTONLY, NINTERP=, EXPTIME=,
;             /LONGCADENCE]
;
; INPUTS:
;
;   PRIORFILE  - An ASCII file with each line containing 2-6 white
;                space delimited columns:
;
;                NAME VALUE WIDTH LOWERLIMIT UPPERLIMIT STARTING_VALUE
;                
;                NAME: must match a parameter.label. If in an array
;                (e.g., of planets), add "_i", where "i" is the
;                zero-indexed index into the array.
;
;                VALUE: If less than 6 columns are supplied, this is
;                the starting point of the fit. If this is the last
;                value specified, this only changes the starting point
;                and no prior is applied. If WIDTH is specified, this
;                is the mean of the gaussian prior. 
; 
;                WIDTH: Standard deviation of a Gaussian prior. If
;                negative, it is ignored (no Gaussian prior). If zero,
;                the parameter is fixed to VALUE.
; 
;                LOWERLIMIT: Lower bound of a uniform prior. If -Inf,
;                it is ignored.
;
;                UPPERLIMIT: Upper bound of a uniform prior. If Inf,
;                it is ignored.
;  
;                STARTING_VALUE: If specified, the fit will start here
;                regardless of the Gaussian prior
; 
;                See
;                $EXOFAST_PATH/parnames.README for a detailed
;                explanation of this file and $EXOFAST_PATH/examples/
;                for templates to use.
;  PREFIX      - Each of the output files will have this string as a
;                prefix. You can specify an absolute or relative
;                path. If the directories do not exist, they will be
;                created. The default is 'fitresults/planet.' to put
;                all outputs in a fitresults subdirectory.
;
; DATA FILE INPUTS
;  RVPATH      - A string specifying the path to the RV data
;                file(s). Wildcards are allowed for multiple files
;                (e.g., '*.rv'). Each file must have at least 3
;                columns:
;                   1) Time (BJD_TDB -- See Eastman et al., 2010)
;                   2) RV (m/s)
;                   3) err (m/s)
;
;                If additional columns are supplied, they will be
;                detrended against.
;
;                The filenames of the RV data files *should* adhere
;                to a certain format:
;
;                planet_name.telescope.whateveryouwant.
;
;                "telescope" is used as the legend text in figures
;                and as header entries in the table.
;
;                NOTE 1: The units must be as specified, or the fit
;                will be wrong or fail. 
;                NOTE 2: Other input time stamps will likely not
;                break the code, but can introduce errors in the
;                reported times of more than a minute. The output
;                table will display BJD_TDB but makes no attempt to
;                convert input times.  See
;                http://adsabs.harvard.edu/abs/2010PASP..122..935E
;                for an explanation of times
;                NOTE 3: If RVPATH is omitted, no RV model will be
;                generated.
;
;  TRANPATH    - A string specifying the path to the transit data
;                file(s). Wildcards are allowed for multiple files
;                (e.g., '*.dat'). Each file must have at least 3
;                columns:
;
;                  1) Time (BJD_TDB -- See Eastman et al., 2010)
;                  2) Normalized flux
;                  3) err
;                  4) Detrend parameter 1
;                  ....
;                  N+3) Detrend parameter N
;       
;                An optional header specifies which parameters
;                should be multiplicatively detrended and which
;                should be additively detrended. Without a header,
;                all detrending parameters are additively detrended.
; 
;                The first character in the header line must be
;                "#". It must have whitespace delimited columns
;                corresponding to the data. Each multiplicative
;                detrending variable is preceeded with "M" (case
;                sensitive). For example, the header line
;
;                #BJD_TDB FLUX ERR MAIRMASS SKY myval
;
;                Will fit a multiplicative detrending to the 4th
;                column (MAIRMASS) and an additive detrending to the
;                5th and 6th columns (SKY and myval).
;
;                The corresponding model will be:
;
;                model = (transit + C0*SKY + C1*myval)*(F0 + M0*MAIRMASS)
;
;                Where C0, C1, and M0 are fitted detrending
;                parameters, F0 is the normalization, and SKY,
;                motion, and MAIRMASS are vectors for each data point
;                to detrend against.
;
;                The names of the files describing the transits
;                *must* adhere to a certain format:
;                nYYYYMMDD.filtername.telescope.whateveryouwant. 
;
;                nYYYYMMDD -- The UTC date of mid transit, where YYYY is
;                the year, MM is the month, and DD is the day. This
;                is only necessary if the data has a single
;                transit. This is used to label the transits in the
;                output plot.
; 
;                filtername -- The name of the observed
;                filter. If NOCLARET is not set, the filter must be
;                among the following list:
;
;                Johnson/Cousins: 'U','B','V','R','I','J','H','K'
;                Sloan: 'Sloanu','Sloang','Sloanr','Sloani','Sloanz'
;                Kepler: 'Kepler'
;                CoRoT: 'CoRoT'
;                TESS: 'TESS'
;                Spitzer: 'Spit36','Spit45','Spit58','Spit80'
;                Stromgren: 'u','b','v','y'
;
;                This is used to define the limb darkening for the
;                transit, dilution, and secondary depth
;                parameters. 
;                
;                NOTE: If NOCLARET is set, you may specify anything
;                here. If a number where "." is replaced with "_", it
;                will be interpreted as a wavelength in microns. If
;                TDELTAVS is also set, a plot of depth vs wavelength
;                will be created for transmission spectroscopy.
;
;                telescope -- a description of the telescope used for
;                the observations. Anything is allowed, but all
;                observations observed with the same telescope should
;                have the same name. This is used in the output plot
;                and color codes the TTV plot.
;
;                whateveryouwant -- any string you want to
;                include for it to make sense to you. This is not
;                used by the code.
;
;                So a transit taken on UTC 2017-01-27 with MINERVA in
;                the V band would be "n20170127.V.MINERVA.dat"
;
;                NOTE 2: The units must be as specified, or the fit
;                will be wrong or fail.
;                NOTE 3: Other input time stamps will likely not
;                break the code, but can introduce errors in the
;                reported times of more than a minute. The output
;                table will display BJD_TDB but makes no attempt to
;                convert input times.  See
;                http://adsabs.harvard.edu/abs/2010PASP..122..935E
;                for an explanation of times
;                NOTE 4: If omitted, just the RV data will bit fit
;
;   *** Astrometry is not currently suppported or thoroughly tested. Future updates are likely to change the expected formats ***
;  ASTROMPATH  - A string specifying the path to the astrometry data
;                file(s). Wildcards are allowed for multiple files
;                (e.g., '*.dat'). Each file must have 5 or 8 columns:
;
;                   1) Time (BJD_TDB -- See Eastman et al., 2010)
;                   2) RA (ICRS deg)
;                   3) DEC (ICRS deg)
;                   4) ra uncertainty (mas)
;                   5) dec uncertainty (mas)
;                   6) X position of observatory (AU from Solar System Barycenter)
;                   7) Y position of observatory (AU from Solar System Barycenter)
;                   8) Z position of observatory (AU from Solar System Barycenter)
;
;                The filename must adhere to a certain format:
;                epoch.location.whateveryouwant                
; 
;                epoch -- the reference epoch of the
;                observations in BJD_TDB (e.g., J2000 is
;                2451545).
;                location -- the location of the observatory. Only
;                "Earth" (geocenter), "Gaia" and "Hipparcos" are
;                understood, in which case the observatory position
;                (columns 6-8) will be calculated if not supplied.
;                whateveryouwant -- any string you want to include
;                for it to make sense to you. This is not used by the
;                code.
;
;                If supplied, EXOFASTv2 will fit the following
;                additional parameters: RA, Dec, PMRA, PMDEC, and
;                parallax, plus Omega for each planet. The
;                inclination will span the range from 0 to 180 (-1 <
;                cosi < 1).
;
;  DTPATH      - The (optional) path to the Doppler Tomography fits
;                file(s). If supplied, the code will fit a vsini and
;                macro turbulence of the star, as well as a projected
;                spin-orbit alignment for each planet. This must be
;                used in conjunction with the FITDT, an NPLANETS array
;                specifying which planets should have their DT signal
;                modeled.
;
;                Each fits file is a 2D array describing all DT
;                observations during a single transit with extensions
;                describing the axes of the array. The pixel values of
;                the 2D array contain the fractional flux decrement at
;                a given BJD_TDB (spectrum) and velocity (pixel value
;                of the CCF). The first extension should specify the
;                BJD_TDB corresponding to each Y pixel and the second
;                extension should specify the velocity corresponding
;                to each X pixel.
;
;                Such a file can be generated given the 2D array of
;                fractional flux decrements (DT), a time array
;                (BJD_TDB) and a velocity array (VEL, in km/s) like this:
; 
;                writefits,'nYYYYMMDD.pname.instrument.resolution.fits',DT
;                writefits,'nYYYYMMDD.pname.instrument.resolution.fits',BJD_TDB, /append
;                writefits,'nYYYYMMDD.pname.instrument.resolution.fits',VEL, /append
;
;                The names of these files *must* adhere to a certain format:
;                nYYYYMMDD.pname.instrument.resolution.whateveryouwant.fits 
;
;                nYYYYMMDD -- The UTC date of mid transit. This is
;                used to label the transits in the output plot.
; 
;                pname -- The name of the planet, including the
;                         letter for the planet order (e.g., KELT-24b)
;
;                instrument -- the name of the instrument used for
;                the observations. Anything is allowed, but all
;                observations observed with the same telescope should
;                have the same name. This is used in the labels.
;
;                resolution -- The R value of the instrument. This is
;                required to accurately model the DT signal and must
;                be correct.
;
;                whateveryouwant - any string you want to include (or
;                nothing) for it to make sense to you (e.g, target
;                name). This is not used by the code.
;
;                So a transit taken on UTC 2017-01-27 with the TRES
;                spectrograph (R=44000) would be called
;                "n20170127.TRES.44000.fits"
; 
;                NOTE 0: When supplied this way, the DT data does not
;                constrain vsini. A gaussian prior must be given in
;                the prior file.
;
;                NOTE 1: NAXIS1 must equal the number of velocities and
;                NAXIS2 must equal the number of times.
;
;                NOTE 2: Not using BJD_TDB time stamps will likely
;                not break the code, but can introduce errors in the
;                reported times of more than a minute or cause
;                internal inconsistencies if different data sets use
;                different timestamps. The output table will display
;                BJD_TDB but makes no attempt to convert or reconcile
;                input times between data sets.
;
;                See
;                http://adsabs.harvard.edu/abs/2010PASP..122..935E 
;                for an explanation of times
;
;                NOTE 3: If DTPATH is omitted, Doppler Tomography is
;                not included in the global fit
;
; SED MODEL INPUTS:
;  FLUXFILE    - An ASCII file with each line containing 3-5 white
;                space delimited columns: 
;
;                FILTER MAG UNCERTAINTY CATALOG_UNCERTAINTY STARNDX
;
;                FILTER: The filter name. See
;                $EXOFAST_PATH/sed/mist/filternames2.txt for
;                a complete list of allowed filters.
;
;                MAG: The apparent magnitude in the specified filter
;
;                UNCERTAINTY: The uncertainty in the apparent magnitude
;
;                CATALOG_UNCERTAINTY: Often catalog errors ignore
;                systematic error. This is the catalog error for 
;                reference, but is ignored. 
;
;                STARNDX: The index of the star(s) this magnitude
;                corresponds to. Multiple stars can be specified as a
;                comma delimited list (e.g., "0,1"). Differential
;                magnitudes can be specified with a "-" (e.g.,
;                "2-0,1"). This is ignored when the FLUXFILE model.
;
;                NOTE: a parallax prior should be specified or the SED
;                model will simply determine a photometric parallax.
;
;                $EXOFAST_PATH/mkticsed.pro will automatically
;                generate this file using queryVizier to query trusted
;                catalogs with appropriate systematic error floors
;                based on the object's name or coordinates. Exercise
;                caution that the correct object was selected for each
;                catalog.
; 
;                The FLUXFILE SED model assumes top hat filter
;                profiles to sum NEXTGEN stellar atmospheres. It does
;                not support Gaia bands, multiple stars, differential
;                photometry, or spectrophotometry. It is no longer
;                recommended for general use but is left for
;                comparison.
;
;  MISTSEDFILE - Same as FLUXFILE, but uses the MIST bolometric
;                correction tables based on the C3K atmospheric
;                grids. Supports multiple stars and differential
;                photometry, and a much wider range of filters, but
;                not spectrophotometry.
;
;  SEDFILE     - Same as FLUXFILE, but uses detailed filter profiles
;                to determine the model flux. The user can supply
;                their own filter profiles or generate any from SVO
;                using
;                $EXOFAST_PATH/sed/filtercurves/getfilter.pro. This
;                model supports multiple stars, differential
;                photometry, and spectrophotometry.
;
;  SPECPHOTPATH- The path to spectrophotometry files. Can only be used
;                in conjunction with the SEDFILE SED model.
;
;  NOAVPRIOR   - By default, when multiple stars are fit, the
;                exinction (AV) is assumed to be a monotonically
;                increasing function of distance. Set this keyword to
;                allow any AV regardless of distance.
; 
;  FBOLSEDFLOOR- The SED model's constraint on the bolometric flux is
;                typically dominated by systematic uncertainty. This
;                sets that floor. Default is 0.024 (2.4%). See
;                Tayar+2022.
;
;  TEFFSEDFLOOR- The SED model's constraint on the effective
;                temperature is typically dominated by systematic
;                uncertainty. This sets that floor. Default is 0.02
;                (2.0%). See Tayar+2022.
;
;  FEHSEDFLOOR - The SED model's constraint on the metalicity
;                occasionally exceeds the systematic uncertainty. This
;                sets that floor. Default is none. If set, a new
;                FEHSED is fit.
;
;  ONED        - For the FLUXFILE and SEDFILE models, the reddened
;                NEXTGEN atmospheres are trilinearly interpolated in
;                3D in Teff, logg, and [Fe/H]. If set, the nearest
;                gridpoints in logg and [Fe/H] are used (typically
;                ~0.5 dex), and the interpolation is only done in
;                Teff. Use of this keyword is not recommended.
; 
; EVOLUTIONARY MODEL INPUTS:
;
;  NOMIST     - An NSTARS boolean array specifying not to use the MIST
;               evolutionary tracks see
;               http://waps.cfa.harvard.edu/MIST/ And please cite
;               http://adsabs.harvard.edu/abs/2016ApJS..222....8D
;               http://adsabs.harvard.edu/abs/2016ApJ...823..102C
;
;  PARSEC     - An NSTARS boolean array specifying to use the PARSEC
;               evolutionary tracks. This should
;               be accompanied by the NOMIST keyword (but is not
;               enforced).
;
;  YY         - An NSTARS boolean array specifying to use the YY
;               evolutionary tracks see
;               http://adsabs.harvard.edu/abs/2001ApJS..136..417Y to
;               constrain the mass/radius of the star. YY should not
;               be used for low-mass stars (~< 0.5 msun). This should
;               be accompanied by the NOMIST keyword (but is not
;               enforced).
;
;  TORRES     - An NSTARS boolean array specifying to use the Torres
;               relations to constrain the mass and radius of the
;               star. This may be useful to investigate potential
;               systematics and should be accompanied by the NOMIST
;               keyword (but is not enforced). The Torres relations
;               are not applicable for low mass (< 0.6 msun) stars.
;
;  MANNRAD    - An NSTARS boolean array specifying to use the Mann+
;               2015 Ks-Rstar relations to constrain the radius of the
;               star. A prior on the apparent Ks band magnitude (e.g.,
;               from 2MASS) and parallax should be given, and should
;               be accompanied by the NOMIST keyword. The Mann
;               relations are not applicable to high mass (> 0.7 msun)
;               stars. When the transit-based constraint on rhostar
;               exceeds ~5%, use of this relation is not recommended.
;
;  MANNMASS   - An NSTARS boolean array specifying to use the Mann+
;               2019 Ks-Mstar relations to constrain the mass of the
;               star. A prior on the apparent Ks-band magnitude (e.g.,
;               from 2MASS) and parallax should be given, and should
;               be accompanied by the NOMIST keyword. The Mann
;               relations are not applicable to high mass (> 0.7 msun)
;               stars.
;
;  TEFFEMFLOOR- When using MIST, PARSEC, or YY, the evolutionary
;               model's constraint on the effective temperature is
;               typically dominated by systematic uncertainty. This
;               sets that floor. The default is a function of stellar
;               mass and evolutionary model. See Tayar+2022.
;
;  FEHEMFLOOR - When using MIST, PARSEC, or YY, the evolutionary
;               model's constraint on the metallicity is typically
;               dominated by systematic uncertainty. This sets that
;               floor. Default is 0.08 dex.
;
;  RSTAREMFLOOR
;            - When using MIST, PARSEC, or YY, the evolutionary
;              model's constraint on the stellar radius is typically
;              dominated by systematic uncertainty. This sets that
;              floor. The default is a function of stellar mass and
;              evolutionary model. See Tayar+2022.
;
;  AGEEMFLOOR- The evolutionary model's constraint on the stellar age
;              is typically dominated by systematic uncertainty
;              (maybe?). This sets that floor. The default is a
;              function of stellar mass and evolutionary model.
;
; BEER MODEL INPUTS:
;
;  FITTHERMAL- An NTRANSITS boolean array specifying which transits to
;              fit thermal emission for. This is what you want to set
;              to model an isolated secondary eclipse. When set, the
;              transit will be modeled with a baseline of 1 between t2
;              and t3 of the secondary eclipse and 1 + thermal
;              emission (in PPM) out of eclipse.
;
;              NOTE: earlier versions of this array specified the band
;              name and all similar bands were linked together. This
;              functionality can be recreated by linking the THERMAL
;              parameters in the prior file for common bands.
;
;  FITELLIP  - A string array specifying which bands to fit an
;              ellipsoidal variation amplitude for. Note: this does
;              *not* feedback into the global model through F&M 2011,
;              eq 2
;              NOTE: this should be (has been?) converted like FITTHERMAL
;
;  FITREFLECT- A string array specifying which bands to fit reflected
;              light for. Set this along with thermal if you're
;              fitting a full phase curve. It will be modeled as a
;              sinusoid with the orbital period, a minimum at the
;              primary transit, and a fitted amplitude (in PPM).
;
;  FITPHASE  - An NPLANETS boolean array that specifies if a phase
;              shift to the sinusoidal reflection signal should be
;              fit. This appears to be broken. Use with caution.
;
;  FITBEAM   - An NPLANETS boolean array that specifies which planets
;              to model a doppler beaming signal. If this beaming
;              amplitude is fit directly, it will *not* constrain the
;              planetary mass (set DERIVEBEAM instead).
;
;  DERIVEBEAM- An NPLANET boolean array that specifies which planets
;              to model a doppler beaming signal. If this beaming
;              amplitude is derived (from K), it will constrain the
;              planetary mass through Faigler & Mazeh, 2011, eq 1 (set
;              FITBEAM instead to model the beaming but break this
;              connection).
;
; STAR INPUTS:
;
;  NSTARS    - The number of stars to model. Default=1. Must be larger
;              than 0.
; 
;  STARNDX   - An NPLANETS long array that specifies the index of the
;              star each planet orbits. The default is 0 for all.
; 
;  SEDDEBLEND- An NTRANSITSxNSTARS boolean array specifying which
;              transits are blended with which stars. These will
;              automatically be deblended according to the SED models
;              of all stars specified here, assuming an error of 10%.
;              NOTE: TESS lightcurves from SPOC are already deblended
;              with known companions (in TICv8.1). They can be
;              undeblended with the -u option to
;              $EXOFAST_PATH/getdata.py.
;
;  FITDILUTE -A string array specifying which bands to fit a dilution
;             term for. It will fit the fractional contribution from
;             the companion. This should be set if the star is blended
;             with a neighbor and you expect color-dependent depth
;             variations in the transit or photocenter variations in
;             the astrometry.
;
;             If you just have the flux ratio, X, solve this system of
;             equations for Flux_companion:
;         
;             Flux_companion/Flux_primary = X
;             Flux_companion + Flux_primary = 1
;             Flux_companion = X/(X+1)
;
;             Note 1: For transits, this is likely to be degenerate
;             with F0 (transit normalization). You probably need to
;             apply a prior (DILUTE_#, where # corresponds to the band
;             index) from some external information to constrain it
;             (e.g., an SED fit). 
; 
;             Note 2: For astrometry, this is degenerate with the
;             semi-major axis/mass. You must supply a prior on the
;             dilution (e.g., from AO) in at least one band if one of
;             those is not independently determined.
;
;             Note 3: This dilution does not impact the SED fit. The
;             SED broadband photometry is assumed to be for the
;             single, primary component.
;
; PLANET INPUTS:
;
;  NPLANETS  - The number of planets you wish to fit to the
;              data. Default is 1. If 0 is set, it will just fit the
;              star.
;
;  FITTRAN   - An NPLANETS boolean array that specifies which planets
;              should be fit with a transit model.  If TRANPATH is
;              specified, default is bytarr(nplanets) + 1B.  If
;              TRANPATH is not specified, default is bytarr(nplanets).
;              At least one of FITRV and FITTRAN should be true for
;              each planet
;
;  FITRV     - An NPLANETS boolean array that specifies which planets
;              should be fit with an RV model.  If RVPATH is
;              specified, default is bytarr(nplanets) + 1B.  If RVPATH
;              is not specified, default is bytarr(nplanets).  At
;              least one of FITRV and FITTRAN should be true for each
;              planet
;
;  ROSSITER  - An NPLANETS boolean array that specifies which planets
;              should fit the Rossiter-McLaughlin effect to the RV
;              data using the Ohta approximation.  See
;              https://ui.adsabs.harvard.edu/abs/2005ApJ...622.1118O/abstract
;              If a run was dedicated to RM, it should be separated
;              into its own RV file to fit a different zero point and
;              jitter parameter to it.  Default is bytarr(nplanets)
;
;  FITDT    - An NPLANETS boolean array that specifies which planets
;             should be fit using a Doppler Tomography model.  Default
;             is bytarr(nplanets)
;
;  CIRCULAR - An NPLANETS boolean array that specifies which planets
;             should be fixed to be circular (1) or left with
;             eccentricity free (0). Default is all eccentric planets
;             bytarr(nplanets)
;
;  TIDES    - If set, when (1-Rstar/a-rp/a) < e < (1-3*Rstar/a), we
;             set the eccentricity to zero, presuming that the tidal
;             circularization timescale is much much smaller than the
;             age of the system.
;
;  ALLOWORBITCROSSING
;           - By default, models that have projected orbits within
;             each others' hill spheres are rejected, presuming
;             they are unstable and constraining eccentricity. This
;             may be too strict when mutual inclinations are
;             significant (it would exclude Pluto/Neptune). Set this
;             keyword to disable this constraint and allow orbits to
;             cross one another.
;
;  CHEN     - An NPLANETS boolean array that specifies which planets
;             should have the Chen & Kipping, 2017 mass-radius
;             relation applied. By default CHEN = FITRV xor
;             FITTRAN. That is, only apply the mass-radius prior when
;             RV is not fit (to derive the planet mass) or when a
;             transit is not fit (to derive the planet radius). If the
;             defaults have been overridden and FITRV and CHEN are
;             both false for a given planet, the RV semi-amplitude
;             (planet mass) and all derived parameters will not be
;             quoted. Multi-planet systems will be constrained not to
;             cross orbits, but the Hill Sphere will be set to zero.
;             If FITTRAN and CHEN are both false for a given planet,
;             the planetary radius and all derived parameters will not
;             be quoted.
;
;  I180     - An NPLANETS boolean array specifying which planets'
;             inclination should be allowed to be between 0 and 180
;             instead of the default of 0 to 90. Note that a typical
;             transiting planet has a perfect degeneracy between i and
;             180-i which is likely to cause convergence
;             problems. Additional information (e.g., from astrometry
;             or mutual eclipses) must be used for this keyword to
;             have the desired impact. Even in the case of mutual
;             eclipses (which is currently not supported), at least
;             one planet must be arbitrarily constrained from 0 to 90.
;
; RV INPUTS
;
;  FITSLOPE - If set, it will fit a linear trend to the RV data.
;
;  FITQUAD  - If set, it will fit a quadratic trend to the RV
;             data. FITSLOPE is automatically set if FITQUAD is set.
;
;  RVEPOCH  - This is the BJD_TDB of a reference epoch for RV linear
;             and quadratic trends. Default is the midpoint of all RV
;             observations.
;
;             RV = planet + gamma + SLOPE*(t-RVEPOCH) + QUAD*(t-RVEPOCH)^2
; 
; TRANSIT INPUTS
;
;  NOCLARET - An NTRANSITFILES boolean array that, if set, will ignore
;             the Claret & Bloeman limb darkening tables and just fit
;             the limb darkening. This should be specified for
;             low-mass stars where these tables are unreliable or for
;             high SNR transits where they can constrain the LD
;             directly.
;
;  TTVS     - An NTRANSITFILESxNPLANETS boolean array, specifying which
;             transit files should have a timing offset applied
;             corresponding to which planets. 
;             NOTE: if the same transit applies to more than one
;             planet, the same offset will be applied to both. Such
;             transits should probably be removed from the analysis.
;
;  TIVS     - An NTRANSITFILESxNPLANETS boolean array, specifying which
;             transit files should have an inclination offset applied
;             corresponding to which planets.  
;             NOTE: if the same transit applies to more than one
;             planet, the same offset will be applied to both. Such
;             transits should probably be removed from the analysis.
;
;  TDELTAVS - An NTRANSITFILESxNPLANETS boolean array, specifying which
;             transit files should have a depth offset applied
;             corresponding to which planets.  NOTE: if the same
;             transit applies to more than one planet, the same offset
;             will be applied to both. Such transits should probably
;             be removed from the analysis.
;
;  LONGCADENCE
;           - An NTRANSITFILES boolean array, which sets/overwrites
;             EXPTIME=29.425 and NINTERP=10 to handle Kepler long
;             cadence data. This is measurably different from TESS
;             long cadence data (EXPTIME=30, NINTERP=10).
; 
;  EXPTIME  - An NTRANSITFILES array of exposure times, in minutes,
;             for each transit. The model data point will be an
;             average of NINTERP samples over a duration of EXPTIME
;             minutes centered on the input time. Generally, a
;             sampling at least every 3 minutes is recommended. For a
;             dicussion on binning, see:
;             http://adsabs.harvard.edu/abs/2010MNRAS.408.1758K
;
;  NINTERP  - An NTRANSITFILES long array for the number of points to
;             average throughout the EXPTIME. If set, the each model
;             data point will be an average of NINTERP samples over a
;             duration of EXPTIME minutes centered on the input
;             time. Generally, a sampling at least every 3 minutes is
;             recommended.
;
;  REJECTFLATMODEL
;           - An NTRANSITFILES byte array specifying which transits
;             must have a signal. If model without a signal is
;             generated and REJECTFLATMODEL==1B, the model will be
;             rejected a priori. If a flat model is accepted when a
;             transit is expected, there is an infinite volume of
;             parameter space with equal likelihood, which breaks the
;             Markov condition, and it will never converge. Note the
;             default allows flat models because 1) the user may wish
;             to model out of transit constraints on the ephemeris and
;             2) rejecting flat models a priori may give the user an
;             undue confidence in the signficance of the signal. This
;             is likely required for low SNR transits or when using
;             parallel tempering (see NTEMP and TF).
;
; NOPRIMARY - An NPLANETS byte array specifying which planets are
;             allowed to not have a primary transit. Normally, b>1+p
;             models are rejected because the transit fit is
;             unconstrained and poorly behaved, but in the rare cases
;             a planet only has a secondary and not a primary, this is
;             required. You probably also want to specify
;             REQUIRESECONDARY=NOPRIMARY.
;
; REQUIRESECONDARY - 
;             An NPLANETS byte array specifying which planets are
;             required to have a secondary eclipse. Models with bs>1+p
;             will be rejected. When fitting a secondary, this could
;             bias non-significant detections. When only fitting a
;             secondary, this may be required to keep the fit
;             reasonably bounded.
;
;  FITSPLINE- An NTRANSITFILES byte array specifying which transits
;             should be flatted with Andrew Vanderburg's
;             keplerspline. This should only be used for long baseline
;             data like Kepler, K2, or TESS.
;
;  SPLINESPACE
;           - An NTRANSITFILES array specifying the knot spacing
;             between adjacent points in the spline, in days. To be
;             used in conjunction with FITSPLINE. Default is 0.75
;             days. This spacing should be large compared to the
;             transit duration (>~ 3x), or there is significant risk
;             of introducing a strong covariance between the spline
;             and the transit parameters, biasing the inferred values,
;             inflating the uncertainties, and increasing convergence
;             times. It should be short compared to the variability
;             you wish to remove.
;
;  FITRAMP  - An NTRANSITFILES byte array specifying which transits
;             should be fit with an exponential ramp (common for JWST
;             and Spitzer LCs). When set, the corresponding lightcurve
;             will be multiplied by 1+A*exp((time[0]-time)/tau), where
;             A and tau (unique to each LC) are fitted parameters
;             reported alongside other transit parameters and time is
;             the user-supplied time from the transit file. A is
;             typically negative.
;
;  FITWAVELET
;           - An NTRANSITFILES array specifying which to fit Carter's
;             wavelet error to handle correlated errors. See
;             https://ui.adsabs.harvard.edu/abs/2009ApJ...704...51C/abstract
;
; REPARAMETERIZATION INPUTS:
;
;  FITLOGMP - By default, planet mass is parameterized as MP. Set this
;             keyword to parameterize it as LOGMP. This may impose a
;             more realistic prior, but is problematic for low SNR RVs
;             where the small-planet volume is infinite. It also
;             excludes negative masses which biases the mass high.
;
;  NOVCVE   - When only transits are fit, eccentricity is
;             parameterized as VCVE. Set this to keep the usual
;             sqrt(e)cos(omega) sqrt(e)sin(omega) parameterization.
;
;  NOCHORD  - When only transits are fit, cosi is parameterized as
;             chord. Set this to keep cosi.
;
;  FITSIGN  - When only transits are fit, e is parameterized as VCVE,
;             but solving for e requires a quadratic. By default, a
;             sign is fit to choose the solution. Set this to zeros to
;             use the degenerate parameter L to determine the sign. 
;             See Eastman, 2023 (https://arxiv.org/abs/2309.14410)
; 
;  FITTT    - TT is the minimum projected separation, the intuitive
;             defintion of the transit time, but one that must be
;             numerically computed. By default, we fit for the time of
;             conjunction, TC. For eccentric, inclined orbits, this
;             can introduce significant non-linear covariances and it
;             may be faster to take the hit to derive TT.
;
;  EARTH    - If set, the output units of Mp and Rp are in Earth
;             units, not Jupiter units.
;
; PLOTTING INPUTS
;  TRANSITRANGE
;           - A six element array that defines the
;             [XMIN,XMAX,YMIN,YMAX,O-C_YMIN,O-CYMAX] plotting limits
;             for the transit figure. Any value set to !values.d_nan
;             will use the default value. This allows the user a level
;             of customization of the plot, but was designed to be
;             used in conjunction with mkfitgif to ensure the plot
;             ranges are constant throughout the fit.
;
;  RVRANGE
;           - A six element array that defines the
;             [XMIN,XMAX,YMIN,YMAX,O-C_YMIN,O-CYMAX] plotting limits
;             for the RV figure. Any value set to !values.d_nan will
;             use the default value. This allows the user a level of
;             customization of the plot, but was designed to be used
;             in conjunction with mkfitgif to ensure the plot ranges
;             are constant throughout the fit.
;
;  SEDRANGE
;           - A six element array that defines the
;             [XMIN,XMAX,YMIN,YMAX,O-C_YMIN,O-CYMAX] plotting limits
;             for the SED figure. Any value set to !values.d_nan will
;             use the default value. This allows the user a level of
;             customization of the plot, but was designed to be used
;             in conjunction with mkfitgif to ensure the plot ranges
;             are constant throughout the fit.
;
;  EMRANGE
;           - A four element array that defines the
;             [XMIN,XMAX,YMIN,YMAX] plotting limits of the
;             evolutionary model (MIST, PARSEC, YY) figure. Any value
;             set to !values.d_nan will use the default value. This
;             allows the user a level of customization of the plot,
;             but was designed to be used in conjunction with mkfitgif
;             to ensure the plot ranges are constant throughout the
;             fit.
;
; DEBUGGING INPUTS
;  DEBUG    - When this keyword is set, plots will be displayed to the
;             screen for each model. This is intended to help identify
;             problematic parameters that need additional constraints,
;             though it's often too slow to be helpful. Generally,
;             VERBOSE serves this purpose better. It can also
;             help build intuition about how the different stages of
;             EXOFASTv2 work. 
;
;  VERBOSE  - When this keyword is set, more detailed information is
;             printed to the screen and logfile at each step,
;             including specific reasons why a model was rejected
;             (which parameters were out of bounds), the values of all
;             fitted parameters, and the chi^2 components for each
;             fit. This is often very helpful in identifying
;             problematic parameters that need additional constraints.
;
;  DELAY    - If set, the chi^2 function runs a dummy for loop to
;             count to the DELAY value specified, arbitrarily
;             inflating the runtime of the chi^2 computation. This is
;             designed to benchmark the performance of multi-threading
;             and is not intended for general use.
;
; MCMC INPUTS:
;  
;  NOTE: The MCMC portion of the fit will automatically stop when one
;  of the following conditions are met:
;
;           GR < MAXGR and TZ > MINTZ at 5 consecutive evaluations
;               (and DONTSTOP is not set)
;           STEPS > MAXSTEPS
;           TIME > MAXTIME
;           !STOPNOW=1
;
;  MAXSTEPS - The maximum number of steps to take in the MCMC
;             chain. Note that a 32-bit installation of IDL cannot
;             allocate more than 260 million double-precision numbers,
;             and redundant copies of each parameter are
;             required. Even a 64-bit installation may require very
;             slow disk swapping. A very large number for MAXSTEPS
;             will cause memory management problems. The default
;             chooses the number that will use 1 GB of RAM, and larger
;             values are strongly discouraged (increase NTHIN instead
;             if the chains are not well-mixed). Decrease MAXSTEPS for
;             preliminary runs.
;
;             NOTE1: If the MCMC chains run to MAXSTEPS, doubling this
;             value will double the runtime. Short test runs
;             (MAXSTEPS=100) are strongly encouraged before running
;             week-long fits.
;
;  NTHIN    - If set, only every NTHINth element will be kept. High
;             values typically don't degrade the resultant fit because
;             there is the correlation length between adjacent steps
;             is long. This has the advantage of improved memory
;             management, reduced storage, and faster generation of
;             the final plots. Default is roughly set to pass the
;             mixing criteria, but this is highly uncertain. Revising
;             this parameter based on preliminary fits is strongly
;             encouraged.
;
;             NOTE: Only saved links in the chain count toward
;             MAXSTEPS, so if the MCMC chains run to MAXSTEPS,
;             doubling this value will double the runtime.
;
;  MAXTIME  - The maximum runtime for the DEMC portion of the code, in
;             seconds. If set, the fit will begin wrapping up after
;             MAXTIME seconds. This is intended for non-interactive
;             use with a time constraint (e.g., running on a super
;             computer). Note this only applies to the MCMC portion of
;             the code, and it finishes the loops over NTEMPS and
;             NTHIN before it will check the time. The time AMOEBA
;             takes or the time to synthesize the results is not
;             counted, therefore, it would be wise to set this at
;             least 1800 seconds (30 minutes) less than any hard
;             limit.
;
;  MAXGR    - The maximum Gelman Rubin statistic that is considered
;             well-mixed (default=1.01). 
;
;  MINTZ    - The minimum number of independent draws that is
;             considered well-mixed (default=1000). The fit will
;             automatically stop when GR < MAXGR and TZ > MINTZ at 5
;             consecutive evaluations.
;
;  DONTSTOP - A keyword that, when set, the fit will take MAXSTEPS
;             steps regardless of GR and TZ statistics. This is
;             statistically unnecessary, but sometimes useful to
;             generate prettier output plots (smoother contours and
;             PDFs).
;
;  !STOPNOW - This global keyword can be set during a fit to stop a
;             fit early and have it summarize the results. This is
;             helpful for prelimary or problematic runs. During a run
;             that you want to stop, type:
;
;                 control + c 
;                 !stopnow=1 
;                 .con
;
;             NOTE: This is not possible with license-free use.
;
;  NTEMPS   - This enables parallel tempering and sets the number of
;             chains to run at different temperatures, uniformly
;             distributed between 1 and TF. The more temperatures, the
;             easier it is to make swaps between temperatures, but the
;             longer it can take. The default is 1 (no parallel
;             tempering). 8 is a good value, but more may be required
;             for high dimensional fits. Note that only the T=1 chain
;             is kept by default. Parallel tempering is better at
;             finding the global minimum in a rough likelihood surface
;             and sampling widely spaced, multi-modal distributions
;             (often the case for MIST models). Nominally, NTEMPS more
;             steps are taken, requiring NTEMPS times longer to
;             run. However, many models are rejected as out of bounds,
;             and it may lead to better mixing and therefore faster
;             convergence, so the effect on runtime is system
;             dependent and not obvious. You are strongly encouraged
;             to use NTEMPS=8 for initial fits in order to optimize
;             the starting location of the next fit.
;
;  TF       - The temperature of the hottest chain, when NTEMPS >
;             1. The higher the temperature, the more parameter space
;             it will probe, but the chains will be less likely to
;             swap. Default is 200, which can robustly find minima
;             separated by 50 sigma. A value of 1 will result in no
;             temperature difference between chains.
;
;  KEEPHOT  - By default, the hot chains (and their corresponding
;             chi^2) from parallel tempering (see NTEMPS and TF) are
;             dynamically discarded to save memory. If this keyword is
;             set and NTEMPS > 1, the steps from the hot chains
;             (HOTPARS) and the corresponding chi2 (HOTCHI2) will be
;             stored in memory during the fit and saved in the output
;             IDL file. This is not useful to determine the posterior,
;             and should only be used for debugging (e.g., determining
;             the extent of parameter space explored), as it
;             significantly increases the memory footprint.
;
;  RANDOMFUNC
;           - A string specifying the name of the random number
;             generator to use. This generator must be able to return
;             1,2 or 3 dimensional uniform or normal random
;             deviates. 
;             
;             When using IDL versions < v8.2, IDL's built-in RANDOMU
;             is based on numerical recipe's RAN1, and is not good
;             enough for MCMC. In this case, we use 'EXOFAST_RANDOM'
;             (which is based on Numerical Recipe's RAN3 and is slow
;             but robust).
;
;             When using IDL versions >= v8.2, we use IDL's
;             built-in RANDOMU.
;
;  SEED     - The seed to the random number generator used for the
;             MCMC fits. Be sure not to mix seeds from different
;             generators. The default is -systime(/seconds). This is
;             helpful if the user wishes to generate the exact same
;             sequence of "random" numbers, the proposed steps will
;             quickly diverge if anything else is changed, but this
;             can be particularly useful make bugs reproducible.
;
;  STRETCH  - By default, EXOFASTv2 uses a differential evolution
;             method to generate the next proposed set of
;             parameters. Set this keyword to use the Affine Invariant
;             "stretch" step instead (popularized by emcee). In our
;             tests, STRETCH is inferior to DEMC and is not generally
;             recommended.
;
;  NTHREADS - The number of threads to use during the MCMC portion of
;             the fit. For short preliminary runs, often the overhead
;             of initializing the threads exceeds its benefit. But for
;             longer runs, especially for complex models (that take >~
;             0.005 seconds), the scaling is near perfect up to the
;             number of cores on your machine (4 cores => NTHREAD=4 =>
;             4x faster). The default is !CPU.HW_NCPU, meaning it will
;             use all available cores on your machine for a dramatic
;             speed up. For shared computers/servers, this is likely
;             an inappropriate use and counterproductive for you
;             (requesting more CPU power than is available adds
;             overhead without benefit).
;
; GENERAL INPUTS:
;
;  SKIPTT   - A keyword that skips the expensive numerical computation
;             of TT for each step at the end. This is rarely useful
;             unless doing TTVs.
;
;  USERNOTE - A string that is printed to the top of the log file,
;             intended to help the user keep track of fits. 
;
;  MKSUMMARYPG
;           - A keyword that generates a quick look summary page of
;             all MCMC plots and chain outputs at the end of the
;             run. This spawns several commands to the terminal,
;             including gs, ps2pdf, grep, and convert, so
;             probably only works for a fraction of linux users, no
;             windows users, and maybe some mac users.
;
;  NOCOVAR  - A keyword that skips the generation of the covariance
;             corner plot. This is somewhat expensive, large, and ugly for
;             high dimensional fits. GDL users cannot generate this
;             anyway. Users interested in better corner plots should
;             look at $EXOFAST_PATH/remake_corner.pro.
;
;  PLOTONLY - Will stop after generating START models. Useful for
;             setting up an initial fit.
;
;  BESTONLY - Will stop after generating AMOEBA models. Useful for
;             setting up an initial fit.
;
;  BADSTART - A return value that is true if the starting model was
;             bad. Useful to check for/handle failures during batch
;             fits.
;  
; OUTPUTS:
;
;   Each of the output files will be preceeded by PREFIX (defined
;   above).
;
;   mcmc.idl   - An IDL save file that contains the stellar structure
;                with the full chains for each parameter, including
;                derived parameters, the chi2 at each link, and the
;                index of the burn-in period, called MCMCSS (MCMC
;                stellar structure).
;   chain.ps   - A postscript plot of each parameter as a function of
;                chain link, 8 to a page. Each color is a separate
;                chain.
;   pdf.ps     - A postscript plot of each posterior distribution
;                function, 8 to a page. Each color is a separate
;                chain. The thick black line is the average of all
;                chains.
;   covar.ps   - A postscript plot of the covariances between each
;                fitted parameter in a corner plot.
;   median.tex - The LaTeX source code for a deluxe table of the
;                median values and 68% confidence interval, rounded to
;                two sig figs in the uncertainty.
;   median.csv - A machine readable table of the median values and 68%
;                confidence interval, rounded to two sig figs in the
;                uncertainty.
;   log        - A file containing everything printed to screen during
;                the fit for later review.
;   
;   Each of the following outputs are created at 3 distinct phases in
;   the fit. 
;
;   1) Files preceeded with "start." are the starting values of the
;   fit. Often times it is necessary to tweak these by hand until the
;   starting values are a pretty good match to the data.
;
;   2) Files preceeded with "amoeba." are the models after the amoeba
;   minimization is complete. If this model is not a good fit to the
;   data, the fit is likely to fail and you should refine your
;   starting values, priors, and constraints. At a minimum, the MCMC
;   is likely to take significantly longer to converge if this is not
;   a good fit.
;
;   3) Files preceeded with "mcmc." are the final plots meant for
;   publication. They represent the best fit among all links of all
;   chains.
;
;   transit.ps - A multi-page postscript plot of the transit
;                model. The first page has each transit (file) stacked
;                on the same plot, offset by a constant. Subsequent
;                pages are the detrended, phase-folded transits for
;                each planet.
;   rv.ps      - A multi-page postscript plot of the RV model. The
;                first page is the unphased model with gamma and any
;                fitted linear or quadratic terms subtracted from both
;                the data and model. If more than one instrument is
;                used, a legend is produced in the corner. Subsequent
;                pages are each planet, phased to its period, with all
;                other planets, gamma, and any fitted linear and
;                quadratic terms subtracted.
;   mist.eps   - A plot of the star with its MIST isochrone
;                overplotted. The black point is the best-fit value,
;                the red point is the corresponding model value. Only
;                generated if the /NOMIST keyword is not set.
;   yy.eps     - A plot of the star with its YY isochrone
;                overplotted. The black point is the best-fit value,
;                the red point is the corresponding model value. Only
;                generated if the /YY keyword is not set.
;   sed.eps    - A plot of the broadband photometry and best fit
;                SED. Only generated if FLUXFILE is given.
; 
;  We recognize that in many cases, the outputs plots are not
;  publication quality. The following text outputs are intended to
;  allow the user to generate their own, more complicated or
;  customized plots without having to recreate the model or
;  detrending. These files are placed in a modelfiles subdirectory within
;  the output directory.
;
;  residuals.telescope_[i].txt --
;      A text file containing the residuals to the best-fit model RV
;      for the ith telescope.
;   model.telescope_[i].txt --
;      A text file containing the best-fit model (with gamma, slope,
;      and quad subtracted) RVs for each time given for telescope.
;      One for each telescope is generated.
;   model.transit_[i].txt - 
;      The model of the ith transit at each input data point,
;      including detrending.
;   residuals.transit_[i].txt -- 
;      The residuals of the ith transit file at each data point.
;   detrendedmodel.transit_[i].planet_[j].txt -- 
;      The detrended model of the ith transit of the jth planet minus
;      1 at each input data point. Sum all transit[i].planet_*.txt and add
;      one for the combined model. Useful for separating overlapping
;      transits.
;   prettymodel.transit_[i].planet_[j].txt -- 
;      Same as detrendedmodel.transit_[i].planet_[j].txt, but sampled
;      every minute from the first to the last data point. This is
;      intended to make a pretty model plot if there are gaps in the
;      data.
;
;   NOTE 1: The detrended data can by derived by adding the residuals
;      to the detrended model

;   NOTE 2: To extract a single page out of a multi-page PS file, use:
;
;         psselect -p# input.ps output.ps 
;
;      where # is the (1-indexed) page number to extract.
;
; COMMON BLOCKS:
;   CHI2_BLOCK:
;     SS      - A structure that describes the entire stellar system
;
; EXAMPLES:
;   See the $EXOFAST_PATH/examples directory for several full examples.
;
; MODIFICATION HISTORY
; 
;  2017/01 -- Complete rewrite of exofast.pro. More general (fits
;             multiple planets, mulitple bands, multiple instrumental
;             offsets). Now easily extensible.
;  2023/12 -- Missed a lot of updates here (see git
;             history). Documentation cleanup.
;  2024/01 -- Add FITRAMP (primarily for JWST transits)
;-
pro exofastv2, priorfile=priorfile, $
               prefix=prefix,$
               ;; data file inputs
               rvpath=rvpath, tranpath=tranpath, $
               astrompath=astrompath, dtpath=dtpath, $
               ;; SED model inputs
               fluxfile=fluxfile,mistsedfile=mistsedfile,$
               sedfile=sedfile,specphotpath=specphotpath,$
               noavprior=noavprior,$
               fbolsedfloor=fbolsedfloor,teffsedfloor=teffsedfloor,$
               fehsedfloor=fehsedfloor, oned=oned,$
               ;; evolutionary model inputs
               yy=yy, nomist=nomist, parsec=parsec, $
               torres=torres, mannrad=mannrad,mannmass=mannmass, $
               teffemfloor=teffemfloor, fehemfloor=fehemfloor, $
               rstaremfloor=rstaremfloor,ageemfloor=ageemfloor,$
               ;; BEER model inputs
               fitthermal=fitthermal, fitellip=fitellip, $
               fitreflect=fitreflect, fitphase=fitphase, $
               fitbeam=fitbeam, derivebeam=derivebeam, $
               ;; star inputs
               nstars=nstars, starndx=starndx, $
               seddeblend=seddeblend, fitdilute=fitdilute, $
               ;; planet inputs
               nplanets=nplanets, $
               fittran=fittran, fitrv=fitrv, $
               rossiter=rossiter, fitdt=fitdt, $
               circular=circular, tides=tides, $ 
               alloworbitcrossing=alloworbitcrossing, $
               chen=chen, i180=i180, $
               ;; RV inputs
               fitslope=fitslope, fitquad=fitquad, rvepoch=rvepoch,$
               ;; transit inputs
               noclaret=noclaret, $
               ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs,$
               longcadence=longcadence, exptime=exptime, ninterp=ninterp, $
               rejectflatmodel=rejectflatmodel,$
               noprimary=noprimary, requiresecondary=requiresecondary,$
               fitspline=fitspline, splinespace=splinespace, $
               fitramp=fitramp, fitwavelet=fitwavelet, $
               ;; reparameterization inputs
               fitlogmp=fitlogmp,$
               novcve=novcve, nochord=nochord, fitsign=fitsign, $
               fittt=fittt, earth=earth, $
               ;; plotting inputs
               transitrange=transitrange,rvrange=rvrange,$
               sedrange=sedrange,emrange=emrange, $
               ;; debugging inputs
               debug=debug, verbose=verbose, delay=delay, $
               ;; MCMC inputs
               maxsteps=maxsteps, nthin=nthin, maxtime=maxtime, $
               maxgr=maxgr, mintz=mintz, $
               dontstop=dontstop, $
               ntemps=ntemps, tf=tf, keephot=keephot, $
               randomfunc=randomfunc, seed=seed,$
               stretch=stretch, $
               nthreads=nthreads, $              
               ;; General inputs
               skiptt=skiptt, $
               usernote=usernote, $
               mksummarypg=mksummarypg,$
               nocovar=nocovar, $
               plotonly=plotonly, bestonly=bestonly, $
               badstart=badstart
               
;; this is the stellar system structure
COMMON chi2_block, ss

;; if a virtual machine or runtime license, read the arguments from args.txt
if lmgr(/vm) or lmgr(/runtime) then begin
   ;; IDL_IDLBridge disabled in the virtual machine
   ;; no multi-threading without a license :(
   nthreads = 1L 

   par = command_line_args(count=numargs)
   if numargs eq 1 then begin
      argfile = par[0]
   endif else argfile = 'args.txt'

   if not file_test(argfile) then message, argfile + ', containing desired arguments to EXOFASTv2, does not exist'
   readargs, argfile, priorfile=priorfile, $
             prefix=prefix,$
             rvpath=rvpath, tranpath=tranpath, $
             astrompath=astrompath, dtpath=dtpath, $
             fluxfile=fluxfile,mistsedfile=mistsedfile,$
             sedfile=sedfile,specphotpath=specphotpath,$
             fbolsedfloor=fbolsedfloor,teffsedfloor=teffsedfloor, $
             fehsedfloor=fehsedfloor, oned=oned,$
             yy=yy, nomist=nomist, parsec=parsec, $
             torres=torres, mannrad=mannrad, mannmass=mannmass, $
             teffemfloor=teffemfloor, fehemfloor=fehemfloor, $
             rstaremfloor=rstaremfloor,ageemfloor=ageemfloor,$
             fitthermal=fitthermal, fitellip=fitellip, $
             fitreflect=fitreflect, fitphase=fitphase, $
             fitbeam=fitbeam, derivebeam=derivebeam, $
             nstars=nstars,starndx=starndx, $
             seddeblend=seddeblend,fitdilute=fitdilute, $
             nplanets=nplanets, $
             fittran=fittran, fitrv=fitrv, $
             rossiter=rossiter, fitdt=fitdt, $
             circular=circular, tides=tides, $ 
             alloworbitcrossing=alloworbitcrossing, $
             chen=chen, i180=i180, $
             fitslope=fitslope, fitquad=fitquad, rvepoch=rvepoch, $
             noclaret=noclaret, $
             ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs, $
             longcadence=longcadence, exptime=exptime, ninterp=ninterp, $
             rejectflatmodel=rejectflatmodel,$
             noprimary=noprimary, requiresecondary=requiresecondary,$
             fitspline=fitspline, splinespace=splinespace, $
             fitramp=fitramp, fitwavelet=fitwavelet, $              
             fitlogmp=fitlogmp,$
             novcve=novcve, nochord=nochord, fitsign=fitsign, $
             fittt=fittt, earth=earth, $             
             transitrange=transitrange,rvrange=rvrange,$
             sedrange=sedrange,emrange=emrange, $
             debug=debug, verbose=verbose, delay=delay,$
             maxsteps=maxsteps, nthin=nthin, maxtime=maxtime, $
             maxgr=maxgr, mintz=mintz, $
             dontstop=dontstop, $              
             ntemps=ntemps,tf=tf,keephot=keephot,$
             randomfunc=randomfunc, seed=seed,$
             stretch=stretch,$
             skiptt=skiptt, $
             usernote=usernote,$
             mksummarypg=mksummarypg,$
             nocovar=nocovar,$
             plotonly=plotonly, bestonly=bestonly,$
             logname=logname

endif

;; default prefix for all output files (filename without extension)
if n_elements(prefix) eq 0 then prefix = 'fitresults/planet.'
basename = file_basename(prefix)

;; output to log file and the screen
logname = prefix + 'log'
file_delete, logname, /allow_nonexistent
if n_elements(usernote) ne 0 then printandlog, usernote, logname

;; name of the chi square function
chi2func = 'exofast_chi2v2'

;; compile all routines now to keep output legible 
;; resolve_all doesn't interpret execute; it's also broken prior to IDL v6.4(?)
defsysv, '!GDL', exists=runninggdl  

;; default to NCORES threads, if we're running a full copy of IDL
if runninggdl or lmgr(/runtime) or lmgr(/vm) then begin
   nthreads=1
endif else if n_elements(nthreads) ne 1 then begin
   nthreads = !cpu.hw_ncpu
endif ;; else use the user's input   

if double(!version.release) ge 6.4d0 and ~lmgr(/vm) and ~lmgr(/runtime) and ~runninggdl then $
   resolve_all, resolve_either=[chi2func,'exofast_random','ramp_func'],skip_routines=['cggreek'],/cont,/quiet

;; output to log file too
logname = prefix + 'log'
file_delete, logname, /allow_nonexistent


exofast_path = getenv('EXOFAST_PATH')
;if (strpos(exofast_path,'//') ne -1) or 
if (strpos(exofast_path,'EXOFASTv2') eq -1) then begin
   printandlog, "ERROR: EXOFAST_PATH (" + exofast_path + ") not set properly",logname
   printandlog, "Typically, EXOFAST_PATH is something like ${HOME}/idl/EXOFASTv2/", logname
;   printandlog, "NOTE: '//' is not allowed as it is inconsistenly resolved", logname
   return
endif

;; if the directory doesn't exist, make it
dirname = file_dirname(prefix)
if dirname ne '.' then begin
   if ~file_test(dirname,/directory) then begin
      file_mkdir, dirname
   endif
endif

modeldirname = dirname + path_sep() + 'modelfiles'
if ~file_test(modeldirname,/directory) then begin
   file_mkdir, modeldirname
endif

;; some error checking on NTHREADS (multi-threading)
if nthreads gt !cpu.hw_ncpu then begin
   printandlog, "WARNING: Using more threads (" + strtrim(nthreads,2) + ") than physical cores (" + strtrim(!cpu.hw_ncpu,2) + "); this is likely to be inefficient", logname
;endif else if nthreads eq !cpu.hw_ncpu then begin
;   printandlog, "WARNING: Doing fit with " + strtrim(nthreads,2) + " threads. This may impact performance of other tasks (See NTHREADS argument).", logname
endif else begin
   printandlog, "Doing fit with " + strtrim(nthreads,2) + ' threads', logname
endelse

;; insert the commit id into the log
spawn, 'git -C $EXOFAST_PATH rev-parse HEAD', output, stderr
if stderr[0] eq '' then printandlog, "Using EXOFASTv2 commit " + output[0], logname

;; create the master structure 
;; keyword inheritance would be helpful here, but it might break multi-threading
ss = mkss(priorfile=priorfile, $
          prefix=prefix,$
          ;; data file inputs
          rvpath=rvpath, tranpath=tranpath, $
          astrompath=astrompath, dtpath=dtpath, $
          ;; SED model inputs
          fluxfile=fluxfile, mistsedfile=mistsedfile, $
          sedfile=sedfile, specphotpath=specphotpath,$
          noavprior=noavprior,$
          fbolsedfloor=fbolsedfloor,teffsedfloor=teffsedfloor,$
          fehsedfloor=fehsedfloor, oned=oned,$
          ;; evolutionary model inputs
          yy=yy, nomist=nomist, parsec=parsec, $
          torres=torres, mannrad=mannrad, mannmass=mannmass,$
          teffemfloor=teffemfloor, fehemfloor=fehemfloor, $
          rstaremfloor=rstaremfloor,ageemfloor=ageemfloor,$
          ;; BEER model inputs
          fitthermal=fitthermal, fitellip=fitellip, $
          fitreflect=fitreflect, fitphase=fitphase,$
          fitbeam=fitbeam, derivebeam=derivebeam, $
          ;; star inputs
          nstars=nstars, starndx=starndx, $
          seddeblend=seddeblend, fitdilute=fitdilute, $
          ;; planet inputs
          nplanets=nplanets, $
          fittran=fittran,fitrv=fitrv,$
          rossiter=rossiter, fitdt=fitdt,$ 
          circular=circular, tides=tides, $
          alloworbitcrossing=alloworbitcrossing,$
          chen=chen, i180=i180,$
          ;; RV inputs
          fitslope=fitslope, fitquad=fitquad, rvepoch=rvepoch, $
          ;; transit inputs
          noclaret=noclaret,$
          ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs, $
          longcadence=longcadence, exptime=exptime, ninterp=ninterp, $
          rejectflatmodel=rejectflatmodel,$
          noprimary=noprimary, requiresecondary=requiresecondary,$
          fitspline=fitspline, splinespace=splinespace, $
          fitramp=fitramp, fitwavelet=fitwavelet, $            
          ;; reparameterization inputs
          fitlogmp=fitlogmp,$
          novcve=novcve, nochord=nochord, fitsign=fitsign, $
          fittt=fittt, earth=earth, $
          ;; plotting inputs
          transitrange=transitrange,rvrange=rvrange,$
          sedrange=sedrange,emrange=emrange, $
          ;; debugging inputs
          debug=debug, verbose=verbose, delay=delay, $
          ;; internal inputs
          chi2func=chi2func, $
          logname=logname)

if (size(ss))[2] ne 8 then begin
   badstart=1
   return
endif else badstart=0

npars = ss.npars
nfit = n_elements((*(ss.tofit))[0,*])

;; this is where threads were previously initialized

if n_elements(mintz) eq 0 then mintz = 1000d0
if n_elements(maxgr) eq 0 then maxgr = 1.01d0

;; by default, use 1 GB of RAM and set NTHIN so it converges
if n_elements(maxsteps) eq 0 then begin
   maxsteps = round(1024d0^3/(double(ss.nchains)*npars*8d0))
   
   ;; set to 10x the limit expected from an idealized Gaussian case of N parameters
   if n_elements(nthin) eq 0 then begin
      a = [173881.7852504589d0,-3934.6300584134d0,512.5700554749d0]
      nthin = ceil((a[0]+a[1]*nfit+a[2]*nfit^2)/(2d0*nfit)/maxsteps*10d0) > 1
      printandlog, 'NTHIN set to ' + strtrim(nthin,2) + ' for the fit to converge.', logname
      printandlog, 'This is an extremely rough estimate and varies wildly depending on the ',logname
      printandlog, 'details of the fit and the accuracy of the starting conditions.',logname
      printandlog, 'You should monitor the progress of the fit and increase NTHIN if necessary.', logname
      printandlog, '', logname
   endif
endif 

memrequired = double(ss.nchains)*double(maxsteps)*npars*8d0/(1024d0^3)
printandlog, 'Fit will require ' + strtrim(memrequired,2) + ' GB of RAM for the final structure', logname
if memrequired gt 2d0 then begin
   printandlog, 'WARNING: this likely exceeds your available RAM and may crash after the end of a very long run. You likely want to reduce MAXSTEPS and increase NTHIN by the same factor. If you would like to proceed anyway, type ".con" to continue', logname
   if ~lmgr(/vm) then stop
endif
printandlog, '', logname

if n_elements(nthin) eq 0 then nthin = 1L
if nthin lt 1L then nthin=1L

printandlog, 'MAXSTEPS set to ' + strtrim(maxsteps,2), logname
printandlog, 'NTHIN set to ' + strtrim(nthin,2), logname

pars = str2pars(ss,scale=scale,name=name, angular=angular)

;; plot the data + starting guess
modelfile = prefix + 'start'
ss.verbose=1B
bestchi2 = call_function(chi2func, pars, psname=modelfile)
ss.verbose = keyword_set(verbose)
if ~finite(bestchi2) then begin
   printandlog, 'Starting model is out of bounds; cannot recover. You must change the starting parameter(s) via the prior file.', logname
   printandlog, 'Re-running starting model with /VERBOSE flag to identify the parameter', logname
   ss.verbose=1B
   bestchi2 = call_function(chi2func, pars, psname=modelfile)
   printandlog, 'Starting model is out of bounds; cannot recover. You must change the starting parameter(s) via the prior file.', logname
   badstart=1
   return
endif else begin
   printandlog, 'The loglike of the starting model was ' + strtrim(-bestchi2/2d0,2), logname
endelse
if keyword_set(plotonly) then return

;; do it again for accurate timing 
;; after loading all the files into cache, not including plotting
t0 = systime(/seconds)
bestchi2 = call_function(chi2func, pars)
modeltime = systime(/seconds)-t0

;; these are the starting values for all step parameters
printandlog, 'These are the starting values for all fitted parameters', logname
printandlog, 'and the amoeba stepping scale, which is roughly the range', logname
printandlog, 'of parameter space it will explore around the starting value', logname
printandlog, 'and is equal to 3x any Gaussian width. When priors are not ', logname
printandlog, 'specified in ' + priorfile + ', a default guess is used.', logname
printandlog, 'The parameter number is useful for tracking down unconstrained',logname
printandlog, 'parameters', logname
printandlog,'',logname 
printandlog, '**************************************************************', logname
printandlog, '*** IT IS WISE TO MAKE SURE THESE AGREE WITH YOUR          ***', logname
printandlog, '*** EXPECTATION FROM THE PRIORFILE. IF NOT, YOUR           ***', logname
printandlog, '*** STARTING PRIORS MAY NOT BE TRANSLATED CORRECTLY INTO   ***', logname
printandlog, '*** THE FITTED PARAMETERIZATION BY                         ***', logname
printandlog, '*** $EXOFAST_PATH/pars2step.pro. NOT ALL PARAMETER         ***', logname
printandlog, '*** COMBINATIONS ARE ALLOWED/SUPPORTED.                    ***', logname
printandlog, '*** WHEN IN DOUBT, SET PRIORS/CHANGE STARTING VALUES       ***', logname
printandlog, '*** DIRECTLY IN THESE FITTED PARAMETERS.                   ***', logname
printandlog, '**************************************************************', logname
printandlog, '', logname
printandlog, 'Par #      Par Name    Par Value       Amoeba Scale', logname
for i=0, n_elements(name)-1 do printandlog, string(i, name[i], pars[i], scale[i], format='(i3,x,a15,x,f14.6,x,f14.6)'), logname
printandlog, '', logname

if ss.debug and ~lmgr(/vm) then begin
   printandlog, 'program halted to give you time to inspect the priors. Type ".con" to continue', logname
   stop
end

nmax = 1d5
printandlog, 'It takes ' + strtrim(modeltime,2) + ' seconds to calculate a single model', logname
printandlog, 'Beginning AMOEBA fit; this may take up to ' + string(modeltime*nmax/60d0,format='(f0.1)') + ' minutes if it takes the maximum allowed steps (' + strtrim(nmax,2) + ')', logname

;; do the AMOEBA fit
ss.amoeba = 1B
ss.delay =0
best = exofast_amoeba(1d-5,function_name=chi2func,p0=pars,scale=scale,nmax=nmax)

ss.delay = delay
if best[0] eq -1 then begin
   printandlog, 'ERROR: Could not find best combined fit; adjust your starting values and try again. You may want to set the /DEBUG keyword.', logname
   return
endif
printandlog, 'Finished AMOEBA fit', logname
save, best, filename=prefix + 'amoeba.idl'

;; update the parameter array with the chen-derived logks/rp (is this necessary?)
;best = str2pars(ss,scale=scale,name=name) 

;; try again?
;printandlog, 'restarting AMOEBA with chen enabled', logname
;printandlog, call_function(chi2func,best,modelrv=modelrv,modelflux=modelflux, psname=prefix + 'model'), logname
;best = exofast_amoeba(1d-8,function_name=chi2func,p0=best,scale=scale,nmax=nmax)
;printandlog, call_function(chi2func,best,modelrv=modelrv,modelflux=modelflux, psname=prefix + 'model'), logname

;; output the best-fit model fluxes/rvs
bestchi2 = call_function(chi2func,best,modelrv=modelrv,modelflux=modelflux, psname=prefix + 'amoeba')
printandlog, 'The best loglike found by AMOEBA was ' + strtrim(-bestchi2/2d0,2), logname
printandlog, 'It should only be compared against the loglike of the same model with different starting points', logname

;; initialize the threads
if nthreads gt 1 then begin
   printandlog, 'Initializing ' + strtrim(nthreads,2) + ' threads', logname

   ;; load the stellar structure into the common block for each thread
   ;; (can't pass structures between threads, so we have to create it in each thread)
   thread_array = replicate(create_struct('obridge',obj_new("IDL_IDLBridge", output=''),$
                                          'j',-1L, 'm',-1L, 'k', -1L, $
                                          'newpars',dblarr(nfit),'status',0B, 'fac', 1d0,$
                                          'start', systime(/seconds)),nthreads)
   
   ;; get the current working directory (make sure all threads start there)
   cd, current=cwd

   for i=0L, nthreads-1 do begin
      ;; replicate copies the IDLBridge by reference, not by value!! reinitialize here
      thread_array[i].obridge = obj_new("IDL_IDLBridge", output='')
      
      ;; can't share a structure directly with a thread
      ;; must share components and create the structure within the thread
      ;; share every variable in memory to make it future proof
      help, output=helpoutput
      for j=0L, n_elements(helpoutput)-2 do begin
         ;; done with variables, we're done
         if helpoutput[j] eq 'Compiled Procedures:' then break
         
         ;; undefined variable; skip it
         if strpos(helpoutput[j],'UNDEFINED') eq 16 then continue

         entries = strsplit(helpoutput[j],/extract)

         ;; either a long variable name that spans two lines, or not a variable
         if strpos(helpoutput[j],'=') ne 26 then begin

            ;; not a variable; skip it
            if strpos(helpoutput[j+1],'=') ne 26 then continue

            ;; undefined variable, skip it
            if strpos(helpoutput[j+1],'UNDEFINED') eq 16 then continue

            ;; if it's not a lone (long) variable name, skip it
            if n_elements(entries) ne 1 then continue

            ;; this line is the name of a long variable name
            ;; skip the next line, which is its value
            j++

         endif else if n_elements(entries) gt 4 then begin
            ;; catch N-dimentional arrays and strings with spaces
            if strpos(helpoutput[j],'Array') ne 28 and strpos(helpoutput[j],'STRING') ne 16 then continue
         endif else if n_elements(entries) ne 4 then continue            
         ;; declare it in the thread
         ;; EXECUTE is ok here, since we can't use VM with IDLBridge anyway
         varname = entries[0]
         junk = execute("thread_array[i].obridge->setvar,'" + varname + "'," + varname)
;         print, 'setting ' + varname + ' to '
;         junk = execute('print, ' + strtrim(varname,2))
      endfor

      ;; disable NaN warnings inside each thread
      thread_array[i].obridge->setvar,'!except',0

      ;; make sure each thread is run from the current working directory
      thread_array[i].obridge->setvar,'cwd',cwd
      thread_array[i].obridge->execute,'cd, cwd'
      
      ;; compile all the codes in each thread so compilation messages don't pollute the screen
      if double(!version.release) ge 6.4d0 and ~lmgr(/vm) and ~lmgr(/runtime) and ~runninggdl then $
         thread_array[i].obridge->execute, "resolve_all, resolve_either=[chi2func,'exofast_random','ramp_func'], resolve_procedure=['exofastv2'],skip_routines=['cggreek'],/cont,/quiet"
         
      ;; create the stellar stucture within each thread
      thread_array[i].obridge->execute,$
         'ss = mkss(priorfile=priorfile, prefix=prefix,'+$
         'rvpath=rvpath, tranpath=tranpath,'+$
         'astrompath=astrompath, dtpath=dtpath,'+$
         'fluxfile=fluxfile, mistsedfile=mistsedfile,'+$
         'sedfile=sedfile,specphotpath=specphotpath,'+$
         'noavprior=noavprior,'+$
         'fbolsedfloor=fbolsedfloor,teffsedfloor=teffsedfloor,'+$
         'fehsedfloor=fehsedfloor, oned=oned,'+$
         'yy=yy, nomist=nomist, parsec=parsec,'+ $
         'torres=torres, mannrad=mannrad, mannmass=mannmass,'+$         
         'teffemfloor=teffemfloor, fehemfloor=fehemfloor,'+$
         'rstaremfloor=rstaremfloor,ageemfloor=ageemfloor,'+$
         'fitthermal=fitthermal, fitellip=fitellip,'+ $
         'fitreflect=fitreflect,fitphase=fitphase,'+ $
         'fitbeam=fitbeam, derivebeam=derivebeam,'+ $
         'nstars=nstars,starndx=starndx,'+ $         
         'seddeblend=seddeblend, fitdilute=fitdilute,'+$
         'nplanets=nplanets,'+$
         'fittran=fittran, fitrv=fitrv,'+$
         'rossiter=rossiter,fitdt=fitdt,'+$
         'circular=circular,tides=tides,'+$
         'alloworbitcrossing=alloworbitcrossing,'+$
         'chen=chen, i180=i180,'+$
         'fitslope=fitslope, fitquad=fitquad,rvepoch=rvepoch,'+$
         'noclaret=noclaret,'+$
         'ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs,'+$
         'longcadence=longcadence,exptime=exptime,ninterp=ninterp,'+$ 
         'rejectflatmodel=rejectflatmodel,'+$
         'noprimary=noprimary, requiresecondary=requiresecondary,'+$
         'fitspline=fitspline, splinespace=splinespace,'+$
         'fitramp=fitramp, fitwavelet=fitwavelet,'+$
         'fitlogmp=fitlogmp,'+$
         'novcve=novcve, nochord=nochord, fitsign=fitsign,'+$
         'fittt=fittt, earth=earth,'+$
         'transitrange=transitrange,rvrange=rvrange,'+$
         'sedrange=sedrange,emrange=emrange,'+$
         'debug=debug, verbose=verbose,delay=delay,'+$
         '/silent,'+$
         'chi2func=chi2func,'+$
         'logname=logname)'
   endfor
endif

;; do the MCMC fit
if not keyword_set(bestonly) then begin

   exofast_demcpt_multi, best, chi2func, pars, chi2=chi2,$
                         nthin=nthin,maxsteps=maxsteps, maxtime=maxtime, $
                         ntemps=ntemps, tf=tf, dontstop=dontstop, $
                         burnndx=burnndx, goodchains=goodchains, seed=seed, randomfunc=randomfunc, $
                         gelmanrubin=gelmanrubin, tz=tz, maxgr=maxgr, mintz=mintz, $
                         stretch=stretch, logname=logname, angular=angular, $
                         keephot=keephot, hotpars=hotpars, hotchi2=hotchi2, thread_array=thread_array
   
   if pars[0] eq -1 then begin
      printandlog, 'MCMC Failed to find a stepping scale. This usually means one or more parameters are unconstrained by the data or priors.', logname
   endif

   bad = where(tz lt mintz or gelmanrubin gt maxgr,nbad)
   if bad[0] ne -1 then begin
      printandlog, 'WARNING: The Gelman-Rubin statistic indicates ' + $
                   'the following parameters are not well-mixed', logname
      printandlog, '    Parameter   Rz     Tz', logname
      for i=0, nbad-1 do printandlog, string(name[bad[i]], gelmanrubin[bad[i]],tz[bad[i]], format='(a13,x,2(f0.4,x))'), logname
   endif
   printandlog, 'Synthesizing results; for long chains and/or many fitted parameters, this may take up to 15 minutes', logname

   ;; combine all chains
   sz = size(pars)
   npars = sz[1]
   nsteps = sz[2]
   nchains = sz[3]
   pars = reform(pars,npars,nsteps*nchains)
   chi2 = reform(chi2,nsteps*nchains)
   minchi2 = min(chi2,bestndx)
   
   printandlog, 'The best loglike found by MCMC was ' + strtrim(-minchi2/2d0,2), logname
   printandlog, 'It should only be compared against the loglike of the same model with different starting points', logname
   printandlog, '', logname
   printandlog, 'Use BIC and AIC to compare different models', logname
   bic = nfit*alog(ss.ndata) + minchi2
   aic = 2d0*nfit + minchi2
   printandlog, 'NDATA = ' + strtrim(ss.ndata,2),logname
   printandlog, 'NFIT = ' + strtrim(nfit,2),logname
   printandlog, 'BIC = ' + strtrim(bic,2),logname
   printandlog, 'AIC = ' + strtrim(aic,2),logname
   if minchi2 lt bestchi2 then begin
      printandlog, 'WARNING: MCMC found a better model that AMOEBA.', logname
      printandlog, 'Using mkprior to refine your starting values and refitting', logname
      printandlog, 'may result in a faster and more robust answer.', logname
      printandlog, 'Look at the chain plot before you trust this fit.', logname
   endif
endif else begin
   pars = reform(best[tofit],n_elements(tofit),1)
   bestndx = 0
endelse

;; generate the model fit from the best MCMC values, not AMOEBA
bestamoeba = best
best = pars[*,bestndx]
modelfile = prefix + 'mcmc'
ss.verbose = 1
bestchi2 = call_function(chi2func,best,psname=modelfile, $
                         modelrv=modelrv, modelflux=modelflux)
ss.verbose = keyword_set(verbose)

;; make a new stellar system structure with only fitted and derived
;; parameters, populated by the pars array
;mcmcss = mcmc2str(pars, ss)
mcmcss = mkss(priorfile=priorfile, $
              prefix=prefix,$
              ;; data file inputs
              rvpath=rvpath, tranpath=tranpath, $
              astrompath=astrompath, dtpath=dtpath, $
              ;; SED model inputs
              fluxfile=fluxfile, mistsedfile=mistsedfile, $
              sedfile=sedfile, specphotpath=specphotpath,$
              noavprior=noavprior,$
              fbolsedfloor=fbolsedfloor,teffsedfloor=teffsedfloor,$
              fehsedfloor=fehsedfloor, oned=oned,$
              ;; evolutionary model inputs
              yy=yy, nomist=nomist, parsec=parsec, $
              torres=torres, mannrad=mannrad, mannmass=mannmass,$
              teffemfloor=teffemfloor, fehemfloor=fehemfloor, $
              rstaremfloor=rstaremfloor,ageemfloor=ageemfloor,$
              ;; BEER model inputs
              fitthermal=fitthermal, fitellip=fitellip, $
              fitreflect=fitreflect, fitphase=fitphase,$
              fitbeam=fitbeam, derivebeam=derivebeam, $
              ;; star inputs
              nstars=nstars, starndx=starndx, $
              seddeblend=seddeblend, fitdilute=fitdilute, $
              ;; planet inputs
              nplanets=nplanets, $
              fittran=fittran,fitrv=fitrv,$
              rossiter=rossiter, fitdt=fitdt,$ 
              circular=circular, tides=tides, $
              alloworbitcrossing=alloworbitcrossing,$
              chen=chen, i180=i180,$
              ;; RV inputs
              fitslope=fitslope, fitquad=fitquad, rvepoch=rvepoch, $
              ;; transit inputs
              noclaret=noclaret,$
              ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs, $
              longcadence=longcadence, exptime=exptime, ninterp=ninterp, $
              rejectflatmodel=rejectflatmodel,$
              noprimary=noprimary, requiresecondary=requiresecondary,$
              fitspline=fitspline, splinespace=splinespace, $
              fitramp=fitramp, fitwavelet=fitwavelet, $            
              ;; reparameterization inputs
              fitlogmp=fitlogmp,$
              novcve=novcve, nochord=nochord, fitsign=fitsign, $
              fittt=fittt, earth=earth, $
              ;; plotting inputs
              transitrange=transitrange,rvrange=rvrange,$
              sedrange=sedrange,emrange=emrange, $
              ;; debugging inputs
              debug=debug, verbose=verbose, delay=delay, $
              ;; internal inputs
              nvalues=nsteps*nchains,$ 
              /silent, $
              chi2func=chi2func, $
              logname=logname, $
              best=best)

if (size(mcmcss))[2] ne 8 then return

mcmcss.nchains = nchains
mcmcss.burnndx = burnndx
*(mcmcss.goodchains) = goodchains
*(mcmcss.chi2) = chi2

pars2str, pars, mcmcss

;; populate residuals for the mcmcss file
for i=0L, mcmcss.ntran-1 do $
   (*mcmcss.transit[i].transitptrs).residuals = (*ss.transit[i].transitptrs).residuals
for i=0L, mcmcss.ntel-1 do $
   (*mcmcss.telescope[i].rvptrs).residuals = (*ss.telescope[i].rvptrs).residuals
   
;; derive all parameters
derivepars, mcmcss, logname=logname

spawn, 'git -C $EXOFAST_PATH rev-parse --short HEAD', output, stderr
if output[0] ne '' then versiontxt = ", created using EXOFASTv2 commit number " + output[0] $
else versiontxt = ''

;; output filenames
label = "tab:" + basename
caption = "Median values and 68\% confidence interval for " + basename + versiontxt
parfile = prefix + 'pdf.ps'
covarfile = prefix + 'covar.ps'
chainfile = prefix + 'chain.ps'
texfile = prefix + 'median.tex'
csvfile = prefix + 'median.csv'

exofast_plotdist_corner, mcmcss, pdfname=parfile, covarname=covarfile,nocovar=nocovar,logname=logname, csvfile=csvfile
;; save the chains for additional analysis 
;; wait until here because many things are populated in exofast_plotdist_corner
idlfile = prefix + 'mcmc.idl'

;; make a new prior file to start at the best fit found here
priorfileparts = strsplit(priorfile,'.',/extract,/preserve_null)
suffix = priorfileparts[n_elements(priorfileparts)-1]
if valid_num(suffix) then priorfileparts[n_elements(priorfileparts)-1] = strtrim(suffix+1L,2) $
else priorfileparts = [priorfileparts,'2']
priorfile2 = file_dirname(prefix) + path_sep() + file_basename(strjoin(priorfileparts,'.'))
mkprior2, mcmcss=mcmcss, priorfilename=priorfile2

;; GDL compatibility
if runninggdl then begin
   if keyword_set(keephot) and ntemps gt 1 then cmsave, mcmcss, hotpars, hotchi2, filename=idlfile $
   else cmsave, mcmcss, filename=idlfile
endif else begin
   if keyword_set(keephot) and ntemps gt 1 then save, mcmcss, hotpars, hotchi2, filename=idlfile $
   else save, mcmcss, filename=idlfile
endelse

exofast_latextab2, mcmcss, caption=caption, label=label,texfile=texfile
exofast_plotchains, mcmcss, chainfile=chainfile, logname=logname

if (total(mcmcss.ttvs) gt 0 or ~keyword_set(skiptt)) and mcmcss.ntran gt 0 then begin
   printandlog, 'The fit is done and can be interrupted without losing any results', logname
   printandlog, 'Now generating a table of the numerically solved times of ', logname
   printandlog, 'minimum projected separation, depth, and impact parameters for',logname
   printandlog, 'each transit file. This may take a while, but can be done at any',logname
   printandlog, 'point with the following command:', logname
   printandlog, "junk = exofast_gettt(filename='" + idlfile + "', filebase='" + prefix + "')", logname
   junk = exofast_gettt(mcmcss, filebase=prefix)

   if total(mcmcss.ttvs) gt 1 then begin
      ;; generate an O-C diagram for each planet
      readcol, prefix + 'transits.csv',label,planet,epoch,time,hierr,loerr, format='a,a,l,d,d,d', delimiter=',', comment='#',/silent
      err = (double(hierr) + double(loerr))/2d0
      telescope = label
      for i=0L, n_elements(label)-1 do telescope[i] = (strsplit(label[i],' UT ',/regex,/extract))[0]

      sorted = sort(planet)
      uniqplanets = planet[sorted[uniq(planet[sorted])]]
      for i=0L, n_elements(uniqplanets)-1 do begin
         match = where(planet eq uniqplanets[i]); and mcmcss.ttvs[*,i])
         omc, time[match], err[match], telescope=telescope[match], epsname=prefix + uniqplanets[i] + '.ttv.eps', $
              period=median(mcmcss.planet[i].period.value), t0=median(mcmcss.planet[i].tc.value),logname=logname
      endfor
   endif
   
   ;; generate a plot of the tdeltavs 
   ;; *** I think this breaks for multiple planets***
   if total(mcmcss.tdeltavs) gt 1 then begin
      plot_tdeltav, prefix + 'transits.csv'
   endif

endif

;; this makes a quick-look summary page, but requires gs, ps2pdf, grep,
;; convert and likely only works on linux
if keyword_set(mksummarypg) then begin
   mksummaryframe,idlfile=idlfile,base=prefix
endif

end
