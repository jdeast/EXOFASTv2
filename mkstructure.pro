;; This creates a structure that describes observations of a stellar
;; system
;; It describes an arbitrary number of planets, observed bands, rv
;; telescopes (a new zero point for each), and observed transits. 
;; This structure is the input to everything in EXOFASTv2 and allows
;; for much easier manipulation/generalization of fitting planetary
;; systems.
 

pro mkstructure, nplanets=nplanets, bands=bands, circular=circular, chen=chen,$
                    fitslope=fitslope, fitquad=fitquad, ttv=ttv, tdv=tdv, $
                    rm=rm, eprior4=eprior4, fittran=fittran, fitrv=fitrv

nplanets = 3
ntel = 3
ntran = 10
nbands = 4

;; each parameter
parameter = create_struct('value',0d0,$
                          'prior',0d0,$
                          'priorwidth',!values.d_infinity,$
                          'unit','',$
                          'description','',$
                          'latex','',$
                          'label','',$
                          'cgs',1d0,$         ;; conversion to cgs
                          'fit',0B,$          ;; is it fit?
                          'link',ptr_new(1),$ ;; a pointer to a linked variable
                          'derive',1B)        ;; is it derived?

mstar = parameter
mstar.unit = 'msun'
mstar.description = 'Mass'
mstar.latex = 'M_{*}'
mstar.label = 'mstar'
mstar.cgs = 1.9891d33
mstar.fit = 1
mstar.derive=0

logg = parameter
logg.unit = 'cgs'
logg.description = 'Surface gravity'
logg.latex = '\log{g}'
logg.label = 'logg'
logg.fit = 1
logg.derive=0

teff = parameter
teff.unit = 'K'
teff.description = 'Effective Temperature'
teff.latex = 'T_{eff}'
teff.label = 'teff'
teff.fit=1
teff.derive=0

feh = parameter
feh.description = 'Metalicity'
feh.latex = '[Fe/H]'
feh.label = 'feh'
feh.fit=1
feh.derive=0

vsini = parameter
vsini.unit = 'm/s'
vsini.description = 'Projection rotational velocity'
vsini.latex = 'v\sin{i}'
vsini.label = 'vsini'
vsini.cgs = 100d0
if keyword_set(rm) then vsini.fit = 1
vsini.derive = 0

distance = parameter
distance.unit = 'pc'
distance.description = 'Distance'
distance.latex = 'd'
distance.label = 'd'
distance.cgs = 3.08567758d18
distance.derive = 0 ;; not supported

slope = parameter
slope.unit = 'm/s/day'
slope.description = 'RV slope'
slope.latex = '\dot{\gamma}'
slope.label = 'slope'
slope.cgs = 100d0/86400d0
if keyword_set(fitslope) then slope.fit = 1
slope.derive=0

quad = parameter
quad.unit = 'm/s/day^2'
quad.description = 'RV quad'
quad.latex = '\ddot{\gamma}'
quad.label = 'quad'
quad.cgs = 100d0/86400d0^2
if keyword_set(fitquad) then quad.fit=1
quad.derive=0

rstar = parameter
rstar.unit = '\rsun'
rstar.description = 'Radius'
rstar.latex = 'R_{*}'
rstar.label = 'radius'
rstar.cgs = 6.955d10

age = parameter
age.unit = 'Gyr'
age.description = 'Age'
age.latex = 'Age'
age.label = 'age'
age.cgs = 3600d0*24d0*365.242d0*1d9

lstar = parameter
lstar.unit = '\lsun'
lstar.description = 'Lstar'
lstar.latex = 'Lstar'
lstar.label = 'lstar'
lstar.cgs = 3600d0*24d0*365.242d0*1d9

rhostar = parameter
rhostar.unit = 'cgs'
rhostar.description = 'Density'
rhostar.latex = '\rho_*'
rhostar.label = 'rhostar'

tc = parameter
tc.unit = '\bjdtdb'
tc.description = 'Time of Transit'
tc.latex = 'T_C'
tc.label = 'tc'
tc.cgs = 86400d0
tc.fit = 1
tc.derive=0

logp = parameter
logp.unit = ''
logp.description = 'Log of Period'
logp.latex = 'log P'
logp.label = 'logp'
logp.cgs = !values.d_nan
logp.fit = 1
logp.derive=0

qesinw = parameter
qesinw.unit = ''
qesinw.description = 'qesinw'
qesinw.latex = 'e^{1/4} sin{\omega_*}'
qesinw.label = 'qesinw'
qesinw.cgs = !values.d_nan

qecosw = parameter
qecosw.unit = ''
qecosw.description = 'qecosw'
qecosw.latex = 'e^{1/4} cos{\omega_*}'
qecosw.label = 'qecosw'
qecosw.cgs = !values.d_nan

sesinw = parameter
sesinw.unit = ''
sesinw.description = 'sesinw'
sesinw.latex = '\sqrt{e} sin{\omega_*}'
sesinw.label = 'sesinw'
sesinw.cgs = !values.d_nan

secosw = parameter
secosw.unit = ''
secosw.description = 'secosw'
secosw.latex = '\sqrt{e} cos{\omega_*}'
secosw.label = 'secosw'
secosw.cgs = !values.d_nan

if keyword_set(eprior4) then begin
   ;; step in qesinw (prior favoring lower e)
   qesinw.fit = 1
   qecosw.fit = 1
   qesinw.derive = 0
   qecosw.derive = 0 
endif else begin
   ;; step in sqrt(e)sinw (uniform prior)
   sesinw.fit = 1
   secosw.fit = 1
   sesinw.derive = 0
   secosw.derive = 0
endelse

esinw = parameter
esinw.unit = ''
esinw.description = 'esinw'
esinw.latex = 'esin{\omega_*}'
esinw.label = 'esinw'
esinw.cgs = !values.d_nan

ecosw = parameter
ecosw.unit = ''
ecosw.description = 'ecosw'
ecosw.latex = 'ecos{\omega_*}'
ecosw.label = 'ecosw'
ecosw.cgs = !values.d_nan

e = parameter
e.unit = ''
e.description = 'Eccentricity'
e.latex = 'e'
e.label = 'e'
e.cgs = !values.d_nan

omega = parameter
omega.unit = 'Radians'
omega.description = 'Argument of Periastron'
omega.latex = '\omega_*'
omega.label = 'omega'
omega.cgs = 1d0

omegadeg = parameter
omegadeg.unit = 'degrees'
omegadeg.description = 'Argument of Periastron'
omegadeg.latex = '\omega_*'
omegadeg.label = 'omegadeg'
omegadeg.cgs = !dpi/180d0

p = parameter
p.unit = ''
p.description = 'Radius of planet in stellar radii'
p.latex = 'R_P/R_*'
p.label = 'p'
p.fit = 1

ar = parameter
ar.unit = ''
ar.description = 'Radius of planet in stellar radii'
ar.latex = 'a/R_*'
ar.label = 'ar'

cosi = parameter
cosi.unit = ''
cosi.description = 'Cos of inclination'
cosi.latex = 'cos i'
cosi.label = 'cosi'
cosi.fit = 1

inc = parameter
inc.unit = 'Radians'
inc.description = 'Inclination'
inc.latex = 'i'
inc.label = 'inc'

incdeg = parameter
incdeg.unit = 'degrees'
incdeg.description = 'Inclination'
incdeg.latex = 'i'
incdeg.label = 'incdeg'
incdeg.cgs = !dpi/180d0


;; for each star
star = create_struct('mass',mstar,$
                     'logg',logg,$
                     'teff',teff,$
                     'feh',feh,$
                     'vsini',vsini,$
                     'distance',distance,$
                     'slope',slope,$
                     'quad',quad,$
                     'ellip',0d0,$
                     'ra',0d0,$
                     'dec',0d0,$
                     'omega',0d0,$
                     'pm',[0d0,0d0],$
                     'radius',rstar,$
                     'age',age,$
                     'luminosity',lstar,$
                     'density',rhostar)
            
;; for each planet
planet = create_struct('tc',tc,$
                       'logP',logp,$
                       'logK',0d0,$
                       'cosi',cosi,$
                       'p',p,$
                       'f0',0d0,$
                       'lambda',0d0,$
                       'beam',0d0,$
                       'e',e,$
                       'omega',omega,$
                       'omegadeg',omegadeg,$
                       'ecosw',ecosw,$
                       'esinw',esinw,$
                       'secosw',secosw,$
                       'sesinw',sesinw,$
                       'qecosw',qecosw,$
                       'qesinw',qesinw,$
                       'tp',0d0,$
                       'ts',0d0,$
                       'ta',0d0,$
                       'td',0d0,$
                       'bs',0d0,$
                       'tfwhms',0d0,$
                       'taus',0d0,$
                       't14s',0d0,$
                       'ps',0d0,$ ;; A priori non-grazing eclipse prob
                       'psg',0d0,$ ;; A priori eclipse prob
                       'K',0d0,$
                       'period',0d0,$
                       'ar',0d0,$
                       'a',0d0,$
                       'mass',0d0,$
                       'radius',0d0,$
                       'density',0d0,$
                       'logg',0d0,$
                       'teq',0d0,$
                       'safronov',0d0,$
                       'fave',0d0,$
                       'msini',0d0,$
                       'q',0d0,$
                       'i',0d0,$
                       'b',0d0,$
                       'delta',0d0,$
                       'tfwhm',0d0,$
                       'tau',0d0,$
                       't14',0d0,$
                       'pt',0d0,$
                       'ptg',0d0,$
                       'depth',0d0,$
                       'dr',0d0)
                       
;; for each wavelength
band = create_struct('u1',0d0,$ ;; linear limb darkening
                     'u2',0d0,$ ;; quadratic limb darkening
                     'u3',0d0,$ ;; 1st non-linear limb darkening
                     'u4',0d0,$ ;; 2nd non-linear limb darkening
                     'thermal',0d0,$   ;; thermal emission
                     'dilution',0d0,$  ;; dilution
                     'reflect',0d0,$   ;; reflection
                     'mag',0d0,$
                     'name','')

;; for each telescope
telescope = create_struct('gamma',0d0,$
                          'rvptrs',ptr_new(1),$
                          'name','')

;; for each transit
transit = create_struct('ttv',0d0,$ ;; Transit Timing Variation
                        'tbv',0d0,$ ;; Transit b variation
                        'tdv',0d0,$ ;; Transit depth variation
                        'transitptrs',ptr_new(1),$ ;; Data
                        'band','',$
                        'exptime',0d0,$
                        'ninterp',0d0,$
                        'name','')
                        
;; for each system
ss = create_struct('star',star,$
                   'planet',replicate(planet,nplanets),$
                   'band',replicate(band,nbands),$
                   'telescope',replicate(telescope,ntel),$
                   'transit',replicate(transit,ntran),$
                   'tofit',ptr_new(1),$
                   'priors',ptr_new(1))

for i=0, nplanets-1 do begin
   ;; circular orbit, don't fit e or omega
   if circular[i] then begin
      ss.planet[i].qesinw.derive = 0
      ss.planet[i].qecosw.derive = 0
      ss.planet[i].sesinw.derive = 0
      ss.planet[i].secosw.derive = 0
      ss.planet[i].qesinw.fit = 0
      ss.planet[i].qecosw.fit = 0
      ss.planet[i].sesinw.fit = 0
      ss.planet[i].secosw.fit = 0
   endif 
endfor

;; this assumes a certain structure of the parameters... is that ok?
tofit = [-1,-1,-1]
priors = [-1,-1,-1]
for i=0, n_tags(ss)-1 do begin
   for j=0, n_elements(ss.(i))-1 do begin
      for k=0, n_tags(ss.(i)[j])-1 do begin
         print, n_tags(ss.(i)[j].(k))
         if n_tags(ss.(i)[j].(k)) ne 0 then begin
            if tag_exist(ss.(i)[j].(k),'fit') then begin
               if ss.(i)[j].(k).fit then tofit = [[tofit],[i,j,k]]
               if finite(ss.(i)[j].(k).priorwidth) then priors=[[priors],[i,j,k]]
            endif
         endif
      endfor
   endfor
endfor
tofit = tofit[*,1:*]
priors = priors[*,1:*]

*(ss.tofit) = tofit
*(ss.priors) = priors

;; populate the structure with the parameter array
for i=0, n_elements(pars)-1 do begin
   tofit = *(ss.tofit)
   ss.(tofit[i,0])[tofit[i,1]].(tofit[i,2]).value = pars[i]
endfor

;; this function populates the priors in the structure
function mkprior, str, priorname, priorval, priorwidth

  ;; for each input prior
  for j=0, n_elements(priorname)-1 do begin
     
     ;; determine the subscript
     tmp = strsplit(priorname[i],'_',/extract)
     if n_elements(tmp) eq 2 then priornum = long(tmp[1]) $
     else priornum = 0
     priorlabel = tmp[0]
     
     ;; look for it in the structure
     for i=0, n_tags(str)-1 do begin
        found = 0
        if (n_elements(str.(i))-1) ge priornum then begin
           for k=0, n_tags(str.(i)[priornum])-1 do begin
              if tag_exist(str.(i)[priornum],priorlabel,index=ndx) then begin
                 ;; found it!
                 str.(i)[priornum].(ndx).prior = priorval[j]
                 str.(i)[priornum].(ndx).priorwidth = priorwidth[j]
                 break
              endif
           endfor
        endif
        ;; didn't find it, warn user
        if not found then message, 'WARNING: No parameter name ' + $
                                   priorname[j] + '; not applying prior'
     endfor

  endfor ;; each input prior

end

stop

end
