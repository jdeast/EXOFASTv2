function multifast_dopptom, pars

COMMON dopptom_block, useDoppTom, dopptomptrs

tc        = pars[1]       ;; transit center time
Period    = 10^(pars[2])  ;; Period of orbit
sqrtecosw = pars[3]       ;; ecosw
sqrtesinw = pars[4]       ;; esinw
cosi      = pars[6]       ;; cosi
p         = pars[7]       ;; rp/rstar
ar        = 10^(pars[8])
logg      = pars[9]       ;; density of the Primary (log(rho_star)/(rho_sun))
teff      = pars[10]      ;; effective temperature of the Primary
feh       = pars[11]      ;; metallicity of the Primary
vsini     = pars[13]/1000.      ;; vsini of the star
lambda    = pars[14]      ;; sky-projected spin-orbit angle
macturb   = pars[15]/1000.

if macturb lt 0d0 then return, !values.d_infinity ; invalid line width
if vsini lt 0d0 then return, !values.d_infinity ; invalid vsini
chisqr = 0.0d0;
convolve_limit = 5. ;; The gaussian broadening has to be less than this factor larger than vsini for it to try and include the rotation kernel. E.g., gaussian broadening ('GaussTerm') = 10km/s, vsini needs to be at least 2km/s. Otherwise just use a gaussian for the shadow.

inc = acos(cosi)
b = ar * cosi
e = sqrtecosw^2 + sqrtesinw^2
if e eq 0d0 then omega = !dpi/2d0 else omega = atan(sqrtesinw, sqrtecosw)

for d = 0, n_elements(dopptomptrs) - 1 do begin
    DTccf2d = (*(dopptomptrs[d])).ccf2d
    DTbjd = (*(dopptomptrs[d])).bjd
    DTvel = (*(dopptomptrs[d])).vel
    DTstepsize = (*(dopptomptrs[d])).stepsize 
    DTrms = (*(dopptomptrs[d])).rms
    DTchisqr0 = (*(dopptomptrs[d])).chisqr0
    DTRspec = (*(dopptomptrs[d])).Rspec
    
    velsini = DTvel/ss.star.vsini.value
    stepsize = DTstepsize/ss.star.vsini.value
    meanstepsize = mean(stepsize)
    
    IndepVels = (3d5/DTRspec) / (meanstepsize*ss.star.vsini.value) ; accounts for the fact that we're supersampling within the actual spectral res. of the spectrograph, to adjust chi^2 later
    
    GaussTerm = sqrt(ss.star.macturb.value^2. + (3d5/DTRspec)^2.) ; so we can add in both the broadening from the instrument's spectral resolution, and the macturb parameter
    
    GaussRel = GaussTerm/ss.star.vsini.value
    
    nper = round((mean(DTbjd)-ss.planet[i].tc.value) / ss.planet[i]/period.value)
    localtc = ss.planet[i].tc.value + ss.planet[i].Period.value*nper
    phase = (DTbjd-localtc) / ss.planet[i].Period.value
    
    ;; init for planet shadow model
    velwidth = p
    c1 = 2d0 / (!dpi*velwidth)
    
    xp = ss.planet[i].ar.value * sin(phase*2d0*!dpi)
    zp = ar * cos(phase*2d0*!dpi) * cosi
    up = xp*cos(ss.planet[i].lambda.value) + zp*sin(ss.planet[i].lambda.value)
    
    ;; Calculate Collier-Cameron's "beta" as the actual transit LC in V
    ;; This presumes we're using an ~optical spectrograph for the observations
    beta = phase*0d0 + 1d0
    
    ldcoeffs = quadld(ss.star.logg.value, ss.star.teff.vale, ss.star.feh.value, 'V')
    u1 = ldcoeffs[0]
    u2 = ldcoeffs[1]
    if ((not finite(u1)) or (not finite(u2))) then  return, !values.d_infinity
    
    z = exofast_getb(DTbjd, i=inc, a=ar, tperiastron=localtc, period=Period, e=e,omega=omega,z=depth)
    
    primary = where(depth gt 0,complement=secondary)     
    if primary[0] ne -1 then begin
       exofast_occultquad, z[primary], u1, u2, p, mu1
       beta[primary] = mu1
    endif
    
    beta = 1d0-beta ; to just get the amount of light blocked, rather than the ligthcurve
    
    DTmodel = DTccf2d * 0d0 + median(DTccf2d)


    relevantVels = where((velsini gt -1.5d0) and (velsini lt 1.5d0)) ; we only care about these, to speed things up

    good = where(beta gt 0d0 and up lt 1d0 and (velsini gt -1.5d0) and (velsini lt 1.5d0))
    if GaussRel lt convolve_limit then begin ; rotation is significant
       G = velsini[good]*0.0d0
       c2 = ((velsini[good]-up[good])/velwidth)^2.
       validlocs = where(c2 lt 1.0d0)
       if validlocs[0] ne -1 then G[validlocs] = c1*sqrt(1-c2[validlocs])
       normG = 1d0 / total(G*stepsize)
       if normG ne !values.d_infinity then begin ; so we don't waste time on zeros
          rotprofile = G * normG
          unnormalized = multifast_gaus_convol(velsini[good], rotprofile, GaussRel)
          normalization = 1 / total(unnormalized*stepsize)
          DTmodel[good,i] += (beta[i] * normalization * unnormalized)
       endif
    endif else begin ; gaussian line width dominates
       expterm = (up[good]-velsini[relevantVels])^2/(2d0*GaussRel^2)
       unnormalized = exp(expterm)
       normalization = 1d0 / total(unnormalized*stepsize)
       DTmodel[relevantVels,i] += (beta[i] * normalization * unnormalized)
    endelse

    nphase = n_elements(phase)
    for i = 0, nphase-1 do begin
    	if beta[i] gt 0.0d0 then begin ; only do this if the planet is on the star
    		if up[i] lt 1.0d0 then begin ; also only if the DT thinks it's on the star
    			if GaussRel lt convolve_limit then begin ; rotation is significant
    				G = velsini[relevantVels]*0.0d0
    				c2 = ((velsini[relevantVels]-up[i])/velwidth)^2.
    				validlocs = where(c2 lt 1.0d0)
    				;; Need to check if we have any valid locations
    				if validlocs[0] ne -1 then G[validlocs] = c1*sqrt(1-c2[validlocs])
    				normG = 1 / total(G*stepsize)
    				if normG ne !values.d_infinity then begin ; so we don't waste time on zeros
    					rotprofile = G * normG
    					unnormalized = multifast_gaus_convol(velsini[relevantVels], rotprofile, GaussRel)
    					normalization = 1 / total(unnormalized*stepsize)
    					DTmodel[relevantVels,i] += (beta[i] * normalization * unnormalized)
    				endif
    			endif else begin ; gaussian line width dominates
    				expterm = -1.*(velsini[relevantVels]-up[i])^2./(2*GaussRel^2.)
    				unnormalized = exp(expterm)
    				normalization = 1 / total(unnormalized*stepsize)
    				DTmodel[relevantVels,i] += (beta[i] * normalization * unnormalized)
    			endelse
    		endif
    	endif
    endfor
    
    DTresid = DTccf2d - DTmodel
    ;chisqr0 = DTchisqr0 / IndepVels
    chisqr += total(sqrt(DTresid^2./DTrms^2.))/IndepVels; - chisqr0
    (*(dopptomptrs[d])).model = DTmodel
    (*(dopptomptrs[d])).resid = DTresid
endfor

return, chisqr

end
