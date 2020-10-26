;; examples from https://arxiv.org/pdf/2001.04973.pdf
;istar, 9.8d0, 2d0, 10d0, 2d0
;istar, 4d0, 2d0, 10d0, 2d0
;istar, 1d0, -1d0, 10d0, 2d0, vsini_max=5d0
;istar, 1d0, -1d0, 10d0, 2d0, vsini_max=30d0


pro istar, vsini, sigma_vsini, vrot, sigma_vrot, maxsteps=maxsteps, nthin=nthin, $
           vsini_min=vsini_min, vsini_max=vsini_max,$
           vrot_min=vrot_min, vrot_max=vrot_max, dontstop=dontstop

if n_elements(maxsteps) eq 0 then maxsteps = 1d6
if n_elements(vsini_min) eq 0 then vsini_min = 0d0
if n_elements(vsini_max) eq 0 then vsini_max = !values.d_infinity
if n_elements(vrot_min) eq 0 then vrot_min = 0d0
if n_elements(vrot_max) eq 0 then vrot_max = !values.d_infinity

common chi2_block, priors


priors = [[vsini, sigma_vsini, vsini_min, vsini_max],$
          [vrot, sigma_vrot, vrot_min, vrot_max]]

cosi = cos(asin(vsini/vrot))
best = [vrot, cosi]

chi2func = 'istar_chi2'

exofast_demcpt_multi, best, chi2func, pars, chi2=chi2, nthin=nthin, maxsteps=maxsteps, dontstop=dontstop

burnndx = getburnndx(chi2,goodchains=goodchains)

vrot = pars[0,burnndx:*,goodchains]
cosi = pars[1,burnndx:*,goodchains]

hist = histogram(cosi, locations=x, nbin=100)
plot, x, hist/total(hist), psym=10, xtitle='cosi', ytitle='Normalized Probability'



end
