function istar_chi2, pars, determinant=determinant, derived=derived

common chi2_block, priors
vsini_prior = priors[0,0]
sigma_vsini = priors[1,0]
vsini_min = priors[2,0]
vsini_max = priors[3,0]

vrot_prior = priors[0,1]
sigma_vrot = priors[1,1]
vrot_min = priors[2,1]
vrot_max = priors[3,1]


determinant=1

vrot = pars[0] 
cosi = pars[1]

istar = acos(cosi)
vsini = vrot*sin(istar)

if vrot lt 0d0 then return, !values.d_infinity
if cosi lt 0d0 then return, !values.d_infinity
if cosi gt 1d0 then return, !values.d_infinity

if vrot lt vrot_min then return, !values.d_infinity
if vrot gt vrot_max then return, !values.d_infinity
if vsini lt vsini_min then return, !values.d_infinity
if vsini gt vsini_max then return, !values.d_infinity


chi2 = 0d0
if sigma_vsini gt 0d0 then chi2 += ((vsini - vsini_prior)/sigma_vsini)^2
if sigma_vrot gt 0d0 then chi2 += ((vrot - vrot_prior)/sigma_vrot)^2

;print, ((vsini - vsini_prior)/sigma_vsini)^2, ((vrot - vrot_prior)/sigma_vrot)^2, pars[0], pars[1], chi2


return, chi2

end
