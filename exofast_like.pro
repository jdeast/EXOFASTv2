;+
; NAME:
;       EXOFAST_LIKE
; PURPOSE:
;       Calculates the log likelihood for the given residuals and
;       errors. A simpler (and faster) alternative to the wavelet
;       analysis of Carter & Winn 2009.
;
; CALLING SEQUENCE:
;       result = exofast_like(residuals, sigma_r, sigma_w [,/chi2])
;
; INPUTS:
;   RESIDUALS - The residuals of a fit (data-model)
;   SIGMA_R   - Red noise amplitude parameter
;   SIGMA_W   - White noise amplitude parameter
;
; OPTIONAL KEYWORDS:
;   CHI2    - If set, will return the "equivalent chi2". The result
;             (-2d0*loglike) is not actually the chi^2, but (a hack),
;             designed to be the required output of CHI2FUNC called by
;             EXOFAST_DEMC such that it uses the likelihood rather
;             than the chi^2. It is equivalent because the difference
;             between this and the true chi^2 is a constant that
;             cancels when comparing two models (and is ultimately a
;             deficiency in using the chi^2 rather than the
;             likelihood). This allows us to fit a red noise term
;             (sigma_r) while maintaining backward compatibility with
;             EXOFAST_DEMC.
;  TRUECHI2 - If set, will return the actual chi2.
;
; OUTPUTS:
;    RESULT - alog(likelihood), as defined in Eastman et al., 2013, eqn
;             1.
;
; EXAMPLES:
;
;    ;; return the effective chi2 (for use with EXOFAST_DEMC) with an 
;    ;; error scaling term but no red noise
;    chi2 = exofast_like(residuals, 0d0, err*scale,/chi2)
;
; REVISION HISTORY:
;       2013/11 - Written, Jason Eastman, LCOGT
;-   

function exofast_like, residuals, var_r, sigma_w0, chi2=chi2, truechi2=truechi2

  ;; if sigma_w0 is only one elements, replicate it
  if n_elements(sigma_w0) eq 1 then $
     sigma_w = replicate(sigma_w0,n_elements(residuals)) $
  else sigma_w = sigma_w0

  good = where(finite(sigma_w))
  if good[0] eq -1 then begin
     if keyword_set(chi2) then return, !values.d_infinity
     return, -!values.d_infinity
  endif

  chisq = total(residuals[good]^2/(sigma_w[good]^2 + var_r))
  if keyword_set(truechi2) then return, chisq

  ;; Calculate the natural log of the likelihood
  loglike = -0.5d0*(total(alog(sqrt(2d0*!dpi*(sigma_w[good]^2+var_r)))) + chisq)

  if ~finite(loglike) then begin
     if keyword_set(chi2) then return, !values.d_infinity
     return, -!values.d_infinity
  endif    

  ;; This is not actually the chi2; see above for explanation
  if keyword_set(chi2) then return, -2d0*loglike
     
  return, loglike

end
