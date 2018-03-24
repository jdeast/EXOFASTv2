;+
; PURPOSE:
;  This function calculates the Planck function I_nu or
;  I_lambda.
;    I_nu     = 2h nu^3 c^-2 ( exp(h nu / kT) -1)^-1
;    I_lambda = 2h c^2 lambda^-5 (exp ( hc / lambda kt) -1)^-1
;
; CATEGORY:
;  Astrophysics
;
; CALLING SEQUENCE:
;  result = BLACKBODY(temp, x, [/freq, /wave, /cgs, /mks])
; 
; INPUTS:
;  temp: The temperature, in K. Scalar or vector. If it is a scalar,
;        the same temperature is assumed for each value of x.
;     x: By default, the wagelength in mks units. If /cgs is set, then
;     this number is in cgs units. If /freq is set, this number is the
;     frequency. 
;
; KEYWORD PARAMETERS:
;  wave: Set to receive the input as a wavelength, and output I_lambda
;  freq: Set to receive the input as a frequency, and output I_nu
;        (default) 
;  cgs:  Set to receive the input and return the output in cgs units
;  mks:  Set to receive the input and return the output in mks units
;        (default) 
;
; OUTPUTS:
;  I_nu or, if /wave is set, I_lambda
;
; MODIFICATION HISTORY:
;  May 2009: Written by Chris Beaumont
;-
function blackbody, temp, x, freq = freq, wave = wave, cgs = cgs, mks = mks
  compile_opt idl2
  on_error, 2

  ;- check inputs
  if n_params() ne 2 then begin
     print, 'blackbody calling sequence:'
     print, 'result = blackbody(temp, x, [/freq, /wave, /cgs, /mks])'
     print, '         default: /freq, /mks'
     return, !values.f_nan
  endif

  ntemp = n_elements(temp)
  nx = n_elements(x)
  
  if ntemp ne nx && ntemp ne 1 then $
     message, 'temp must be a scalar, or temp and x must have the same number of elements'

  doCGS = keyword_set(cgs)
  if doCGS && keyword_set(mks) then $
     message, 'do not set both /cgs and /mks!'

  doWAVE = keyword_set(wave)
  if doWAVE && keyword_set(freq) then $
     message, 'do not set both /freq and /wave!'

  ;- constants
  if (doCGS) then begin
     h = 6.626d-27
     k = 1.38065d-16
     c = 2.99792d10
  endif else begin
     h = 6.626d-34
     k = 1.38065d-23
     c = 2.99792d8
  endelse

  nu = doWave ? c / x : x
  
  ;- I_nu
  result = 2 * h * nu^3 * c^(-2) * 1 / (exp(h * nu / (k * temp)) - 1)
  
  ;- convert to I_lambda if desired
  if (doWAVE) then result = result * nu^2 / c

  return, result
end
  
