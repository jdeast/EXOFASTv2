function exofast_wavelike,x,sigma_r,sigma_w,zeropad=zeropad

;+
; NAME:
;       WAVELETLIKE
; PURPOSE:
;       Calculates the (log) likelihood for a given vector to be
;       described for parameters sigma_r and sigma_w as 1/f noise plus
;       white noise as detailed by Carter & Winn (2009)
;
; CALLING SEQUENCE:
;       result = waveletlike(x, sigma_r, sigma_w)
;
; INPUTS:
;       X - Data vector, must be a power of two
;       SIGMA_R - Red noise amplitude parameter (not equal to
;                 1/f-component RMS)
;       SIGMA_W - White noise amplitude parameter (approximately equal
;                 to white-component RMS)
;
; OUTPUTS:
;       RESULT - log(likelihood) as defined in Eqn. (32) of Carter &
;                Winn (2009)
;
; NOTES:
;       Refer to the paper by Carter & Winn (2009) for theory and
;       details.  In the notation of that paper, gamma=1 for this
;       algorithm.
;
; REVISION HISTORY:
;       Written,    J. Carter               September, 2009
;       Removed inner for loop (~7x speed gain). JDE, Dec 2012
;-   

gamma=1d0

d=x
els = n_elements(d)
pow = ceil(alog(els)/alog(2.))
if (2^pow ne els and keyword_set(zeropad)) then begin
   diff = 2.^pow-els
   left = floor(diff/2.)
   right = diff-left
   if diff gt 1 then x = [replicate(0d,left),d,replicate(0d,right)] $
   else x = [d,replicate(0d,right)] 
endif else x = d

J = alog(n_elements(x))/alog(2)
if (abs(J-fix(J)) ne 0) then message,'Data length must be a power of two'
J = fix(J)

info = WV_FN_DAUBECHIES(2,wavelet,scaling,ioff,joff)
wv = wv_dwt(x,wavelet,scaling,ioff,joff)

sm2 = sigma_r^2*(gamma eq 1 ? 1d0/(2d0*alog(2d0)) : 2-2^gamma)+sigma_w^2
sum = wv[0]^2/sm2+alog(2d0*!DPI*sm2)

for i=1,J do begin
   sm2 = sigma_r^2*2^(-gamma*i)+sigma_w^2
   sum += total(wv[2^(i-1):2^i-1]^2)/sm2 + 2^(i-1)*alog(2d0*!DPI*sm2)
endfor

x=d

return,-sum/2d0

end
