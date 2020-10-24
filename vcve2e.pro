function vcve2e, vcve0, omega=omega0, lsinw=lsinw, lcosw=lcosw, sign=sign

;; don't pass back the modified values
vcve = vcve0

;; derive omega
if n_elements(omega0) eq 0 then begin
   if n_elements(lsinw) eq 0 or n_elements(lcosw) eq 0 then begin
      message, 'must specify either omega or lsinw and lcosw'
   endif
   omega = atan(lsinw,lcosw)
endif else omega = omega0

;; allow multiple omegas for 1 vcve
if n_elements(vcve0) eq 1 then begin
   if n_elements(omega) ne 1 then begin
      vcve = replicate(vcve,n_elements(omega))
   endif
endif

;; allow multiple vcves for 1 omega
if n_elements(omega) eq 1 then begin
   if n_elements(vcve) ne 1 then begin
      omega = replicate(omega,n_elements(vcve))
   endif
endif

;; allow multiple vcve/omega with 1 sign
if n_elements(sign) eq 1 then begin
   if n_elements(vcve) ne 1 then begin
      sign = replicate(sign,n_elements(vcve))
   endif
endif

;; make sure the inputs are reasonable
npoints = n_elements(vcve)
if n_elements(omega) ne npoints then $
   message, 'vcve and omega must have the same number of elements'

;; solve the quadratic for e
a = vcve^2*sin(omega)^2 + 1d0
b = 2d0*vcve^2*sin(omega)
c = vcve^2-1d0
epos = (-b + sqrt(b^2 - 4d0*a*c))/(2d0*a) 
eneg = (-b - sqrt(b^2 - 4d0*a*c))/(2d0*a) 

;; pick the desired solution
e = dblarr(npoints)
if n_elements(sign) ne npoints then begin
   if n_elements(lsinw) eq 0 or n_elements(lcosw) eq 0 then begin
      ;; no sign specified, either solution is allowed
      ;; pick the one that is physical (if there is one)
      ;; if they're both good, default to the lower eccentricity
      useneg = where(finite(eneg) and $
                     eneg ge 0d0 and $
                     eneg lt 1d0 and $
                     eneg lt epos, complement=usepos)
      
      ;; return the sign that will return this value
      sign = bytarr(npoints)
      if useneg[0] ne -1 then sign[useneg] = 1B
      if npoints eq 1 then sign = sign[0]
   endif else begin
      ;; we use L to choose the solution
      ;; This choice reduces gaps in parameter space
      useneg = where((lsinw^2 + lcosw^2) ge 0.5d0, complement=usepos) 
   endelse
endif else begin
   ;; we use the sign to choose the solution
   useneg = where(floor(sign),complement=usepos)
endelse

if usepos[0] ne -1 then e[usepos] = epos[usepos]
if useneg[0] ne -1 then e[useneg] = eneg[useneg]

;; debugging statement
;print, 'vcve = ' + strtrim(vcve,2) + ', omega = ' + strtrim(omega,2) +  ',epos = ' + strtrim(epos,2) + ', eneg = ' + strtrim(eneg,2) + ', sign = ' + strtrim(long(sign),2)

;; if scalar input, scalar output
if npoints eq 1 then return, e[0]
return, e

end
