;+
; NAME:
;   EXOFAST_MOD
;
; PURPOSE: 
;   Return postive-only mod.
; 
;   IDL's built-in mod:
;     result = a mod b 
;   returns -b < result < +b
;
;   EXOFAST_MOD:
;     result = exofast_mod(a,b)
;   returns 0 <= result < b 
;
; CALLING SEQUENCE:
;   result = exofast_mod(a,b)
; 
; OPTIONAL KEYWORDS:
;   NEGATIVE: if set, returns:
;     result = exofast_mod(a,b,/negative)
;   returns -b/2 <= result < b/2
; INPUTS:
;
;  A - Number
;  B - Divisor
;
; EXAMPLES:
;
;  ;; compute the phase of a planet
;  phase = exofast_mod(time,period)/period
;
;  ;; return times folded and centered around the transit
;  phasetime = exofast_mod(time-tc,period,/negative)
; MODIFICATION HISTORY
; 
;  2023/02/09 -- JDE, Written
;-
function exofast_mod, a, b, negative=negative

  if keyword_set(negative) then begin
     ;; -b/2 <= result < +b/2
     a0 = a+b/2d0
     return, a0-floor(a0/b)*b - b/2d0
  endif else begin
     ;; 0 <= result < b     
     return, a-floor(a/b)*b  
  endelse

  ;; significantly faster than double mod version (for n_elements(a) >~ 20)
;  return, ((a mod b) + b) mod b


ntrials = 100d0

alltimemod1 = 0d0
alltimemod2 = 0d0

for i=0L, ntrials-1 do begin

   t0 = systime(/seconds)
   mod1 = ((a mod b) + b) mod b
   t1 = systime(/seconds)
   mod2 = a-floor(a/b)*b
   t2 = systime(/seconds)
   mod3 = a mod b
   t3 = systime(/seconds)
;   mod4 = ((a mod b) + b) mod b
;   t4 = systime(/seconds)

   timemod1 = t1-t0
   timemod2 = t2-t1
   timemod3 = t3-t2
;   timemod4 = t4-t3
   
   alltimemod1 += timemod1
   alltimemod2 += timemod2

;   print, max(abs(mod1-mod2))
;   print, 'mod2 took ' + strtrim(timemod2/timemod3,2) + ' times longer than built-in mod'
;   print, 'mod1 took ' + strtrim(timemod4/timemod3,2) + ' times longer than built-in mod'
   
endfor

print, 'double mod took ' + strtrim(alltimemod1/alltimemod2,2) + ' times longer than floor mod'

return, a-floor(a/b)*b ;; faster

end
