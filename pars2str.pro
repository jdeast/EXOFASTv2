;+
; NAME:
;   PARS2STR
;
; PURPOSE: 
;   Translates a parameter array (for use with EXOFAST_DEMC) into the
;   stellar system structure (for use with EXOFASTv2)
;
; CALLING SEQUENCE:
;    pars2str(pars,str)
;-

pro pars2str, pars, str, best=best

  tofit = *(str.tofit)
  npars = n_elements(tofit[0,*])

  ;; just a sanity check to make sure the mapping is correct
  errmsg = "ERROR: There is a serious bug in the mapping between " + $
           "the parameter steps and the data structure! Please notify " + $
           "jason.eastman@cfa.harvard.edu and include everything required " + $
           "to run your fit. Do not trust any results from a run that shows this message."

  if (size(pars))[0] eq 1 then begin
     ;; populate the structure with the parameter array
     for i=0, npars-1 do begin
        if tofit[3,i] eq -1 then begin
           if not str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).fit then message, errmsg
           if keyword_set(best) then str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).best = pars[i] $
           else str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).value = pars[i]
        endif else begin
           if not  (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].fit then message, errmsg
           if keyword_set(best) then (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].best = pars[i] $
           else (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].value = pars[i]
        endelse
     endfor
  endif else begin
     for i=0, npars-1 do begin
        if tofit[3,i] eq -1 then begin
           if not str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).fit then message, errmsg
           if keyword_set(best) then str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).best = transpose(pars[i,*]) $
           else str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).value = transpose(pars[i,*])
        endif else begin
           if not  (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].fit then message, errmsg
           if keyword_set(best) then (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].best = transpose(pars[i,*]) $
           else (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].value = transpose(pars[i,*])
        endelse
     endfor
  endelse     

end
