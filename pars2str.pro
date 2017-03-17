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

pro pars2str, pars, str

  tofit = *(str.tofit)
  npars = n_elements(tofit[0,*])

  if (size(pars))[0] eq 1 then begin
     ;; populate the structure with the parameter array
     for i=0, npars-1 do begin
        if tofit[3,i] eq -1 then begin
           str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).value = pars[i]
        endif else begin
           (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].value = pars[i]
        endelse
     endfor
  endif else begin
     for i=0, npars-1 do begin
        if tofit[3,i] eq -1 then begin
           str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).value = transpose(pars[i,*])
        endif else begin
           (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].value = transpose(pars[i,*])
        endelse
     endfor
  endelse     

end
