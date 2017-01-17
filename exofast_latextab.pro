;+
; NAME:
;   EXOFAST_LATEXTAB
;
; PURPOSE:
;   Prints the latex source to create a delux table of parameter values and
;   their errors. The errors are automatically rounded to two
;   significant digits and the values are rounded to the last
;   significant digit of the error.
;
;   This is intended to reduce typos and alleviate the tedium of making
;   tables.
;
; CALLING SEQUENCE:
;   exofast_latextab, array, ['TEXFILE', PARNAMES=,
;                     UNITS=,TITLE=,LABEL=, CAPTION=,SIDELABLES=,ORDER]
;
; INPUTS: 
;   ARRAY    - A 2 or 3 column by NPARS rows array with the value,
;              +error, -error (e.g., the output from EXOFAST_PLOTDIST). If
;              only two columns are specified, all errors are assumed
;              to be symmetric. 
;
; OPTIONAL INPUTS:
;   TEXFILE    - If specified, the program will write to this named file
;                instead of the screen
;   PARNAMES   - If specified, the table will include the parameter
;                names as specified by this keyword. Otherwise they will
;                be left blank. Parameter names will be bracketed by
;                "$" so they are written in math mode.
;   UNITS      - If specified, the table will include this
;                discriptor/units next to the parameters. Otherwise they will
;                be left blank.
;   TITLE      - The title of the table
;   LABEL      - The label for the table
;   CAPTION    - The caption of the table
;   ORDER      - An array of indices that specify the order of
;                outputs. Default is the input order. If SIDELABLES
;                are used, -1 should be inserted in each place the
;                next label should be printed.
;   SIDELABELS - An array of side labels to print. The number of
;                elements must match the number of -1s in the ORDER
;                array (this is where each lable will be inserted).
; OUTPUTS:
;   Prints the latex source code to generate a table of values.
;
; Modification History:
;  2011/03/10 - Written by Jason Eastman (OSU)
;-

pro exofast_latextab, array, texfile, parnames=parnames, units=units, $
                      title=title, label=label, caption=caption, $
                      sidelabels=sidelabels, order=order

sz = size(array)
if sz[0] ne 2 then message, 'ERROR: array must be two dimensional'

;; by default, print the parameters in order
if n_elements(order) eq 0 then order = indgen(sz[2])

;;  check if SIDELABELS is properly formatted
if n_elements(sidelabels) ne 0 then begin
   if n_elements(where(order eq -1)) ne n_elements(sidelabels) then $
      message, "ERROR: 'ORDER' and 'SIDELABELS' do not match: " + $
               "see documentation for proper syntax"
   prefix = '~~~$' ;; indent rows other than the label
endif else prefix = '$'

;; don't modify input parameter names
parnames0 = parnames

;; open the texfile (or print to screen instead)
if n_elements(texfile) ne 0 then begin
    openw, lun, texfile, /get_lun
endif else lun = -1 ;; print to screen

;; number of columns, duplicate errors if necessary
if sz[1] eq 2 then begin
    array = [array,transpose([transpose(array[1,*])])]
endif else if sz[1] ne 3 then message, 'ERROR: array must have 2 or 3 columns'

;; must specify parnames
if n_elements(parnames) ne sz[2] then message, 'ERROR: the dimensions of PARNAMES must match ARRAY[0,*]'
parnames = prefix + parnames + '$'

;; format units
for i=0, sz[2]-1 do begin
   if n_elements(units) eq sz[2] then begin
      if units[i] ne '' then parnames[i]+='\dotfill &'+units[i]+'\dotfill' $
      else parnames[i]+='\dotfill & \dotfill'
   endif else parnames[i]+='\dotfill  & \dotfill'
endfor

if n_elements(parnames) ne sz[2] then parnames = strarr(sz[2])

nchar = max(strlen(parnames))
nlabel = 0

; preamble stuff
printf, lun, '\documentclass{emulateapj}'
printf, lun, '\newcommand{\bjdtdb}{\ensuremath{\rm {BJD_{TDB}}}}'
printf, lun, '\newcommand{\feh}{\ensuremath{\left[{\rm Fe}/{\rm H}\right]}}'
printf, lun, '\newcommand{\teff}{\ensuremath{T_{\rm eff}}}'
printf, lun, '\newcommand{\ecosw}{\ensuremath{e\cos{\omega_*}}}'
printf, lun, '\newcommand{\esinw}{\ensuremath{e\sin{\omega_*}}}'
printf, lun, '\newcommand{\msun}{\ensuremath{\,M_\Sun}}'
printf, lun, '\newcommand{\rsun}{\ensuremath{\,R_\Sun}}'
printf, lun, '\newcommand{\lsun}{\ensuremath{\,L_\Sun}}'
printf, lun, '\newcommand{\mj}{\ensuremath{\,M_{\rm J}}}'
printf, lun, '\newcommand{\rj}{\ensuremath{\,R_{\rm J}}}'
printf, lun, '\newcommand{\me}{\ensuremath{\,M_{\rm E}}}'
printf, lun, '\newcommand{\re}{\ensuremath{\,R_{\rm E}}}'
printf, lun, '\newcommand{\fave}{\langle F \rangle}'
printf, lun, '\newcommand{\fluxcgs}{10$^9$ erg s$^{-1}$ cm$^{-2}$}'
printf, lun, '\usepackage{apjfonts}'
printf, lun, '\begin{document}'

; begin table
printf, lun, '\begin{deluxetable}{lcc}'
if n_elements(caption) ne 0 then printf, lun, '\tablecaption{' + caption + '}'
printf, lun, '\tablehead{\colhead{~~~Parameter} & \colhead{Units} & \colhead{Value}}'
printf, lun, '\startdata'

npar = n_elements(order)
for i=0, npar-1 do begin

   ;; print the side label
   if order[i] eq -1 then begin
      printf, lun, '\sidehead{' + sidelabels[nlabel] + '}'
      nlabel++
   endif else begin

      ;; round the high error to 2 sig figs
      if array[1,order[i]] eq !values.d_infinity then begin
         errhi = '\infty'
         exphi = !values.d_infinity
      endif else begin
         exphi=fix(alog10(array[1,order[i]]))
         if (array[1,order[i]] lt 1d0) then exphi=exphi-1
         roundhi=round(array[1,order[i]]/10.d0^(exphi-1d0),/L64)*10.d0^(exphi-1d0)
         if (roundhi gt 10) then errhi = strtrim(round(roundhi,/L64),2) $
         else errhi = string(roundhi,format='(f255.'+strtrim(1-exphi,2)+')')
      endelse

      ;; round the low error to 2 sig figs
      if array[2,order[i]] eq !values.d_infinity then begin
         errlo = '\infty'
         explo = !values.d_infinity
      endif else begin
         explo=fix(alog10(array[2,order[i]]))
         if (array[2,order[i]] lt 1d0) then explo=explo-1
         roundlo=round(array[2,order[i]]/10.d0^(explo-1d0),/L64)*10.d0^(explo-1d0)
         if (roundlo gt 10) then errlo = strtrim(round(roundlo,/L64),2) $
         else errlo = string(roundlo,format='(f255.'+strtrim(1-explo,2)+')')
      endelse

      ;; round the value to the greater number of sig figs
      ndec = long(1 - (exphi < explo))
      if ndec eq 0 then value = string(array[0,order[i]],format='(i255)') $
      else if ndec lt 0 then $
         value=round(round(array[0,order[i]]/10.d0^(-ndec),/L64)*10.d0^(-ndec),/L64) $
      else value = string(array[0,order[i]],format='(f255.'+strtrim(ndec,2)+')')
      
      ;; print a row in the table
      if errlo eq errhi then begin
         if errlo eq 0.d0 then begin
            ;; no error (value)
            if i eq npar-1 then format='("$",a,"$")' $
            else format='("$",a,"$ \\")'
            row = string(array[0,order[i]], format=format)
         endif else begin
            ;; same error (value +/- err)
            if i eq npar-1 then format='("$",a," \pm ",a,"$")' $
            else format='("$",a," \pm ",a,"$ \\")'
            row = string(value,errlo,format=format)
         endelse
      endif else begin
         ;; asymmetric error (value_{-errlo}^{+errhi})

         if i eq npar-1 then format='("$",a,"_{-",a,"}^{+",a,"}$")' $
         else format='("$",a,"_{-",a,"}^{+",a,"}$ \\")'
         row = string(value, errlo, errhi, format=format)
      endelse
      
      printf,lun,string(parnames[order[i]],format='(a'+strtrim(nchar,2)+')')+$
             ' & ' + strcompress(row,/remove_all)
   endelse
   
endfor

; finish up
printf, lun, '\enddata'
if n_elements(label) ne 0 then printf, lun, '\label{' + label + '}'
printf, lun, '\end{deluxetable}'
printf, lun, '\end{document}'

if n_elements(texfile) ne 0 then free_lun, lun

parnames = parnames0

end


