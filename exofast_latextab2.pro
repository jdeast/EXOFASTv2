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

pro exofast_latextab2, ss, title=title, label=label, caption=caption, texfile=texfile

;; open the texfile (or print to screen instead)
if n_elements(texfile) ne 0 then begin
    openw, lun, texfile, /get_lun
endif else lun = -1 ;; print to screen

; preamble stuff makes it an "includable" file in your tex document
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

;; number of columns
maxvalues = 0
for i=0, n_tags(ss)-1 do begin
   nvalues = n_elements(ss.(i)[*])
   if nvalues gt maxvalues then maxvalues=nvalues
endfor
printf, lun, '\begin{deluxetable}{lc' + strjoin(replicate('c',maxvalues)) + '}'

if n_elements(caption) ne 0 then printf, lun, '\tablecaption{' + caption + '}'
printf, lun, '\tablehead{\colhead{~~~Parameter} & \colhead{Units} & \multicolumn{' + strtrim(maxvalues,2) + '}{c}{Values}}'
printf, lun, '\startdata'

;; format for the line in the table
;; $symbol$ ... description (units) ... $value^{+upper}_{-lower}$\\
for i=0, n_tags(ss)-1 do begin
   
   nvalues = n_elements(ss.(i)[*])

   ;; the Side labels
   if n_tags(ss.(i)[0]) ne 0 then begin
      if tag_exist(ss.(i)[0], 'rootlabel') then $
         printf, lun, ss.(i)[0].rootlabel, format='("\sidehead{",a,"}")'
;      printf, lun, ss.(i)[0].rootlabel, ss.(i)[*].label, format='("\sidehead{",a,"&",' + strtrim(nvalues,2) + '("&",a),"}")'
;; how to do this (latex problem).. separate tables?
   endif

   ;; for each tag
   for k=0, n_tags(ss.(i)[0])-1 do begin
      
      if n_tags(ss.(i)[0].(k)) ne 0 then begin
         
         if tag_exist(ss.(i)[0].(k),'derive') then begin
            
            if (where(ss.(i)[*].(k).derive))[0] ne -1 then begin
               
               ;; format the units
               if ss.(i)[0].(k).unit eq '' then unit = '' $
               else unit = "(" + ss.(i)[0].(k).unit + ")"
               
               ;; before all the values
               header = string(ss.(i)[0].(k).latex,ss.(i)[0].(k).description,unit,$
                               format='("~~~~$",a,"$\dotfill &",a,x,a,"\dotfill ")')
               
               ;; put each value of the array in a new column
               values = ''
               for j=0, nvalues-1 do begin
                  
                  if ss.(i)[j].(k).derive then begin
                     ;; format the errors
                     if ss.(i)[j].(k).upper eq ss.(i)[j].(k).lower then begin
                        ;; fixed value, no errors
                        if ss.(i)[j].(k).upper eq '0.00' then error = '' $ 
                        else error = '\pm' + ss.(i)[j].(k).upper ;; symmetric errors, use +/-
                     endif else begin
                        ;; assymetric errors, use ^+_-
                        error = "^{+" + ss.(i)[j].(k).upper + "}_{-" +  $
                                ss.(i)[j].(k).lower + "}"
                     endelse
                     values = values + '&$' + ss.(i)[j].(k).medvalue + error + '$'
                  endif else values += '&--'
               endfor
      
               ;; print the line of the latex table
               printf, lun, header + values + '\\'

            endif
         endif
      endif
   endfor
endfor

; finish the table
printf, lun, '\enddata'
if n_elements(label) ne 0 then printf, lun, '\label{' + label + '}'
printf, lun, '\end{deluxetable}'
printf, lun, '\end{document}'

if n_elements(texfile) ne 0 then free_lun, lun

end


