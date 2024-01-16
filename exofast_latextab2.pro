;+
; NAME:
;   EXOFAST_LATEXTAB2
;
; PURPOSE:
;   Prints the latex source to create a deluxe table of parameter values and
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

pro exofast_latextab2, ss, title=title, label=label, caption=caption, texfile=texfile, emulateapj=emulateapj

;; open the texfile (or print to screen instead)
if n_elements(texfile) ne 0 then begin
    openw, lun, texfile, /get_lun
endif else lun = -1 ;; print to screen

; preamble stuff makes it an "includable" file in your tex document
if keyword_set(emulateapj) then printf, lun, '\documentclass{emulateapj}' $
else printf, lun, '\documentclass{aastex62}'

printf, lun, '\providecommand{\bjdtdb}{\ensuremath{\rm {BJD_{TDB}}}}'
printf, lun, '\providecommand{\tjdtdb}{\ensuremath{\rm {TJD_{TDB}}}}'
printf, lun, '\providecommand{\feh}{\ensuremath{\left[{\rm Fe}/{\rm H}\right]}}'
printf, lun, '\providecommand{\teff}{\ensuremath{T_{\rm eff}}}'
printf, lun, '\providecommand{\teq}{\ensuremath{T_{\rm eq}}}'
printf, lun, '\providecommand{\ecosw}{\ensuremath{e\cos{\omega_*}}}'
printf, lun, '\providecommand{\esinw}{\ensuremath{e\sin{\omega_*}}}'
printf, lun, '\providecommand{\msun}{\ensuremath{\,M_\Sun}}'
printf, lun, '\providecommand{\rsun}{\ensuremath{\,R_\Sun}}'
printf, lun, '\providecommand{\lsun}{\ensuremath{\,L_\Sun}}'
printf, lun, '\providecommand{\mj}{\ensuremath{\,M_{\rm J}}}'
printf, lun, '\providecommand{\rj}{\ensuremath{\,R_{\rm J}}}'
printf, lun, '\providecommand{\me}{\ensuremath{\,M_{\rm E}}}'
printf, lun, '\providecommand{\re}{\ensuremath{\,R_{\rm E}}}'
printf, lun, '\providecommand{\fave}{\langle F \rangle}'
printf, lun, '\providecommand{\fluxcgs}{10$^9$ erg s$^{-1}$ cm$^{-2}$}'
printf, lun, '\usepackage{apjfonts}'
printf, lun, '\begin{document}'
if keyword_set(emulateapj) then printf, lun, '\LongTables' $
else printf, lun, '\startlongtable'

;; number of columns
maxvalues = 0
for i=0, n_tags(ss)-1 do begin
   nvalues = n_elements(ss.(i)[*])
   if nvalues gt maxvalues then maxvalues=nvalues
endfor
printf, lun, '\begin{deluxetable*}{lc' + strjoin(replicate('c',maxvalues)) + '}'

if n_elements(caption) ne 0 then printf, lun, '\tablecaption{' + caption + '}'
printf, lun, '\tablehead{\colhead{~~~Parameter} & \colhead{Description} & \multicolumn{' + strtrim(maxvalues,2) + '}{c}{Values}}'
printf, lun, '\startdata'

if tag_exist(ss,'rvepoch') then rvepoch = ss.rvepoch $
else rvepoch = 0d0

;; format for the line in the table
;; $symbol$ ... description (units) ... $value^{+upper}_{-lower}$\\
masternotes = ['Uses measured mass and estimated radius from \citet{Chen:2017}',$                                ;; 0
               'Uses measured radius and estimated mass from \citet{Chen:2017}',$                                ;; 1
               'Reference epoch = ' + string(rvepoch,format='(f0.6)'),$                                          ;; 2
               'This value ignores the systematic error and is for reference only',$                             ;; 3 (Rstar,SED, Teff,SED)
               "Time of conjunction is commonly reported as the ``transit time''",$                              ;; 4 
               "Time of minimum projected separation is a more correct ``transit time''",$                       ;; 5
               'At the epoch that minimizes the covariance between $T_C$ and Period',$                           ;; 6
               'In RV-only fits, we marginalize over a uniform cosi prior',$                                     ;; 7
               'Estimated eclipse depth assumes blackbodies',$                                                   ;; 8
               'Assumes no albedo and perfect redistribution',$                                                  ;; 9
               "Corresponds to static points in a star's evolutionary history. See \S2 in \citet{Dotter:2016}",$ ;; 10
               'The metallicity of the star at birth',$                                                          ;; 11
               'See \citet{Eastman:2023b} for a detailed description',$                                          ;; 12
               "\tjdtdb is the target's barycentric frame and corrects for light travel time",$                  ;; 13 
               "Use this to predict future transit times",$                                                      ;; 14
               "Use this to model TTVs, e",$                                                                     ;; 15
               "This If you are modeling TTVs, this is the transit time you want",$                              ;; 16
               "Time of superior conjunction is commonly reported as the ``occultation time''",$                 ;; 17
               'See Table 3 in \citet{Eastman:2019} for a detailed description of all parameters']               ;; 18

;; these are the names of parameters that are influenced by the Chen &
;; Kipping mass estimate (given a radius)
notendx = [['rstarsed','3'],$
           ['teffsed','3'],$
           ['tco','4'],$ ;; bjdtdb
           ['tc','4,13'],$ ;; tjdtdb
           ['tt','5,13,15'],$ ;; tjdtdb 
           ['t0','5,6,14'],$ ;; bjdtdb
           ['tso','4'],$ ;; bjdtdb
           ['ts','4,13,'],$ ;; tjdtdb
           ['te','5,13,15'],$ ;; tjdtdb
           ['te0','5,6,14'],$ ;; bjdtdb
           ['eclipsedepth36','8'],$
           ['eclipsedepth45','8'],$
           ['teq','9'],$
           ['eep','10'],$
           ['initfeh','11']]

if ss.ntran eq 0 then begin
   notendx = [[notendx],$
              ['mp','7'],$
              ['mpearth','7'],$
              ['mpsun','7'],$
              ['q','7'],$
              ['rhop','7'],$
              ['loggp','7'],$
              ['cosi','7'],$
              ['ideg','7'],$
              ['i','7']]
endif

;; these notes are only neccesary if a slope is fit
if ss.star[0].slope.fit then begin
   notendx = [[notendx],$
              ['gamma','2'],$
              ['slope','2'],$
              ['quad','2']]
endif

;; these notes are only neccesary if vcve is fit
if ss.planet[0].vcve.fit then begin
   notendx = [[notendx],$
              ['vcve','12']]
endif

;; these notes are only neccesary if the Chen & Kipping relation is used
junk = where(ss.chen,nusechen)
if nusechen gt 0 then begin
   if ss.ntel eq 0 then begin
      notendx = [[notendx],$
                 ['rhop','1'],$
                 ['loggp','1'],$
                 ['k','1'],$
                 ['logk','1'],$
                 ['logk','1'],$
                 ['mp','1'],$
                 ['mpearth','1'],$
                 ['mpsun','1'],$
                 ['msini','1'],$
                 ['msiniearth','1'],$
                 ['q','1']]

   endif else if ss.ntran eq 0 then begin
      notendx = [[notendx],$
                 ['p','0'],$
                 ['rp','0'],$
                 ['rpsun','0'],$
                 ['rpearth','0'],$
                 ['rhop','0'],$
                 ['loggp','0']]
   endif
endif

notes = ['']
nnotes = 0L           

for i=0, n_tags(ss)-1 do begin
   
   nvalues = n_elements(ss.(i)[*])

   ;; the Side labels
   if n_tags(ss.(i)[0]) ne 0 then begin
      if tag_exist(ss.(i)[0], 'rootlabel') then begin
;         printf, lun, ss.(i)[0].rootlabel, format='("\sidehead{",a,"}")'
         point_lun, -lun, pos
         printf, lun, '\smallskip\\\multicolumn{2}{l}{' + ss.(i)[0].rootlabel + '}' + strjoin('&' + ss.(i).label) + '\smallskip\\' 
;        if ss.(i)[0].label ne '' then printf, lun, '&' + strjoin('&' + ss.(i).label) + '\\'
;      printf, lun, ss.(i)[0].rootlabel, ss.(i)[*].label, format='("\sidehead{",a,"&",' + strtrim(nvalues,2) + '("&",a),"}")'
;; how to do this (latex problem)... separate tables?
      endif
   endif

   npars = 0L
   ;; for each tag
   for k=0, n_tags(ss.(i)[0])-1 do begin      

      ;; this captures the detrending variables
      if (size(ss.(i)[0].(k)))[1] eq 10 then begin ;; if it's a pointer
         if ptr_valid(ss.(i)[0].(k)) then begin ;; if it's not empty
            for l=0L, n_tags(*(ss.(i)[0].(k)))-1 do begin ;; loop through each tag
               if (size((*(ss.(i)[0].(k))).(l)))[2] eq 8 then begin ;; if it's an array of structures

                  ;; now I need to loop through all transits (j => columns) and
                  ;; all detrending parameters for each transit (m =>
                  ;; rows) to map out the table
                  maxdetrend = lonarr(nvalues)-1
                  for j=0L, nvalues-1 do begin
                     maxdetrend[j] = n_elements((*(ss.(i)[j].(k))).(l))
                  endfor
                  ;; the index here is for the transit with the most
                  ;; detrending variables
                  nrows = max(maxdetrend,ndx)

                  for m=0L, nrows-1 do begin ;; loop through each structure

                     ;; format the line
                     if tag_exist((*(ss.(i)[ndx].(k))).(l)[m],'derive') then begin ;; if it's a parameter
                           
                        ;; format the units
                        if (*(ss.(i)[ndx].(k))).(l)[m].unit eq '' then unit = '' $
                        else unit = "(" + (*(ss.(i)[ndx].(k))).(l)[m].unit + ")"
                        
                        ;; before all the values
                        header = string((*(ss.(i)[ndx].(k))).(l)[m].latex,$
                                        (*(ss.(i)[ndx].(k))).(l)[m].description,unit,$
                                        format='("~~~~$",a,"$\dotfill &",a,x,a,"\dotfill ")')
                     endif

                     ;; put each value of the array in a new column
                     values = ''
                     ngood = 0L
                     for j=0, nvalues-1 do begin
                        if n_elements((*(ss.(i)[j].(k))).(l)) gt (m) then begin
                           if (*(ss.(i)[j].(k))).(l)[m].derive then begin
                              ;; format the errors
                              if (*(ss.(i)[j].(k))).(l)[m].upper eq (*(ss.(i)[j].(k))).(l)[m].lower then begin
                                 ;; fixed value, no errors
                                 if (*(ss.(i)[j].(k))).(l)[m].upper eq '0.00' then error = '' $ 
                                 else error = '\pm' + (*(ss.(i)[j].(k))).(l)[m].upper ;; symmetric errors, use +/-
                              endif else begin
                                 ;; asymmetric errors, use ^+_-
                                 error = "^{+" + (*(ss.(i)[j].(k))).(l)[m].upper + "}_{-" +  $
                                         (*(ss.(i)[j].(k))).(l)[m].lower + "}"
                              endelse
                              scinote = (*(ss.(i)[j].(k))).(l)[m].scinote
                              if (*(ss.(i)[j].(k))).(l)[m].medvalue eq '0.00' and error eq '' then values += '&--' $
                              else values = values + '&$' + (*(ss.(i)[j].(k))).(l)[m].medvalue + error + scinote + '$'
                              ngood++
                           endif else values += '&--'
                        endif else values += '&--'
                     endfor
                     if ngood gt 0 then begin
                        ;; print the line of the latex table
                        printf, lun, header + values + '\\'
                        npars++
                     endif
                  endfor
               endif
            endfor
         endif
      endif else if n_tags(ss.(i)[0].(k)) ne 0 then begin
         ;; this captures everything else
         if tag_exist(ss.(i)[0].(k),'derive') then begin            
            if (where(ss.(i)[*].(k).derive))[0] ne -1 then begin
               
               ;; format the units
               if ss.(i)[0].(k).unit eq '' then unit = '' $
               else unit = "(" + ss.(i)[0].(k).unit + ")"

               footnote = ''
               ;; add appropriate tablenotes
               match = where(notendx[0,*] eq ss.(i)[0].(k).label)
               if match[0] ne -1 then begin
                  splitnotes = long(strsplit(notendx[1,match[0]],',',/extract))
                  notestr = lonarr(n_elements(splitnotes))
                  for ii=0L, n_elements(splitnotes)-1 do begin
                     matchexisting = (where(notes eq masternotes[splitnotes[ii]]))[0]
                     
                     if matchexisting eq -1 then begin
                        nnotes++
                        notes = [notes,masternotes[splitnotes[ii]]]
                        notestr[ii] = nnotes
                     endif else notestr[ii] = matchexisting
                        
                  endfor
                  sorted = sort(notestr)
                  footnote += '$^{' + strjoin(strtrim(notestr[sorted],2),',') + '}$'
               endif
               
               ;; before all the values
               header = string(ss.(i)[0].(k).latex,ss.(i)[0].(k).description+footnote,unit,$
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
                        ;; asymmetric errors, use ^+_-
                        error = "^{+" + ss.(i)[j].(k).upper + "}_{-" +  $
                                ss.(i)[j].(k).lower + "}"
                     endelse
                     scinote = ss.(i)[j].(k).scinote
                     values = values + '&$' + ss.(i)[j].(k).medvalue + error + scinote + '$'
                  endif else values += '&--'

               endfor
      
               ;; print the line of the latex table
               printf, lun, header + values + '\\'
               npars++

            endif
         endif
      endif
   endfor
   ;; if no parameters in this section, rewind the file to erase the header
   if npars eq 0L then point_lun, lun, pos $
   else point_lun, -lun, pos
      
endfor

; finish the table
printf, lun, '\enddata'
if n_elements(label) ne 0 then printf, lun, '\label{' + label + '}'
printf, lun, '\tablenotetext{}{' + masternotes[18] + '}'
if nnotes gt 0 then begin
   for i=1L, nnotes do begin
      printf, lun, '\tablenotetext{' + strtrim(i,2) + '}{' + notes[i] + '}'
   endfor
endif
printf, lun, '\end{deluxetable*}'

;; create the bibliography
;; use $EXOFAST_PATH/References.bib
printf, lun, '\bibliographystyle{apj}'
printf, lun, '\bibliography{References}'

printf, lun, '\end{document}'
defsysv, '!GDL', exists=runninggdl
if not runninggdl then truncate_lun, lun

if n_elements(texfile) ne 0 then free_lun, lun

end


