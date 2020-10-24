pro combinetex2, file1, file2, outfile

files = [file1,file2]

for j=0, n_elements(files)-1 do begin

   openr, lun, files[j], /get_lun
   nlines = file_lines(files[j])
   lines = strarr(nlines)
   names = strarr(nlines)

   tmp = ''
   for i=0, nlines-1 do begin
      readf, lun, tmp
      lines[i] = tmp
      names[i] = (strsplit(tmp,'&',/extract))[0]
   endfor
   free_lun, lun

   if j eq 0 then begin
      lines1 = lines 
      names1 = names
   endif else begin
      lines2 = lines
      names2 = names
   endelse

endfor

;; latex name, variablename, value
varnames = [['                               ~~~$M_{*}$\dotfill ',    'mstar',''],$
            ['                             ~~~$R_{*}$\dotfill ',      'rstar',''],$
            ['                         ~~~$L_{*}$\dotfill ',          'lstar',''],$
            ['                             ~~~$\rho_*$\dotfill ',     'rhostar',''],$
            ['                                    ~~~$Age$\dotfill ', 'age',''],$
            ['                  ~~~$\log{g_*}$\dotfill ',             'logg',''],$
            ['~~~~$T_{\rm eff}$\dotfill','teff',''],$
            ['~~~~$[{\rm Fe/H}]$\dotfill','feh',''],$
            ['~~~~$[{\rm Fe/H}]_{0}$\dotfill','initfeh',''],$            
            ['                  ~~~$\teff$\dotfill *',                'teff',''],$
            ['                                 ~~~$\feh$\dotfill ',   'feh',''],$
            ['             ~~~$v\sin{I_*}$\dotfill ',                 'vsini',''],$
            ['           ~~~$\lambda$\dotfill ',                      'lambda',''],$
            ['                                  ~~~$d$\dotfill ',     'dist',''],$
            ['                                   ~~~$e$\dotfill ',    'eccen',''],$
            ['        ~~~$\omega_*$\dotfill ',                        'omega',''],$
            ['                                  ~~~$P$\dotfill ',     'period',''],$
            ['                           ~~~$a$\dotfill ',            'a',''],$
            ['                                 ~~~$M_{P}$\dotfill ',  'mp',''],$
            ['                               ~~~$R_{P}$\dotfill ',    'rp',''],$
            ['                           ~~~$\rho_{P}$\dotfill ',     'rhop',''],$
            ['                      ~~~$\log{g_{P}}$\dotfill ',       'loggp',''],$
            ['               ~~~$T_{eq}$\dotfill ',                   'teq',''],$
            ['                           ~~~$\Theta$\dotfill ',       'theta',''],$
            ['                   ~~~$\fave$\dotfill ',                'fave',''],$
            ['       ~~~$T_C$\dotfill ',                              'tc',''],$
            ['               ~~~$T_{P}$\dotfill ',                    'tp',''],$
            ['                        ~~~$K$\dotfill ',               'k',''],$
            ['                           ~~~$K_R$\dotfill ',        'kr',''],$
            ['                    ~~~$M_P\sin{i}$\dotfill ', 'msini',''],$
            ['                           ~~~$M_{P}/M_{*}$\dotfill ','massratio',''],$
            ['                       ~~~$u$\dotfill ',          'u',''],$
            ['                              ~~~$\gamma_{EXPERT}$\dotfill ', 'gamma',''],$
            ['\ecosw',     'ecosw',''],$
            ['\esinw',     'esinw',''],$
            ['M_P\sin{i}', 'msini',''],$
            ['f(m1,m2)',   'massfunc',''],$
            ['R_{P}/R_{*}','rprs',''],$
            ['a/R_*',      'ar',''],$
            ['i',          'inc',''],$
            ['b',          'b',''],$
            ['\delta',     'delta',''],$
            ['T_{FWHM}',   'tfwhm',''],$
            ['\tau',       'tau',''],$
            ['T_{14}',     't14',''],$
            ['P_{T}',      'pt',''],$
            ['P_{T,G}',    'ptg',''],$
            ['TTV_{?}',    'ttv',''],$
            ['F_{?,?}',    'f',''],$
            ['C_{?,?}',    'c',''],$
            ['T_{C,?}',    'tc',''],$
            ['u_{1?}',     'u1',''],$
            ['u_{2?}',     'u2',''],$
            ['T_{S}',      'ts',''],$
            ['b_{S}',      'bs',''],$
            ['T_{S,FWHM}', 'tsfwhm',''],$
            ['\tau_{S}',   'taus',''],$
            ['P_{S}',      'ps',''],$
            ['P_{S,G}',    'psg','']]

values = strarr(n_elements(varnames))

varnames = strarr(4,n_elements(lines1))

openw, outlun, outfile, /get_lun
deffile = 'def.tex'
openw, deflun, deffile, /get_lun



for i=0, n_elements(lines1)-1 do begin

   match = where(names2 eq names1[i],nmatch)
   entries1 = strsplit(lines1[i],'&',/extract)

   if n_elements(entries1) eq 1 then begin
      if entries1 eq '\begin{deluxetable}{lcc}' then begin
         printf, outlun, '\begin{deluxetable*}{lccc}'
      endif else if entries1 eq '\end{deluxetable}' then begin
         printf, outlun, '\end{deluxetable*}'
      endif else if entries1 eq '\end{document}' then begin
         break
      endif else begin
         printf, outlun, lines1[i]
         if nmatch eq 1 then names2[match] = ''
      endelse
   endif else begin
      if strpos(lines1[i],'multicolumn') ne -1 then begin
         ;; skip it
      endif else if entries1[0] eq '\tablehead{\colhead{~~~Parameter} ' then begin
         printf, outlun, '\tablehead{\colhead{~~~Parameter} & \colhead{Units} & \colhead{Value 1} & \colhead{Value 2}}'
      endif else begin

         entries1 = strsplit(lines1[i],'&',/extract)
         entries1[2] = (strsplit(entries1[2],'\\\\',/extract,/regex))[0]

         if nmatch eq 0 then begin
            entries2[2] = '--'
         endif else begin
            entries2 = strsplit(lines2[match],'&',/extract)
            entries2[2] = (strsplit(entries2[2],'\\\\',/extract,/regex))[0]
            names2[match] = ''
         endelse
         

         ;; strip the latexname and make it a legal variable name
         varname1 = (strsplit(entries1[0],'$',/extract))[1]
;         if strpos(varname1,'^') ne -1 then stop
         varname1 = strjoin(strsplit(varname1,'*',/extract),'star')
         varname1 = strjoin(strsplit(varname1,'\^2',/extract,/preserve_null,/regex),'sq')
         varname1 = strjoin(strsplit(varname1,'\^\{2\}',/extract,/preserve_null,/regex),'sq')
         varname1 = strjoin(strsplit(varname1,' ',/extract,/preserve_null),'')
         varname1 = strjoin(strsplit(varname1,'.',/extract,/preserve_null),'point')
         varname1 = strjoin(strsplit(varname1,'1',/extract,/preserve_null),'one')
         varname1 = strjoin(strsplit(varname1,'2',/extract,/preserve_null),'two')
         varname1 = strjoin(strsplit(varname1,'3',/extract,/preserve_null),'three')
         varname1 = strjoin(strsplit(varname1,'4',/extract,/preserve_null),'four')
         varname1 = strjoin(strsplit(varname1,'5',/extract,/preserve_null),'five')
         varname1 = strjoin(strsplit(varname1,'6',/extract,/preserve_null),'six')
         varname1 = strjoin(strsplit(varname1,'7',/extract,/preserve_null),'seven')
         varname1 = strjoin(strsplit(varname1,'8',/extract,/preserve_null),'eight')
         varname1 = strjoin(strsplit(varname1,'9',/extract,/preserve_null),'nine')
         varname1 = strjoin(strsplit(varname1,'0',/extract,/preserve_null),'zero')
         varname1 = '\'  + strlowcase(strjoin(strsplit(varname1,'\{_}()/,',/extract),'')) + 'val'
         
         varname2 = varname1 + 'two'
         varname1 = varname1 + 'one'
         
         varnames[0,i] = varname1
         varnames[1,i] = varname2
         
         printf, outlun, strtrim(entries1[0:1],2), varname1, varname2, format='(3(a50," &"),a50,"\\")'

         val1 = '\ensuremath{' + (strsplit(strtrim(entries1[2],2),'$',/extract))[0] + '}'
         val2 = '\ensuremath{' + (strsplit(strtrim(entries2[2],2),'$',/extract))[0] + '}'

         printf, deflun, '\newcommand{' + varname1 + '}{' + val1 + '}'
;         printf, deflun, '\newcommand{' + varname1 + '}{' + strtrim(entries1[2],2) + '}'

         printf, deflun, '\newcommand{' + varname2 + '}{' + val2 + '}'
;         printf, deflun, '\newcommand{' + varname2 + '}{' + strtrim(entries2[2],2) + '}'
         
      endelse
   endelse

endfor
free_lun, outlun
free_lun, deflun

end
