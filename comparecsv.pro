pro comparecsv, csv1name, csv2name, threshhold=threshhold, sigdiff=sigdiff

if n_elements(threshhold) eq 0 then threshhold = 0.5d0

readcol, csv1name, parname1, value1, uerr1, lerr1, format='a,d,d,d', delimiter=',', comment='#',/silent
readcol, csv2name, parname2, value2, uerr2, lerr2, format='a,d,d,d', delimiter=',', comment='#',/silent

sigdiff = value1*0d0
matches = value1*0L -1L

for i=0L, n_elements(parname1)-1 do begin
   match = (where(parname2 eq parname1[i]))[0]
   if match eq -1 then continue
   matches[i] = match

   if value1[i] gt value2[match] then begin
      sigdiff[i] = (value1[i]-value2[match])/sqrt(lerr1[i]^2+uerr2[match]^2)
   endif else begin
      sigdiff[i] = (value1[i]-value2[match])/sqrt(uerr1[i]^2+lerr2[match]^2)
   endelse

endfor

bad = where(abs(sigdiff) gt threshhold)

if bad[0] ne -1 then forprint, parname1[bad], value1[bad], value2[matches[bad]], sigdiff[bad], format='(a25, x, f14.6, x, f14.6, x, f14.6)'

end

