pro extractpgno, filename, primary=primary, secondary=secondary, $
                 stacked_primary=stacked_primary, $
                 stacked_secondary=stacked_secondary

openr, lun, filename, /get_lun
line = ''
pageno = 0
primary = []
secondary = []
stacked_primary = []
stacked_secondary = []
while not eof(lun) do begin
   readf, lun, line
   ;; this marks a new page
   if strpos(line, "BeginPageSetup") ne -1 then begin
      pageno = pageno+1
      isprimary = 0
      issecondary = 0
      stacked=0
   endif

   if strpos(line, '(S)') ne -1 then issecondary = 1
   if strpos(line, '(C)') ne -1 then isprimary = 1

   ;; stacked plot
   if strpos(line, '(Normalized flux + Constant)') ne -1 then begin
      if issecondary then stacked_secondary = [stacked_secondary,pageno]
      if isprimary then stacked_primary = [stacked_primary,pageno]
   endif else if strpos(line, '(O-C)') ne -1 then begin
      if issecondary then secondary = [secondary,pageno]
      if isprimary then primary = [primary,pageno]
   endif      

endwhile
free_lun, lun

end
