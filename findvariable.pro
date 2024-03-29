function findvariable, structure, varname, depth=depth, silent=silent, logname=logname, count=count

tagname = 'fit' ;; only parameter structures have the 'fit' tag
mapline = [-1,-1,-1,-1,-1,-1]

;if strpos(varname,'_') eq -1 then varnames = [strupcase(varname),'0'] $
;else 
;stop

varnames = strupcase(strsplit(varname,'_',/extract))
if n_elements(varnames) eq 1 then varnames = [varnames,'0']

ntags = n_tags(structure)
for i=0L, ntags-1 do begin

   ;; if it's not a pointer or structure, skip it   
   sz = size((structure.(i)))
   if sz[1] ne 10 and sz[2] ne 8 then continue 

   for j=0L, n_elements(structure.(i))-1 do begin

      sz = size((structure.(i))[j])

      ;; if it's a pointer, recursively search for more parameters
      if sz[1] eq 10 then begin
         if (structure.(i))[j] eq !null then continue
         if n_elements(depth) eq 0 then newdepth=[i,j]
         if n_elements(depth) eq 2 then newdepth=[depth,i,j]
         mapline = findvariable(*((structure.(i))[j]), varname, depth=newdepth, /silent, count=count)
         if mapline[0] ne -1 then return, mapline
         continue
      endif

      ;; otherwise, if it's not a structure, skip it
      if sz[2] ne 8 then continue
      
      ;; it's a structure, but not a parameter structure,
      ;; recursively search for parameters
      if ~tag_exist((structure.(i))[j], tagname, /top_level) then begin
         if n_elements(depth) eq 0 then newdepth=[i,j]
         if n_elements(depth) eq 2 then newdepth=[depth,i,j]
         mapline = findvariable((structure.(i))[j], varname, depth=newdepth, /silent, count=count)
         if mapline[0] ne -1 then return, mapline
         continue ;; didn't find it in stucture.(i))[j], keep searching the primary structure
      endif

      ;if varnames[0] eq 'C1' and strupcase((structure.(i))[j].label) eq 'C1' then stop

      ;; its label doesn't match varname, keep searching
      if strupcase((structure.(i))[j].label) ne varnames[0] then continue
     
      ;; add to the depth of the map
      if n_elements(depth) eq 0 then depth = [i,j] $
      else depth = [depth,i,j]

      if n_elements(depth) ge 4 and depth[1] eq long(varnames[1]) then begin
         ;; hack for variables in pointers... 
         ;; there must be a more elegant way
         mapline[0:n_elements(depth)-1] = depth
         return, mapline
      endif else if count eq long(varnames[1]) then begin
         mapline[0:n_elements(depth)-1] = depth
         return, mapline
      endif else count++

   endfor
endfor

if mapline[0] eq -1 then begin
   if ~keyword_set(silent) then begin
      printandlog, 'WARNING: ' + varname + ' not found!', logname
   endif
endif

return, mapline

end
