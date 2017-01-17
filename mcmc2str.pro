;; copies str while stripping all parameters that don't have
;; fit=1 or derived=1, and expand the value tag to an NSTEPS x NCHAINS
;; array 


pro mcmc2str, pars, str

newpar

for i=0, n_tags(str)-1 do begin

   if n_tags(str.(i)) ne 0 then begin
      ;; create new object structures
      for k=0, n_tags(str.(i))-1 do begin
;      if n_tags(str.(i))
         if str.(i)[0].(k).fit or str.(i)[0].(k).derive then begin
            newstr = create_struct(newstr, str.(i)[0].(k))
         help, str.(i)[0].(k)

stop
      endfor

   endif else begin
      ;; copy metadata directly

   endelse

endfor

stop
newpar = 

for i=0, n_tags(str)-1 do begin
   for j=0, n_elements(str.(i))-1 do begin
      for k=0, n_tags(str.(i)[j])-1 do begin

         if n_tags(str.(i)[j].(k)) eq 0 then begin
            
            if tag_exist(str.(i)[j].(k),'fit') then begin
               
               if str.(i)[j].(k).fit then tofit = [[tofit],[i,j,k]]
               
            endif else begin
               ;; copy metadata directly
               newstr = create_struct(newstr, newpar)

            endelse
         endif else begin
            ;; copy metadata directly

         endelse

      endfor
   endfor
endfor

  tofit = *(str.tofit)
  npars = n_elements(tofit[0,*])

  for i=0, npars-1 do $
     str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).value = transpose(pars[i,*])

  names = ['a','b']
  values = [1,1]
  newstr = create_struct(names,values)

end
