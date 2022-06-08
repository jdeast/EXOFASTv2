function find_by_tag, str, tag


;; for each planet 
;; populate the planet fitting parameters
;; planetary labels, a bit optimistic...
plabels = ['b','c','d','e','f','g','h','i','j','k','l','m','n',$
           'o','p','q','r','s','t','u','v','w','x','y','z'] 


;; determine the subscript
tmp = strsplit(tag,'_',/extract)
if n_elements(tmp) eq 2 then tagnum = tmp[1] $
else tagnum = '0'
taglabel = tmp[0]

;; allow labels 'b','c','d', etc.
numndx = (where(plabels eq strtrim(tagnum,2)))[0]
if numndx ne -1 then tagnum = numndx
tagnum = long(tagnum)

;; look for the name in the structure
for i=0L, n_tags(str)-1 do begin
   if (n_elements(str.(i))-1) ge tagnum then begin
      for k=0, n_tags(str.(i)[tagnum])-1 do begin
         
         ;; this captures the detrending variables
         if (size(str.(i)[tagnum].(k)))[1] eq 10 then begin ;; if it's a pointer
            if ptr_valid(str.(i)[tagnum].(k)) then begin
               for l=0L, n_tags(*(str.(i)[tagnum].(k)))-1 do begin

                  ;; if not a structure, skip it
                  if (size((*(str.(i)[tagnum].(k))).(l)))[2] ne 8 then continue
                     
                  valid = 0
                  if strpos(strupcase(taglabel),'C') eq 0 then begin
                     if strmid((*(str.(i)[tagnum].(k))).(l)[0].label,0,1) eq 'C' then begin ;; if it's the additive variable
                        if (*(str.(i)[tagnum].(k))).nadd gt 0 then begin
                           detrendnum = long((strsplit(strupcase(taglabel),'C',/extract))[0])
                           if detrendnum lt (*(str.(i)[tagnum].(k))).nadd then valid = 1
                        endif
                     endif
                  endif else if strpos(strupcase(taglabel),'M') eq 0 then begin
                     if strmid((*(str.(i)[tagnum].(k))).(l)[0].label,0,1) eq 'M' then begin
                        if (*(str.(i)[tagnum].(k))).nmult gt 0 then begin
                           detrendnum = long((strsplit(strupcase(taglabel),'M',/extract))[0])
                           if detrendnum lt (*(str.(i)[tagnum].(k))).nmult then valid = 1
                        endif
                     endif
                  endif
                  
                  if valid then begin
                     if (*(str.(i)[tagnum].(k))).(l)[detrendnum].label eq strupcase(taglabel) then return, (*(str.(i)[tagnum].(k))).(l)[detrendnum]
                  endif
               
               endfor
            endif          
         endif else if n_tags(str.(i)[tagnum].(k)) ne 0 then begin
            ;; and this captures everything else
            if tag_exist(str.(i)[tagnum],taglabel,index=ndx) then return, str.(i)[tagnum].(ndx)
         endif         
      endfor
   endif
endfor

return, -1

end              
  
