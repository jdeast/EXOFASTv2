function bitmask2arr, mask, nbits=nbits

if n_elements(nbits) eq 0 then begin
   if mask eq 0 then nbits=1 $
   else nbits = ceil(alog(mask)/alog(2d0))+1
endif

arr = bytarr(nbits)

mask0 = mask ;; don't overwrite input
while mask0 gt 0 do begin

   bit = floor(alog(mask0)/alog(2))
   arr[bit] = 1B
   mask0 -= 2^bit

endwhile

;; reverse function:
;; print, total(2^where(arr))

return, arr

end
