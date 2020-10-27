function bitmask2arr, mask

nbits = ceil(alog(mask)/alog(2d0))+1
arr = bytarr(nbits)

mask0 = mask ;; don't overwrite input
while mask0 gt 0 do begin

   bit = floor(alog(mask0)/alog(2))
   arr[bit] = 1B
   mask0 -= 2^bit

endwhile

;; reverse function:
print, total(2^arr)

return, arr

end
