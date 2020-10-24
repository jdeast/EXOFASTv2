function str2num, string

;; this converts a string to a number by multiplying the numeric
;; values together (a = 1, b=2, c=3, etc)

ndx = ['-1','1','2','3','4','5','6','7','8','9',$
       '0','a','b','c','d','e','f','g','h','i',$
       'j','k','l','m','n','o','p','q','r','s',$
       't','u','v','w','x','y','z']

num = 1UL
for i=0L, strlen(string)-1 do begin
   match = (where(ndx eq strlowcase(strmid(string,i,1))))[0]
   if match ne -1 then num *= match
endfor

return, num

end
