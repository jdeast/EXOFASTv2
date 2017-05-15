pro gridpars

files = file_search('sed/nextgenfin/nextgenfin2/*.idl',count=nfiles)


teff = []
logg = []
feh = []
alpha = []

for i=0L, nfiles-1 do begin
   
   filename = file_basename(files[i],'.NextGen.spec.idl')
   len = strlen(filename)
   teff = [teff,double(strmid(filename,3,len-15))*100d0]
   logg = [logg,double(strmid(filename,len-12,4))]
   feh = [feh,double(strmid(filename,len-8,4))]
   alpha = [alpha,double(strmid(filename,len-4,4))]

endfor

tefforig = teff
loggorig = logg
fehorig = feh
alphaorig = alpha

teff = teff[sort(teff)]
teffs = teff[uniq(teff)]

logg = logg[sort(logg)]
loggs = logg[uniq(logg)]

feh = feh[sort(feh)]
fehs = feh[uniq(feh)]

alpha = alpha[sort(alpha)]
alphas = alpha[uniq(alpha)]

nteffs = n_elements(teffs)
nloggs = n_elements(loggs)
nfehs = n_elements(fehs)
nalphas = n_elements(alphas)
grid = bytarr(nteffs,nloggs,nfehs,nalphas)
for i=0,nteffs-1L do begin
   for j=0, nloggs-1L do begin
      for k=0, nfehs-1L do begin
         for l=0, nalphas-1L do begin
            match = where(tefforig eq teffs[i] and loggorig eq loggs[j] and $
                          fehorig eq fehs[k] and alphaorig eq alphas[l])
            if match[0] ne -1 then begin
               grid[i,j,k,l] = 1B
               oplot, teffs[i],loggs[j]
            endif
         endfor
      endfor
   endfor
endfor

stop
end
