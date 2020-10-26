pro testcorr,lat,lon,update=update

jd = julday(3,20,2009,11,44);,julday(6,21,2009,5,45),$
;      julday(9,22,2009,21,18),julday(12,21,2009,17,47)]

jd = julday(12,21,2009,12,00)
jd = julday(6,21,2009,12,00)


;jd = jd + 365.242199d0*dindgen(5)/4

jd = systime(/julian,/utc) + findgen(10000)/9999.d0*365.25
ra = 270.d0
dec = 0.d0
;obsname = 'winer'


tbase = 2400000.5;floor(jd)
jd -= tbase

taifile = find_with_def('tai-utc.dat','ASTRO_DATA')
clock_corr = tai_utc(jd + tbase, filename=taifile) + 32.184d0
jd_tt = jd + clock_corr/86400.d0
bjd = barycen(jd_tt, ra, dec,/time_diff) + clock_corr


npoints = 1
;lat = -23.44d0
;lon = 0.d0;dindgen(npoints)/(npoints-1.d0)*360

t = utc2tcb(jd, ra, dec, lat, lon, obsname=obsname, tbase=tbase,/time_diff,update=update)
print, string(bjd - t, format='(f30.20)')

stop
x = r_obs[0]
y = r_obs[1]
xrange=[x-10000,x+10000]
yrange=[y-10000,y+10000]
plot, [0],[0],/iso, xrange=xrange,yrange=yrange

for i=0, npoints-1 do begin
    t = utc2tcb(jd, ra, dec, lat, lon[i], obsname=obsname, tbase=tbase,/time_diff,r_obs=r_obs)
    oplot, [r_obs[0]],[r_obs[1]],psym=1

endfor







;stop





stop

end
