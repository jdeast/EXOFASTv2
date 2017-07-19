pro updatetime

;; big correction (~1 second/yr)
spawn, 'wget --ftp-user=anonymous -NP $ASTRO_DATA ftp://maia.usno.navy.mil/ser7/tai-utc.dat'

;; small corrections (ms)
spawn, 'wget --ftp-user=anonymous -NP $ASTRO_DATA ftp://maia.usno.navy.mil/ser7/finals.all'
file_copy, getenv('ASTRO_DATA') + '/finals.all', getenv('ASTRO_DATA') + '/iers_final_a.dat', /overwrite

;; smallest correction (us)

;; did they stop updating this?? Is this obsolete?
caldat, systime(/utc,/julian)-365, month, day, year
yr = strmid(strtrim(year,2),2,2)
filename = 'TTBIPM.'

spawn, 'wget --ftp-user=anonymous -NP $ASTRO_DATA ftp://tai.bipm.org/TFG/TT\(BIPM\)/' + filename + yr
spawn, 'cp ' + getenv('ASTRO_DATA') + '/' + filename + yr + ' ' + getenv('ASTRO_DATA') + '/bipmfile'

end
