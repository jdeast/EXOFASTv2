;; (x,y,z) = (0,0,0) is the center of the star
;; x -- An NPlanets x NTimes array of X positions of each planet, in stellar radii
;; y -- An NPlanets x NTimes array of Y positions of each planet, in stellar radii
;; z -- An NPlanets x NTimes array of Z positions of each planet, in stellar radii
;; p -- An NPlanets array of the radius of each planet, in stellar radii
;; u1 -- linear limb darkening parameter
;; u2 -- quadratic limb darkening parameter

function mutual_eclipse, x, y, z, p, u1, u2

sz = size(x)
nplanets = sz[1]
ntimes = sz[2]
flux = dblarr(ntimes) + 1d0

for i=0, nplanets-1 do begin
   

endfor

return, flux

end
