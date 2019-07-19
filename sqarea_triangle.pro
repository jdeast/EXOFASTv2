pro sqarea_triangle,z0,p0,sqarea
; Computes sixteen times the square of the area of a triangle
; with sides 1, z0 and p0 using Kahan method (Goldberg 1991).
sqarea = dblarr(n_elements(z0))
; There are six cases to consider
pz1 = where((p0 le z0) and (z0 le 1))
if pz1[0] ne -1 then sqarea[pz1] = (p0+(z0[pz1]+1d0))*(1d0-(p0-z0[pz1]))*(1d0+(p0-z0[pz1]))*(p0+(z0[pz1]-1d0))
zp1 = where((z0 le p0) and (p0 le 1))
if zp1[0] ne -1 then sqarea[zp1] = (z0[zp1]+(p0+1d0))*(1d0-(z0[zp1]-p0))*(1d0+(z0[zp1]-p0))*(z0[zp1]+(p0-1d0))
p1z = where((p0 le 1d0) and (1d0 le z0))
if p1z[0] ne -1 then sqarea[p1z] = (p0+(1d0+z0[p1z]))*(z0[p1z]-(p0-1d0))*(z0[p1z]+(p0-1d0))*(p0+(1d0-z0[p1z]))
z1p = where((z0 le 1d0) and (1d0 le p0))
if z1p[0] ne -1 then sqarea[z1p] = (z0[z1p]+(1d0+p0))*(p0-(z0[z1p]-1d0))*(p0+(z0[z1p]-1d0))*(z0[z1p]+(1d0-p0))
onepz = where((1d0 le p0) and (p0 le z0))
if onepz[0] ne -1 then sqarea[onepz] = (1d0+(p0+z0[onepz]))*(z0[onepz]-(1d0-p0))*(z0[onepz]+(1d0-p0))*(1d0+(p0-z0[onepz]))
onezp = where((1d0 le z0) and (z0 le p0))
if onezp[0] ne -1 then sqarea[onezp] = (1d0+(z0[onezp]+p0))*(p0-(1d0-z0[onezp]))*(p0+(1d0-z0[onezp]))*(1d0+(z0[onezp]-p0))
return
end
