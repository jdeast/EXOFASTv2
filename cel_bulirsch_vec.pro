; More efficient "Vector" version of cel:
pro cel_bulirsch_vec,k2,kc,p,a1,a2,a3,b1,b2,b3,f1,f2,f3
; This assumes first value of a and b uses p, the rest have p=1.
ca = sqrt(k2*2.2d-16)
; Avoid undefined k2=1 case:
indx = where((k2 eq 1.0) or (kc eq 0.0))
if indx[0] ne -1 then kc[indx] = 2.22d-16
; Initialize values:
ee = kc & m=0.0*kc+1.0
pos = where(p ge 0.0,complement=neg)
pinv = 0.0*k2
if pos[0] ne -1 then begin & $
  p[pos] = sqrt(p[pos]) & pinv[pos] = 1d0/p[pos] & b1[pos] *= pinv[pos] &$
endif 
if neg[0] ne -1 then begin &$
  q=k2[neg] & g=1.0-p[neg] & f = g-k2[neg] &$
  q *= (b1[neg]-a1[neg]*p[neg]) & ginv = 1d0/g   &$
  p[neg]=sqrt(f*ginv) & a1[neg]=(a1[neg]-b1[neg])*ginv & pinv[neg] = 1d0/p[neg] &$
  b1[neg] = -q*ginv^2*pinv[neg]+a1[neg]*p[neg] &$
endif
; Compute recursion:
f1=a1
; First compute the first integral with p:
a1 += b1*pinv & g=ee*pinv & b1 += f1*g & b1 += b1 & p +=g
; Next, compute the remainder with p = 1:
p1 = 1d0 & g1=ee
f2 = a2 & f3 = a3
a2 += b2 & b2 += f2*g1 & b2 += b2
a3 += b3 & b3 += f3*g1 & b3 += b3
p1 +=g1
g = m
m += kc
iter = 0 & itmax = 50
while max(abs(g-kc) gt g*ca) and (iter lt itmax) do begin
  kc = sqrt(ee)
  kc += kc
  ee = kc*m
  f1 = a1 & f2=a2 & f3=a3
  pinv = 1d0/p
  pinv1 = 1d0/p1
  a1 += b1*pinv
  a2 += b2*pinv1
  a3 += b3*pinv1
  g = ee*pinv
  g1= ee*pinv1
  b1 += f1*g
  b2 += f2*g1
  b3 += f3*g1
  b1 += b1
  b2 += b2
  b3 += b3
  p  += g
  p1 += g1
  g  = m
  m += kc
  iter +=1
endwhile
f1 = 0.5*!dpi*(a1*m+b1)/(m*(m+p))
f2 = 0.5*!dpi*(a2*m+b2)/(m*(m+p1))
f3 = 0.5*!dpi*(a3*m+b3)/(m*(m+p1))
return
end
