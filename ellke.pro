pro ellke,k,ek,kk
;+
; NAME:
;   ELLKE
;
; PURPOSE: 
;   Computes Hasting's polynomial approximation for the complete 
;   elliptic integral of the first (ek) and second (kk) kind. Combines
;   the calculation of both so as not to duplicate the expensive
;   calculation of alog10(1-k^2).
;
; CALLING SEQUENCE:
;    ellke, k, ek, kk
;
; INPUTS:
;
;    k - The elliptic modulus.
;
; OUTPUTS:
;
;    ek - The elliptic integral of the first kind
;    kk - The elliptic integral of the second kind
;
; MODIFICATION HISTORY
; 
;  2009/04/06 -- Written by Jason Eastman (Ohio State University)
;-

m1=1.d0-k^2
logm1 = alog(m1)

a1=0.44325141463d0
a2=0.06260601220d0
a3=0.04757383546d0
a4=0.01736506451d0
b1=0.24998368310d0
b2=0.09200180037d0
b3=0.04069697526d0
b4=0.00526449639d0
ee1=1.d0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*(-logm1)
ek = ee1+ee2
         
a0=1.38629436112d0
a1=0.09666344259d0
a2=0.03590092383d0
a3=0.03742563713d0
a4=0.01451196212d0
b0=0.5d0
b1=0.12498593597d0
b2=0.06880248576d0
b3=0.03328355346d0
b4=0.00441787012d0
ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*logm1
kk = ek1-ek2

end
