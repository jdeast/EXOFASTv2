pro exofast_occultquad_cel,z0,u1,u2,p0,muo1,mu0,d=d
;+
; NAME:
;   EXOFAST_OCCULTQUAD_CEL
;
; PURPOSE: 
;   This routine computes the lightcurve for occultation of a
;   quadratically limb-darkened source without microlensing.  Please
;   cite Mandel & Agol (2002) and Eastman et al., (2013) if you make use
;   of this routine in your research.  Please report errors or bugs to
;   agol@astro.washington.edu and jason.eastman@cfa.harvard.edu
;
; Limb darkening has the form:
;  I(r)=[1-u1*(1-sqrt(1-(r/rs)^2))-u2*(1-sqrt(1-(r/rs)^2))^2]/(1-u1/3-u2/6)/pi
; 
; CALLING SEQUENCE:
;    exofast_occultquad_cel, z0, u1, u2, p0, muo1, mu0, d=d
;
; INPUTS:
;
;    z0 - impact parameter in units of rs
;    u1 - linear    limb-darkening coefficient (gamma_1 in paper)
;    u2 - quadratic limb-darkening coefficient (gamma_2 in paper)
;    p0 - occulting star size in units of rs
;
; OUTPUTS:
;
;    muo1 - fraction of flux at each z0 for a limb-darkened source
;
; OPTIONAL OUTPUTS:
;
;    mu0  - fraction of flux at each z0 for a uniform source
;    d    - The coefficients required to analytically calculate the
;           limb darkening parameters (see Eastman et al, 2013). For
;           backward compatibility, u1 and u2 are required, but not
;           necessary if d is used.
;
; EXAMPLES:  
;
;; Calculate the same geometric transit with two different sets of
;; limb darkening coefficients
; p = 0.1d0
; b = 0.5d0
; x = (dindgen(300)/299d0 - 0.5d0)*2d0
; z = sqrt(x^2 + b^2)
; u1 = 0.25d0
; u1 = 0.75d0
; exofast_occultquad, z, u1, u2, p, muo1, muo, d=d
;
; MODIFICATION HISTORY
; 
;  2002 -- Eric Agol 2002
;
;  2009/04/06 -- Eric Agol (University of Washington) and 
;                Jason Eastman (Ohio State University)
;    fixed bugs for p > 1 (e.g., secondary eclipses)
;    used double precision constants throughout
;    replaced rj,rc,rf with ellpic_bulirsch
;      faster, more robust integration for specialized case of
;      elliptic integral of the third kind
;    vectorized
;    more efficient case handling
;    combined ellk and ellec into ellke for speed
;    200x speed improvement over previous IDL code in typical case
;    allow negative values of p (anti-transit) to avoid Lucy-Sweeney like bias
;  2018/12/12 -- Now uses more stable integrations, 25% faster
;-

nz=n_elements(z0)
lambdad=dblarr(nz)
etad=dblarr(nz)
lambdae=dblarr(nz)

;; ensure double precision for equalities
;; don't modify inputs
p = double(abs(p0))
z = double(z0)

x1=(p-z)^2
x2=(p+z)^2
x3=p^2-z^2

;; Case 1 - the star is unocculted:
;; only consider points with z lt 1+p
;; exit if there is no planet (p <= 0)
notusedyet = where(z lt (1.d0 + p) and p gt 0.d0)
if notusedyet[0] eq -1 then goto, final

;; Case 11 - the  source is completely occulted:
if p ge 1.d0 then begin
    occulted = where(z[notusedyet] le p-1.d0,complement=notused2)
    if occulted[0] ne -1 then begin
        ndxuse = notusedyet[occulted]
        etad[ndxuse] = 0.5d0 ;; corrected typo in paper
        lambdae[ndxuse] = 1.d0
        ;; lambdad = 0 already
    endif
    if notused2[0] eq -1 then goto, final
    notusedyet = notusedyet[notused2]
endif

;; Case 2, 7, 8 - ingress/egress (uniform disk only)
inegressuni = where(z[notusedyet] ge abs(1.d0-p) and z[notusedyet] lt 1.d0+p)
if inegressuni[0] ne -1 then begin
    ndxuse = notusedyet[inegressuni]
    sqarea_triangle,z[ndxuse],p,sqarea
    kite_area2 = sqrt(sqarea)
    kap1 = atan(kite_area2,(1d0-p)*(p+1d0)+z[ndxuse]^2)
    kap0 = atan(kite_area2,(p-1d0)*(p+1d0)+z[ndxuse]^2)
    lambdae[ndxuse] = (p^2*kap0+kap1 - 0.5d0*kite_area2)/!dpi
    ;; eta_1
    etad[ndxuse] = 1.d0/(2.d0*!dpi)*(kap1+p^2*(p^2+2.d0*z[ndxuse]^2)*$
                      kap0-0.25d0*(1.d0+5.d0*p^2+z[ndxuse]^2)*kite_area2)
endif

;; Case 5, 6, 7 - the edge of planet lies at origin of star
ocltor = where(z[notusedyet] eq p, complement=notused3)
if ocltor[0] ne -1 then begin
    ndxuse = notusedyet[ocltor]
    if p lt 0.5d0 then begin
        ;; Case 5
        q=2.d0*p  ; corrected typo in paper (2k -> 2p)
        ellke, q, Ek, Kk
        ;; lambda_4
        lambdad[ndxuse] = 1.d0/3.d0+2.d0/9.d0/!dpi*$
          (4.d0*(2.d0*p^2-1.d0)*Ek+(1.d0-4.d0*p^2)*Kk)
        ;; eta_2
        etad[ndxuse] = 3.d0*p^4/2.d0; p^2/2.d0*(p^2+2.d0*z[ndxuse]^2)        
        lambdae[ndxuse] = p^2 ;; uniform disk
    endif else if p gt 0.5d0 then begin
        ;; Case 7
        q=0.5d0/p ; corrected typo in paper (1/2k -> 1/2p)
        ellke, q, Ek, Kk
        ;; lambda_3
        lambdad[ndxuse] = 1.d0/3.d0+16.d0*p/9.d0/!dpi*(2.d0*p^2-1.d0)*$
          Ek-(32.d0*p^4-20.d0*p^2+3.d0)/9.d0/!dpi/p*Kk
        ;; etad = eta_1 already
    endif else begin
        ;; Case 6
        lambdad[ndxuse] = 1.d0/3.d0-4.d0/!dpi/9.d0
        etad[ndxuse] = 3.d0/32.d0
    endelse
endif
if notused3[0] eq -1 then goto, final
notusedyet = notusedyet[notused3]

;; Case 3, 4, 9, 10 - planet completely inside star
inside = where(p lt 1.d0 and z[notusedyet] le (1.d0-p), complement=notused5)
if inside[0] ne -1 then begin
    ndxuse = notusedyet[inside]

    ;; eta_2
    etad[ndxuse] = p^2/2.d0*(p^2+2.d0*z[ndxuse]^2)

    ;; uniform disk
    lambdae[ndxuse] = p^2

    ;; Case 4 - edge of planet hits edge of star
    edge = where(z[ndxuse] eq 1.d0-p, complement=notused6)
    if edge[0] ne -1 then begin
        ;; lambda_5
        lambdad[ndxuse[edge]] = 2.d0/3.d0/!dpi*acos(1.d0-2.d0*p)-$
          4.d0/9.d0/!dpi*sqrt(p*(1.d0-p))*$
          (3.d0+2.d0*p-8.d0*p^2)-2.d0/3.d0*(p gt 0.5d0)
        if notused6[0] eq -1 then goto, final
        ndxuse = ndxuse[notused6]
    endif

    ;; Case 10 - origin of planet hits origin of star
    origin = where(z[ndxuse] eq 0, complement=notused7)
    if origin[0] ne -1 then begin
        ;; lambda_6
        lambdad[ndxuse[origin]] = -2.d0/3.d0*(1.d0-p^2)^1.5d0
        if notused7[0] eq -1 then goto, final
        ndxuse = ndxuse[notused7]
    endif

    onembpr2 = (1-z[ndxuse]-p)*(1+z[ndxuse]+p) & onembmr2=(p-z[ndxuse]+1)*(1-p+z[ndxuse]) & fourbr = 4*z[ndxuse]*p & fourbrinv = 1d0/fourbr
    k2 = onembpr2*fourbrinv+1
    onembmr2inv = 1d0/onembmr2 & k2inv = 1d0/k2 & kc2 =onembpr2*onembmr2inv & kc = sqrt(kc2)
    bmrdbpr = (z[ndxuse]-p)/(z[ndxuse]+p)
    mu = 3*bmrdbpr*onembmr2inv
    p_bulirsch = bmrdbpr^2*onembpr2*onembmr2inv
    cel_bulirsch_vec,k2inv,kc,p_bulirsch,1+mu,1d0,1d0,p_bulirsch+mu,kc2,0d0,Piofk,Eofk,Em1mKdm
    lambdad[ndxuse] = 2*sqrt(onembmr2)*(onembpr2*Piofk -(4-7*p^2-z[ndxuse]^2)*Eofk)/(9d0*!dpi)
 endif


;; Case 2, Case 8 - ingress/egress (with limb darkening)
inegress = notused5
if inegress[0] ne - 1 then begin
    ndxuse = notusedyet[inegress]
    onembpr2 = (1-z[ndxuse]-p)*(1+z[ndxuse]+p) & onembmr2=(p-z[ndxuse]+1)*(1-p+z[ndxuse]) & fourbr = 4*z[ndxuse]*p & fourbrinv = 1d0/fourbr
    k2 = onembpr2*fourbrinv+1
    kc2 = -onembpr2*fourbrinv & kc = sqrt(kc2)
    cel_bulirsch_vec,k2,kc,(z[ndxuse]-p)^2*kc2,0d0,1d0,1d0,3*kc2*(z[ndxuse]-p)*(z[ndxuse]+p),kc2,0d0,Piofk,Eofk,Em1mKdm 
    lambdad[ndxuse]  = onembmr2*(Piofk+ (-3+6*p^2+2*z[ndxuse]*p)*Em1mKdm-fourbr*Eofk)/(9*!dpi*sqrt(z[ndxuse]*p))
endif

final:
omega=1.d0-u1/3.d0-u2/6.d0

;; avoid Lutz-Kelker bias (negative values of p0 allowed)
if p0 gt 0 then begin
    ;; limb darkened flux
    muo1 =1.d0-((1.d0-u1-2.d0*u2)*lambdae+(u1+2.d0*u2)*$
                (lambdad+2.d0/3.d0*(p gt z))+u2*etad)/omega
    ;; uniform disk
    mu0=1.d0-lambdae

    ;; coeffs for quadratic limb darkening fit
    if arg_present(d) then $
      d = transpose([[1.d0-lambdae],$
                     [2d0/3d0*(lambdae - (p gt z)) - lambdad],$
                     [lambdae/2d0 - etad]])
endif else begin
    ;; limb darkened flux
    muo1 =1.d0+((1.d0-u1-2.d0*u2)*lambdae+(u1+2.d0*u2)*$
                (lambdad+2.d0/3.d0*(p gt z))+u2*etad)/omega
    ;; uniform disk
    mu0=1.d0+lambdae

    ;; coeffs for quadratic limb darkening fit
    if arg_present(d) then $
      d = transpose([[1.d0+lambdae],$
                     [2d0/3d0*((p gt z) - lambdae) + lambdad],$
                     [etad - lambdae/2d0]])
endelse

end
