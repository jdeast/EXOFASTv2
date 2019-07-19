pro exofast_occultquad,z0,u1,u2,p0,muo1,mu0,d=d
;+
; NAME:
;   EXOFAST_OCCULTQUAD
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
;    occultquad, z0, u1, u2, p0, muo1, mu0, d=d
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
;    allow negative values of p (anti-transit) to avoid Lucy-Sweeney
;    like bias
; 2018/12/12 -- Catch rare bug that hangs the code when planet and
;               stellar limb are within 10^-13 of one another.
;-

nz=n_elements(z0)
lambdad=dblarr(nz)
etad=dblarr(nz)
lambdae=dblarr(nz)

;; ensure double precision for equalities
;; don't modify inputs
p = double(abs(p0))
z = double(z0)

;; tolerance for double precision equalities
;; special case integrations
tol = 1d-14
bad = where(abs(p-z) lt tol)
if bad[0] ne -1 then z[bad] = p
bad = where(abs((p-1.d0)-z) lt tol)
if bad[0] ne -1 then z[bad] = p-1.d0
bad = where(abs((1.d0-p)-z) lt tol)
if bad[0] ne -1 then z[bad] = 1.d0-p
bad = where(z lt tol)
if bad[0] ne -1 then z[bad] = 0.d0

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
inegressuni = where(z[notusedyet] ge abs(1.d0-p) and z[notusedyet] lt 1.d0+p, complement=notused6)
if inegressuni[0] ne -1 then begin
    ndxuse = notusedyet[inegressuni]

    kap1 = acos(-1.d0 > ((1.d0-p^2+z[ndxuse]^2)/2.d0/z[ndxuse]) < 1.d0)
    kap0 = acos(-1.d0 > ((p^2+z[ndxuse]^2-1.d0)/2.d0/p/z[ndxuse]) < 1.d0)

    lambdae[ndxuse] = (p^2*kap0+kap1 - 0.5d0*sqrt(4.d0*z[ndxuse]^2-$
                              (1.d0+z[ndxuse]^2-p^2)^2  > 0.d0))/!dpi
    ;; eta_1
    etad[ndxuse] = 1.d0/2.d0/!dpi*(kap1+p^2*(p^2+2.d0*z[ndxuse]^2)*$
                      kap0-(1.d0+5.d0*p^2+z[ndxuse]^2)/4.d0*$
                      sqrt((1.d0-x1[ndxuse])*(x2[ndxuse]-1.d0)))

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

;; Case 2, Case 8 - ingress/egress (with limb darkening)
q=sqrt((1.d0-x1[notusedyet])/(x2[notusedyet]-x1[notusedyet]))
inegress = where((z[notusedyet] gt 0.5d0+abs(p-0.5d0) and z[notusedyet] lt 1.d0+p) or $
                 (p gt 0.5d0 and z[notusedyet] gt abs(1.d0-p) and z[notusedyet] lt p) and (q ne 1d0), complement=notused4)
if inegress[0] ne - 1 then begin
    ndxuse = notusedyet[inegress]
    q=sqrt((1.d0-x1[ndxuse])/(x2[ndxuse]-x1[ndxuse]))
    ellke, q, Ek, Kk
    n=1.d0/x1[ndxuse]-1.d0
    ;; lambda_1: 
    lambdad[ndxuse]=2.d0/9.d0/!dpi/sqrt(x2[ndxuse]-x1[ndxuse])*$
      (((1.d0-x2[ndxuse])*(2.d0*x2[ndxuse]+x1[ndxuse]-3.d0)-3.d0*$
        x3[ndxuse]*(x2[ndxuse]-2.d0))*Kk+(x2[ndxuse]-x1[ndxuse])*$
       (z[ndxuse]^2+7.d0*p^2-4.d0)*Ek-3.d0*x3[ndxuse]/x1[ndxuse]*$
       ellpic_bulirsch(n,q))

endif
if notused4[0] eq -1 then goto, final
notusedyet = notusedyet[notused4]

;; Case 3, 4, 9, 10 - planet completely inside star
inside = where(p lt 1.d0 and z[notusedyet] le (1.d0-p), complement=notused5)
if inside[0] ne -1 then begin
    ndxuse = notusedyet[inside]

    ;; eta_2
    etad[ndxuse] = p^2/2.d0*(p^2+2.d0*z[ndxuse]^2)

    ;; uniform disk
    lambdae[ndxuse] = p^2

    ;; Case 4 - edge of planet hits edge of star
    q = sqrt((x2[ndxuse]-x1[ndxuse])/(1.d0-x1[ndxuse]))
    edge = where(z[ndxuse] eq 1.d0-p or q eq 1d0, complement=notused6)
    if edge[0] ne -1 then begin
        ;; lambda_5
        lambdad[ndxuse[edge]] = 2.d0/3.d0/!dpi*acos(1.d0-2.d0*p)-$
          4.d0/9.d0/!dpi*sqrt(p*(1.d0-p))*$
          (3.d0+2.d0*p-8.d0*p^2)-2.d0/3.d0*(p gt 0.5d0)
        if notused6[0] eq -1 then goto, final
        ndxuse = ndxuse[notused6]
    endif

    ;; Case 10 - origin of planet hits origin of star
    origin = where(z[ndxuse] eq 0d0, complement=notused7)
    if origin[0] ne -1 then begin
        ;; lambda_6
        lambdad[ndxuse[origin]] = -2.d0/3.d0*(1.d0-p^2)^1.5d0
        if notused7[0] eq -1 then goto, final
        ndxuse = ndxuse[notused7]
    endif
   
    q=sqrt((x2[ndxuse]-x1[ndxuse])/(1.d0-x1[ndxuse]))
    n=x2[ndxuse]/x1[ndxuse]-1.d0
    ellke, q, Ek, Kk

    ;; Case 3, Case 9 - anywhere in between
    ;; lambda_2
    lambdad[ndxuse] = 2.d0/9.d0/!dpi/sqrt(1.d0-x1[ndxuse])*$
      ((1.d0-5.d0*z[ndxuse]^2+p^2+x3[ndxuse]^2)*Kk+(1.d0-x1[ndxuse])*$
       (z[ndxuse]^2+7.d0*p^2-4.d0)*Ek-3.d0*x3[ndxuse]/x1[ndxuse]*$
       ellpic_bulirsch(n,q))
 endif
;; if there are still unused elements, there's a bug in the code
;; (please report it)
if notused5[0] ne -1 then begin
    print, 'Undefined case -- please report to jason.eastman@cfa.harvard.edu'
    message, "ERROR: the following values of z didn't fit into a case:" +$
      strtrim(z[notusedyet[notused5]],2)
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
