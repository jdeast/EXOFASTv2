;+
; NAME:	    
;
;       GAUS_CONVOL
; 
; PURPOSE:
; 
;	Convolve a function with a gaussian.
; 
; CALLING SEQUENCE:
; 
;	Result=GAUS_CONVOL(X,Y,SIG)
;
; INPUTS:
;
;	X - vector of independent variables.
; 
;	Y - vector of dependent variables.
;
;	SIG - width of gaussian, in same units as X.
;
; KEYWORDS:
;
;       NSIGMA - A factor that determines how many points to use for
;                the gaussian kernel. The exact number of points used
;                depends on the value of SIG and the range of X values
;                as well, but it will be roughly equal to 2*NSIGMA*SIG.
;                Default=2.5.
;
; PROCEDURE:
; 
;       CONVOL is used to convolve y with a gaussian whose width is
;       sig, in x axis units.  The gaussian is defined as Gaussian =
;       1./sqrt(2.*!pi*sigma)*exp(-.5*(x/sigma)^2)
; 
; MODIFICATION HISTORY:
; 
;   David L. Windt, Bell Labs, July 1990.
;   
;   November 1990: added factor of 2 in gaussian argument.
;   
;   February 1997: made sure gaussian kernel has at least 2 points.
;   
;   windt@bell-labs.com
;
;   16-SEP-1998  - Introduced the keyword NSIGMA and changed the
;                  calculation of the number of points
;                  n_pts. NSIGMA=5.  Made conversion parameter DOUBLE
;                  precision and corrected the calculation of this
;                  parameter.  Note, width in xvar is only approximate
;                  (and always a little larger than +- NSIGMA).  -
;                  
;   Roger Dejus (dejus@aps.anl.gov)
;
function gaus_convol,x,y,sig,nsigma=nsigma

on_error,2
c=check_math(0,1)
if n_params() ne 3 then message,'usage: result=gaus_convol(x,y,sig)'

if n_elements(nsigma) ne 1 then nsigma=2.5      ; approximate width in
                                                ; units of sigma
nsigma2= NSIGMA*2
n      = n_elements(x)
conv   = (max(double(x))-min(double(x)))/(n -1) ; conversion, units/point
n_pts  = ceil(nsigma2*sig/conv)                 ; number of points
n_pts  = (n_pts > 2) < (n-2)                    ; restrictions on n_pts
if (n_pts/2*2) eq n_pts then n_pts = n_pts +1   ; make odd number of points
xvar = (dindgen(n_pts)/(n_pts-1)-0.5d0)*n_pts   ; approx. - NSIGMA < x < +NSIGMA

gaus=exp(-.5*(xvar/(sig/conv))^2)               ; gaussian of width sig.
return,convol(y,gaus,/center)/total(gaus)       ; do convolution.
end
