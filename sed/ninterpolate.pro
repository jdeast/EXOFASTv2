;+
; NAME:
;
;    NINTERPOLATE
;
; PURPOSE:
; 
;    Perform n-dimensional linear interpolation on arrays of
;    arbitrary dimensionality (ala the bilinear and trilinear
;    interpolation of INTERPOLATE).
;
; CALLING SEQUENCE:
;
;    value=ninterpolate(data,point)
;
; INPUTS:
;
;    DATA: The n-dimensional array to interpolate.
; 
;    POINT: A vector of length n specifying the (single) point for
;       which the interpolated value is desired.
;
; OUTPUTS:
;
;    The interpolated value at the specified point.
;
; RESTRICTIONS:
;
;    Unlike the INTERPOLATE command, only one point can be
;    interpolated at a time.  Ill-behaved on points outside the array
;    boundary.  Requires IDL >v5.6.
;
; EXAMPLE:
;
;  r=randomu(sd,40,50,60,70)
;  v=ninterpolate(r,[22.4,18.6,52.2,60.4])
;
; MODIFICATION HISTORY:
;
;    Mon Jul 21 12:33:30 2003, J.D. Smith <jdsmith@as.arizona.edu>
;		Written.
;==========================================================================
; Copyright (C) 2003, J.D. Smith
;-
function ninterpolate,data,point
  n=n_elements(point) & two_n=2L^n & reb=[n,two_n]
  if n ne size(data,/N_DIMENSIONS) then $
     message,'Point must specify 1 coordinate for each dimension'
  if n eq 1 then return, interpolate(data,point[0])
  inds=(rebin(indgen(1,two_n),reb) AND rebin([ishft(1,indgen(n))],reb)) gt 0
  base=floor(point) & f=float(point)-base
  f=[(1.-f),f]
  multi_to_single=[1,product((size(data,/DIMENSIONS))[0:n-2],/CUMULATIVE)]
  data_ind=rebin(base,reb)+inds
  data_ind=long(total(data_ind*rebin(multi_to_single,reb),1))
  return, total(product(f[rebin(indgen(n),reb)+inds*n],1)*data[data_ind])
end 

