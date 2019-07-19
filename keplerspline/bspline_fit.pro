;+
; NAME:
;   bspline_fit
;
; PURPOSE:
;   Calculate a B-spline in the least-squares sense 
;     based on two variables: x which is sorted and spans a large range
;				  where bkpts are required
;  		and 	      y which can be described with a low order
;				  polynomial	
;
; CALLING SEQUENCE:
;   error_code = bspline_fit(x, y, invvar, sset, $
;    fullbkpt=, nord=, x2=, npoly=, yfit=)
;
; INPUTS:
;   x          - Independent variable
;   y          - Dependent variable
;   invvar     - Inverse variance of Y
;   sset       - Structure to be returned with all fit parameters
;                (if not set, then it is created)
;
; OPTIONAL KEYWORDS:
;   fullbkpt   - Pass fullbkpt to seed structure; if not set, then
;                this is generated with CREATE_BSPLINESET()
;   nord       - Order of b-splines; default to 4 (cubic)
;   x2         - Orthogonal dependent variable
;   npoly      - Order of x2 polynomial fit; default to the value in the
;                SSET structure, or to 1.
;
; OUTPUTS:
;   error_code - Non-negative numbers indicate ill-conditioned bkpts
;                 0 is good
;                -1 is dropped breakpoints, try again
;                -2 is failure, should abort
;   sset       - Structure with all fit parameters
;
; OPTIONAL OUTPUTS:
;   yfit       - Evaluation of b-spline at X (and X2)
;
; COMMENTS:
;   This code replaces efcmn and efc2d calls in the slatec library.
;
; BUGS:
;   Do we need to sort X for this routine???
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   bspline_action()
;   bspline_maskpoints()
;   bspline_valu()
;   cholesky_band()
;   cholesky_solve()
;   create_bsplineset()
;
; REVISION HISTORY:
;   20-Aug-2000 Written by Scott Burles, FNAL
;-
;------------------------------------------------------------------------------
function bspline_fit, xdata, ydata, invvar, sset, fullbkpt=fullbkpt, $
 x2=x2, npoly=npoly, nord=nord, yfit=yfit

   if NOT keyword_set(nord) then nord = 4L

   if (size(sset, /tname) NE 'STRUCT') then $
    sset = create_bsplineset(fullbkpt, nord, npoly=npoly) 

   if ((where(tag_names(sset) EQ 'NPOLY'))[0] NE -1) then npoly = sset.npoly
   if (NOT keyword_set(npoly)) then npoly = 1L
    
   goodbk = where(sset.bkmask[nord:*] NE 0, nbkpt)

   nord = sset.nord

   if (nbkpt LT nord) then begin
      if (arg_present(yfit)) then yfit = fltarr(n_elements(ydata))
      return, -2L
   endif

   nn = nbkpt 
   nfull = nn * npoly
   bw = npoly * nord   ; this is the bandwidth

   ;  The next line is REQUIRED to fill a1

   a1 = bspline_action(xdata, sset, x2=x2, lower=lower, upper=upper)

   a2 = a1 * (invvar # replicate(1,bw))

   alpha = dblarr(bw,nfull+bw)
   beta = dblarr(nfull+bw)

   bi = lindgen(bw)
   bo = lindgen(bw)
   for i=1L, bw-1 do bi = [bi, lindgen(bw-i)+(bw+1)*i]
   for i=1L, bw-1 do bo = [bo, lindgen(bw-i)+bw*i]

   for i=0L, nn-nord do begin

      itop = i * npoly
      ibottom = (itop < nfull + bw) - 1
       
      ict = upper[i] - lower[i] + 1
  
      if (ict GT 0) then begin

         work = a2[lower[i]:upper[i],*] ## transpose(a1[lower[i]:upper[i],*])
         wb =  ydata[lower[i]:upper[i]] # a2[lower[i]:upper[i],*] 

         alpha[bo+itop*bw] = alpha[bo+itop*bw] + work[bi]
         beta[itop:ibottom] = beta[itop:ibottom] + wb

      endif
   endfor

;-----------------------------------------------------------------------------
;
;     This was used to test the full matrix inversion, just to compare with
;       band matrix inversion
;
;
;    afull = fltarr(nfull+bw, nx)
;    afull[0:bw-1,*] = transpose(a1)
;    for i=0,nx-1 do afull[*,i] = shift(afull[*,i],(indx[i]-min(indx))*npoly)
;    afull = afull[0:nfull-1, *]
;    jj = afull # transpose(afull)
;    choldc, jj, p
;    res = cholsol(jj, p, transpose(ydata # transpose(afull)))
;
;-----------------------------------------------------------------------------

   ; Drop break points where minimal influence is located

   min_influence = 1.0e-10 * total(invvar) / nfull

   ; This call to cholesky_band operates on alpha and changes contents
 
 a = alpha
   errb = cholesky_band(alpha, mininf=min_influence) 

   if (errb[0] NE -1) then begin 
      if (arg_present(yfit)) then $
       yfit = bspline_valu(xdata, sset, x2=x2, action=a1, upper=upper, $
        lower=lower)
      return, bspline_maskpoints(sset, errb, npoly)
   endif
 
   ; this changes beta to contain the solution

   errs = cholesky_solve(alpha, beta)   
   if (errs[0] NE -1) then begin
      if (arg_present(yfit)) then $
       yfit = bspline_valu(xdata, sset, x2=x2, action=a1, upper=upper, $
        lower=lower)
      return, bspline_maskpoints(sset, errs, npoly)
   endif

   sc = size(sset.coeff)
   if (sc[0] EQ 2) then begin
      sset.icoeff[*,goodbk] = reform(alpha[0,lindgen(nfull)],npoly,nn)
      sset.coeff[*,goodbk] = reform(beta[lindgen(nfull)], npoly, nn)
   endif else begin
      sset.icoeff[goodbk] = alpha[0,lindgen(nfull)]
      sset.coeff[goodbk] = beta[lindgen(nfull)]
   endelse

   if (arg_present(yfit)) then $
    yfit = bspline_valu(xdata, sset, x2=x2, action=a1, upper=upper, $
     lower=lower)

   return, 0L
end
;------------------------------------------------------------------------------
