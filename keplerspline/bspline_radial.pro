;+
; NAME:
;   bspline_radial
;
; PURPOSE:
;   Calculate a Radial B-spline in the least squares sense with rejection
;
; CALLING SEQUENCE:
;   sset = bspline_radial( )
;
; INPUTS:
;   r          - independent radial coordinate
;   theta      - independent angular coordinate (in radians)
;   data       - Data values
;
; OPTIONAL KEYWORDS:
;   invvar     - Inverse variance of data; if not set, then all data points
;                 have weight of 1.0
;   ntheta     - Vector of multipole identifiers
;                 Each entry will correspond to an independent angular mode
;                 Negative elements represent sin(abs(ntheta)*theta)
;                 Positive elements represent cos(ntheta*theta)
;                 Default is [0,-2,2]
;   rbkpt      - This is analogous to bkpt keyword in bspline_iterfit
;                This sets the breakpoint positions in the radial dimension
;                 and does not include the boundary breakpoints.
;                One can choose to use bkspace or everyn to also set rbkpt.
;   _EXTRA     - Keywords for BSPLINE_BKPTS() and/or DJS_REJECT().
;
; OUTPUTS:
;   sset       - Structure describing spline fit.
;                Return 0 if too few good (INVVAR NE 0) points are passed.
;
; OPTIONAL OUTPUTS:
;   outmask    - Output mask, set =1 for good points, =0 for bad points.
;   yfit       - B-spline fit evaluated at each data point.
;
; COMMENTS:
;
;   The calling is analogous to bspline_iterfit, so hopefully you're
;   familiar with the 1d bsplines.  This fits a function in polar coordinates,
;   and r, theta, and data must be supplied.  The number of free parameters
;   is given by (Nbkpt + Nord - 1) * N_multipole.
;     Nbkpt       = n_elements(bkpt)
;     Nord        = the order of the bspline (default is 4)
;     N_multipole = n_elements(ntheta)  
;
;   See Bolton et al. 2005 for an official description of radial bsplines
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   bspline_bkpts()
;   bspline_workit()
;   bspline_action()
;   create_bsplineset()
;   djs_reject()
;
; REVISION HISTORY:
;   22-Feb-2005  Written by S. Burles & A. Bolton, MIT
;-
;------------------------------------------------------------------------------
function bspline_radial, r, theta, data, invvar=invvar, _EXTRA=EXTRA, $
                        rbkpt=rbkpt, ntheta=ntheta, yfit=yfit, outmask=outmask

      if N_elements(ntheta) EQ 0 then ntheta=[0,-2,2] 
      if NOT keyword_set(nord) then nord = 4L
      if (n_elements(maxiter) EQ 0) then maxiter = 25L
      if NOT keyword_set(upper) then upper=5.0
      if NOT keyword_set(lower) then lower=5.0


      nr = n_elements(r)
      nt = n_elements(theta)
      if (nr NE nt) OR (nt EQ 0) OR (nr EQ 0) then return, 0

      nd = n_elements(data)
      if (nr NE nd) OR (nd EQ 0) then return, 0

      ni = n_elements(invvar)
      if ni EQ 0 then invvar=r*0 + 1.0

      rs = sort(r)
      maskwork = (invvar[rs] GT 0)
      mask_r = where(maskwork, ngood)
      if ngood LT 10 then return, 0

       
      if NOT keyword_set(fullbkpt) then $
        fullbkpt = bspline_bkpts(r[rs[mask_r]], nord=nord, bkpt=rbkpt, $
                     _EXTRA=EXTRA)

      npoly = n_elements(ntheta)
      sset = create_bsplineset(fullbkpt, nord, npoly=npoly)
      sset.funcname = 'Radial B-Spline'
      ncoeff = n_elements(sset.coeff[0,*])
      sset = struct_addtags(sset, $
         { ntheta : ntheta, reff : fltarr(ncoeff), fweight: fltarr(ncoeff) } )
      sset.xmax = 2*!Pi

  
      if (ngood LT nord) then begin
         print, 'Number of good data points fewer the NORD'
         return, sset
      endif

      ;----------
      ; Iterate spline fit

      iiter = 0
      error = 0
      tempin = 0
      reduced_chi = 0
      if NOT keyword_set(buff) then buff = ''

      while (((error[0] NE 0) OR (keyword_set(qdone) EQ 0)) $
        AND iiter LE maxiter) do begin
 

        ngood = total(maskwork)
        goodbk = where(sset.bkmask NE 0)

        if (ngood LE 1 OR goodbk[0] EQ -1) then begin
           sset.coeff = 0
           iiter = maxiter + 1; End iterations
        endif else begin
 
        bf1 = bspline_action(r[rs], sset, lower=loweraction, upper=upperaction)

        dummy1 = fltarr(ncoeff)
        dummy2 = fltarr(ncoeff)
        ll = lindgen(nord)

        for j=0, n_elements(loweraction)-1 do begin
           l = loweraction[j]
           u = upperaction[j]
           if u GE l then begin
              dummy1[ll+j] = dummy1[ll+j] + total(bf1[l:u,*] * $
                       (r[rs[l:u]] # replicate(1,nord)),1)
              dummy2[ll+j] = dummy2[ll+j] + total(bf1[l:u,*],1)
           endif
        endfor

        sset.reff = dummy1 / (dummy2 + (dummy2 EQ 0)) * (dummy2 GT 0)  
        sset.fweight = dummy2 

        filler = replicate(1,nord)
        action = dblarr(nr, npoly*nord)

        for j=0, n_elements(sset.ntheta)-1 do begin
          if sset.ntheta[j] EQ 0 then action[*,lindgen(nord)*npoly+j] = bf1 $
          else begin
            if sset.ntheta[j] LT 0 then ft =sin(-1.0*sset.ntheta[j]*theta[rs]) $
            else ft = cos(sset.ntheta[j]*theta[rs])
            action[*,lindgen(nord)*npoly+j] = bf1 * (ft # filler)
          endelse
        endfor


        if total(finite(action) EQ 0) GT 0 then begin
         print, '!! Infinities in action matrix,messed up!!!'
            return, 0
        endif

        error = bspline_workit(r[rs], data[rs], invvar[rs]*maskwork, $
              action, sset, alpha=alpha, $
              lower=loweraction, upper=upperaction, $
              yfit=yfittemp, covariance=covariance)

        endelse

        iiter = iiter + 1

        if (error[0] EQ -2L) then begin
           ; All break points have been dropped.
           return, sset
        endif else if (error[0] EQ 0) then begin
           ; Iterate the fit -- next rejection iteration.
            reduced_chi = total((data[rs] - yfittemp)^2 * $
                                (invvar[rs]*maskwork))/ $
                                (ngood - 2*(n_elements(goodbk)+nord) - 1)


            relative_factor = 1.0
            if keyword_set(relative) then $
                relative_factor = sqrt(reduced_chi) > 1.0

           qdone = djs_reject(data[rs], yfittemp, invvar=invvar[rs], $
                            inmask=tempin, $
                            outmask=maskwork, upper=upper*relative_factor, $
                            lower=lower*relative_factor, _EXTRA=EXTRA)
          tempin = maskwork


           print, format='(i4, f7.3,i5)', iiter, reduced_chi, $
                                               total(maskwork EQ 0)
          print, format='(a, $)', buff

       endif
   endwhile
   ;----------
   ;  Re-sort the output arrays OUTMASK and YFIT to agree with the input data.

   outmask = byte(r * 0)
   outmask[rs] = maskwork
   
   yfit = r*0.  ; do this to have same dimensions as r
   yfit[rs] = yfittemp

   return, sset
end

