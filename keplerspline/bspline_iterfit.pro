;+
; NAME:
;   bspline_iterfit
;
; PURPOSE:
;   Calculate a B-spline in the least squares sense with rejection
;
; CALLING SEQUENCE:
;   sset = bspline_iterfit( )
;
; INPUTS:
;   xdata      - Data x values
;   ydata      - Data y values
;
; OPTIONAL KEYWORDS:
;   invvar     - Inverse variance of y; if not set, then set to be
;                consistent with the standard deviation.  This only matters
;                if rejection is being done.
;   nord       - Order for spline fit; default to 4.
;   x2         - 2nd dependent variable for 2-D spline fitting.
;   npoly      - Polynomial order to fit over 2nd variable (X2); default to 2.
;   xmin       - Normalization minimum for X2; default to MIN(XDATA).
;   xmax       - Normalization maximum for X2; default to MAX(XDATA).
;   oldset     - If set, then use values of FULLBKPT, NORD, XMIN, XMAX, NPOLY
;                from this structure.
;   funcname   - If OLDSET is not specified and this is a 2-D B-spline,
;                then the function for the second variable may be passed.
;                The default is 'legendre' in the call to CREATE_BSPLINESET().
;   maxiter    - Maximum number of rejection iterations; default to 10;
;                set to 0 to disable rejection.
;   upper      - Upper rejection threshhold; default to 5 sigma.
;   lower      - Lower rejection threshhold; default to 5 sigma.
;   requiren   - force at at least (requiren) data points between each set of bkpts.
;   _EXTRA     - Keywords for BSPLINE_BKPTS(), BSPLINE_FIT and/or DJS_REJECT().
;
; OUTPUTS:
;   sset       - Structure describing spline fit.
;                Return 0 if too few good (INVVAR NE 0) points are passed.
;
; OPTIONAL OUTPUTS:
;   outmask    - Output mask, set =1 for good points, =0 for bad points.
;   fullbkpt   - If OLDSET is not specified, then the break points are
;                chosen with a call to BSPLINE_BKPTS() which can be returned
;                with this keyword.
;   yfit       - B-spline fit evaluated at each data point.
;
; COMMENTS:
;   Data points can be masked either by setting their weights to zero
;   (INVVAR[]=0), or by using INMASK and setting bad elements to zero.
;   INMASK is passed to DJS_REJECT().
;
;   If OLDSET is used, then the output structure SSET will be a structure
;   with the same name as OLDSET.  This will allow the two structures to
;   be concatented, i.e.
;     > junk = [oldset, sset]
;
;   Although I'm not sure how to treat data points which fall outside
;   minmax(bkpt), now I will set them equal to minmax with invvar = 0
;
; EXAMPLES:
;   Construct a random function, and fit a B-spline to it without
;   any rejection:
;     IDL> x = findgen(1000)
;     IDL> y = smooth(randomu(1234,1000),10)
;     IDL> sset = bspline_iterfit(x,y,nord=3,maxiter=0,bkspace=10)
;     IDL> yfit = bspline_valu(x,sset)
;
; PROCEDURES CALLED:
;   bspline_bkpts()
;   bspline_fit()
;   create_bsplineset()
;   djs_reject()
;
; REVISION HISTORY:
;   05-Sep-2000  Written by D. Schlegel & S. Burles
;-
;------------------------------------------------------------------------------
function bspline_iterfit, xdata, ydata, invvar=invvar, nord=nord, $
 x2=x2, npoly=npoly, xmin=xmin, xmax=xmax, yfit=yfit, $
 bkpt=bkpt, oldset=oldset, maxiter=maxiter, upper=upper, lower=lower, $
 requiren=requiren, $
 outmask=outmask, fullbkpt=fullbkpt, funcname=funcname, _EXTRA=EXTRA

   if (n_params() LT 2) then begin
      print, 'Syntax -  sset = bspline_iterfit( )'
      return, 0
   endif

   ;----------
   ; Check dimensions of inputs

   nx = n_elements(xdata)
   if (n_elements(ydata) NE nx) then $
    message, 'Dimensions of XDATA and YDATA do not agree'

   if (NOT keyword_set(nord)) then nord = 4L
   if n_elements(upper) EQ 0 then upper = 5
   if n_elements(lower) EQ 0 then lower = 5

   if (keyword_set(invvar)) then begin
      if (n_elements(invvar) NE nx) then $
       message, 'Dimensions of XDATA and INVVAR do not agree'
   endif 

   if (keyword_set(x2)) then begin
      if (n_elements(x2) NE nx) then $
       message, 'Dimensions of X and X2 do not agree'
      if (NOT keyword_set(npoly)) then npoly = 2L
   endif
   if (n_elements(maxiter) EQ 0) then maxiter = 10

   yfit = 0 * ydata ; Default return values

   if (NOT keyword_set(invvar)) then begin
      var = variance(ydata, /double)
      if (var EQ 0) then var = 1
      invvar = 0.0 * ydata + 1.0/var
   endif

   if (n_elements(invvar) EQ 1) then outmask = 1B $
    else outmask = make_array(size=size(invvar), /byte) + 1B

   xsort = sort(xdata)
   maskwork = (outmask * (invvar GT 0))[xsort]
   these = where(maskwork, nthese)
 
   ;----------
   ; Determine the break points and create output structure

   if (keyword_set(oldset)) then begin
      sset = oldset
      sset.bkmask = 1
      sset.coeff = 0
      tags = tag_names(oldset)
      if ((where(tags EQ 'XMIN'))[0] NE -1 AND NOT keyword_set(x2)) then $
       message, 'X2 must be set to be consistent with OLDSET'

   endif else begin

      if (nthese EQ 0) then begin
         print, 'No valid data points'
         fullbkpt = 0
         return, 0
      endif

      if NOT keyword_set(fullbkpt) then $
       fullbkpt = bspline_bkpts(xdata[xsort[these]], nord=nord, bkpt=bkpt, $
                                 _EXTRA=EXTRA,/silent)

      sset = create_bsplineset(fullbkpt, nord, npoly=npoly) 

      if (nthese LT nord) then begin
         print, 'Number of good data points fewer the NORD'
         return, sset
      endif

      ;----------
      ; Condition the X2 dependent variable by the XMIN, XMAX values.
      ; This will typically put X2NORM in the domain [-1,1].

      if keyword_set(x2) then begin
         if n_elements(xmin) NE 1 then xmin = min(x2)
         if n_elements(xmax) NE 1 then xmax = max(x2)
         if (xmin EQ xmax) then xmax = xmin + 1
         sset.xmin = xmin
         sset.xmax = xmax

         if keyword_set(funcname) then sset.funcname=funcname
      endif

   endelse


   ;----------
   ; It's okay now if the data fall outside breakpoint regions, the
   ; fit is just set to zero outside.

   ;----------
   ; Sort the data so that X is in ascending order.

   xwork = xdata[xsort]
   ywork = ydata[xsort]
   invwork = invvar[xsort]
   if (keyword_set(x2)) then x2work = x2[xsort]


   ;----------
   ; Iterate spline fit

   iiter = 0
   error = 0

   while (((error[0] NE 0) OR (keyword_set(qdone) EQ 0)) $
    AND iiter LE maxiter) do begin

      ngood = total(maskwork)
      goodbk = where(sset.bkmask NE 0, ngb)

      if (ngood LE 1 OR goodbk[0] EQ -1) then begin
         sset.coeff = 0
         iiter = maxiter + 1; End iterations
      endif else begin

        if keyword_set(requiren) then begin

          ; Locate where there are two break points in a row with no good
          ; data points in between, and drop (mask) one of those break points.
          ; The first break point is kept.
           i = 0L
           while(xwork[i] LT sset.fullbkpt[goodbk[nord]] AND i LT nx-1) do i = i+ 1

           ct = 0L
           for ileft=nord, ngb-nord do begin
             while(xwork[i] GE sset.fullbkpt[goodbk[ileft]] AND $
                   xwork[i] LT sset.fullbkpt[goodbk[ileft+1]] AND i LT nx-1) do begin
                ct = ct + (invwork[i] * maskwork[i] GT 0)
                i=i+1
             endwhile
             if ct GE requiren then ct = 0L else sset.bkmask[goodbk[ileft]] = 0
           endfor

         endif

        ; Do the fit.  Return values for ERROR are as follows:
        ;    0 if fit is good
        ;   -1 if all break points are masked
        ;   -2 if everything is screwed
        error = bspline_fit(xwork, ywork, invwork*maskwork, sset, $
         x2=x2work, yfit=yfit, nord=nord, _EXTRA=EXTRA)
      endelse

      iiter = iiter + 1

      inmask = maskwork

      if (error[0] EQ -2L) then begin
         ; All break points have been dropped.
         return, sset
      endif else if (error[0] EQ 0) then begin
         ; Iterate the fit -- next rejection iteration.
         qdone = djs_reject(ywork, yfit, invvar=invwork, inmask=inmask, $
          outmask=maskwork, upper=upper, lower=lower, _EXTRA=EXTRA)
      endif

   endwhile

   ;----------
   ; Re-sort the output arrays OUTMASK and YFIT to agree with the input data.

   outmask[xsort] = maskwork

   temp = yfit
   yfit[xsort] = temp

   return, sset
end
;------------------------------------------------------------------------------
