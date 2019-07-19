;+
; NAME:
;   bspline_psf_action
;
; PURPOSE:
;   Calculate the action matrix for a Radial B-spline with PSF convolution
;
; CALLING SEQUENCE:
;   sset = bspline_psf_action( )
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
;   16-Aug-2005  Written by S. Burles & A. Bolton, MIT
;-
;------------------------------------------------------------------------------
function bspline_psf_action, r, theta, psf=psf, sset=sset, $
                          rbkpt=rbkpt, ntheta=ntheta, _EXTRA=EXTRA

      if N_elements(ntheta) EQ 0 then ntheta=[0,-2,2] 
      if NOT keyword_set(nord) then nord = 4L

      nr = n_elements(r)
      nt = n_elements(theta)
      if (nr NE nt) OR (nt EQ 0) OR (nr EQ 0) then return, 0

      if (size(r))[0] NE 2 then begin
        print, 'r should be 2-dimensional'
        return, 0
      endif

      rs = sort(r)

      if NOT keyword_set(fullbkpt) then begin
        if keyword_set(rbkpt) then $
          fullbkpt = bspline_bkpts(r[rs], nord=nord, bkpt=rbkpt,_EXTRA=EXTRA) $
        else begin
          print, 'Need to provide either rbkpt or fullbkpt'
          return, 0
        endelse
      endif

      npoly = n_elements(ntheta)
      sset = create_bsplineset(fullbkpt, nord, npoly=npoly)
      sset.funcname = 'Radial B-Spline PSF'
      ncoeff = n_elements(sset.coeff[0,*])
      sset = struct_addtags(sset, $
         { ntheta : ntheta, reff : fltarr(ncoeff), fweight: fltarr(ncoeff) } )
      sset.xmax = 2*!Pi

      if keyword_set(psf) then $
          sset = struct_addtags(sset, { psf : psf} )

      npar = ncoeff*npoly 
      final_action = dblarr(nr, npar)

      bw = npoly*nord
      action = dblarr(nr, bw)
      filler = replicate(1,nord)
      ll = lindgen(nord)


      bf1 = bspline_action(r[rs], sset, lower=loweraction, upper=upperaction)

      for j=0L, npoly-1L do begin
        if sset.ntheta[j] EQ 0 then action[*,ll*npoly+j] = bf1 $
        else begin
            if sset.ntheta[j] LT 0 then $
              ft =sin(-1.0*sset.ntheta[j]*theta[rs]) $
            else ft = cos(sset.ntheta[j]*theta[rs])
            action[*,ll*npoly+j] = bf1 * (ft # filler)
        endelse
      endfor


      for j=0, n_elements(loweraction)-1 do begin
        l = loweraction[j]
        u = upperaction[j]
        if u GE l then $
          final_action[rs[l:u],j*npoly:bw+j*npoly-1] = action[l:u,*]
      endfor

 
      ncol = (size(r))[1]
      nrow = (size(r))[2]
      final_action = reform(temporary(final_action),ncol,nrow, npar)

      if NOT keyword_set(psf) then return, final_action

;      final_action[*,*,0] = convolve(final_action[*,*,0], psf, ft_psf=ft_psf)
;      for i=1,npar-1 do $
;        final_action[*,*,i] = convolve(final_action[*,*,i], psf, ft_psf=ft_psf)

      final_action[*,*,0] = convol(final_action[*,*,0], psf)
      for i=1,npar-1 do $
        final_action[*,*,i] = convol(final_action[*,*,i], psf)

      return, final_action
end

