function create_bsplineset, fullbkpt, nord, npoly=npoly

      numbkpt = n_elements(fullbkpt)
      numcoeff = numbkpt - nord

      if (NOT keyword_set(npoly)) then begin
         sset = $
         { fullbkpt: fullbkpt             , $
           bkmask  : bytarr(numbkpt) + 1  , $
           nord    : long(nord)           , $
           coeff   : fltarr(numcoeff)     , $
           icoeff  : fltarr(numcoeff) }
      endif else begin
         sset = $
         { fullbkpt: fullbkpt            , $
           bkmask  : bytarr(numbkpt) + 1 , $
           nord    : long(nord)          , $
           xmin    : 0.0                 , $
           xmax    : 1.0                 , $
           funcname: 'legendre'          , $
           npoly   : long(npoly)         , $
           coeff   : fltarr(npoly,numcoeff) , $
           icoeff  : fltarr(npoly,numcoeff) }
      endelse

      return, sset
end

