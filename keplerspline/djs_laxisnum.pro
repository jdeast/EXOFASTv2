;-----------------------------------------------------------------------
;+
; NAME:
;   djs_laxisnum
;
; PURPOSE:
;   Return a longword integer array with the specified dimensions.
;   Each element of the array is set equal to its index number in
;   the specified axis.
;
; CALLING SEQUENCE:
;   result = djs_laxisnum( dimens, [ iaxis= ] )
;
; INPUT:
;   dimens:     Vector of the dimensions for the result.
;               Only up to 3 dimensions can be specified.
;   iaxis:      Axis number to use for indexing RESULT.  The first dimension
;               is axis number 0, the second 1, etc.  Default to 0.
;
; OUTPUTS:
;   result:     Output array
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   15-Jun-2001  Written by D. Schlegel, Princeton
;-
;-----------------------------------------------------------------------
function djs_laxisnum, dimens, iaxis=iaxis

   ; Need one parameter
   if N_params() LT 1 then begin
      print, 'Syntax - result = djs_laxisnum( dimens, [iaxis= ] )'
      return, -1
   endif

   if (NOT keyword_set(iaxis)) then iaxis = 0

   ndimen = N_elements(dimens)
   naxis = long(dimens) ; convert to type LONG

   result = make_array(dimension=naxis, /long)

   case ndimen of
   1 : $
      begin
         result[*] = 0
      end
   2 : $
      begin
         case iaxis of
            0: for ii=0, naxis[0]-1 do result[ii,*] = ii
            1: for ii=0, naxis[1]-1 do result[*,ii] = ii
         endcase
      end
   3 : $
      begin
         case iaxis of
            0: for ii=0, naxis[0]-1 do result[ii,*,*] = ii
            1: for ii=0, naxis[1]-1 do result[*,ii,*] = ii
            2: for ii=0, naxis[2]-1 do result[*,*,ii] = ii
         endcase
      end
   else : $
      begin
         print, ndimen, ' dimensions not supported'
         result = -1
      end
   endcase

   return, result
end 
;-----------------------------------------------------------------------
