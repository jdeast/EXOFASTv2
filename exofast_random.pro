;+
; NAME:
;        EXOFAST_RANDOM
;
; PURPOSE: 
;
;        EXOFAST_RANDOM is a wrapper for pg_ran, the IDL translation
;        of ran from "Numerical Recipes Third Edition", page 342 by
;        Paolo C. Grigis, to make it a drop in replacement of RANDOMU
;        for normal or uniform deviates.
;        http://hea-www.cfa.harvard.edu/~pgrigis/idl_stuff/pg_ran.pro
;
;        While statistically, it is a far superior routine to IDL's
;        built-in, ran1 based, RANDOMU, it is also substantially (120x)
;        slower, as random number generators require serial
;        operations, which IDL does poorly.
;
;        If sequences larger than ~10^7 are required, IDL's built-in
;        routine is insufficient.
;
;        Call with a NEGATIVE or UNDEFINED SEED TO INITIALIZE
;        thereafter, do not alter seed between successive deviates in
;        a sequence.
;
; CALLING SEQUENCE:
;
;        Result = EXOFAST_RANDOM(Seed [, D1, ..., Dn,/normal])
;
; INPUTS:
;         Seed   - A named variable containing the seed value for
;                  random number generation. Seed is updated
;                  once for each random number generated. The initial
;                  value of Seed should be set to different NEGATIVE
;                  values in order to obtain different random
;                  sequences. An undefined variable name can be
;                  specified and it will begin the seed with
;                  -systime(/seconds). It is an error not to specify a
;                  seed.
; OPTIONAL INPUTS:
;        Di      - The dimensions of the result. The dimension
;                  parameters can be any scalar expression. Up to eight
;                  dimensions can be specified. If no dimensions are
;                  specified, RAN2 returns a scalar result
; OPTIONAL KEYWORDS:
;        NORMAL  - Set this keyword to generate random normal deviates
;                  which uses an IDL translation of the Numerical
;                  Recipes routine gasdev with the uniform deviates
;                  generated from pg_ran.
;
; OUTPUTS:
;        This function returns uniform random deviates between 0.0 and
;        1.0 (exclusive of the endpoint values) or normal deviates if
;        /NORMAL is set, with array dimensions given by the input
;        parameters, D1, ..., Dn.
;
; MODIFICATION HISTORY:
;        Written by: Jason Eastman, OSU 3/2/2009 
;          Based on Numerical Recipes Third Edition page 342
;-

function exofast_random, seed, D1, D2, D3, D4, D5, D6, D7, D8, NORMAL=normal

Ds = 1L
NP   = N_PARAMS()

SWITCH NP of
;    9    : Ds*=long(D8)
    8    : Ds*=long(D7)
    7    : Ds*=long(D6)
    6    : Ds*=long(D5)
    5    : Ds*=long(D4)
    4    : Ds*=long(D3)
    3    : Ds*=long(D2)
    2    : Ds*=long(D1)
    1    : break
    else : message,'Must have 1-8 parameters, Seed [, D1,... D7]'
ENDSWITCH

;; initialize the seed
if n_elements(seed) eq 0 then seed = systime(/seconds)

if keyword_set(normal) then begin
    ;; a translation of the
    ;; gasdev algorithm from Numerical Recipes
    nran = Ceil(Ds/2.)
    ndx = lindgen(nran)
    v1 = dblarr(nran)
    v2 = dblarr(nran)
    rsq = dblarr(nran)
    
    while nran ne 0 do begin
        v1[ndx]=2.d0*pg_ran(seed,nran,/double)-1.d0
        v2[ndx]=2.d0*pg_ran(seed,nran,/double)-1.d0
        rsq[ndx]=v1[ndx]^2 + v2[ndx]^2
        ndx = where(rsq gt 1 or rsq eq 0,nran)
    endwhile
    
    fac=sqrt(-2.d0*alog(rsq)/rsq)
    r = ([v1*fac,v2*fac])(0:Ds-1) 
endif else begin
    r = pg_ran(seed,Ds,/double)
endelse

;; return the desired dimensions
CASE NP of
    1    : return, r[0]
    2    : return, r
    3    : return, reform(r,D1,D2)
    4    : return, reform(r,D1,D2,D3)
    5    : return, reform(r,D1,D2,D3,D4)
    6    : return, reform(r,D1,D2,D3,D4,D5)
    7    : return, reform(r,D1,D2,D3,D4,D5,D6)
    8    : return, reform(r,D1,D2,D3,D4,D5,D6,D7)
;    9    : return, reform(r,D1,D2,D3,D4,D5,D6,D7,D8)
ENDCASE

end
