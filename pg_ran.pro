;+
; NAME:
;
;  pg_ran
;
; PURPOSE:
;
; generates uniform random numbers
;
; CATEGORY:
;
; numeric utils
;
; CALLING SEQUENCE:
;
; data=pg_ran(seed,n)
;
; INPUTS:
;
; n: number of uniform random deviates to produce
;
; OPTIONAL INPUTS:
;
; seed: seed for the random number generator. Input is either a scalar integer
;       or an array of 3 UL64 integers. The internal seed of the routine is the
;       array of 3 UL64 integers.
;
; KEYWORD PARAMETERS:
;
; double: if set, output is double precision instead of ulong64
;
; OUTPUTS:
;
; data: a set of ulong64 integers, uniformly randomly distributed
;       between 0ULL and 2ULL^64-1
;       (if keyword /double is set - doubles between 0d and 1d will
;       be returned instead)
;
; OPTIONAL OUTPUTS:
;
; seed: the seed to be used to get the next number in the sequence.
;       It's an array of 3 L64 integers.
;
; none
;
; COMMON BLOCKS:
;
; none
;
; SIDE EFFECTS:
;
; none
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
; IDL implementation of the routine "ran" of numerical
; recipes 3d edition (page 342). This is a high quality
; random generator with a period of ~ 3.138x10^57.
;
; EXAMPLE:
;
;
; seed=17
; print,pg_ran(seed,10)
;    269952321389814056
;   7477734313819993120
;  16294976781531816119
;  17039904789424739738
;   4945048831639962635
;   1565409385732501729
;   7095006703038622919
;  13927236388846696772
;    150171266583137103
;   8092874815854888167
; seed=17
; print,pg_ran(seed,10,/double),format='(d17.15)'
;     0.014634144665917
;     0.405368789415652
;     0.883352461356882
;     0.923735089582023
;     0.268071634315548
;     0.084861012842019
;     0.384621084061685
;     0.754996997475338
;     0.008140800673717
;     0.438715622850155
;     
;
; MODIFICATION HISTORY:
;
; 2008 Paolo C. Grigis written
; 07-AUG-2009 Paolo C. Grigis updated documentation
; 12-NOV-2009 PG added /double keyword 
;-

;this is the inner loop of the function, works only
;with long64 integers and updates the seed
FUNCTION pg_ran_64,u,v,w

u=u*2862933555777941757ULL+7046029254386353087ULL

v XOR= ishft(v,-17)
v XOR= ishft(v, 31)
v XOR= ishft(v, -8)

w=4294957665ULL*( w AND  4294967295ULL) +ishft(w,-32)

x= u XOR ishft(u,21)

x XOR= ishft(x,-35)
x XOR= ishft(x,  4)

RETURN, (x+v) XOR w

END

;main routine
FUNCTION pg_ran,seed,n,double=double

;clean up seed
IF n_elements(seed) EQ 3 THEN BEGIN 
   u=seed[0]
   v=seed[1]
   w=seed[2]
ENDIF $ 
ELSE BEGIN 

   IF n_elements(seed) EQ 0 THEN j=1ULL ELSE j=ulong64(seed)

   v=4101842887655102017ULL
   w=1ULL

   ;initialize 
   u=j XOR v
   dummy=pg_ran_64(u,v,w)
   v=u
   dummy=pg_ran_64(u,v,w)
   w=v
   dummy=pg_ran_64(u,v,w)
   
ENDELSE 

;allocates array for result
randomdev=ulon64arr(n)

;loops over array
FOR i=0L,n-1 DO BEGIN 
   randomdev[i]=pg_ran_64(u,v,w)
ENDFOR

IF keyword_set(double) THEN randomdev=double(randomdev*5.42101086242752217d-20)

seed=[u,v,w]

RETURN,randomdev

END







