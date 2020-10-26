;; Written by Thomas G Beatty ~2015
;; Modified to integrate into EXOFASTv2, 2017/06/28 JDE

function exofast_readdt, filename, lambdaRange

basename = file_basename(filename)

;; First index is velocity, second index is BJD
ccf2d = readfits(filename, /silent)
bjd = readfits(filename, exten_no=1, /silent)
vel = readfits(filename, exten_no=2, /silent)

if n_elements(lambdarange) eq 0 then lambdarange=0

model = ccf2d * 0.0 + median(ccf2d)
resid = ccf2d - model

telescope = (strsplit(basename,'.',/extract))(2)
rspec = (strsplit(basename,'.',/extract))(3)
if rspec eq 0d0 then message, 'ERROR: the DT filename does not properly encode the instrument resolution. See documentation for DTPATH in exofastv2.pro for details.'

night = strmid(basename,1,4)+'-'+strmid(basename,5,2)+'-'+strmid(basename,7,2)
label = 'UT ' + night + ' ' + telescope
pletters = ['b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
pname = (strsplit(basename,'.',/extract))(1)
pletter = strmid(pname,strlen(pname)-1,1)
planetndx = (where(pletters eq pletter))[0] 
if planetndx eq -1 then planetndx = 0
  
;stepsize = vel*0.
;nvels = n_elements(vel)
;for i=0,(nvels-2) do begin
;	stepsize[i] = vel[i+1]-vel[i]
;endfor
;meanstepsize = mean(stepsize) ;; bug! (includes zero from initialization)
;stepsize[-1] = meanstepsize

;; ~100x faster
stepsize = shift(vel,-1)-vel
stepsize[n_elements(stepsize)-1] = mean(stepsize[0:n_elements(stepsize)-2])

;; determine rough uncertainties for the observations
rms0 = stddev(ccf2d)
errscale = sqrt(total(((ccf2d-median(ccf2d))/rms0)^2) / (n_elements(ccf2d)-3d0))
rms = rms0 * errscale

chisqr0 = total(((ccf2d-median(ccf2d))/rms)^2)

dopptom=create_struct('ccf2d',ccf2d,'bjd',bjd,'vel',vel,'stepsize',stepsize,'rms',rms, $
                      'chisqr0',chisqr0,'Rspec',Rspec,'model',model,'resid',resid, $
                      'lambdaRange',lambdaRange,'label',label,'tel',telescope,'night',night, 'planetndx',planetndx)

return, dopptom

end
