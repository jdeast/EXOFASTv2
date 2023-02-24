pro mkfitgif, idlfile, skip=skip, chain=chain, maxstep=maxstep, redoframe=redoframe, redocorner=redocorner

;; this is the (1 element) structure used by exofast_chi2v2
COMMON chi2_block, ss

restore, idlfile
ss = mcmcss2ss(mcmcss)

;; fix the plot ranges so the animations look better
ss.sedrange=[0.1,30,-12,-9,-3,3]
ss.transitrange=[-2,2,0.984,1.0055,-0.006,0.006]
ss.rvrange = [0,1,-100,100,-10,10]
ss.emrange = [5800,4200,5,3]

;; wasp-76
ss.sedrange = [0.1,30,-13,-8,-5,5]
ss.transitrange =[-4,4,0.984,1.0055,-0.006,0.006]
ss.emrange = [7100,4400,5,3]
ss.rvrange = [0,1,-150,150,-50,50]
transitpage = 5

;; by default, do the whole chain in 100 steps
nsteps = mcmcss.nsteps/mcmcss.nchains
nchains = mcmcss.nchains
if n_elements(skip) eq 0 then skip = floor(nsteps/100d0)
if n_elements(chain) eq 0 then chain = 0L
if n_elements(maxstep) eq 0 then maxstep=nsteps
if maxstep gt nsteps then maxstep=nsteps

;; format for each frame; include enough digits for up to maxsteps
ndigits = ceil(alog10(maxstep))
format = '(i0' + strtrim(ndigits,2) + ')'

dirname = file_dirname(mcmcss.prefix) + path_sep() + 'gif' + path_sep()
if not file_test(dirname) then file_mkdir, dirname
gifname = dirname + file_basename(mcmcss.prefix) + 'gif'

rgb = colortable(39,ncolors=nchains)
xsize=500
ysize=200
set_plot, 'Z',/copy
device, set_resolution=[xsize,ysize], set_pixel_depth=24, set_font='Times',/tt_font
chi2 = reform(*mcmcss.chi2,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains)
loglike = -chi2/2d0
;loglike = -reform(*mcmcss.chi2,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains)/2d0

ymin = min(loglike,max=ymax)
if ymin lt 0 then sign = ' + ' else sign = ' - '
latex = textoidl('\Delta log(Like)') ; + sign + strtrim(abs(ymax),2)
white = 'ffffff'x
white = 256L^2*255L + 256L*255 + 255
black = '000000'x  



for i=0L, maxstep-1, skip do begin

   base = dirname + file_basename(mcmcss.prefix) + 'frame.' + string(i,format=format)

   print, 'creating frame ' + base + '.exofastv2.png'

   if file_test(base + '.exofastv2.png') and not keyword_set(redoframe) then continue

   if i gt 3 then begin
      burnndx = getburnndx(chi2[0:i,*],goodchains=goodchains,/silent)
      burnndx = (i-3) < burnndx
   endif else begin
      burnndx = 0
      goodchains = lindgen(nchains)
   endelse


   ;; Likelihood chain plot
   yminplot = min(loglike[burnndx:i,goodchains],max=ymaxplot)
   yminplot = -3000+ymax
   ymaxplot = 0d0+ymax
   plot, [0],[0], xrange=[0,i],yrange=[yminplot,ymaxplot]-ymax, $
         ytitle=latex,xtitle='Chain Link number',/xs,background=white,color=black
   for j=0, n_elements(goodchains)-1 do begin
      color = rgb[goodchains[j],0]*256L^2 + rgb[goodchains[j],1]*256L + rgb[goodchains[j],2]
      oplot, loglike[0:i,goodchains[j]]-ymax,color=color
   endfor
   oplot, loglike[0:i,chain]-ymax,color=black,thick=5

;   plot, loglike[0:i,0]-ymax, ytitle=latex, $
;         xtitle='Chain link number',xrange=[0,i],/xs,background=white,color=black
   oplot, [burnndx,burnndx],[-9d9,9d9],color=black

;stop

;   print, burnndx

   chain_png_name = base + '.chi2chain.png'
   write_png, chain_png_name, tvrd(/true)

   if keyword_set(redocorner) or not file_test(base+'.star_corner.ps') then begin
      ;; make a corner plot for the stellar parameters
      sample = burnndx + lindgen((i-burnndx) > 1)
      remake_corner, idlfile, tags=['mstar','rstar','teff','av','age'     ],$
                     sample=sample,/keepburn, psname=base+'.star_corner.ps'
      ;; make a corner plot for the planet parameters
      remake_corner, idlfile, tags=['mp'   ,   'rp', 'b'  ,'e' ,'omegadeg'],$
                     sample=sample,/keepburn, psname=base+'.planet_corner.ps'
   endif
   
   mksummaryframe, mcmcss, i + chain*nsteps, base=base, transitpage=transitpage
endfor

spawn, 'convert ' + dirname + '*.exofastv2.png ' + gifname
display=1
if keyword_set(display) then spawn, 'firefox ' + gifname + ' &'

stop

end
