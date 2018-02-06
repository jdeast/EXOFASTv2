pro mksimdata

G = 2942.71377d0 ;; R_sun^3/(m_sun*day^2), Torres 2010
au = 215.094177d0
period0 = 7.8457d0
period1 = 1312.6
teff = 5768.1665
feh = 0
gamma0 = -314
gamma1 = 654
tc0 = julday(1,1,2013,0,0,0)
tc1 = julday(6,1,2013,0,0,0)
inc = !dpi/2d0
k0 = 40
k1 = 400
p = 0.1
mstar = 1d0
age = 5d0
e0 = 0d0
omega0 = !dpi/2d0
e1 = 0d0
omega1 = !dpi/2d0


dummy = massradius_yy3(mstar, feh, age, teff, yyteff=yyteff, yyrstar=rstar)
logg = alog10(27443.4141*mstar/rstar^2)
coeffs = quadld(logg, teff, feh,'Sloanz')
u1 = coeffs[0]
u2 = coeffs[1]
f0 = 1d0

mp = ktom2(K0, e0, inc, period0, mstar)

arsun = (G*(mstar+mp)*period0^2/(4d0*!dpi^2))^(1d0/3d0) ;; rsun
a = arsun/AU ;; au
ar = arsun/rstar ;; rstar

phase=exofast_getphase(e0,omega0,/pri)
tp0 = tc0 - period0*phase

phase=exofast_getphase(e1,omega1,/pri)
tp1 = tc1 - period1*phase


ntransit = 300
tranerr = dblarr(ntransit) + 0.001 
trantimes = tc0 + ((dindgen(ntransit)/(ntransit-1)) - 0.5d0)*0.45

modelflux = exofast_tran(trantimes, inc, ar, tp0, period0, e0, omega0, p, u1, u2, f0, rstar=rstar/au, tc=tc) 

exofast_forprint, trantimes, modelflux+randomn(seed,ntransit)*tranerr, tranerr, format='(3(f0.6,x))',/nocomment,textout='n20130101.z.test.dat'
exofast_forprint, trantimes, modelflux+randomn(seed,ntransit)*tranerr, tranerr, format='(3(f0.6,x))',/nocomment,textout='n20130101.z.test2.dat'
exofast_forprint, trantimes, modelflux+randomn(seed,ntransit)*tranerr, tranerr, format='(3(f0.6,x))',/nocomment,textout='n20130101.z.test3.dat'
exofast_forprint, trantimes, modelflux+randomn(seed,ntransit)*tranerr, tranerr, format='(3(f0.6,x))',/nocomment,textout='n20130101.z.test4.dat'
exofast_forprint, trantimes, modelflux+randomn(seed,ntransit)*tranerr, tranerr, format='(3(f0.6,x))',/nocomment,textout='n20130101.z.test5.dat'
exofast_forprint, trantimes, modelflux+randomn(seed,ntransit)*tranerr, tranerr, format='(3(f0.6,x))',/nocomment,textout='n20130101.z.test6.dat'

plot, trantimes-tc0, modelflux, psym=1,/ys,/xs

nrv = 50
rverr = 3d0 + dblarr(nrv)
rvtimes = tc1 + 2*period1*((dindgen(nrv)/(nrv-1)) - 0.5d0)


modelrv = dblarr(nrv)
modelrv += exofast_rv(rvtimes,tp0,period0,0d0,K0,e0,omega0)
modelrv += exofast_rv(rvtimes,tp1,period1,0d0,K1,e1,omega1)
exofast_forprint, rvtimes, modelrv+randomn(seed,nrv)*rverr+gamma0, rverr, format='(3(f0.6,x))',/nocomment,textout='test.SIM.rv'
exofast_forprint, rvtimes, modelrv+randomn(seed,nrv)*rverr+gamma1, rverr, format='(3(f0.6,x))',/nocomment,textout='test.SIM2.rv'

plot, rvtimes,modelrv, psym=1

stop
end
