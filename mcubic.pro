;; 1D Interpolation by Steffen, 1990
;; Guarantees a monotonic behavior of the interpolating function. The
;; result is a smooth curve with continuous first order derivatives
;; that passes through any given set of data points without spurious
;; oscillations.
;; http://adsabs.harvard.edu/full/1990A%26A...239..443S
;; Taken from gitlab
function mcubic, yinput, xinput, xout
;
    nout=n_elements(xout)
    nx  =n_elements(xinput)
    ny  =n_elements(yinput)
    if (nx eq ny) then begin
       n=nx
    endif else begin
        print, 'xinput and yinput must have the same number of elements!'
        stop
    endelse
;
;    print, 'xinput: ',xinput
;    print, 'yinput: ',yinput
;    print, 'xout:',xout
; --- we assume that xinput is a sequence at least 4 elements,
;     composed of unique, monotonically increasing or decraesing values ---
    if (n lt 4) then begin
       print, 'xinput, yinput must have at least 4 elements!'
       stop
    endif
;    
    if (xinput[0] lt xinput[n-1]) then begin
       xin=xinput
       yin=yinput
    endif else begin
       xin=reverse(xinput)
       yin=reverse(yinput)
    endelse
;
;  --- Compute x-step h, linear slope s and slope of parabola p ---
    h= xin[1:n-1] - xin[0:n-2]
    s=(yin[1:n-1] - yin[0:n-2])/h
    p=(s[0:n-3]*h[1:n-2]+s[1:n-2]*h[0:n-3])/(h[0:n-3]+h[1:n-2])
    q0=signum(s[0:n-3])+signum(s[1:n-2])
    q1=[[ abs(s[0:n-3])], [abs(s[1:n-2])], [0.5*abs(p)] ]
    q2=q0*MIN(q1,dim=2)
    y1=dblarr(n)
    y1[1:n-2]=q2
;    print, 'h:',h
;    print, 's:',s
;    print, 'p:',p
;    print, 'q0:',q0
;    print, 'q1:',q1
;    print, 'q2:',q2,q2.dtype
;    print, 'y1:',y1,y1.dtype
;
; --- Boundary handling d2y/dx2=0 [bcfac=1.5] ---
; --- y1(0)  =1.5*s(0)  -0.5*y1(1)
; --- y1(n-1)=1.5*s(n-2)-0.5*y1(n-2)
    yleft=1.5*s[0]   - 0.5*y1[1]
    yrght=1.5*s[n-2] - 0.5*y1[n-2]
    y1[0]  =yleft
    y1[n-1]=yrght
; --- Compute coefficients for cubic polynomial (inner points) ---
    a=(y1[0:n-2]+y1[1:n-1]-2.0*s)/h^2
    b=(3.0*s-2.0*y1[0:n-2]-y1[1:n-1])/h
    c=y1[ 0:n-2]
    d=yin[0:n-2]
;    print, 'a:',a
;    print, 'b:',b
;    print, 'c:',c
;    print, 'd:',d
;
; --- Linear extrapolation at upper boundary ---
    a = [a,0.0]
    b = [b,0.0]
    c = [c,y1[ n-1]]
    d = [d,yin[n-1]]
;    print, 'a:',a
;    print, 'b:',b
;    print, 'c:',c
;    print, 'd:',d
;
; --- Linear extrapolation at lower boundary ---
    a = [0.0,a]
    b = [0.0,b]
    c = [y1[ 0],c]
    d = [yin[0],d]
;    print, 'a:',a
;    print, 'b:',b
;    print, 'c:',c
;    print, 'd:',d
;
; --- Determine interval indices and relative x value ---
    iout=dblarr(nout)
    xrel=dblarr(nout)
    for j=0, nout-1 do begin
        if (xout[j] ge xin[0] and xout[j] le xin[n-1]) then begin
            for i=1, n-1 do begin
                if (xin[i] ge xout[j]) then begin
                    iout[j]=i
                    xrel[j]=xout[j]-xin[i-1]
                    break
                endif
            endfor
        endif else begin    
            if (xout[j] lt xin[0]) then begin
                iout[j]=0
                xrel[j]=xout[j]-xin[0]
            endif
            if (xout[j] gt xin[n-1]) then begin
                iout[j]=n
                xrel[j]=xout[j]-xin[n-1]
            endif
        endelse
    endfor
;
;    print, 'iout:',iout
;    print, 'xrel:',xrel
; --- Interpolation ---
    yout=((a[iout]*xrel+b[iout])*xrel+c[iout])*xrel+d[iout]
;
    return, yout
;    
end                            ; mcubic
    
