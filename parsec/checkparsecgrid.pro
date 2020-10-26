pro checkparsecgrid

;; check for allowed [Fe/H]
;; z grid points
allowedz = [0.0001d0,0.0002d0,0.0005d0,0.0010d0,0.0020d0,$
            0.0040d0,0.0060d0,0.0080d0,0.0140d0,0.0170d0,$
            0.0100d0,0.0200d0,0.0300d0,0.0400d0,0.0600d0]
allowedy = [0.2490d0,0.2490d0,0.2490d0,0.2500d0,0.2520d0,$
            0.2560d0,0.2590d0,0.2630d0,0.2730d0,0.2790d0,$
            0.2670d0,0.2840d0,0.3020d0,0.3210d0,0.3560d0]
allowedx = 1d0-allowedy-allowedz
allowedmh = alog10((allowedz/allowedx)/0.0207d0)
allowedfeh = allowedmh ;; assumes alpha = 0

;; mass grid points
allowedmass = [0.10d0,0.12d0,0.14d0,0.16d0,0.20d0,0.25d0,0.30d0,0.35d0,$
               0.40d0,0.45d0,0.50d0,0.55d0,0.60d0,0.65d0,0.70d0,0.75d0,$
               0.80d0,0.85d0,0.90d0,0.95d0,1.00d0,1.05d0,1.10d0,$
               1.20d0,1.25d0,1.30d0,1.35d0,1.40d0,1.45d0,1.50d0,1.55d0,$
               1.60d0,1.65d0,1.70d0,1.75d0,1.80d0,1.85d0,1.95d0,2.00d0,$
               2.05d0,2.10d0,2.15d0,2.20d0,2.30d0,2.40d0,2.60d0,2.80d0,$
               3.00d0,3.20d0,3.40d0,3.60d0,3.80d0,4.00d0,4.20d0,4.40d0,$
;               4.60d0,4.80d0,5.00d0,5.20d0,5.40d0,       5.80d0,6.00d0,$
               4.60d0,4.80d0,5.00d0,5.20d0,5.40d0,5.60d0,5.80d0,6.00d0,$ ;; 5.6 is missing a grid point
;                      6.40d0,7.00d0,8.00d0,       10.0d0,       12.0d0,$ ;; 6.2, 9.0, 11 are missing a grid point
               6.20d0,6.40d0,7.00d0,8.00d0,9.00d0,10.0d0,11.0d0,12.0d0,$
;               14.0d0,16.0d0,18.0d0,20.0d0,       30.0d0,35.0d0,40.0d0,$ ;; 24 missing a grid point
               14.0d0,16.0d0,18.0d0,20.0d0,24.0d0,30.0d0,35.0d0,40.0d0,$

;                      50.0d0,55.0d0,60.0d0,       70.0d0,              $ ;; 45, 65, 75, 80 are missing grid points
              45.0d0,50.0d0,55.0d0,60.0d0,65.0d0,70.0d0,75.0d0,80.0d0,$
               90.0d0,100d0,120d0,150d0,200d0,250d0,300d0,350d0]

for i=0L, n_elements(allowedfeh)-1 do begin
   for j=0L, n_elements(allowedmass)-1 do begin
      eep = readeep_parsec(allowedmass[j],allowedfeh[i])
;      print, i,j, allowedmass[j], allowedfeh[i]
   endfor
endfor

end

