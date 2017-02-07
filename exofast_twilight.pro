;+
; NAME:
;   exofast_twilight
; PURPOSE:
;   Calculate the time of twilight (sunrise/sunset) to ~2 minute
;   accuracy (as compared to sunpos/eq2hor, which claims 2.2"/1"
;   accuracy, respectively).
;
;   Algorithm from http://williams.best.vwh.net/sunrise_sunset_algorithm.htm
;   Almanac for Computers, 1990
;   published by Nautical Almanac Office
;   United States Naval Observatory
;   Washington, DC 20392
;
; INPUTS:
;   JD        - The Julian date of days to calculate the twilight, in
;               UT. Scalar or vector.
;   OBSNAME   - Name of the observatory (calls
;               observatory.pro). Either OBSNAME or
;               LATITUDE/LONGITUDE must be specified.
;   LATITUDE  - Latitude of the observatory, in degrees.
;   LONGITUDE - Longitude of the observatory, in degrees West.
;
; OPTIONAL INPUTS:
;   HORIZON   - The Sun altitude that is considered twilight, in
;               degrees.
;               Civil Twilight = -6 degrees.
;               Nautical Twilight = -12 degrees. (Default)
;               Astronomical Twilight = -18 degrees.
;   TZ        - The timezone of the observatory, in hours west of
;               Greenwich. Ignored unless /LOCAL is set.
;
; OPTIONAL KEYWORDS:
;   SUNSET    - Return the time of sunset instead of sunrise.
;   LOCAL     - Return the result in local time. TZ must be specified
;               if this keyword is set.
;               NOTE: This does not take into account DST.
;
; RESULT:
;   Returns the julian date of sunrise (or sunset, if /SUNSET is
;   set). If the sun is always down for the particular date, the
;   result will be -1. If the sun is always up for the particular
;   date, the result will be -2.
;
; REVISION HISTORY:
;   2011/01/27 - Written: Jason Eastman - The Ohio State University
;-

function exofast_twilight, jd, sunset=sunset, horizon=horizon, $
                           latitude=latitude, longitude=longitude, $
                           obsname=obsname, tz=tz, local=local

if n_elements(horizon) eq 0 then horizon=-12

if keyword_set(obsname) then begin
    observatory, obsname, obs
    latitude = obs.latitude
    longitude = obs.longitude
    tz = obs.tz
endif else begin
    if n_elements(latitude) eq 0 or $
      n_elements(longitude) eq 0 then $
      message, 'Must specify either obsname or latitude/longitude'
    
    if n_elements(tz) eq 0 and keyword_set(local) then $
      message, 'Must specify TZ if /LOCAL is set'

endelse

zenith = (90.d0 - horizon)*!dpi/180.d0

;; day of year
caldat, jd, month, day, year
doy = julday(month, day, year,0,0,0) - julday(1,0,year,12,0,0)

;; convert the longitude to hour value
lnghour = -longitude/15.d0
if keyword_set(sunset) then begin
    t = doy + ((18.d0 - lngHour) / 24.d0)
endif else t = doy + ((6.d0 - lngHour) / 24.d0)

;; Sun's mean anomaly (radians)
M = (0.9856d0*t - 3.289d0)*!dpi/180.d0

;; Sun's true longitude (radians)
L = M + (1.916d0*sin(M) + 0.02d0*sin(2.d0*M) + 282.634d0)*!dpi/180.d0
L = (L + 2.d0*!dpi) mod (2.d0*!dpi)

;; calculate the Sun's right ascension (Hours)
RA = atan(0.91764d0*tan(L))
RA = ((RA + 2.d0*!dpi) mod (2.d0*!dpi))*12.d0/!dpi

;; right ascension value needs to be in the same quadrant as L
Lquadrant  = (floor(2.d0*L/!dpi))*6d0
RAquadrant = (floor(RA/6.d0))*6d0
RA += (Lquadrant - RAquadrant)

;; calculate the Sun's declination
sinDec = 0.39782d0*sin(L)
cosDec = cos(asin(sinDec))

;; calculate the Sun's local hour angle
cosH = (cos(zenith) - sinDec*sin(latitude*!dpi/180.d0) ) / $
  ( cosDec*cos(latitude*!dpi/180.d0) )

;; the sun never rises on this location (on the specified date)
norise = where(cosH gt 1)
noset = where(cosH lt -1)

;; finish calculating H and convert into hours
if keyword_set(sunset) then begin
    H = acos(cosH)*12.d0/!dpi
endif else H = 24.d0 - acos(cosH)*12.d0/!dpi

;;calculate local mean time of rising/setting
time = H + RA - 0.06571d0*t - 6.622d0

;; adjust back to UTC
UT = Time - lngHour
UT = ((UT + 24) mod 24)
twilightjd = julday(month, day, year, UT)

;; make sure its the right time
if keyword_set(sunset) then begin
    toolow = where(twilightjd - floor(twilightjd) gt 0.5d)
    if toolow[0] ne -1 then twilightjd[toolow]++
endif

;; adjust to local time
if keyword_set(LOCAL) then twilightjd -= tz/24.d0

;; special cases
if norise[0] ne -1 then twilightjd[norise] = -1
if noset[0] ne -1 then twilightjd[noset] = -2

return, twilightjd

end
