;+
; NAME:
;   UPDATETIME
;
; NOTE: READ THE DOCUMENTATION -- ITS USE IS EASY TO GET WRONG!!!!
;
; PURPOSE:
;   This function updates three files required to maintain the most
;   precise time. These files are used by EXOFAST and related codes.
;
;   1) The leap second file, updates whenever there is a new leap
;      second. Typically every ~6 months. It is requires for ~1 second
;      precision.
;
;   2) The IERS file keeps track of meausured changes in the
;      Earth's rotation and is required for ms precision.
;
;   3) The BIPM file tracks the difference between TT(TAI) and
;      TT(BIPM) and is required for microsecond precison. Note that
;      this can only be done precisely months after the time in
;      question.
;
; NOTE:
;   There are no requirements that these files are maintained in a
;   backward compatible way. Changes to the source file formats or
;   file locations may occur, which could break this code. Care should
;   be taken to ensure this function is operating as intended. We will
;   update it as required, but may be slow to identify such
;   changes. Please email me (jason.eastman@cfa.harvard.edu) if the
;   most up to date version does not work.
;
; CALLING SEQUENCE:
;   updatetime
;
; NOTE: The files this code grabs update regularly. This function
; should be run at least monthly. The easiest way to ensure this is to
; set it up as a cronjob by adding a line like this in your crontab:
;
;   0 0 15 * * . /h/onion0/.bashrc ; idl -e "updatetime"
;
; ****** BE SURE TO TEST THIS! IT IS TRICKIER THAN IT SEEMS ******
;
; The first part
;   0 0 15 * * 
; will run it at midnight on the 15th of every month. You may wish to
; randomize the day and time to distribute the load on their FTP
; server. The second part
;   . /h/onion0/.bashrc 
; sources a bashrc file that defines the required environment
; variables ($ASTRO_DATA, $IDL_PATH, and $PATH). Just because it works
; in your terminal does not mean it will work inside crontab. You must
; change it so it sources your profile file. If your crontab shell is
; different than your usual shell, you may need to write a new profile
; file for your crontab. The last part
;   idl -e "updatetime"
; runs this program. idl must be in your $PATH for this to work (or
; you can specify the full path to IDL in the call). 
;
; DEPENDENCIES:
;   You must have the $ASTRO_DATA environment variable set. This is
;   the directory where the files will be copied.
;
; REVISION HISTORY:
;   2021/04/13 - Updated and documented, Jason Eastman (CfA)
;-
pro updatetime, forceupdate=forceupdate

path = getenv("ASTRO_DATA") + path_sep()

now = systime(/julian,/utc)

;; check to see when it was last updated, 
;; skip update if done in the past month
updatefile = path + 'exofast_lastupdate'
if query_ascii(updatefile) then begin ;; if a valid ascii file
   lastupdated = (double((read_ascii(updatefile)).field1))[0]
   if ~keyword_set(forceupdate) and (now - lastupdated) lt 0.9d0 then begin
      print, 'Last updated ' + strtrim(now - lastupdated,2) + ' days ago; skipping update'
      print, 'Use /FORCEUPDATE to update anyway'
      return
   endif
endif

caldat, now, month, day, year, hour, minute, second
fmt = '(i04,"-",i02,"-",i02,x,i02,":",i02,":",f06.3)'
nowstr = string(year,month,day,hour,minute,second,format=fmt)
bipm_filename = 'TTBIPM.' + strtrim(year-1,2)

;; big correction (~1 second/yr)
;; Craig Markwardt's TAI_UTC was written around the USNO
;; format, which is no longer available
;; Recreate it based on the file maintained and hosted by IERS
spawn, 'wget -NP $ASTRO_DATA https://hpiers.obspm.fr/iers/bul/bulc/ntp/leap-seconds.list'
readcol, path + 'leap-seconds.list', ntp_time, dtai, format='d,d', comment='#',/silent
info = file_info(path+'leap-seconds.list')
caldat, julday(1,1,1970,0,0,info.mtime), leapmo, leapday, leapyr, leaphr, leapmin, leapsec
leap_update_str = string(leapyr,leapmo,leapday,leaphr,leapmin,leapsec,format=fmt)
jd = julday(1,1,1900,0,0,ntp_time)
caldat, jd, month, day, year
monthstrs = ['','JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
taifile = path + 'tai-utc.dat'

exofast_printf, lun, ' 1961 JAN  1 =JD 2437300.5  TAI-UTC=   1.4228180 S + (MJD - 37300.) X 0.001296 S ', filename=taifile, /new
exofast_printf, lun, ' 1961 AUG  1 =JD 2437512.5  TAI-UTC=   1.3728180 S + (MJD - 37300.) X 0.001296 S ', filename=taifile
exofast_printf, lun, ' 1962 JAN  1 =JD 2437665.5  TAI-UTC=   1.8458580 S + (MJD - 37665.) X 0.0011232S ', filename=taifile
exofast_printf, lun, ' 1963 NOV  1 =JD 2438334.5  TAI-UTC=   1.9458580 S + (MJD - 37665.) X 0.0011232S ', filename=taifile
exofast_printf, lun, ' 1964 JAN  1 =JD 2438395.5  TAI-UTC=   3.2401300 S + (MJD - 38761.) X 0.001296 S ', filename=taifile
exofast_printf, lun, ' 1964 APR  1 =JD 2438486.5  TAI-UTC=   3.3401300 S + (MJD - 38761.) X 0.001296 S ', filename=taifile
exofast_printf, lun, ' 1964 SEP  1 =JD 2438639.5  TAI-UTC=   3.4401300 S + (MJD - 38761.) X 0.001296 S ', filename=taifile
exofast_printf, lun, ' 1965 JAN  1 =JD 2438761.5  TAI-UTC=   3.5401300 S + (MJD - 38761.) X 0.001296 S ', filename=taifile
exofast_printf, lun, ' 1965 MAR  1 =JD 2438820.5  TAI-UTC=   3.6401300 S + (MJD - 38761.) X 0.001296 S ', filename=taifile
exofast_printf, lun, ' 1965 JUL  1 =JD 2438942.5  TAI-UTC=   3.7401300 S + (MJD - 38761.) X 0.001296 S ', filename=taifile
exofast_printf, lun, ' 1965 SEP  1 =JD 2439004.5  TAI-UTC=   3.8401300 S + (MJD - 38761.) X 0.001296 S ', filename=taifile
exofast_printf, lun, ' 1966 JAN  1 =JD 2439126.5  TAI-UTC=   4.3131700 S + (MJD - 39126.) X 0.002592 S ', filename=taifile
exofast_printf, lun, ' 1968 FEB  1 =JD 2439887.5  TAI-UTC=   4.2131700 S + (MJD - 39126.) X 0.002592 S ', filename=taifile
fmt2 = '(i5,x,a3,x,i2," =JD ", f9.1,"  TAI-UTC=",f6.1,"       S + (MJD - 41317.) X 0.0      S ")'
for i=0L, n_elements(year)-1 do $
   exofast_printf, lun, string(year[i], monthstrs[month[i]], day[i], jd[i], dtai[i],format=fmt2), filename=taifile



;exofast_printf, lun, '; Source: $EXOFAST_PATH/updatetime.pro', filename=taifile,/new
;exofast_printf, lun, '; Updated: ' + nowstr, filename=taifile
;exofast_printf, lun, '; Based on $ASTRO_DATA/leap-seconds.list updated on ' + leap_update_str, filename=taifile
;exofast_printf, lun, '; Leap Seconds Table - used by CDF', filename=taifile
;exofast_printf, lun, '; Update it when a leap second(s) is added.', filename=taifile
;exofast_printf, lun, "; Comment lines starts with ';' at column 1.", filename=taifile
;exofast_printf, lun, '; Year Month Day Leap Seconds      Drift', filename=taifile
;exofast_printf, lun, '  1960   1    1    1.4178180  37300.0  0.001296', filename=taifile
;exofast_printf, lun, '  1961   1    1    1.4228180  37300.0  0.001296', filename=taifile
;exofast_printf, lun, '  1961   8    1    1.3728180  37300.0  0.001296', filename=taifile
;exofast_printf, lun, '  1962   1    1    1.8458580  37665.0  0.0011232', filename=taifile
;exofast_printf, lun, '  1963  11    1    1.9458580  37665.0  0.0011232', filename=taifile
;exofast_printf, lun, '  1964   1    1    3.2401300  38761.0  0.001296', filename=taifile
;exofast_printf, lun, '  1964   4    1    3.3401300  38761.0  0.001296', filename=taifile
;exofast_printf, lun, '  1964   9    1    3.4401300  38761.0  0.001296', filename=taifile
;exofast_printf, lun, '  1965   1    1    3.5401300  38761.0  0.001296', filename=taifile
;exofast_printf, lun, '  1965   3    1    3.6401300  38761.0  0.001296', filename=taifile
;exofast_printf, lun, '  1965   7    1    3.7401300  38761.0  0.001296', filename=taifile
;exofast_printf, lun, '  1965   9    1    3.8401300  38761.0  0.001296', filename=taifile
;exofast_printf, lun, '  1966   1    1    4.3131700  39126.0  0.002592', filename=taifile
;exofast_printf, lun, '  1968   2    1    4.2131700  39126.0  0.002592', filename=taifile
;for i=0L, n_elements(year)-1 do $
;   exofast_printf, lun, string(year[i], month[i], day[i], dtai[i], '0.0','0.0', $
;                       format='(i6,i4,i5,f7.1,a15,a5)'), filename=taifile

;; small corrections (ms)
spawn, 'wget --ftp-user=anonymous -NP ' + path + ' https://datacenter.iers.org/data/latestVersion/finals.all.iau2000.txt'
file_copy, path + 'finals.all.iau2000.txt', path + 'iers_final_a.dat', /overwrite

;; smallest correction from TT(TAI) to TT(BIPM) (microseconds)
;; only available months after the time in question
bipm_url = 'https://webtai.bipm.org/ftp/pub/tai/ttbipm/'
spawn, 'wget --ftp-user=anonymous -NP ' + path + ' ' + bipm_url + bipm_filename
file_copy, path + bipm_filename, path + 'bipmfile', /overwrite

;; update the "last updated" file to now
exofast_printf, lun, strtrim(now,2), filename=updatefile, /new

end
