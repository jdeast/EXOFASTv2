;+
; NAME:
;   splog
;
; PURPOSE:
;   Logging routine for writing to standard output and/or a log file.
;
; CALLING SEQUENCE:
;   splog, v1, v2 ..., [, _EXTRA=extra, /noname, prelog=, filename=, $
;    /append, /close, /no_stdout ]
;
; INPUTS:
;   v1, v2 ... - The expressions to be passed to PRINT or PRINTF
;
; OPTIONAL KEYWORDS:
;   _EXTRA     - Any keywords for PRINT or PRINTF, such as FORMAT
;   noname     - If set, then suppress printing name of calling routine
;   prelog     - If set, then print this string immediately after the
;                name of the calling routine on each line, i.e. 'b1'
;   filename   - If specified, then open this file for text output
;   append     - If set at the same time as FILENAME, then append to this file;
;                default is to overwrite file
;   close      - If set, then close the output file
;   no_stdout  - If set, then do not print anything to the standard output.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The output is formatted just like with the IDL PRINT command, except
;   that extraneous whitespace is removed from non-STRING elements unless
;   a FORMAT keyword is used.
;
; EXAMPLES:
;   Open a file for text output, write to it, then close it:
;     splog, filename='test.log'
;     splog, 'This is a test', 123.4, format='(a,f8.2)'
;     splog, /close
;   Alternatively, this can all be done on one line
;     splog, filename='test.log', /close, $
;      'This is a test', 123.4, format='(a,f8.2)'
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-Nov-1999  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
; Trim extra whitespace (multiple blanks) from any string conversion
; iff the FORMAT keyword is not specified, and this isn't already a string.
; Any strings that are passed will not have whitespace removed.
function splog_trim, v

   common com_splog_trim, qcompress

   if (n_elements(v) NE 0) then $
    vt = (qcompress AND size(v,/tname) NE 'STRING') $
     ? strcompress(string(v)) : v $
   else $
    vt = '' ; unused return value

   return, vt
end

;------------------------------------------------------------------------------
pro splog, noname=noname, prelog=prelog, $
 filename=filename, append=append, close=close, no_stdout=no_stdout, $
 v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, _EXTRA=extra

   ; Declare LOGLUN in a common block so that it is remembered between calls.
   common com_splog, loglun, fullprelog
   common com_splog_trim, qcompress

   if (keyword_set(filename)) then begin
      ; First close a file if one is already open
      if (keyword_set(loglun)) then splog, /close

      ; Now open the file
      get_lun, loglun
      openw, loglun, filename, append=append
   endif

   if (n_elements(prelog) EQ 1) then fullprelog = prelog

   ; Determine the name of the calling routine
   help, calls=calls

;   fname = (str_sep(calls[1], ' '))[0] + ': '
   fname = (strsplit(calls[1], ' ',/extract))[0] + ': '

   ; Add spaces for depth of routine
   for i=0,n_elements(calls)-4 do fname = ' ' + fname 

   ;----------
   ; If there is no FORMAT keyword specified, then compress the output string
   ; (remove extraneous whitespace).

   qcompress = 1
   if (keyword_set(extra)) then begin
      tags = tag_names(extra)
      if ( (where(tags EQ 'FORMAT'))[0] NE -1) then qcompress = 0
   endif

   vt1 = splog_trim(v1)
   vt2 = splog_trim(v2)
   vt3 = splog_trim(v3)
   vt4 = splog_trim(v4)
   vt5 = splog_trim(v5)
   vt6 = splog_trim(v6)
   vt7 = splog_trim(v7)
   vt8 = splog_trim(v8)
   vt9 = splog_trim(v9)
   vt10 = splog_trim(v10)
   vt11 = splog_trim(v11)
   vt12 = splog_trim(v12)

   ;----------
   ; Construct the output text string

   nv = n_params()
;   if (nv GT 0 OR keyword_set(extra)) then begin
   if (nv GT 0) then begin
      case nv of
;      0: textstring = string('', _EXTRA=extra) ; Does not work! IDL bug?
      1: textstring = string(vt1, _EXTRA=extra)
      2: textstring = string(vt1, vt2, _EXTRA=extra)
      3: textstring = string(vt1, vt2, vt3, _EXTRA=extra)
      4: textstring = string(vt1, vt2, vt3, vt4, _EXTRA=extra)
      5: textstring = string(vt1, vt2, vt3, vt4, vt5, _EXTRA=extra)
      6: textstring = string(vt1, vt2, vt3, vt4, vt5, vt6, _EXTRA=extra)
      7: textstring = string(vt1, vt2, vt3, vt4, vt5, vt6, vt7, _EXTRA=extra)
      8: textstring = string(vt1, vt2, vt3, vt4, vt5, vt6, vt7, vt8, _EXTRA=extra)
      9: textstring = string(vt1, vt2, vt3, vt4, vt5, vt6, vt7, vt8, vt9, _EXTRA=extra)
      10: textstring = string(vt1, vt2, vt3, vt4, vt5, vt6, vt7, vt8, vt9, vt10, _EXTRA=extra)
      11: textstring = string(vt1, vt2, vt3, vt4, vt5, vt6, vt7, vt8, vt9, vt10, vt11, _EXTRA=extra)
      else: textstring = string(vt1, vt2, vt3, vt4, vt5, vt6, vt7, vt8, vt9, vt10, vt11, vt12, _EXTRA=extra)
      endcase
   endif

   ;----------
   ; Write to standard out

   if (NOT keyword_set(no_stdout) $
    AND (nv GT 0 OR keyword_set(extra))) then begin
      if (NOT keyword_set(noname)) then print, fname, format='(a,$)'
      if (keyword_set(fullprelog)) then $
       print, fullprelog+': ', format='(a,$)'
      if (nv EQ 0) then print, _EXTRA=extra $
       else print, textstring
      IF NOT (LMGR(/DEMO)) THEN flush,-1
  endif

   ;----------
   ; Write to the output file (if it exists)

   if ((nv GT 0 OR keyword_set(extra)) AND keyword_set(loglun)) then begin
      if (NOT keyword_set(noname)) then printf, loglun, fname, format='(a,$)'
      if (keyword_set(fullprelog)) then $
       printf, loglun, fullprelog+': ', format='(a,$)'
      if (nv EQ 0) then printf, loglun, _EXTRA=extra $
       else printf, loglun, textstring
      flush, loglun
   endif

   ;----------
   ; Close the output file

   if (keyword_set(close) AND keyword_set(loglun)) then begin
      close, loglun
      free_lun, loglun
      loglun = 0
      fullprelog = ''
   endif

   return
end
;------------------------------------------------------------------------------
