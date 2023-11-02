;+
; NAME:
;   READIPAC
;
; PURPOSE:
;   Reads the http://exoplanets.org/csv-files/exoplanets.csv database into an
;   IDL structure for easy plotting and manipulation.
;
;   See http://exoplanets.org/help/common/data for the definitions
;   of each tag.
;
; OPTIONAL KEYWORDS:
;   UPDATE - By default, if exoplanets.csv exists, it will not
;            retrieve a new copy. Set this keyword to wget the csv
;            file and overwrite exoplanets.csv. It will only overwrite
;            it if the copy on the server is newer. This program is
;            about 25% slower when this keyword is specified
;            unnecessarily and 2x slower when it actually needs to
;            update.
;
; CALLING SEQUENCE:
;   result = readexo([/update])
;
; EXAMPLE - Print all data for HAT-P-13 b:
;
;   planetdata = readexo(/update)
;   match = where(planetdata.name eq 'HAT-P-13 b')
;   names = tag_names(planetdata)
;   if match[0] ne -1 then for i=0, n_tags(planetdata)-1 do $
;     print, names[i], planetdata.(i)[match], format='(a15,a100)'
;
; DEPENDENCIES:
;   wget
;
; ENVIRONMENT VARIABLES:
;   EXOFAST_PATH  - The environment variable specifying the directory
;                   in which to save/look for the CSV file. If not set,
;                   the program will use the current working directory.
;
; NOTES: 
;   All data types are strings for robustness and generality (no
;   update to this program is required if new rows or columns are
;   added to the database). Converting to numeric values is trivial,
;   but be sure to use the appropriate precision (eg, to convert ra
;   from hours to degrees, do radeg = planetdata.ra*15.d0, not radeg =
;   planetdata.ra*15).
;
;   Be aware that the order of tags or even the tagnames may change if
;   the database changes.
;
; REVISION HISTORY
;   Created: 2011/01/19 Jason Eastman, OSU
;            2012/06/15 Updated for change to exoplanets.org.
;            2013/01/14 Updated for change to URL.
;-

function readipac, update=update, default_only=default_only

;; update the file if desired or not present
path = getenv('EXOFAST_PATH')
if path eq '' then path = './'
filename = path + 'planets.csv'
url = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+ps&format=csv'

if not file_test(filename) or keyword_set(update) then begin
    spawn, 'wget -P ' + path + $
      ' "' + url + '"',output,err,exit_status=status
    if status ne 0 then message, strjoin(err,string(10B))
    wgetname = 'sync?query=select+*+from+ps&format=csv'
    spawn, 'mv "' + path + wgetname + '" ' + filename
endif

;; read the header line
openr, table, filename, /get_lun
line = ''
repeat begin
   readf, table, line
endrep until strpos(line,'#') ne 0

;; parse the headers into tag names
;; See http://exoplanetarchive.ipac.caltech.edu/docs/API_exoplanet_columns.html
;; for the definitions of each tag.
tagnames = strcompress(strsplit(line,',',/extract,/preserve_null),/remove_all)
ntags = n_elements(tagnames)
nentries = file_lines(filename) - 1

;; read the table into a string array
entry=0
array = strarr(ntags,nentries)
while not eof(table) do begin
    readf, table, line
    values = strtrim(strsplit(line,',',/extract,/preserve_null),2)
    if n_elements(values) gt ntags then begin
        message, "WARNING: too many values for row " + strtrim(entry,2) + $
          ", skipping row (check for extraneous commas in CSV file)",/con
    endif else array[0:n_elements(values)-1,entry] = values
    entry++
endwhile
free_lun, table

;; create a structure with the tag names
structure = create_struct(tagnames[ntags-1],transpose(array[ntags-1,*]))

if default_only then begin

stop
   ndx = (where(tagnames eq 'default_flag'))
   selection = where(array[ndx,*] eq '1')
endif else selection = lindgen(n_elements(array[0,*]))

for i=ntags-2,0,-1 do begin
    structure = create_struct(tagnames[i],transpose(array[i,selection]),structure)
endfor

;; return the structure
return, structure

end


 
