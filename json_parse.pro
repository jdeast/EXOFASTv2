; $Id: //depot/InDevelopment/scrums/ENVI_Yukon/idl/idldir/lib/idl_base64.pro#3 $
;
; Copyright (c) 2012-2015, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;----------------------------------------------------------------------------
;+
; :Description:
;    The JSON_PARSE function takes a JSON string and converts it into
;    an IDL variable.
;
; :Returns:
;    If *String* begins with a curly brace "{", then the result is
;    a HASH, containing the unordered name-value pairs from the JSON string.
;    If *String* begins with a square bracket "[", then the result is
;    a LIST, containing the ordered collection of values.
;    If *String* is a scalar string that does not begin with "{" or "[",
;    then it is assumed to be a file name. If the file exists, then the
;    contents of the file are read and parsed.
;
;    *Note*: When converting JSON strings into IDL variables,
;    the following rules are used:
;      "null" becomes !NULL.
;      "false" becomes byte 0, "true" becomes byte 1.
;      Integer values become IDL variables of type LONG64.
;      Floating-point numbers become IDL variables of type DOUBLE.
;      Strings will be converted to IDL strings, and all \ escaped characters
;      will be converted back to normal characters.
;      JSON arrays will become IDL LIST variables.
;      JSON objects will become IDL HASH variables.
;   
;    *Note*: Since the HASH stores its name-value pairs in an arbitrary
;    order, the HASH::Keys() method may return the name-value pairs
;    in a different order compared to the original JSON string.
;   
; :Params:
;    String:
;      *String* must be a valid JSON string containing either a
;      JSON *object* of name-value pairs, or a JSON *array* of values.
;
; :Keywords:
;    None
;
;
; :Author:
;   CT, VIS, Jan 2012. Based on the IDLffJSON, written by Dawn Lenz.
;-
function JSON_Parse, String, DEBUG=debug, _EXTRA=extra

  compile_opt idl2, hidden
  
  if (~KEYWORD_SET(debug)) then begin
    ON_ERROR, 2
    ; Catch errors from our method and fake the call stack.
    CATCH, iErr
    if (iErr ne 0) then begin
      CATCH, /CANCEL
      MESSAGE, !ERROR_STATE.msg
    endif
  endif
  
  ; Make sure the object definition is compiled, in case this
  ; will be included in save files.
  IDLffJSON__Define

  ; See if our input is actually a file name.
  if (ISA(String, 'STRING', /SCALAR)) then begin

    char = STRMID(String, 0, 1)
    if (char ne '"' && char ne '[' && char ne '{') then begin

      if (FILE_TEST(String, /READ) && QUERY_ASCII(String, info)) then begin
        data = STRARR(info.lines)
        OPENR, lun, String, /GET_LUN
        READF, lun, data
        FREE_LUN, lun
        obj = OBJ_NEW("IDLffJSON")
        result = obj->Parse(STRJOIN(data), DEBUG=debug, _STRICT_EXTRA=extra)
        OBJ_DESTROY, obj
        return, result
      endif

    endif

  endif


  obj = OBJ_NEW("IDLffJSON")
  result = obj->Parse(String, DEBUG=debug, _STRICT_EXTRA=extra)
  OBJ_DESTROY, obj
  return, result
end

