; Repacked into one routine - D. Finkbeiner 23 Feb 2000
; Fixed () -> [] bugs here and there - Hogg 15 Jun 2000
; Fixed one more () -> [] bug - Hogg 12 Jul 2000
;
;+
; NAME:
;       TEXTOIDL_STRTRANS
; PURPOSE:
;       Translate all occurences of one substring to another.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       new = textoidl_strtrans(oldstr,from,to,ned)
; INPUTS:
;       oldstr -- string on which to operate.              in
;                 May be an array.
;       from   -- substrings to be translated. May be      in
;                 an array.
;       to     -- what strings in from should be           in
;                 translated to. May be an array.
; KEYWORD PARAMETERS:
;       /HELP  -- Set this to print useful message and 
;                 exit.
; OUTPUTS:
;       new    -- Translated string. Array if oldstr is    out          
;                 an array.
;       ned    -- number of substitutions performed in     out
;                 oldstr.  Array if oldstr is an array.
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
;       - Any of old, from, and to can be arrays.  
;       - from and to must have the same number of elements.
; EXAMPLE:
;       inp='Many*bad!chars+in_here'
;       from=['*','!','+','_']
;       to  =[' ',' ',' ',' ']
;       out = textoidl_strtrans(inp,from,to,ned)
;       Will produce out='Many bad chars in here', and set ned to 4.
; MODIFICATION HISTORY:
;       $Id: textoidl.pro,v 1.9 2007-08-28 17:04:45 blanton Exp $
;       $Log: not supported by cvs2svn $
;       Revision 1.8  2004/07/19 16:48:04  schlegel
;       Replace calls to the obsolete function RSTRPOS() with STRPOS(/REVERSE_SEARCH).
;
;       Revision 1.7  2000/11/20 02:27:24  dfink
;       Added \AA option for Anstroms
;
;       Revision 1.7  2000/11/19 18:25:00  dfink / Jonathan Swift
;       Added \AA option for Anstroms
;
;       Revision 1.6  2000/07/12 14:31:33  hogg
;       fixed another ()/[] bug.
;
;       Revision 1.5  2000/06/15 18:21:23  hogg
;       fixed tiny () -> [] bugs
;
;       Revision 1.4  2000/06/03 17:12:42  dfink
;       Fink's new textoidl - all procedures in one file; avoid name conflict
;
;       Revision 1.3  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.2  1996/05/09 00:22:17  mcraig
;       Sped up significantly by using str_sep to handle the translation.  No longer
;       relies on routines fromother user libraries.
;
;       Revision 1.1  1996/01/31 18:47:37  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: not supported by cvs2svn $
;
; COPYRIGHT:
;  Copyright (C) 1996 The Regents of the University of California, All
;  Rights Reserved.  Written by Matthew W. Craig.
;  See the file COPYRIGHT for restrictions on distrubting this code.
;  This code comes with absolutely NO warranty; see DISCLAIMER for details.
;-
;
FUNCTION textoidl_strtrans, InputString, from, to, ned,  $
                   HELP=Help

; Bomb out to caller if error.
    On_error, 2

; Offer help if we don't have at least InputString, from, and to, or
; if the user asks for it.
    IF (n_params() LT 3) OR keyword_set(help) THEN BEGIN
        offset = '   '
        print, offset+'Translate all occurences of one substring to another.'
        print, offset+'new = textoidl_strtrans(oldstr,from,to,ned)'
        print, offset+'Inputs:'
        print, offset+offset+'oldstr -- string on which to operate.              in'
        print, offset+offset+'          May be an array.'
        print, offset+offset+'from   -- substrings to be translated. May be      in'
        print, offset+offset+'          an array.'
        print, offset+offset+'to     -- what strings in from should be           in'
        print, offset+offset+'          translated to. May be an array.'
        print, offset+'Outputs:'
        print, offset+offset+'new    -- Translated string. Array if oldstr is    out'
        print, offset+offset+'          an array.'
        print, offset+offset+'ned    -- number of substitutions performed in     out'
        print, offset+offset+'          oldstr.  Array if oldstr is an array.'
        print, offset+'Notes:'
        print, offset+offset+'- Any of old, from, and to can be arrays. ' 
        print, offset+offset+'- from and to must have the same number of elements.'
        return, -1
    ENDIF 
    
    strn = InputString

;  Check that From/To have same number of elements.  RETURN if they don't.
    NFrom = n_elements(from)
    NTo = n_elements(to)
    IF (NFrom EQ 0) OR (NTo EQ 0) THEN return, strn
    IF NFrom NE NTo THEN BEGIN
        print,'Error: Number of elements in from/to unequal'
        return,-1
    ENDIF

;  Make sure there are no null strings in From.  RETURN if there are.   
    FromLen = strlen(From)
    IF (total(FromLen EQ 0) GT 0) THEN BEGIN
        print, 'Error: elements of From must have nonzero length.'
        return, -1
    ENDIF 

    NStrings = n_elements(strn)
    ned = lonarr(NStrings)
    tmpned = 0L

; Say strn='a#b#c', from='#' and to='@'.  Then the approach here is to
; first split strn at all occurances of '#', then recombine the pieces
; with '@' inserted instead.  Do this for all elements of strn, and
; all elements of from.
    FOR i = 0L, NStrings-1 DO BEGIN
        ned[i] = 0L
        FOR j=0L, NFrom-1 DO BEGIN
            SepStr = str_sep(strn[i], from[j])
            NSubs = n_elements(SepStr) - 1
            strn[i] = SepStr[0]
            FOR k=1L, NSubs DO strn[i] = strn[i] + To[j] + SepStr[k]
            ned[i] =  ned[i] + NSubs
        ENDFOR 
    ENDFOR

    return, strn
END 



;
;+
; NAME:
;       TEXTOIDL_TABLE
; PURPOSE:
;       Returns a translation table from TeX to IDL.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       table = textoidl_table()
; INPUTS:
;       None
; KEYWORD PARAMETERS:
;       /POSTSCRIPT -- If set, return postscript translation
;                      table rather than vector fonts table.
;                      Default is translations for vector
;                      fonts.
;       /HELP       -- Print help and exit.
; OUTPUTS:
;       table -- a 2D text array.  table(0,*) contains          out
;                the words to be translated away, table(1,*)
;                contains the words to translate them to.   
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
;       To find out what TeX sequences are available, look at
;       table(0,*).
; EXAMPLE:
; MODIFICATION HISTORY:
;       $Id: textoidl.pro,v 1.9 2007-08-28 17:04:45 blanton Exp $
;       $Log: not supported by cvs2svn $
;       Revision 1.8  2004/07/19 16:48:04  schlegel
;       Replace calls to the obsolete function RSTRPOS() with STRPOS(/REVERSE_SEARCH).
;
;       Revision 1.7  2000/11/20 02:27:24  dfink
;       Added \AA option for Anstroms
;
;       Revision 1.6  2000/07/12 14:31:33  hogg
;       fixed another ()/[] bug.
;
;       Revision 1.5  2000/06/15 18:21:23  hogg
;       fixed tiny () -> [] bugs
;
;       Revision 1.4  2000/06/03 17:12:42  dfink
;       Fink's new textoidl - all procedures in one file; avoid name conflict
;
;       Revision 1.7  1996/07/22 23:56:08  mcraig
;       Added \vartheta.
;
;       Revision 1.6  1996/07/12 21:31:42  mcraig
;       Fixed \varphi in vector font, added \circ.
;
;       Revision 1.5  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.4  1996/05/09 00:22:17  mcraig
;       Added command to return to previous font after switching to Greek or
;       symbol font.
;
;       Revision 1.3  1996/02/08 19:49:35  mcraig
;       Removed control sequence \perp because the postscript code for it is '^'.
;
;       Revision 1.2  1996/02/08 18:53:38  mcraig
;       Added translations for PostScript fonts, and added several new TeX
;       control sequences.
;
;       Revision 1.1  1996/01/31 18:47:37  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: not supported by cvs2svn $
;
; COPYRIGHT:
;  Copyright (C) 1996 The Regents of the University of California, All
;  Rights Reserved.  Written by Matthew W. Craig.
;  See the file COPYRIGHT for restrictions on distrubting this code.
;  This code comes with absolutely NO warranty; see DISCLAIMER for details.
;-
;
FUNCTION textoidl_table, POSTSCRIPT=ps, VECTOR=vec,  HELP=Help

; Return to caller if error.
    On_error, 2

; Print help if necessary.
    IF keyword_set(Help)  THEN BEGIN
        offset = '   '
        print, offset+'Returns a translation table from TeX to IDL.'
        print, offset+'table = textoidl_table()'
        print, offset+'Keywords:'
        print, offset+offset+'/POSTSCRIPT -- If set, return postscript translation'
        print, offset+offset+'               table rather than vector fonts table.'
        print, offset+offset+'               Default is translations for vector'
        print, offset+offset+'               fonts.'
        print, offset+offset+'/HELP       -- Print help and exit.'
        print, offset+'Outputs:'
        print, offset+offset+'table -- a 2D text array.  table(0,*) contains          out'
        print, offset+offset+'         the words to be translated away, table(1,*)'
        print, offset+offset+'         contains the words to translate them to.'
        print, offset+'Notes:'
        print, offset+offset+'To find out what TeX sequences are available, look at'
        print, offset+offset+'table(0,*).'
    ENDIF 

    VECFONT=1                   ; index of vector font in translation table
    PSFONT=2			; index of postscript font in trans table
    IF keyword_set(ps) THEN FontSelection=PSFONT ELSE FontSelection=VECFONT

;  Set IDL font sequence needed to switch to Greek letters.
    GreekFont = strarr(3)
    GreekFont(VECFONT) = '!7'
    GreekFont(PSFONT) = '!M'

;  Set IDL font sequence needed to switch to special symbol font.
    SymbolFont = strarr(3)
    SymbolFont(VECFONT) = '!M'
    SymbolFont(PSFONT) = '!M'

;  Set IDL font sequence needed to switch back to initial font.
    PreviousFont = strarr(3)
    PreviousFont(VECFONT) = '!X'
    PreviousFont(PSFONT) = '!X'

;  Set IDL font sequence for roman fonts.
    RomanFont = strarr(3)
    RomanFont(VECFONT) = '!3'
    RomanFont(PSFONT) = '!M'

;lowercase Greek -- 
;    Note there is some trickery involved in getting \varphi
;    to work in the vector fonts, because it is actually
;    a member of the symbol font set, not the Greek font
;    set.  Go figure.  Solution is just to make the vector
;    character a switch to symbol, the proper character from
;    that font, and a switch back out of symbol.  Same comment holds
;    for \vartheta.
;;        TeX SEQUENCE       VECTOR       POSTSCRIPT
    LowercaseGreek = [$
	[ '\alpha',     	'a'     ,     'a'     ],$
	[ '\beta',      	'b'     ,     'b'     ],$
	[ '\gamma',     	'c'     ,     'g'     ],$
	[ '\delta',     	'd'     ,     'd'     ],$
	[ '\epsilon',   	'e'     ,     'e'     ],$
	[ '\zeta',      	'f'     ,     'z'     ],$
	[ '\eta', 		'g'	,     'h'     ],$
	[ '\theta', 		'h'	,     'q'     ],$
	[ '\iota', 		'i'	,     'i'     ],$
	[ '\kappa', 		'j'	,     'k'     ],$
	[ '\lambda', 		'k'	,     'l'     ],$
	[ '\mu', 		'l'	,     'm'     ],$
	[ '\nu', 		'm'	,     'n'     ],$
	[ '\xi', 		'n'	,  '!S !Rx'   ],$
	[ '\pi', 		'p'	,     'p'     ],$
	[ '\rho', 		'q'	,     'r'     ],$
	[ '\sigma', 		'r'	,     's'     ],$
	[ '\tau', 		's'	,     't'     ],$
	[ '\upsilon', 		't'	,     'u'     ],$
	[ '\phi', 		'u'	,     'f'     ],$
	[ '\chi', 		'v'	,     'c'     ],$
	[ '\psi', 		'w'	,     'y'     ],$
	[ '\omega', 		'x'	,     'w'     ],$
	[ '\varpi', 		'p'	,     'v'     ],$
	[ '\varepsilon', 	'e'     ,     'e'     ],$
	[ '\varphi', 	$
            SymbolFont(VECFONT)+'P'+PreviousFont(VECFONT) $
                                        ,     'j'     ],$
	[ '\vartheta', 	$
            SymbolFont(VECFONT)+'t'+PreviousFont(VECFONT) $
                                        ,     'J'     ]$
	                     ]
;Uppercase Greek -- 
;;        TeX SEQUENCE        VECTOR          POSTSCRIPT
    UppercaseGreek = [$
	[ '\Gamma', 		'C'   ,	        'G'         ],$
	[ '\Delta', 		'D'   ,	        'D'         ],$
	[ '\Theta', 		'H'   ,	        'Q'         ],$
	[ '\Lambda', 		'K'   ,	        'L'         ],$
	[ '\Xi', 		'N'   ,      '!S !RX'       ],$
	[ '\Pi', 		'P'   ,	        'P'         ],$
	[ '\Sigma', 		'R'   ,	        'S'         ],$
	[ '\Upsilon', 		'T'   ,  string(byte(161))  ],$
	[ '\Phi', 		'U'   ,	        'F'         ],$
	[ '\Psi', 		'W'   ,	        'Y'         ],$
	[ '\Omega', 		'X'   ,	        'W'         ]$
			   ]
;Roman --
;
;;         TeX SEQUENCE        VECTOR          POSTSCRIPT
    Roman = [$
	[ '\AA', ''+STRING("305B)+''   , ''+STRING("305B)+''   ]$
			   ]

;Special symbols -- 
;  NOTES -- You must leave \infty before \in in the translatation
;           table to avoid having the \in part of \infty translated
;           away. 
;           
;           DO NOT blindly add the control sequence \perp.  Its
;           PostScript code is '^', which leads to thing being
;           interpreted as superscripts which shouldn't be.
;
;;        TeX SEQUENCE        VECTOR          POSTSCRIPT
    Symbols = [$
	[ '\aleph', 		'@'   ,	 string(byte(192))  ],$
	[ '\ast', 		'*'   ,	        '*'         ],$
	[ '\cap', 		'3'   ,	 string(byte(199))  ],$
	[ '\cdot', 		'.'   ,	 string(byte(215))  ],$
	[ '\odot', 		'n'   ,	 '!SO!R!N!I ' + string(183b) + '!X'  ],$
	[ '\cup', 		'1'   ,	 string(byte(200))  ],$
	[ '\exists', 		'E'   ,	        '$'         ],$
	[ '\infty', 		'$'   ,	 string(byte(165))  ],$
	[ '\in', 		'e'   ,	 string(byte(206))  ],$
	[ '\equiv', 		':'   ,	 string(byte(186))  ],$
	[ '\pm', 		'+'   ,	 string(byte(177))  ],$
	[ '\div', 		'/'   ,	 string(byte(184))  ],$
	[ '\subset', 		'0'   ,	 string(byte(204))  ],$
	[ '\superset', 		'2'   ,	 string(byte(201))  ],$
	[ '\leftarrow', 	'4'   ,	 string(byte(172))  ],$
	[ '\downarrow', 	'5'   ,	 string(byte(175))  ],$
	[ '\rightarrow', 	'6'   ,	 string(byte(174))  ],$
	[ '\uparrow', 		'7'   ,	 string(byte(173))  ],$
	[ '\neq', 		'='   ,	 string(byte(185))  ],$
	[ '\propto', 		'?'   ,	 string(byte(181))  ],$
	[ '\sim', 		'A'   ,	 string(byte(126))  ],$
	[ '\partial', 		'D'   ,	 string(byte(182))  ],$
	[ '\nabla', 		'G'   ,	 string(byte(209))  ],$
	[ '\angle', 		'a'   ,	 string(byte(208))  ],$
	[ '\times', 		'X'   ,	 string(byte(180))  ],$
	[ '\geq', 		'b'   ,	 string(byte(179))  ],$
	[ '\leq', 		'l'   ,	 string(byte(163))  ],$
	[ "\'", 		"'"   ,	 string(byte(162))  ],$
	[ '\prime', 		"'"   ,	 string(byte(162))  ],$
	[ '\circ', 		"%"   ,	 string(byte(176))  ]$
	                  ]
    LowercaseGreek(1,*) = $
      GreekFont(FontSelection) $
      + LowercaseGreek(FontSelection,*) $
      + PreviousFont(FontSelection)
    UppercaseGreek(1,*) = $
      GreekFont(FontSelection) +$
      UppercaseGreek(FontSelection,*) $
      + PreviousFont(FontSelection)
    Symbols(1,*) = $
      SymbolFont(FontSelection) $
      + Symbols(FontSelection,*) $
      + PreviousFont(FontSelection)
    Roman(1,*) = $
      RomanFont(FontSelection) $
      + Roman(FontSelection,*) $
      + PreviousFont(FontSelection)



    TranslationTable = [[LowercaseGreek],[UppercaseGreek],[Symbols],[Roman]]
    return,TranslationTable(0:1,*)

END 



;
;+
; NAME:
;       TEXTOIDL_STRTOK
; PURPOSE:
;       Retrieve portion of string up to token.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       new = textoidl_strtok( old, token )
; INPUTS:
;       old   -- String to be split.  Contains text after    in, out
;                token on output.
;       token -- Token to use in splitting old.              in
; KEYWORD PARAMETERS:
;       /TRIM -- set to remove leading blanks from old 
;                before returning.
;       /HELP -- print useful message and exit.
; OUTPUTS:
;       new   -- portion of string up to token.              out
;       old   -- portion of old after token.                 out, in
; COMMON BLOCKS:
; SIDE EFFECTS:
;       Input parameter old is modified.
; NOTES:
;       Token may be one or more characters.
;       If token is not found, returns old and sets old to ''.
; EXAMPLE:
;       If old is 'foo44 bar', then textoidl_strtok( old, '44' ) would return
;       'foo', and upon return, old will be left with ' bar'.  If /TRIM
;       were set, old would be 'bar' on return.
;
;       If old='xyz', then new=textoidl_strtok(old,'a') would return with
;       new='xyz' and old=''.
; THANKS:
;       To D. Linder who wrote GETTOK, part of the goddard library,
;       upon which this is based.
; MODIFICATION HISTORY:
;       $Id: textoidl.pro,v 1.9 2007-08-28 17:04:45 blanton Exp $
;       $Log: not supported by cvs2svn $
;       Revision 1.8  2004/07/19 16:48:04  schlegel
;       Replace calls to the obsolete function RSTRPOS() with STRPOS(/REVERSE_SEARCH).
;
;       Revision 1.7  2000/11/20 02:27:24  dfink
;       Added \AA option for Anstroms
;
;       Revision 1.6  2000/07/12 14:31:33  hogg
;       fixed another ()/[] bug.
;
;       Revision 1.5  2000/06/15 18:21:23  hogg
;       fixed tiny () -> [] bugs
;
;       Revision 1.4  2000/06/03 17:12:42  dfink
;       Fink's new textoidl - all procedures in one file; avoid name conflict
;
;       Revision 1.3  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.2  1996/05/09 00:22:17  mcraig
;       Added built in help.
;
;       Revision 1.1  1996/01/31 18:47:37  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: not supported by cvs2svn $
;
; COPYRIGHT:
;  Copyright (C) 1996 The Regents of the University of California, All
;  Rights Reserved.  Written by Matthew W. Craig.
;  See the file COPYRIGHT for restrictions on distrubting this code.
;  This code comes with absolutely NO warranty; see DISCLAIMER for details.
;-
FUNCTION Textoidl_Strtok, string, token, $
                 TRIM=trim, HELP=Help

; Back to the caller if error occurs.
    On_error, 2

    IF (n_params() NE 2) OR keyword_set(Help) THEN BEGIN 
        offset = '   '
        print, offset+'Retrieve portion of string up to token.'
        print, offset+'new = textoidl_strtok( old, token )'
        print, offset+'Inputs:'
        print, offset+offset+'old   -- String to be split.  Contains text after    in, out'
        print, offset+offset+'         token on output.'
        print, offset+offset+'token -- Token to use in splitting old.              in'
        print, offset+'Keywords:'
        print, offset+offset+'/TRIM -- set to remove leading blanks from old '
        print, offset+offset+'         before returning.'
        print, offset+offset+'/HELP -- print useful message and exit.'
        print, offset+'Outputs:'
        print, offset+offset+'new   -- portion of string up to token.              out'
        print, offset+offset+'old   -- portion of old after token.                 out, in'
        print, offset+'Side effects:'
        print, offset+offset+'Input parameter old is modified.'
        print, offset+'Notes:'
        print, offset+offset+'Token may be one or more characters.'
        print, offset+offset+"If token is not found, returns old and sets old to ''."
        print, offset+'Examples:'
        print, offset+offset+"If old is 'foo44 bar', then textoidl_strtok( old, '44' ) would return'"
        print, offset+offset+"  'foo', and upon return, old will be left with ' bar'.  If /TRIM"
        print, offset+offset+"  were set, old would be 'bar' on return."
;'
        print, offset+offset+"If old='xyz', then new=textoidl_strtok(old,'a') would return with"
        print, offset+offset+"  new='xyz' and old=''."
        return, -1
    ENDIF 

    pos = strpos(string, token)

    IF (pos GE 0) THEN BEGIN
        front = strmid(string, 0, pos) 
        string = strmid(string, pos + strlen(token), strlen(string))
        IF keyword_set(trim) THEN string = strtrim(string, 1)
        return, front
    ENDIF
    
    front = string
    string = ''
    return, front
    
END



;
;+
; NAME:
;       TEXTOIDL_NEXTTOK
; PURPOSE:
;       Find the next occurance of any of a set of characters in a
;       string and return the character which occurs next.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       tok = textoidl_nexttok( strn, tokens )
; INPUTS:
;       strn   -- string to be searched for sub/superscripts    in
;       tokens -- string containing characters to be found.     in
; KEYWORD PARAMETERS:
;       POSITION -- Set to a named variable to get position     out
;                   of next token, or -1 if none found.
;       /HELP    -- Print useful message and exit.
; OUTPUTS:
;       tok    -- Contains the character among tokens which     out
;                 occurs next in strn, or null '' if none found.
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
; EXAMPLE:
;       textoidl_nexttok( 'x^2 + N_j^3', '^_', position=pos ) returns '^' and sets
;       pos to 1.
; MODIFICATION HISTORY:
;       $Id: textoidl.pro,v 1.9 2007-08-28 17:04:45 blanton Exp $
;       $Log: not supported by cvs2svn $
;       Revision 1.8  2004/07/19 16:48:04  schlegel
;       Replace calls to the obsolete function RSTRPOS() with STRPOS(/REVERSE_SEARCH).
;
;       Revision 1.7  2000/11/20 02:27:24  dfink
;       Added \AA option for Anstroms
;
;       Revision 1.6  2000/07/12 14:31:33  hogg
;       fixed another ()/[] bug.
;
;       Revision 1.5  2000/06/15 18:21:23  hogg
;       fixed tiny () -> [] bugs
;
;       Revision 1.4  2000/06/03 17:12:42  dfink
;       Fink's new textoidl - all procedures in one file; avoid name conflict
;
;       Revision 1.3  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.2  1996/05/09 00:22:17  mcraig
;       Generalized so that the next occurence of any of a set of characters will
;       be returned.
;
;       Revision 1.1  1996/01/31 18:41:06  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: not supported by cvs2svn $
;
; COPYRIGHT:
;  Copyright (C) 1996 The Regents of the University of California, All
;  Rights Reserved.  Written by Matthew W. Craig.
;  See the file COPYRIGHT for restrictions on distrubting this code.
;  This code comes with absolutely NO warranty; see DISCLAIMER for details.
;-
FUNCTION textoidl_nexttok, strn, tokens, $
                  POSITION=position, $
                  HELP=Help

;  Return to caller on error.
    On_error, 2

;  Help those in need of it.
    IF (n_params() NE 2) OR keyword_set(Help) THEN BEGIN 
        offset = '   '
        print, offset+'Find the next occurance of any of a set of characters in a'
        print, offset+'string and return the character which occurs next.'
; CALLING SEQUENCE:
        print, offset+'tok = textoidl_nexttok( strn, tokens )'
; INPUTS:
        print, offset+'Inputs:'
        print, offset+offset+'strn   -- string to be searched for sub/superscripts    in'
        print, offset+offset+'tokens -- string containing characters to be found.     in'
; KEYWORD PARAMETERS:
        print, offset+'Keywords:'
        print, offset+offset+'POSITION -- Set to a named variable to get position     out'
        print, offset+offset+'            of next token, or -1 if none found.'
        print, offset+offset+'/HELP    -- Print useful message and exit.'
; OUTPUTS:
        print, offset+'Outputs:'
        print, offset+offset+'tok   -- Contains the character among tokens which      out'
        print, offset+offset+"         occurs next in strn, or null '' if none found."
; EXAMPLE:
        print, offset+'Example:'
        print, offset+offset+"textoidl_nexttok( 'x^2 + N_j^3', '^_', position=pos ) returns '^' and sets"
        print, offset+offset+'pos to 1.'
        return, ''
    ENDIF 

    TmpStr = byte(strn)
    TmpTok = byte(tokens)
    NumToks = n_elements(TmpTok) 

    MatchIdx = 0L
    Matches = 0L
    FOR j=0, NumToks-1 DO BEGIN 
        TmpMatch = where(TmpStr EQ TmpTok(j),  TmpCnt)
        IF (TmpCnt GT 0) THEN BEGIN
            MatchIdx = [MatchIdx, Replicate(j, TmpCnt)]
            Matches = [Matches, TmpMatch]
        ENDIF 
    ENDFOR 

    IF n_elements(MatchIdx) EQ 1 THEN BEGIN 
        Position = -1
        return, ''
    ENDIF 

    MatchIdx = MatchIdx(1:*)
    Matches = Matches(1:*)

    SortInd = sort(Matches)

    Position = Matches(SortInd(0))

    Tok = string(TmpTok(MatchIdx(SortInd(0))))
    
    return, Tok
END



;
;+
; NAME:
;       TEXTOIDL_MATCHDELIM
; PURPOSE:
;        Match open/close delimiters in a string.
; CATEGORY:
;        text/strings
; CALLING SEQUENCE:
;        position = textoidl_matchdelim( strn, [openpos])
; INPUTS:
;        strn        -- a string containing an open                 in
;                       delimiter (e.g. '{') in which you 
;                       want to find the matching closing  
;                       delimiter (e.g. '}')
; KEYWORD PARAMETERS:
;        OPEN_DELIM  -- A single character containing the opening   in
;                       delimiter (e.g. '(').  Default is '{'
;        CLOSE_DELIM -- A single character containing the closing   in
;                       delimiter (e.g. ')').  Default is '}'
; OUTPUTS:
;        position -- returns the position in strn of the            out
;                    closing delimiter, -1 if no closing found.
;        openpos  -- Set to a named variable to receive the         out
;                    position of the first opening delimiter.
;                    Optional.
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
;        - Any pair of (nonidentical) characters can be used as
;          delimiters. 
; EXAMPLE:
;        textoidl_matchdelim('{one{two}}three') returns 9, the character just
;        before 'three'.  
; MODIFICATION HISTORY:
;       $Id: textoidl.pro,v 1.9 2007-08-28 17:04:45 blanton Exp $
;       $Log: not supported by cvs2svn $
;       Revision 1.8  2004/07/19 16:48:04  schlegel
;       Replace calls to the obsolete function RSTRPOS() with STRPOS(/REVERSE_SEARCH).
;
;       Revision 1.7  2000/11/20 02:27:24  dfink
;       Added \AA option for Anstroms
;
;       Revision 1.6  2000/07/12 14:31:33  hogg
;       fixed another ()/[] bug.
;
;       Revision 1.5  2000/06/15 18:21:23  hogg
;       fixed tiny () -> [] bugs
;
;       Revision 1.4  2000/06/03 17:12:42  dfink
;       Fink's new textoidl - all procedures in one file; avoid name conflict
;
;       Revision 1.3  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.2  1996/05/09 00:22:17  mcraig
;       Removed restriction that open delim must be first char.  Added argument
;       to allow for return of position of open delim.
;
;       Revision 1.1  1996/01/31 18:41:06  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: not supported by cvs2svn $
;
; COPYRIGHT:
;  Copyright (C) 1996 The Regents of the University of California, All
;  Rights Reserved.  Written by Matthew W. Craig.
;  See the file COPYRIGHT for restrictions on distrubting this code.
;  This code comes with absolutely NO warranty; see DISCLAIMER for details.
;-
;
FUNCTION Textoidl_Matchdelim, InString, OpenPos, $
                     OPEN_DELIM=OpenDelim, $
                     CLOSE_DELIM=CloseDelim, $
                     HELP=Help

; Return to caller if error.
    On_error, 2

    IF (n_params() LT 1) OR keyword_set(Help) THEN BEGIN
        offset = '   '
        print, offset+'Match open/close delimiters in a string.'
        print, offset+'position = textoidl_matchdelim( strn, [openpos])'
        print, offset+'Inputs:'
        print, offset+offset+'strn        -- a string containing an open                 in'
        print, offset+offset+"               delimiter (e.g. '{') in which you "
        print, offset+offset+'               want to find the matching closing  '
        print, offset+offset+"               delimiter (e.g. '}')"
        print, offset+'Keywords:'
        print, offset+offset+'OPEN_DELIM  -- A single character containing the opening   in'
        print, offset+offset+"               delimiter (e.g. '(').  Default is '{'"
        print, offset+offset+'CLOSE_DELIM -- A single character containing the closing   in'
        print, offset+offset+"               delimiter (e.g. ')').  Default is '}'"
        print, offset+'Outputs:'
        print, offset+offset+'position -- returns the position in strn of the            out'
        print, offset+offset+'            closing delimiter, -1 if no closing found.'
        print, offset+offset+'openpos  -- Set to a named variable to receive the         out'
        print, offset+offset+'            position of the first opening delimiter.'
        print, offset+offset+'            Optional.'
        print, offset+'Example:'
        print, offset+offset+"textoidl_matchdelim('a{one{two}}three') returns 10, the character just"
        print, offset+offset+"  before 'three'.  "
        print, offset+offset+$
          "a=textoidl_matchdelim('aaa[bbb(ccc)]ddd[eee]',f,OP='[',CL=']')"
        print, offset+offset+"  returns a=12 (just before ddd), f=3 "+$
          "(just before bbb)."  
        return, -1
    ENDIF 

; Set default delimiters.
    IF n_elements(OpenDelim) EQ 0 THEN OpenDelim =  '{'
    IF n_elements(CloseDelim) EQ 0 THEN CloseDelim =  '}'

; Make sure InString has more than 1 character.
    length = strlen(InString) 
    IF (length LE 1) THEN return,-1

; Return if no open delimiter
    OpenPos = strpos( InString, OpenDelim )
    IF (OpenPos EQ -1) THEN BEGIN 
        print, 'Error: No opening delimiter'
        return, -1
    ENDIF 
    
; Convert strings to array of integers to speed processing.
    OpenDelim = fix((byte(OpenDelim))(0))
    CloseDelim = fix((byte(CloseDelim))(0))
    TmpStr = fix(byte(strmid( InString, OpenPos, length)))
; Leave the -1* in here.  This forces conversion from BYTE to INTEGER,
; necessary because there are no negative BYTEs.
    TmpStr = (TmpStr EQ OpenDelim) $
              -1*(TmpStr EQ CloseDelim)
    length = n_elements(TmpStr) 

; Initialize count of number of delimiters.  We've found one, the
; first opener.
    BraceCnt = 1
    i=0
    WHILE (BraceCnt GT 0) AND (i LT length-1) DO BEGIN 
        i = i+1
        BraceCnt = BraceCnt + TmpStr(i)
    ENDWHILE 
    
    i = i + OpenPos
    IF (BraceCnt GT 0) THEN i = -1
    return, i
END

	

;
;+
; NAME:
;       TEXTOIDL_STRCNT
; PURPOSE:
;       Count number of occurrences of a substring in a string.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       num = textoidl_strcnt(strn, substring, [pos])
; INPUTS:
;       string    -- The string in which to count occurences.     in
;       substring -- The substring to count occurrences of.       in
;       pos       -- the position at which to begin the search.   [in]
;                    If not supplied, start at beginning of
;                    string.
; KEYWORD PARAMETERS:
;       /HELP     -- Print useful message and return.
; OUTPUTS:
;       num       -- Number of occurances of substring in string. out
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
;       Overlapping occurances are not counted separately.  For
;       example, counting occurances of 'bb' in 'blah bbb' returns one
;       occurance. 
; EXAMPLE:
; MODIFICATION HISTORY:
;       $Id: textoidl.pro,v 1.9 2007-08-28 17:04:45 blanton Exp $
;       $Log: not supported by cvs2svn $
;       Revision 1.8  2004/07/19 16:48:04  schlegel
;       Replace calls to the obsolete function RSTRPOS() with STRPOS(/REVERSE_SEARCH).
;
;       Revision 1.7  2000/11/20 02:27:24  dfink
;       Added \AA option for Anstroms
;
;       Revision 1.6  2000/07/12 14:31:33  hogg
;       fixed another ()/[] bug.
;
;       Revision 1.5  2000/06/15 18:21:23  hogg
;       fixed tiny () -> [] bugs
;
;       Revision 1.4  2000/06/03 17:12:42  dfink
;       Fink's new textoidl - all procedures in one file; avoid name conflict
;
;       Revision 1.3  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.2  1996/05/09 00:22:17  mcraig
;       Added fast processing using BYTE arrays if we are counting occurences of
;       a single character.  Added error handling.
;
;       Revision 1.1  1996/01/31 18:47:37  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: not supported by cvs2svn $
;
; COPYRIGHT:
;  Copyright (C) 1996 The Regents of the University of California, All
;  Rights Reserved.  Written by Matthew W. Craig.
;  See the file COPYRIGHT for restrictions on distrubting this code.
;  This code comes with absolutely NO warranty; see DISCLAIMER for details.
;-
FUNCTION Textoidl_Strcnt, strn, substrn, startpos, $
                 HELP=Help

; Return to caller if error.
    On_error, 2

; Help user, if needed.
    IF (n_params() LT 2) OR keyword_set(Help) THEN BEGIN
        offset = '   '
        print, offset+'Count number of occurrences of a substring in a string.'
        print, offset+'num = textoidl_strcnt(strn, substring, [pos])'
        print, offset+'Inputs:'
        print,offset+offset+'string    -- The string in which to count occurences.     in'
        print,offset+offset+'substring -- The substring to count occurrences of.       in'
        print,offset+offset+'pos       -- the position at which to begin the search.   [in]'
        print,offset+offset+'             If not supplied, start at beginning of'
        print,offset+offset+'             string.'
        print, offset+'Keywords:'
        print,offset+offset+'/HELP     -- Print useful message and return.'
        print, offset+'Outputs:'
        print,offset+offset+'num       -- Number of occurances of substring in string. out'
        return, -1
    ENDIF 

    IF n_params() EQ 2 THEN startpos = 0

;  return if we weren't really given a substring to search for. . .
    IF strlen(substrn) EQ 0 THEN BEGIN 
        print, "Error: Can't count occurances of null string."
        return, -1
    ENDIF 

; . . .or if we were told to start at the end of the string.
    tmpstrn = strmid(strn, startpos, strlen(strn))
    IF strlen(tmpstrn) EQ 0 THEN return, 0

; If looking for occurences of single character, process using BYTE
; array.
    IF strlen(substrn) EQ 1 THEN BEGIN
        tmpstrn = byte(TmpStrn)
        count = n_elements(where(TmpStrn EQ (byte(substrn))(0))) 
    ENDIF ELSE BEGIN 
        count = 0L
        pos = strpos(tmpstrn, substrn, /reverse_search)
        WHILE pos GE 0 DO BEGIN
            count = count + 1
            pos = strpos(tmpstrn, substrn, pos, /reverse_search)
        ENDWHILE
    ENDELSE 

    return, count
END 



;  NOTE to future maintainers:
;   Make sure textoidl_sub_sup stays before textoidl_trans_sub_super.  At least
;   for now, when IDL encounters a function and automatically compiles
;   it, it only compiles the functions in the file up to the named
;   function.  So even if textoidl_sub_sup was declared with
;   FORWARD_FUNCTION in textoidl_trans_sub_super, it would not properly
;   compile. 
;
;+
; SPECIAL NOTE: 
;       The file textoidl_trans_sub_super.pro contains two functions,
;       textoidl_trans_sub_super, and textoidl_sub_sup.  The former is the
;       generic routine for processing TeX sub/superscripts, the
;       latter is used only by textoidl_trans_sub_super and has no general
;       utility.  Hence it lives here.  You will see documentation for
;       textoidl_trans_sub_super second if you use DOC_LIBRARY.
;-
;
;
;+
; NAME:
;       TEXTOIDL_SUB_SUP
; PURPOSE:
;       Return the proper IDL font positioning command for TeX
;       sub/superscripts. 
; CATEGORY:
;       TeXtoIDL
; CALLING SEQUENCE:
;       fnt = textoidl_sub_sup( strn )
; INPUTS:
;       strn      -- Either '^' or '_', the TeX super/subscript       in
;                    characters
; KEYWORD PARAMETERS:
;       /FORCE_UD -- Set this to use !U/!D instead of !E/!I for
;                    sub/superscripts .
;       /HELP     -- Set to print useful message and exit.
; OUTPUTS:
;       fnt       -- Either '!U' or !E' for superscripts,             out
;                    or '!D' or '!I' for subscripts.
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
;       Used only by textoidl_trans_sub_super.  Should be kept in same
;       file. 
; EXAMPLE:
; MODIFICATION HISTORY:
;       $Id: textoidl.pro,v 1.9 2007-08-28 17:04:45 blanton Exp $
;       $Log: not supported by cvs2svn $
;       Revision 1.8  2004/07/19 16:48:04  schlegel
;       Replace calls to the obsolete function RSTRPOS() with STRPOS(/REVERSE_SEARCH).
;
;       Revision 1.7  2000/11/20 02:27:24  dfink
;       Added \AA option for Anstroms
;
;       Revision 1.6  2000/07/12 14:31:33  hogg
;       fixed another ()/[] bug.
;
;       Revision 1.5  2000/06/15 18:21:23  hogg
;       fixed tiny () -> [] bugs
;
;       Revision 1.4  2000/06/03 17:12:42  dfink
;       Fink's new textoidl - all procedures in one file; avoid name conflict
;
;       Revision 1.4  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.3  1996/05/09 00:22:17  mcraig
;       Changed some function calls to reflect changes in those functions, moved
;       some code out of the main loop that didn't need to be there, added
;       documentation.
;
;       Revision 1.1  1996/01/31 18:47:37  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: not supported by cvs2svn $
; COPYRIGHT:
;  Copyright (C) 1996 The Regents of the University of California, All
;  Rights Reserved.  Written by Matthew W. Craig.
;  See the file COPYRIGHT for restrictions on distrubting this code.
;  This code comes with absolutely NO warranty; see DISCLAIMER for details.
;-
FUNCTION textoidl_sub_sup, token,  FORCE_UD = force_ud

; provide help if needed.
    IF (n_params() NE 1) OR keyword_set(Help) THEN BEGIN
        offset = '   '
        print, offset+'Return the proper IDL font positioning command for TeX'
        print, offset+'sub/superscripts. '
        print, offset+'fnt = textoidl_sub_sup( strn )'
        print, offset+'Inputs:'
        print, offset+offset+"strn      -- Either '^' or '_', the TeX super/subscript       in"
        print, offset+offset+'             characters'
        print, offset+'Keywords:'
        print, offset+offset+'/FORCE_UD -- Set this to use !U/!D instead of !E/!I for'
        print, offset+offset+'             sub/superscripts .'
        print, offset+offset+'/HELP     -- Set to print useful message and exit.'
        print, offset+'Outputs:'
        print, offset+offset+"fnt       -- Either '!U' or !E' for superscripts,             out"
        print, offset+offset+"             or '!D' or '!I' for subscripts."
        return, -1
    ENDIF 

    IF keyword_set(force_ud) THEN BEGIN 
        IF (token EQ '^') THEN return, '!U' 
        IF (token EQ '_') THEN return, '!D'
        return, ''
    ENDIF ELSE BEGIN
        IF (token EQ '^') THEN return, '!E' 
        IF (token EQ '_') THEN return, '!I'
        return, ''
    ENDELSE
    
END



;
;+
; NAME:
;       TEXTOIDL_TRANS_SUB_SUPER
; PURPOSE:
;       Translate TeX sub/superscripts to IDL sub/superscripts.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       new = textoidl_trans_sub_super( old )
; INPUTS:
;       old       -- string to be translated from TeX to IDL.   in
; KEYWORD PARAMETERS:
;       /RECURSED -- set if this function is being called 
;                    recursively.                  
;       /HELP     -- Set to print useful message and exit.
; OUTPUTS:
;       new       -- string old converted from TeX to IDL       out
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
;       - For best results, when both a sub and superscript are used,
;         place the shorter of the two first (e.g. 'N^{a}_{bbbb}' is
;         better than 'N_{bbbb}^{a}').
;       - Single character sub/super scripts do not need to be
;         protected by braces.
;       - Sub/superscripts may be nested (e.g. 'N^{N_1^N}').
; EXAMPLE:
;       out = textoidl_trans_sub_super( 'N^2_{big}' )
;       Then out='N!U2!N!Dbig!N' which looks like it should on the
;       display. 
; LIBRARY FUNCTIONS CALLED:
;       textoidl_strtok      -- Text/string (mcraig)
;       textoidl_sub_sup -- contained in this file
; MODIFICATION HISTORY:
;       $Id: textoidl.pro,v 1.9 2007-08-28 17:04:45 blanton Exp $
;       $Log: not supported by cvs2svn $
;       Revision 1.8  2004/07/19 16:48:04  schlegel
;       Replace calls to the obsolete function RSTRPOS() with STRPOS(/REVERSE_SEARCH).
;
;       Revision 1.7  2000/11/20 02:27:24  dfink
;       Added \AA option for Anstroms
;
;       Revision 1.6  2000/07/12 14:31:33  hogg
;       fixed another ()/[] bug.
;
;       Revision 1.5  2000/06/15 18:21:23  hogg
;       fixed tiny () -> [] bugs
;
;       Revision 1.4  2000/06/03 17:12:42  dfink
;       Fink's new textoidl - all procedures in one file; avoid name conflict
;
;       Revision 1.4  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.3  1996/05/09 00:22:17  mcraig
;       Changed some function calls to reflect changes in those functions, moved
;       some code out of the main loop that didn't need to be there, added
;       documentation.
;
;       Revision 1.2  1996/02/08 18:54:20  mcraig
;       Changed default sub/superscript size to be !D/!U rather than !I/!E to
;       improve readability of plat annotations.
;
;       Revision 1.1  1996/01/31 18:47:37  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: not supported by cvs2svn $
;
; COPYRIGHT:
;  Copyright (C) 1996 The Regents of the University of California, All
;  Rights Reserved.  Written by Matthew W. Craig.
;  See the file COPYRIGHT for restrictions on distrubting this code.
;  This code comes with absolutely NO warranty; see DISCLAIMER for details.
;-
FUNCTION textoidl_trans_sub_super, InputString, $
                                   RECURSED=recursed, $
                                   HELP=Help

; allow compilation of recursive function 
  forward_function textoidl_trans_sub_super

; Return to caller if error.
    On_error, 2

; Offer help if needed and/or desired
    IF (n_params() NE 1) OR keyword_set(help) THEN BEGIN
        offset = '   '
        print, offset+'Translate TeX sub/superscripts to IDL sub/superscripts.'
        print, offset+'new = textoidl_trans_sub_super( old )'
        print, offset+'Inputs:'
        print, offset+offset+'old       -- string to be translated from TeX to IDL.   in'
        print, offset+'Keywords:'
        print, offset+offset+'/RECURSED -- set if this function is being called '
        print, offset+offset+'             recursively.                  '
        print, offset+offset+'/HELP     -- Set to print useful message and exit.'
        print, offset+'Outputs:'
        print, offset+offset+'new       -- string old converted from TeX to IDL       out'
        print, offset+'Notes:'
        print, offset+offset+'- For best results, when both a sub and superscript are used,'
        print, offset+offset+"  place the shorter of the two first (e.g. 'N^{a}_{bbbb}' is"
        print, offset+offset+"  better than 'N_{bbbb}^{a}')."
        print, offset+offset+'- Single character sub/super scripts do not need to be'
        print, offset+offset+'  protected by braces.'
        print, offset+offset+"- Sub/superscripts may be nested (e.g. 'N^{N_1^N}')."
        return, -1
    ENDIF 

;  To allow for nested scripts, use !E/!I instead of !U/!D for scripts
;  when called recursively.
    IF (NOT keyword_set(recursed)) THEN $
      ud = 1 $
    ELSE $
      ud = 0

;  Return to the normal level after making sub/superscript unless we
;  are recursed, which indicates we are processing a nested script.
    IF keyword_set(recursed) THEN fontRestore = '' ELSE fontRestore = '!N'

;  Initialize vars for processing scripts.
    SpcByte = (byte(' '))(0)    ;We need the BYTE value for a space below.
    strn = InputString
    pos = 0
    StorePos = ''
    RecallPos = ''
    OldToken =  ''
    LenLastScript = 0

; Grab next sub/superscript.  Token will be either '^' or '_'.
; RETURN if no scripts.
    Token = textoidl_nexttok(strn,  '^_', pos = pos)
    if pos EQ -1 then return, InputString ;nothing to process.

    FntChange =  textoidl_sub_sup(Token)

; Our approach will be to grab the input string up to the next '^' or
; '_', then process the script we've found.
    NewString=textoidl_strtok(strn,Token)

    WHILE  strlen(strn) GT  0 DO  BEGIN
;  Grab first char of sub/superscript.
        Script = strmid(strn, 0, 1)
        EndOfScript = 0         ;Position of end of this script.
        IF (Script EQ '{') THEN BEGIN   ; Scripts of more than 1 char.
            EndOfScript = textoidl_matchdelim(strn)      
            Script = textoidl_trans_sub_super(strmid(strn, 1, EndOfScript-1), /recursed )

        ENDIF 
;     Grab rest of string _after_ the end of the script.        
        strn = strmid(strn, EndOfScript+1, $
                      strlen(strn)-EndOfScript-1)

;     Find the next script and prepare for processing it.
        FntChange = textoidl_sub_sup(Token, FORCE_UD = ud)
        OldToken = Token
        Token = textoidl_nexttok(strn, '^_', POS = pos)

;     If the input is 'n^2_j', we want the '2' to be directly above
;     the 'j', rather than having the 'j' below and to the right of
;     the 2.  In other words, we want the first below, not the second.
;              2               2
;             N               N
;              J                J
;     To accomplish this, we need to save the position at which we
;     begin writing the 2 with a !S, and restore that position with a
;     !R after writing the 2.  The first section in the IF block below
;     handles the 'J' above, the thing after the first script.  We
;     don't care if there is another script following.  We also padd
;     the second script with spaces if it is shorter than the first to
;     make sure that whatever comes out after the scripts starts in
;     the proper place.  The worry is that without the spaces, the
;     input 'N^{looong}_{s} + 1' will end up with the + starting right
;     the 's' ends.
        IF (StorePos EQ '!S') THEN BEGIN
            StorePos = ''
            RecallPos = ''
;     calculate the difference in length between this script and the 
;     previous stacked one, removing font change commands (crudely by
;     guessing that the number of characters this takes is twice the
;     number of exclamation points).  The  + 1 below is a kludge.  I
;     don't know why, but I need one extra space.
            NumSpaces = LenLastScript - (strlen(script) - 2*textoidl_strcnt(Script,'!'))
            NumSpaces = (NumSpaces + 1) > 0
            IF NumSpaces GT 0 THEN $
              Script = Script + string( replicate(SpcByte, NumSpaces) )
        ENDIF ELSE BEGIN
            IF (Token NE OldToken) AND (pos EQ 0) THEN BEGIN
;             The next script immediately folows this one.  Arrange to
;             save the position of the current script so that both begin
;             with the same horizontal position.
                StorePos = '!S'
                RecallPos = '!R'
                LenLastScript = strlen(Script) - 2*textoidl_strcnt(Script,'!')
            ENDIF
        ENDELSE  

;  Continue building the IDL string, adding on our just processed script.
        NewString = NewString + StorePos + FntChange + Script + RecallPos $
          + FontRestore

        IF ( pos NE -1 ) THEN BEGIN     ; more left to process     
            NewString = NewString $
              + textoidl_strtok(strn, Token)   
        ENDIF ELSE BEGIN                ; we are done
            NewString = NewString + strn
            strn = ''
        ENDELSE
    ENDWHILE 
    
    return, NewString
END



;
;+
; NAME:
;       TEXTOIDL
; PURPOSE:
;       Convert a valid TeX string to a valid IDL string for plot labels.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       new = textoidl(old)
; INPUTS:
;       old            -- TeX string to be converted.  Will not be     in
;                         modified.  old may be a string array.
; KEYWORD PARAMETERS:
;       FONT           -- Set to 0 to use hardware font, -1 to use 
;                         vector.  Note that the only hardware font 
;                         supported is PostScript.
;       /TEX_SEQUENCES -- return the available TeX sequences
;       /HELP          -- print out info on use of the function
;                         and exit.
; OUTPUTS:
;       new            -- IDL string corresponding to old.             out
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
;       - Use the procedure SHOWTEX to get a list of the available TeX
;         control sequences.  
;       - The only hardware font for which translation is available is
;         PostScript. 
;       - The only device for which hardware font'
;         translation is available is PostScript.'
;       - The FONT keyword overrides the font selected'
;         by !p.font'
; EXAMPLE:
;       out = TeXtoIDL('\Gamma^2 + 5N_{ed}')
;       The string out may be used in XYOUTS or other IDL text
;       display routines.  It will be an uppercase Gamma, with an
;       exponent of 2, then a plus sign, then an N with the subscript
;       ed.
; MODIFICATION HISTORY:
;       $Id: textoidl.pro,v 1.9 2007-08-28 17:04:45 blanton Exp $
;       $Log: not supported by cvs2svn $
;       Revision 1.8  2004/07/19 16:48:04  schlegel
;       Replace calls to the obsolete function RSTRPOS() with STRPOS(/REVERSE_SEARCH).
;
;       Revision 1.7  2000/11/20 02:27:24  dfink
;       Added \AA option for Anstroms
;
;       Revision 1.6  2000/07/12 14:31:33  hogg
;       fixed another ()/[] bug.
;
;       Revision 1.5  2000/06/15 18:21:23  hogg
;       fixed tiny () -> [] bugs
;
;       Revision 1.4  2000/06/03 17:12:42  dfink
;       Fink's new textoidl - all procedures in one file; avoid name conflict
;
;       Revision 1.4  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.3  1996/05/09 00:22:17  mcraig
;       Added error handling, cleaned up documentation.
;
;       Revision 1.2  1996/02/08 18:52:50  mcraig
;       Added ability to use hardware fonts for PostScript device.
;
;       Revision 1.1  1996/01/31 18:47:37  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: not supported by cvs2svn $
;
; COPYRIGHT:
;  Copyright (C) 1996 The Regents of the University of California, All
;  Rights Reserved.  Written by Matthew W. Craig.
;  See the file COPYRIGHT for restrictions on distrubting this code.
;  This code comes with absolutely NO warranty; see DISCLAIMER for details.
;-
;
FUNCTION Textoidl, InputString, $
                   FONT=fnt, $
                   HELP=hlp, $
                   TEX_SEQUENCES=tex_seq

;  Return to caller if there is an error.
    On_error, 2
;  We begin by deciding on the font.  PostScript = 0 means use vector.
    PostScript = 0
    IF n_elements(fnt) EQ 0 THEN BEGIN     ; get font from !p.font
        IF !p.font NE -1 THEN BEGIN        ; User wants hardware font.
            PostScript=1
        ENDIF
    ENDIF ELSE BEGIN                       ; get font from FONT keyword
        IF fnt NE -1 THEN PostScript = 1
    ENDELSE

;  Bomb out if user wants non-PostScript hardware font.
    IF (PostScript EQ 1) AND (!d.name NE 'PS') THEN BEGIN   
                                              ; Device isn't postscript 
                                              ; and user wants hardware
                                              ; font.  Not good.
        print,'Warning: No translation for device: ',!d.name
        return,InputString               
    ENDIF 
    
    IF keyword_set (tex_seq) THEN BEGIN
        table=textoidl_table()
        return,table(0,*)
    ENDIF 

    IF keyword_set(hlp) OR (n_params() EQ 0) THEN BEGIN
        print, '   Convert a TeX string to an IDL string'
        print, '   new = TeXtoIDL(old)'
        print, '     old = TeX string to translate.                 in'
        print, '     new = resulting IDL string.                    out'
        print, '   Keywords:'
        print, '      FONT       set to -1 to translate for vector fonts '
        print, '                 (DEFAULT) .  Set to 0 to translate for'
        print, '                 hardware font.'
        print, '      /TEX_SEQUENCES -- return the available TeX sequences'
        print, '      /HELP      print this message and exit.'
        print, '   NOTES:  '
        print, '      - Use SHOWTEX to obtain a list of the available'
        print, '        TeX control sequences.'
        print, '      - old may be a string array.  If so, new is too.'
        print, '      - The only device for which hardware font'
        print, '        translation is available is PostScript.'
        print, '      - The FONT keyword overrides the font selected'
        print, '        by !p.font'
        return, -1
    ENDIF
    
; PostScript has been set to 1 if PostScript fonts are desired.
    strn = InputString
    table = textoidl_table(POSTSCRIPT=PostScript)
    
;   Greek sub/superscripts need to be protected by putting braces
;   around them if they are unbraced.  This will have the result the
;   it will be difficult to use \ as a sub/superscript.  Get over it.
    strn =  textoidl_strtrans(strn, '^'+table(0, *), '^{'+table(0, *)+'}')
    strn =  textoidl_strtrans(strn, '_'+table(0, *), '_{'+table(0, *)+'}')

;  First we translate Greek letters and the like.  This makes guessing
;  alignment of sub/superscripts easier, as all special characters will then
;  be one character long.
    strn = textoidl_strtrans(strn, table(0, *), table(1, *))

    FOR i = 0L, n_elements(strn)-1 DO $
      strn[i] = textoidl_trans_sub_super(strn[i]) ; Take care of sub/superscripts

    return,strn
END 
