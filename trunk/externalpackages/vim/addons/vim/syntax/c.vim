" Vim syntax file
" Language:	C
" Maintainer:	Bram Moolenaar <Bram@vim.org>
" Last Change:	1999 Dec 02, 2004 Feb 04

" For version 5.x: Clear all syntax items
" For version 6.x: Quit when a syntax file was already loaded
if version < 600
  syntax clear
elseif exists("b:current_syntax")
  finish
endif
"hi clear

syn match       cName		"\<[a-zA-Z_][a-zA-Z_0-9]*\>"
"syn match       cConstant	"\<[A-Z_][A-Z_0-9]\{2,}[A-Za-z_0-9]*\>[^(:]"me=e-1
"syn match       cConstant	"\<[A-Z_][A-Z_0-9]\{2,}[A-Za-z_0-9]*\>$"
"syn match       cConstant	"\<_[_0-9]\{1,}\>[^(:]"
syn match	cFunction	"\<[a-zA-Z_][a-zA-Z_0-9]*\>[^()]*)("me=e-2
"syn match	cFunction	"\<[a-zA-Z_][a-zA-Z_0-9]*\>\s*)("me=e-2
syn match	cFunction	"\<[a-zA-Z_][a-zA-Z_0-9]*\>\s*("me=e-1
syn match	cBraces		"[{}]"

syn keyword cMC			__DI __EI __asm __set_il
syn keyword cMC			__wait_nop __mul __div __mod __mulu __divu __modu
syn keyword cAnsiFunction	MULU_ DIVU_ MODU_ MUL_ DIV_ MOD_
syn keyword cAnsiFunction	main typeof
syn keyword cAnsiFunction	open close read write lseek dup dup2
syn keyword cAnsiFunction	fcntl ioctl
syn keyword cAnsiFunction	wctrans towctrans towupper
syn keyword cAnsiFunction	towlower wctype iswctype
syn keyword cAnsiFunction	iswxdigit iswupper iswspace
syn keyword cAnsiFunction	iswpunct iswprint iswlower
syn keyword cAnsiFunction	iswgraph iswdigit iswcntrl
syn keyword cAnsiFunction	iswalpha iswalnum wcsrtombs
syn keyword cAnsiFunction	mbsrtowcs wcrtomb mbrtowc
syn keyword cAnsiFunction	mbrlen mbsinit wctob
syn keyword cAnsiFunction	btowc wcsfxtime wcsftime
syn keyword cAnsiFunction	wmemset wmemmove wmemcpy
syn keyword cAnsiFunction	wmemcmp wmemchr wcstok
syn keyword cAnsiFunction	wcsstr wcsspn wcsrchr
syn keyword cAnsiFunction	wcspbrk wcslen wcscspn
syn keyword cAnsiFunction	wcschr wcsxfrm wcsncmp
syn keyword cAnsiFunction	wcscoll wcscmp wcsncat
syn keyword cAnsiFunction	wcscat wcsncpy wcscpy
syn keyword cAnsiFunction	wcstoull wcstoul wcstoll
syn keyword cAnsiFunction	wcstol wcstold wcstof
syn keyword cAnsiFunction	wcstod ungetwc putwchar
syn keyword cAnsiFunction	putwc getwchar getwc
syn keyword cAnsiFunction	fwide fputws fputwc
syn keyword cAnsiFunction	fgetws fgetwc wscanf
syn keyword cAnsiFunction	wprintf vwscanf vwprintf
syn keyword cAnsiFunction	vswscanf vswprintf vfwscanf
syn keyword cAnsiFunction	vfwprintf swscanf swprintf
syn keyword cAnsiFunction	fwscanf fwprintf zonetime
syn keyword cAnsiFunction	strfxtime strftime localtime
syn keyword cAnsiFunction	gmtime ctime asctime
syn keyword cAnsiFunction	time mkxtime mktime
syn keyword cAnsiFunction	difftime clock strlen
syn keyword cAnsiFunction	strerror memset strtok
syn keyword cAnsiFunction	strstr strspn strrchr
syn keyword cAnsiFunction	strpbrk strcspn strchr
syn keyword cAnsiFunction	memchr strxfrm strncmp
syn keyword cAnsiFunction	strcoll strcmp memcmp
syn keyword cAnsiFunction	strncat strcat strncpy
syn keyword cAnsiFunction	strcpy memmove memcpy
syn keyword cAnsiFunction	wcstombs mbstowcs wctomb
syn keyword cAnsiFunction	mbtowc mblen lldiv
syn keyword cAnsiFunction	ldiv div llabs
syn keyword cAnsiFunction	labs abs qsort
syn keyword cAnsiFunction	bsearch system getenv
syn keyword cAnsiFunction	exit atexit abort
syn keyword cAnsiFunction	realloc malloc free
syn keyword cAnsiFunction	calloc srand rand
syn keyword cAnsiFunction	strtoull strtoul strtoll
syn keyword cAnsiFunction	strtol strtold strtof
syn keyword cAnsiFunction	strtod atoll atol
syn keyword cAnsiFunction	atoi atof perror
syn keyword cAnsiFunction	ferror feof clearerr
syn keyword cAnsiFunction	rewind ftell fsetpos
syn keyword cAnsiFunction	fseek fgetpos fwrite
syn keyword cAnsiFunction	fread ungetc puts
syn keyword cAnsiFunction	putchar putc gets
syn keyword cAnsiFunction	getchar getc fputs
syn keyword cAnsiFunction	fputc fgets fgetc
syn keyword cAnsiFunction	vsscanf vsprintf vsnprintf
syn keyword cAnsiFunction	vscanf vprintf vfscanf
syn keyword cAnsiFunction	vfprintf sscanf sprintf
syn keyword cAnsiFunction	snprintf scanf printf
syn keyword cAnsiFunction	fscanf fprintf setvbuf
syn keyword cAnsiFunction	setbuf freopen fopen
syn keyword cAnsiFunction	fflush fclose tmpnam
syn keyword cAnsiFunction	tmpfile rename remove
syn keyword cAnsiFunction	offsetof va_start va_end
syn keyword cAnsiFunction	va_copy va_arg raise signal
syn keyword cAnsiFunction	longjmp setjmp isunordered
syn keyword cAnsiFunction	islessgreater islessequal isless
syn keyword cAnsiFunction	isgreaterequal isgreater fmal
syn keyword cAnsiFunction	fmaf fma fminl
syn keyword cAnsiFunction	fminf fmin fmaxl
syn keyword cAnsiFunction	fmaxf fmax fdiml
syn keyword cAnsiFunction	fdimf fdim nextafterxl
syn keyword cAnsiFunction	nextafterxf nextafterx nextafterl
syn keyword cAnsiFunction	nextafterf nextafter nanl
syn keyword cAnsiFunction	nanf nan copysignl
syn keyword cAnsiFunction	copysignf copysign remquol
syn keyword cAnsiFunction	remquof remquo remainderl
syn keyword cAnsiFunction	remainderf remainder fmodl
syn keyword cAnsiFunction	fmodf fmod truncl
syn keyword cAnsiFunction	truncf trunc llroundl
syn keyword cAnsiFunction	llroundf llround lroundl
syn keyword cAnsiFunction	lroundf lround roundl
syn keyword cAnsiFunction	roundf round llrintl
syn keyword cAnsiFunction	llrintf llrint lrintl
syn keyword cAnsiFunction	lrintf lrint rintl
syn keyword cAnsiFunction	rintf rint nearbyintl
syn keyword cAnsiFunction	nearbyintf nearbyint floorl
syn keyword cAnsiFunction	floorf floor ceill
syn keyword cAnsiFunction	ceilf ceil tgammal
syn keyword cAnsiFunction	tgammaf tgamma lgammal
syn keyword cAnsiFunction	lgammaf lgamma erfcl
syn keyword cAnsiFunction	erfcf erfc erfl
syn keyword cAnsiFunction	erff erf sqrtl
syn keyword cAnsiFunction	sqrtf sqrt powl
syn keyword cAnsiFunction	powf pow hypotl
syn keyword cAnsiFunction	hypotf hypot fabsl
syn keyword cAnsiFunction	fabsf fabs cbrtl
syn keyword cAnsiFunction	cbrtf cbrt scalblnl
syn keyword cAnsiFunction	scalblnf scalbln scalbnl
syn keyword cAnsiFunction	scalbnf scalbn modfl
syn keyword cAnsiFunction	modff modf logbl
syn keyword cAnsiFunction	logbf logb log2l
syn keyword cAnsiFunction	log2f log2 log1pl
syn keyword cAnsiFunction	log1pf log1p log10l
syn keyword cAnsiFunction	log10f log10 logl
syn keyword cAnsiFunction	logf log ldexpl
syn keyword cAnsiFunction	ldexpf ldexp ilogbl
syn keyword cAnsiFunction	ilogbf ilogb frexpl
syn keyword cAnsiFunction	frexpf frexp expm1l
syn keyword cAnsiFunction	expm1f expm1 exp2l
syn keyword cAnsiFunction	exp2f exp2 expl
syn keyword cAnsiFunction	expf exp tanhl
syn keyword cAnsiFunction	tanhf tanh sinhl
syn keyword cAnsiFunction	sinhf sinh coshl
syn keyword cAnsiFunction	coshf cosh atanhl
syn keyword cAnsiFunction	atanhf atanh asinhl
syn keyword cAnsiFunction	asinhf asinh acoshl
syn keyword cAnsiFunction	acoshf acosh tanl
syn keyword cAnsiFunction	tanf tan sinl
syn keyword cAnsiFunction	sinf sin cosl
syn keyword cAnsiFunction	cosf cos atan2l
syn keyword cAnsiFunction	atan2f atan2 atanl
syn keyword cAnsiFunction	atanf atan asinl
syn keyword cAnsiFunction	asinf asin acosl
syn keyword cAnsiFunction	acosf acos signbit
syn keyword cAnsiFunction	isnormal isnan isinf
syn keyword cAnsiFunction	isfinite fpclassify localeconv
syn keyword cAnsiFunction	setlocale wcstoumax wcstoimax
syn keyword cAnsiFunction	strtoumax strtoimax feupdateenv
syn keyword cAnsiFunction	fesetenv feholdexcept fegetenv
syn keyword cAnsiFunction	fesetround fegetround fetestexcept
syn keyword cAnsiFunction	fesetexceptflag feraiseexcept fegetexceptflag
syn keyword cAnsiFunction	feclearexcept toupper tolower
syn keyword cAnsiFunction	isxdigit isupper isspace
syn keyword cAnsiFunction	ispunct isprint islower
syn keyword cAnsiFunction	isgraph isdigit iscntrl
syn keyword cAnsiFunction	isalpha isalnum creall
syn keyword cAnsiFunction	crealf creal cprojl
syn keyword cAnsiFunction	cprojf cproj conjl
syn keyword cAnsiFunction	conjf conj cimagl
syn keyword cAnsiFunction	cimagf cimag cargl
syn keyword cAnsiFunction	cargf carg csqrtl
syn keyword cAnsiFunction	csqrtf csqrt cpowl
syn keyword cAnsiFunction	cpowf cpow cabsl
syn keyword cAnsiFunction	cabsf cabs clogl
syn keyword cAnsiFunction	clogf clog cexpl
syn keyword cAnsiFunction	cexpf cexp ctanhl
syn keyword cAnsiFunction	ctanhf ctanh csinhl
syn keyword cAnsiFunction	csinhf csinh ccoshl
syn keyword cAnsiFunction	ccoshf ccosh catanhl
syn keyword cAnsiFunction	catanhf catanh casinhl
syn keyword cAnsiFunction	casinhf casinh cacoshl
syn keyword cAnsiFunction	cacoshf cacosh ctanl
syn keyword cAnsiFunction	ctanf ctan csinl
syn keyword cAnsiFunction	csinf csin ccosl
syn keyword cAnsiFunction	ccosf ccos catanl
syn keyword cAnsiFunction	catanf catan casinl
syn keyword cAnsiFunction	casinf casin cacosl
syn keyword cAnsiFunction	cacosf cacos assert
syn keyword cAnsiFunction	UINTMAX_C INTMAX_C UINT64_C
syn keyword cAnsiFunction	UINT32_C UINT16_C UINT8_C
syn keyword cAnsiFunction	INT64_C INT32_C INT16_C INT8_C

syn keyword	cMC		__interrupt __nosavereg
syn keyword	cAnsiName	PRId8 PRIi16 PRIo32 PRIu64
syn keyword	cAnsiName	PRId16 PRIi32 PRIo64 PRIuLEAST8
syn keyword	cAnsiName	PRId32 PRIi64 PRIoLEAST8 PRIuLEAST16
syn keyword	cAnsiName	PRId64 PRIiLEAST8 PRIoLEAST16 PRIuLEAST32
syn keyword	cAnsiName	PRIdLEAST8 PRIiLEAST16 PRIoLEAST32 PRIuLEAST64
syn keyword	cAnsiName	PRIdLEAST16 PRIiLEAST32 PRIoLEAST64 PRIuFAST8
syn keyword	cAnsiName	PRIdLEAST32 PRIiLEAST64 PRIoFAST8 PRIuFAST16
syn keyword	cAnsiName	PRIdLEAST64 PRIiFAST8 PRIoFAST16 PRIuFAST32
syn keyword	cAnsiName	PRIdFAST8 PRIiFAST16 PRIoFAST32 PRIuFAST64
syn keyword	cAnsiName	PRIdFAST16 PRIiFAST32 PRIoFAST64 PRIuMAX
syn keyword	cAnsiName	PRIdFAST32 PRIiFAST64 PRIoMAX PRIuPTR
syn keyword	cAnsiName	PRIdFAST64 PRIiMAX PRIoPTR PRIx8
syn keyword	cAnsiName	PRIdMAX PRIiPTR PRIu8 PRIx16
syn keyword	cAnsiName	PRIdPTR PRIo8 PRIu16 PRIx32
syn keyword	cAnsiName	PRIi8 PRIo16 PRIu32 PRIx64

syn keyword	cAnsiName	PRIxLEAST8 SCNd8 SCNiFAST32 SCNuLEAST32
syn keyword	cAnsiName	PRIxLEAST16 SCNd16 SCNiFAST64 SCNuLEAST64
syn keyword	cAnsiName	PRIxLEAST32 SCNd32 SCNiMAX SCNuFAST8
syn keyword	cAnsiName	PRIxLEAST64 SCNd64 SCNiPTR SCNuFAST16
syn keyword	cAnsiName	PRIxFAST8 SCNdLEAST8 SCNo8 SCNuFAST32
syn keyword	cAnsiName	PRIxFAST16 SCNdLEAST16 SCNo16 SCNuFAST64
syn keyword	cAnsiName	PRIxFAST32 SCNdLEAST32 SCNo32 SCNuMAX
syn keyword	cAnsiName	PRIxFAST64 SCNdLEAST64 SCNo64 SCNuPTR
syn keyword	cAnsiName	PRIxMAX SCNdFAST8 SCNoLEAST8 SCNx8
syn keyword	cAnsiName	PRIxPTR SCNdFAST16 SCNoLEAST16 SCNx16
syn keyword	cAnsiName	PRIX8 SCNdFAST32 SCNoLEAST32 SCNx32
syn keyword	cAnsiName	PRIX16 SCNdFAST64 SCNoLEAST64 SCNx64
syn keyword	cAnsiName	PRIX32 SCNdMAX SCNoFAST8 SCNxLEAST8
syn keyword	cAnsiName	PRIX64 SCNdPTR SCNoFAST16 SCNxLEAST16
syn keyword	cAnsiName	PRIXLEAST8 SCNi8 SCNoFAST32 SCNxLEAST32
syn keyword	cAnsiName	PRIXLEAST16 SCNi16 SCNoFAST64 SCNxLEAST64
syn keyword	cAnsiName	PRIXLEAST32 SCNi32 SCNoMAX SCNxFAST8
syn keyword	cAnsiName	PRIXLEAST64 SCNi64 SCNoPTR SCNxFAST16
syn keyword	cAnsiName	PRIXFAST8 SCNiLEAST8 SCNu8 SCNxFAST32
syn keyword	cAnsiName	PRIXFAST16 SCNiLEAST16 SCNu16 SCNxFAST64
syn keyword	cAnsiName	PRIXFAST32 SCNiLEAST32 SCNu32 SCNxMAX
syn keyword	cAnsiName	PRIXFAST64 SCNiLEAST64 SCNu64 SCNxPTR
syn keyword	cAnsiName	PRIXMAX SCNiFAST8 SCNuLEAST8
syn keyword	cAnsiName	PRIXPTR SCNiFAST16 SCNuLEAST16

syn keyword	cAnsiName	errno environ
syn keyword	cAnsiName	stdout stdin stderr

syn keyword	cAnsiName	STDC CX_LIMITED_RANGE
syn keyword	cAnsiName	STDC FENV_ACCESS
syn keyword	cAnsiName	STDC FP_CONTRACT

"syn keyword	cAnsiName	acos sqrt fmod nearbyint
"syn keyword	cAnsiName	asin fabs frexp nextafter
"syn keyword	cAnsiName	atan atan2 tgamma nextafterx
"syn keyword	cAnsiName	acosh cbrt hypot remainder
"syn keyword	cAnsiName	asinh ceil ilogb remquo
"syn keyword	cAnsiName	atanh copysign ldexp rint
"syn keyword	cAnsiName	cos erf lgamma round
"syn keyword	cAnsiName	sin erfc llrint scalbn
"syn keyword	cAnsiName	tan exp2 llround scalbln
"syn keyword	cAnsiName	cosh expm1 log10 trunc
"syn keyword	cAnsiName	sinh fdim log1p carg
"syn keyword	cAnsiName	tanh floor log2 cimag
"syn keyword	cAnsiName	exp fma logb conj
"syn keyword	cAnsiName	log fmax lrint cproj
"syn keyword	cAnsiName	pow fmin lround creal

syn keyword	cAnsiName	and bitor not_eq xor
syn keyword	cAnsiName	and_eq compl or xor_eq
syn keyword	cAnsiName	bitand not or_eq

" A bunch of useful C keywords
"syn keyword	cStatement	goto
syn keyword	cStatement	break return continue asm
syn keyword	cLabel		default
syn keyword	cLabel		case
syn keyword	cConditional	if else switch
syn keyword	cRepeat		while for do

syn keyword	cTodo		contained TODO FIXME XXX
syn match	cTodo		contained "///[A-Z]!*"

" cCommentGroup allows adding matches for special things in comments
syn cluster	cCommentGroup	contains=cTodo

" String and Character constants
" Highlight special characters (those which have a backslash) differently
syn match	cSpecial	display contained "\\\(x\x\+\|\o\{1,3}\|.\|$\)"
if !exists("c_no_utf")
  syn match	cSpecial	display contained "\\\(u\x\{4}\|U\x\{8}\)"
endif
if exists("c_no_cformat")
  syn region	cString		start=+L\="+ skip=+\\\\\|\\"+ end=+"+ contains=cSpecial,@Spell
  " cCppString: same as cString, but ends at end of line
  syn region	cCppString	start=+L\="+ skip=+\\\\\|\\"\|\\$+ excludenl end=+"+ end='$' contains=cSpecial,@Spell
else
  syn match	cFormat		display "%\(\d\+\$\)\=[-+' #0*,]*\(\d*\|\*\|\*\d\+\$\)\(\.\(\d*\|\*\|\*\d\+\$\)\)\=\([hlLjztF]\|ll\|hh\)\=\([bdiuoxXDOUfeEgGcCsSpnAaK]\|\[\^\=.[^]]*\]\)" contained
  syn match	cFormat		display "%%" contained
  syn region	cString		start=+L\="+ skip=+\\\\\|\\"+ end=+"+ contains=cSpecial,cFormat,@Spell
  " cCppString: same as cString, but ends at end of line
  syn region	cCppString	start=+L\="+ skip=+\\\\\|\\"\|\\$+ excludenl end=+"+ end='$' contains=cSpecial,cFormat
  hi link cFormat cSpecial
endif
hi link cCppString cString

syn match	cCharacter	"L\='[^\\]'"
syn match	cCharacter	"L'[^']*'" contains=cSpecial
if exists("c_gnu")
  syn match	cSpecialError	"L\='\\[^'\"?\\abefnrtv]'"
  syn match	cSpecialCharacter "L\='\\['\"?\\abefnrtv]'"
else
  syn match	cSpecialError	"L\='\\[^'\"?\\abfnrtv]'"
  syn match	cSpecialCharacter "L\='\\['\"?\\abfnrtv]'"
endif
syn match	cSpecialCharacter display "L\='\\\o\{1,3}'"
syn match	cSpecialCharacter display "'\\x\x\{1,2}'"
syn match	cSpecialCharacter display "L'\\x\x\+'"

"when wanted, highlight trailing white space
if exists("c_space_errors")
  if !exists("c_no_trail_space_error")
    syn match	cSpaceError	display excludenl "\s\+$"
  endif
  if !exists("c_no_tab_space_error")
    syn match	cSpaceError	display " \+\t"me=e-1
  endif
endif

"catch errors caused by wrong parenthesis and brackets
" also accept <% for {, %> for }, <: for [ and :> for ] (C99)
syn cluster	cParenGroup	contains=cParenError,cIncluded,cSpecial,cCommentSkip,cCommentString,cComment2String,@cCommentGroup,cCommentStartError,cUserCont,cUserLabel,cUserLabel2,cGotoLabel,cBitField,cCommentSkip,cOctalZero,cCppOut,cCppOut2,cCppSkip,cFormat,cNumber,cFloat,cOctal,cOctalError,cNumbersCom
if exists("c_no_bracket_error")
  syn region	cParen		transparent start='(' end=')' contains=ALLBUT,@cParenGroup,cCppParen,cCppString,@Spell
  " cCppParen: same as cParen but ends at end-of-line; used in cDefine
  syn region	cCppParen	transparent start='(' skip='\\$' excludenl end=')' end='$' contained contains=ALLBUT,@cParenGroup,cParen,cString,@Spell
  syn match	cParenError	display ")"
" syn match	cErrInParen	display contained "[{}]\|<%\|%>"
  syn match	cErrInParen	display contained "[]\|<%\|%>"
else
  syn region	cParen		transparent start='(' end=')' contains=ALLBUT,@cParenGroup,cCppParen,cErrInBracket,cCppBracket,cCppString,@Spell
  " cCppParen: same as cParen but ends at end-of-line; used in cDefine
  syn region	cCppParen	transparent start='(' skip='\\$' excludenl end=')' end='$' contained contains=ALLBUT,@cParenGroup,cErrInBracket,cParen,cBracket,cString,@Spell
  syn match	cParenError	display "[\])]"
" syn match	cErrInParen	display contained "[\]{}]\|<%\|%>"
  syn match	cErrInParen	display contained "[\]]\|<%\|%>"
  syn region	cBracket	transparent start='\[\|<::\@!' end=']\|:>' contains=ALLBUT,@cParenGroup,cErrInParen,cCppParen,cCppBracket,cCppString,@Spell
  " cCppBracket: same as cParen but ends at end-of-line; used in cDefine
  syn region	cCppBracket	transparent start='\[\|<::\@!' skip='\\$' excludenl end=']\|:>' end='$' contained contains=ALLBUT,@cParenGroup,cErrInParen,cParen,cBracket,cString,@Spell
  syn match	cErrInBracket	display contained "[);{}]\|<%\|%>"
  "syn region	cBlock		transparent matchgroup=cBraces start='{' end='}' contains=ALLBUT,@cParenGroup,cCppParen,cCppBracket,cCppString,cBraceError,cErrInBracket
  "syn match	cBraceError	"}"
endif

"integer number, or floating point number without a dot and with "f".
syn case ignore
syn match	cNumbers	display transparent "\<\d\|\.\d" contains=cNumber,cFloat,cOctalError,cOctal
" Same, but without octal error (for comments)
syn match	cNumbersCom	display contained transparent "\<\d\|\.\d" contains=cNumber,cFloat,cOctal
syn match	cNumber		display contained "\d\+\(u\=l\{0,2}\|ll\=u\)\>"
"hex number
syn match	cNumber		display contained "0x\x\+\(u\=l\{0,2}\|ll\=u\)\>"
" Flag the first zero of an octal number as something special
syn match	cOctal		display contained "0\o\+\(u\=l\{0,2}\|ll\=u\)\>" contains=cOctalZero
syn match	cOctalZero	display contained "\<0"
syn match	cFloat		display contained "\d\+f"
"floating point number, with dot, optional exponent
syn match	cFloat		display contained "\d\+\.\d*\(e[-+]\=\d\+\)\=[fl]\="
"floating point number, starting with a dot, optional exponent
syn match	cFloat		display contained "\.\d\+\(e[-+]\=\d\+\)\=[fl]\=\>"
"floating point number, without dot, with exponent
syn match	cFloat		display contained "\d\+e[-+]\=\d\+[fl]\=\>"
if !exists("c_no_c99")
  "hexadecimal floating point number, optional leading digits, with dot, with exponent
  syn match	cFloat		display contained "0x\x*\.\x\+p[-+]\=\d\+[fl]\=\>"
  "hexadecimal floating point number, with leading digits, optional dot, with exponent
  syn match	cFloat		display contained "0x\x\+\.\=p[-+]\=\d\+[fl]\=\>"
endif

" flag an octal number with wrong digits
syn match	cOctalError	display contained "0\o*[89]\d*"
syn case match

if exists("xxxc_comment_strings")
  " A comment can contain cString, cCharacter and cNumber.
  " But a "*/" inside a cString in a cComment DOES end the comment!  So we
  " need to use a special type of cString: cCommentString, which also ends on
  " "*/", and sees a "*" at the start of the line as comment again.
  " Unfortunately this doesn't very well work for // type of comments :-(
  syntax match	cCommentSkip	contained "^\s*\*\($\|\s\+\)"
  syntax region cCommentString	contained start=+L\=\\\@<!"+ skip=+\\\\\|\\"+ end=+"+ end=+\*/+me=s-1 contains=cSpecial,cCommentSkip
  syntax region cComment2String	contained start=+L\=\\\@<!"+ skip=+\\\\\|\\"+ end=+"+ end="$" contains=cSpecial
  syntax region  cCommentL	start="//" skip="\\$" end="$" keepend contains=@cCommentGroup,cComment2String,cCharacter,cNumbersCom,cSpaceError,@Spell
  syntax region cComment	matchgroup=cCommentStart start="/\*" end="\*/" contains=@cCommentGroup,cCommentStartError,cCommentString,cCharacter,cNumbersCom,cSpaceError,@Spell
else
  syn region	cCommentL	start="//" skip="\\$" end="$" keepend contains=@cCommentGroup,cSpaceError,@Spell
  syn region	cComment	matchgroup=cCommentStart start="/\*" end="\*/" contains=@cCommentGroup,cCommentStartError,cSpaceError,@Spell
endif
" keep a // comment separately, it terminates a preproc. conditional
syntax match	cCommentError	display "\*/"
syntax match	cCommentStartError display "/\*"me=e-1 contained

syn keyword	cOperator	sizeof
if exists("c_gnu")
  syn keyword	cStatement	__asm__
  syn keyword	cOperator	typeof __real__ __imag__
endif
syn keyword	cType		int long short char void
syn keyword	cType		signed unsigned float double

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"   ISSM special syntax                                                                                             "
"   please try to keep the alphabetical ordering                                                                    "
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"Petsc
syn keyword	cType		Vec Mat SeqVec SeqMat
"other ISSM's syntax
syn keyword	cType		mxArray ErrorException QuadtreeBox
syn keyword	cType		IssmDouble IssmPDouble

"ISSM's objects begin
syn keyword cType BoolInput
syn keyword cType BoolParam
syn keyword cType classes
syn keyword cType Constraint
syn keyword cType Constraints
syn keyword cType Contour
syn keyword cType Contours
syn keyword cType ControlInput
syn keyword cType Covertree
syn keyword cType DatasetInput
syn keyword cType DataSetParam
syn keyword cType Definition
syn keyword cType DependentObject
syn keyword cType DofIndexing
syn keyword cType DoubleArrayInput
syn keyword cType DoubleInput
syn keyword cType DoubleMatArrayParam
syn keyword cType DoubleMatParam
syn keyword cType DoubleParam
syn keyword cType DoubleTransientMatParam
syn keyword cType DoubleVecParam
syn keyword cType Element
syn keyword cType ElementHook
syn keyword cType ElementMatrix
syn keyword cType Elements
syn keyword cType ElementVector
syn keyword cType ExponentialVariogram
syn keyword cType ExternalResult
syn keyword cType FemModel
syn keyword cType FileParam
syn keyword cType Friction
syn keyword cType Gauss
syn keyword cType GaussianVariogram
syn keyword cType gaussobjects
syn keyword cType GaussPenta
syn keyword cType GaussSeg
syn keyword cType GaussTetra
syn keyword cType GaussTria
syn keyword cType GenericExternalResult
syn keyword cType GenericOption
syn keyword cType GenericParam
syn keyword cType GiaDeflectionCoreArgs
syn keyword cType Hook
syn keyword cType Input
syn keyword cType Inputs
syn keyword cType IntInput
syn keyword cType IntMatParam
syn keyword cType IntParam
syn keyword cType IntVecParam
syn keyword cType IoModel
syn keyword cType IssmDirectApplicInterface
syn keyword cType IssmParallelDirectApplicInterface
syn keyword cType krigingobjects
syn keyword cType Load
syn keyword cType Loads
syn keyword cType Masscon
syn keyword cType Massconaxpby
syn keyword cType Massfluxatgate
syn keyword cType Material
syn keyword cType Materials
syn keyword cType Matestar
syn keyword cType Matice
syn keyword cType Matpar
syn keyword cType matrixobjects
syn keyword cType MatrixParam
syn keyword cType Misfit
syn keyword cType Moulin
syn keyword cType Neumannflux
syn keyword cType Nodalvalue
syn keyword cType Node
syn keyword cType Nodes
syn keyword cType Numericalflux
syn keyword cType Observation
syn keyword cType Observations
syn keyword cType Option
syn keyword cType Options
syn keyword cType OptionUtilities
syn keyword cType Param
syn keyword cType Parameters
syn keyword cType Pengrid
syn keyword cType Penpair
syn keyword cType Penta
syn keyword cType PentaInput
syn keyword cType PentaRef
syn keyword cType PowerVariogram
syn keyword cType Profiler
syn keyword cType Quadtree
syn keyword cType Results
syn keyword cType Riftfront
syn keyword cType RiftStruct
syn keyword cType Seg
syn keyword cType SegInput
syn keyword cType Segment
syn keyword cType SegRef
syn keyword cType SpcDynamic
syn keyword cType SpcStatic
syn keyword cType SpcTransient
syn keyword cType SphericalVariogram
syn keyword cType StringArrayParam
syn keyword cType StringParam
syn keyword cType Tetra
syn keyword cType TetraInput
syn keyword cType TetraRef
syn keyword cType TransientInput
syn keyword cType TransientParam
syn keyword cType Tria
syn keyword cType TriaInput
syn keyword cType TriaRef
syn keyword cType Update
syn keyword cType Variogram
syn keyword cType VectorParam
syn keyword cType Vertex
syn keyword cType Vertices
syn keyword cType AdjointBalancethickness2Analysis
syn keyword cType AdjointBalancethicknessAnalysis
syn keyword cType AdjointHorizAnalysis
syn keyword cType Analysis
syn keyword cType Balancethickness2Analysis
syn keyword cType BalancethicknessAnalysis
syn keyword cType BalancethicknessSoftAnalysis
syn keyword cType BalancevelocityAnalysis
syn keyword cType DamageEvolutionAnalysis
syn keyword cType DepthAverageAnalysis
syn keyword cType EnthalpyAnalysis
syn keyword cType EnumToAnalysis
syn keyword cType ExtrapolationAnalysis
syn keyword cType ExtrudeFromBaseAnalysis
syn keyword cType ExtrudeFromTopAnalysis
syn keyword cType FreeSurfaceBaseAnalysis
syn keyword cType FreeSurfaceTopAnalysis
syn keyword cType GiaAnalysis
syn keyword cType HydrologyDCEfficientAnalysis
syn keyword cType HydrologyDCInefficientAnalysis
syn keyword cType HydrologyShreveAnalysis
syn keyword cType HydrologySommersAnalysis
syn keyword cType L2ProjectionBaseAnalysis
syn keyword cType L2ProjectionEPLAnalysis
syn keyword cType LevelsetAnalysis
syn keyword cType LsfReinitializationAnalysis
syn keyword cType MasstransportAnalysis
syn keyword cType MeltingAnalysis
syn keyword cType MeshdeformationAnalysis
syn keyword cType SealevelriseAnalysis
syn keyword cType SmbAnalysis
syn keyword cType SmoothAnalysis
syn keyword cType StressbalanceAnalysis
syn keyword cType StressbalanceSIAAnalysis
syn keyword cType StressbalanceVerticalAnalysis
syn keyword cType ThermalAnalysis
syn keyword cType UzawaPressureAnalysis
"ISSM's objects end
"ISSM's Enums begin
syn keyword cConstant ParametersSTARTEnum
syn keyword cConstant FemModelEnum
syn keyword cConstant FemModelCommEnum
syn keyword cConstant WorldCommEnum
syn keyword cConstant IcecapToEarthCommEnum
syn keyword cConstant NumModelsEnum
syn keyword cConstant ModelIdEnum
syn keyword cConstant EarthIdEnum
syn keyword cConstant AutodiffIsautodiffEnum
syn keyword cConstant AutodiffNumDependentsEnum
syn keyword cConstant AutodiffNumDependentObjectsEnum
syn keyword cConstant AutodiffDependentObjectNamesEnum
syn keyword cConstant AutodiffDependentObjectTypesEnum
syn keyword cConstant AutodiffDependentObjectIndicesEnum
syn keyword cConstant AutodiffDependentObjectsEnum
syn keyword cConstant AutodiffNumIndependentsEnum
syn keyword cConstant AutodiffNumIndependentObjectsEnum
syn keyword cConstant AutodiffIndependentObjectNamesEnum
syn keyword cConstant AutodiffIndependentObjectTypesEnum
syn keyword cConstant AutodiffIndependentObjectsEnum
syn keyword cConstant AutodiffJacobianEnum
syn keyword cConstant AutodiffXpEnum
syn keyword cConstant AutodiffDriverEnum
syn keyword cConstant AutodiffFosForwardIndexEnum
syn keyword cConstant AutodiffFovForwardIndicesEnum
syn keyword cConstant AutodiffFosReverseIndexEnum
syn keyword cConstant AutodiffMassFluxSegmentsPresentEnum
syn keyword cConstant AutodiffKeepEnum
syn keyword cConstant AutodiffObufsizeEnum
syn keyword cConstant AutodiffLbufsizeEnum
syn keyword cConstant AutodiffCbufsizeEnum
syn keyword cConstant AutodiffTbufsizeEnum
syn keyword cConstant AutodiffGcTriggerRatioEnum
syn keyword cConstant AutodiffGcTriggerMaxSizeEnum
syn keyword cConstant BalancethicknessSpcthicknessEnum
syn keyword cConstant BalancethicknessStabilizationEnum
syn keyword cConstant BalancethicknessThickeningRateEnum
syn keyword cConstant BasalforcingsEnum
syn keyword cConstant BasalforcingsGeothermalfluxEnum
syn keyword cConstant BasalforcingsGroundediceMeltingRateEnum
syn keyword cConstant BasalforcingsFloatingiceMeltingRateEnum
syn keyword cConstant BasalforcingsDeepwaterMeltingRateEnum
syn keyword cConstant BasalforcingsDeepwaterElevationEnum
syn keyword cConstant BasalforcingsUpperwaterElevationEnum
syn keyword cConstant BasalforcingsMeltrateFactorEnum
syn keyword cConstant BasalforcingsThresholdThicknessEnum
syn keyword cConstant BasalforcingsUpperdepthMeltEnum
syn keyword cConstant BasalforcingsMantleconductivityEnum
syn keyword cConstant BasalforcingsNusseltEnum
syn keyword cConstant BasalforcingsDtbgEnum
syn keyword cConstant BasalforcingsPlumeradiusEnum
syn keyword cConstant BasalforcingsTopplumedepthEnum
syn keyword cConstant BasalforcingsBottomplumedepthEnum
syn keyword cConstant BasalforcingsPlumexEnum
syn keyword cConstant BasalforcingsPlumeyEnum
syn keyword cConstant BasalforcingsCrustthicknessEnum
syn keyword cConstant BasalforcingsUppercrustthicknessEnum
syn keyword cConstant BasalforcingsUppercrustheatEnum
syn keyword cConstant BasalforcingsLowercrustheatEnum
syn keyword cConstant FloatingMeltRateEnum
syn keyword cConstant LinearFloatingMeltRateEnum
syn keyword cConstant MismipFloatingMeltRateEnum
syn keyword cConstant MantlePlumeGeothermalFluxEnum
syn keyword cConstant BedEnum
syn keyword cConstant BaseEnum
syn keyword cConstant ConstantsGEnum
syn keyword cConstant ConstantsReferencetemperatureEnum
syn keyword cConstant ConstantsYtsEnum
syn keyword cConstant DependentObjectEnum
syn keyword cConstant StressbalanceAbstolEnum
syn keyword cConstant StressbalanceConvergenceNumStepsEnum
syn keyword cConstant StressbalanceIsnewtonEnum
syn keyword cConstant StressbalanceMaxiterEnum
syn keyword cConstant StressbalancePenaltyFactorEnum
syn keyword cConstant StressbalanceReferentialEnum
syn keyword cConstant StressbalanceReltolEnum
syn keyword cConstant StressbalanceNumRequestedOutputsEnum
syn keyword cConstant StressbalanceRequestedOutputsEnum
syn keyword cConstant StressbalanceRestolEnum
syn keyword cConstant StressbalanceRiftPenaltyLockEnum
syn keyword cConstant StressbalanceRiftPenaltyThresholdEnum
syn keyword cConstant StressbalanceShelfDampeningEnum
syn keyword cConstant StressbalanceSpcvxEnum
syn keyword cConstant StressbalanceSpcvyEnum
syn keyword cConstant StressbalanceSpcvzEnum
syn keyword cConstant StressbalanceFSreconditioningEnum
syn keyword cConstant StressbalanceVertexPairingEnum
syn keyword cConstant StressbalanceViscosityOvershootEnum
syn keyword cConstant LoadingforceXEnum
syn keyword cConstant LoadingforceYEnum
syn keyword cConstant LoadingforceZEnum
syn keyword cConstant FlowequationBorderSSAEnum
syn keyword cConstant FlowequationBorderHOEnum
syn keyword cConstant FlowequationBorderFSEnum
syn keyword cConstant FlowequationElementEquationEnum
syn keyword cConstant FlowequationIsSIAEnum
syn keyword cConstant FlowequationIsSSAEnum
syn keyword cConstant FlowequationIsL1L2Enum
syn keyword cConstant FlowequationIsHOEnum
syn keyword cConstant FlowequationIsFSEnum
syn keyword cConstant FlowequationFeSSAEnum
syn keyword cConstant FlowequationFeHOEnum
syn keyword cConstant FlowequationFeFSEnum
syn keyword cConstant FlowequationVertexEquationEnum
syn keyword cConstant FrictionAsEnum
syn keyword cConstant FrictionCoefficientEnum
syn keyword cConstant FrictionCoefficientcoulombEnum
syn keyword cConstant FrictionPEnum
syn keyword cConstant FrictionQEnum
syn keyword cConstant FrictionMEnum
syn keyword cConstant FrictionCEnum
syn keyword cConstant FrictionLawEnum
syn keyword cConstant FrictionGammaEnum
syn keyword cConstant FrictionFEnum
syn keyword cConstant FrictionWaterLayerEnum
syn keyword cConstant FrictionEffectivePressureEnum
syn keyword cConstant FrictionCouplingEnum
syn keyword cConstant GeometryHydrostaticRatioEnum
syn keyword cConstant HydrologyModelEnum
syn keyword cConstant HydrologyshreveEnum
syn keyword cConstant HydrologyshreveSpcwatercolumnEnum
syn keyword cConstant HydrologyshreveStabilizationEnum
syn keyword cConstant HydrologydcEnum
syn keyword cConstant SedimentHeadEnum
syn keyword cConstant SedimentHeadOldEnum
syn keyword cConstant SedimentHeadResidualEnum
syn keyword cConstant EffectivePressureEnum
syn keyword cConstant EplHeadEnum
syn keyword cConstant EplHeadOldEnum
syn keyword cConstant EplHeadSlopeXEnum
syn keyword cConstant EplHeadSlopeYEnum
syn keyword cConstant EplZigZagCounterEnum
syn keyword cConstant HydrologydcMaxIterEnum
syn keyword cConstant HydrologydcRelTolEnum
syn keyword cConstant HydrologydcSpcsedimentHeadEnum
syn keyword cConstant HydrologydcSedimentCompressibilityEnum
syn keyword cConstant HydrologydcSedimentPorosityEnum
syn keyword cConstant HydrologydcSedimentThicknessEnum
syn keyword cConstant HydrologydcSedimentTransmitivityEnum
syn keyword cConstant HydrologydcWaterCompressibilityEnum
syn keyword cConstant HydrologydcSpceplHeadEnum
syn keyword cConstant HydrologydcMaskEplactiveNodeEnum
syn keyword cConstant HydrologydcMaskEplactiveEltEnum
syn keyword cConstant HydrologydcEplCompressibilityEnum
syn keyword cConstant HydrologydcEplPorosityEnum
syn keyword cConstant HydrologydcEplInitialThicknessEnum
syn keyword cConstant HydrologydcEplColapseThicknessEnum
syn keyword cConstant HydrologydcEplMaxThicknessEnum
syn keyword cConstant HydrologydcEplThicknessEnum
syn keyword cConstant HydrologydcEplThicknessOldEnum
syn keyword cConstant HydrologydcEplThickCompEnum
syn keyword cConstant HydrologydcEplConductivityEnum
syn keyword cConstant HydrologydcIsefficientlayerEnum
syn keyword cConstant HydrologydcSedimentlimitFlagEnum
syn keyword cConstant HydrologydcSedimentlimitEnum
syn keyword cConstant HydrologydcTransferFlagEnum
syn keyword cConstant HydrologydcLeakageFactorEnum
syn keyword cConstant HydrologydcPenaltyFactorEnum
syn keyword cConstant HydrologydcPenaltyLockEnum
syn keyword cConstant HydrologydcEplflipLockEnum
syn keyword cConstant HydrologydcBasalMoulinInputEnum
syn keyword cConstant HydrologyLayerEnum
syn keyword cConstant HydrologySedimentEnum
syn keyword cConstant HydrologyEfficientEnum
syn keyword cConstant HydrologySedimentKmaxEnum
syn keyword cConstant HydrologysommersEnum
syn keyword cConstant HydrologyHeadEnum
syn keyword cConstant HydrologyGapHeightEnum
syn keyword cConstant HydrologyBumpSpacingEnum
syn keyword cConstant HydrologyBumpHeightEnum
syn keyword cConstant HydrologyEnglacialInputEnum
syn keyword cConstant HydrologyMoulinInputEnum
syn keyword cConstant HydrologyReynoldsEnum
syn keyword cConstant HydrologyNeumannfluxEnum
syn keyword cConstant HydrologySpcheadEnum
syn keyword cConstant HydrologyConductivityEnum
syn keyword cConstant InversionControlParametersEnum
syn keyword cConstant InversionControlScalingFactorsEnum
syn keyword cConstant InversionCostFunctionThresholdEnum
syn keyword cConstant InversionCostFunctionsCoefficientsEnum
syn keyword cConstant InversionCostFunctionsEnum
syn keyword cConstant InversionGradientScalingEnum
syn keyword cConstant InversionIscontrolEnum
syn keyword cConstant InversionTypeEnum
syn keyword cConstant InversionIncompleteAdjointEnum
syn keyword cConstant InversionMaxParametersEnum
syn keyword cConstant InversionMaxiterPerStepEnum
syn keyword cConstant InversionMaxiterEnum
syn keyword cConstant InversionMaxstepsEnum
syn keyword cConstant InversionFatolEnum
syn keyword cConstant InversionFrtolEnum
syn keyword cConstant InversionGatolEnum
syn keyword cConstant InversionGrtolEnum
syn keyword cConstant InversionGttolEnum
syn keyword cConstant InversionAlgorithmEnum
syn keyword cConstant InversionMinParametersEnum
syn keyword cConstant InversionNstepsEnum
syn keyword cConstant InversionDxminEnum
syn keyword cConstant InversionNumControlParametersEnum
syn keyword cConstant InversionNumCostFunctionsEnum
syn keyword cConstant InversionStepThresholdEnum
syn keyword cConstant InversionThicknessObsEnum
syn keyword cConstant InversionSurfaceObsEnum
syn keyword cConstant InversionVxObsEnum
syn keyword cConstant InversionVyObsEnum
syn keyword cConstant InversionVzObsEnum
syn keyword cConstant MaskIceLevelsetEnum
syn keyword cConstant MaskOceanLevelsetEnum
syn keyword cConstant MaskLandLevelsetEnum
syn keyword cConstant MaterialsBetaEnum
syn keyword cConstant MaterialsHeatcapacityEnum
syn keyword cConstant MaterialsLatentheatEnum
syn keyword cConstant MaterialsMeltingpointEnum
syn keyword cConstant MaterialsMixedLayerCapacityEnum
syn keyword cConstant MaterialsRheologyBEnum
syn keyword cConstant MaterialsRheologyBbarEnum
syn keyword cConstant MaterialsRheologyLawEnum
syn keyword cConstant MaterialsRheologyNEnum
syn keyword cConstant MaterialsRheologyKoEnum
syn keyword cConstant MaterialsRheologyKobarEnum
syn keyword cConstant MaterialsRheologyEcEnum
syn keyword cConstant MaterialsRheologyEcbarEnum
syn keyword cConstant MaterialsRheologyEsEnum
syn keyword cConstant MaterialsRheologyEsbarEnum
syn keyword cConstant DamageIsdamageEnum
syn keyword cConstant DamageDEnum
syn keyword cConstant DamageFEnum
syn keyword cConstant DamageDbarEnum
syn keyword cConstant DamageLawEnum
syn keyword cConstant DamageC1Enum
syn keyword cConstant DamageC2Enum
syn keyword cConstant DamageC3Enum
syn keyword cConstant DamageC4Enum
syn keyword cConstant DamageElementinterpEnum
syn keyword cConstant DamageHealingEnum
syn keyword cConstant DamageStressThresholdEnum
syn keyword cConstant DamageKappaEnum
syn keyword cConstant DamageStabilizationEnum
syn keyword cConstant DamageMaxiterEnum
syn keyword cConstant DamageSpcdamageEnum
syn keyword cConstant DamageMaxDamageEnum
syn keyword cConstant DamageEquivStressEnum
syn keyword cConstant DamageEvolutionNumRequestedOutputsEnum
syn keyword cConstant DamageEvolutionRequestedOutputsEnum
syn keyword cConstant DamageEnum
syn keyword cConstant NewDamageEnum
syn keyword cConstant StressIntensityFactorEnum
syn keyword cConstant CalvingLawEnum
syn keyword cConstant CalvingCalvingrateEnum
syn keyword cConstant CalvingMeltingrateEnum
syn keyword cConstant CalvingLevermannEnum
syn keyword cConstant CalvingDevEnum
syn keyword cConstant CalvingMinthicknessEnum
syn keyword cConstant DefaultCalvingEnum
syn keyword cConstant CalvingRequestedOutputsEnum
syn keyword cConstant CalvinglevermannCoeffEnum
syn keyword cConstant CalvinglevermannMeltingrateEnum
syn keyword cConstant CalvingdevCoeffEnum
syn keyword cConstant CalvingratexEnum
syn keyword cConstant CalvingrateyEnum
syn keyword cConstant CalvingratexAverageEnum
syn keyword cConstant CalvingrateyAverageEnum
syn keyword cConstant StrainRateparallelEnum
syn keyword cConstant StrainRateperpendicularEnum
syn keyword cConstant StrainRateeffectiveEnum
syn keyword cConstant MaterialsRhoIceEnum
syn keyword cConstant MaterialsRhoSeawaterEnum
syn keyword cConstant MaterialsRhoFreshwaterEnum
syn keyword cConstant MaterialsMuWaterEnum
syn keyword cConstant MaterialsThermalExchangeVelocityEnum
syn keyword cConstant MaterialsThermalconductivityEnum
syn keyword cConstant MaterialsTemperateiceconductivityEnum
syn keyword cConstant MaterialsLithosphereShearModulusEnum
syn keyword cConstant MaterialsLithosphereDensityEnum
syn keyword cConstant MaterialsMantleShearModulusEnum
syn keyword cConstant MaterialsMantleDensityEnum
syn keyword cConstant MaterialsEarthDensityEnum
syn keyword cConstant MeshAverageVertexConnectivityEnum
syn keyword cConstant MeshElements2dEnum
syn keyword cConstant MeshElementsEnum
syn keyword cConstant MeshLowerelementsEnum
syn keyword cConstant MeshNumberofelements2dEnum
syn keyword cConstant MeshNumberofelementsEnum
syn keyword cConstant MeshNumberoflayersEnum
syn keyword cConstant MeshNumberofvertices2dEnum
syn keyword cConstant MeshNumberofverticesEnum
syn keyword cConstant MeshUpperelementsEnum
syn keyword cConstant MeshVertexonbaseEnum
syn keyword cConstant MeshVertexonsurfaceEnum
syn keyword cConstant MeshVertexonboundaryEnum
syn keyword cConstant MeshXEnum
syn keyword cConstant MeshYEnum
syn keyword cConstant MeshZEnum
syn keyword cConstant MeshLatEnum
syn keyword cConstant MeshLongEnum
syn keyword cConstant MeshREnum
syn keyword cConstant MeshElementtypeEnum
syn keyword cConstant MeshSegmentsEnum
syn keyword cConstant DomainTypeEnum
syn keyword cConstant DomainDimensionEnum
syn keyword cConstant Domain2DhorizontalEnum
syn keyword cConstant Domain2DverticalEnum
syn keyword cConstant Domain3DEnum
syn keyword cConstant Domain3DsurfaceEnum
syn keyword cConstant MiscellaneousNameEnum
syn keyword cConstant MasstransportHydrostaticAdjustmentEnum
syn keyword cConstant MasstransportIsfreesurfaceEnum
syn keyword cConstant MasstransportMinThicknessEnum
syn keyword cConstant MasstransportPenaltyFactorEnum
syn keyword cConstant MasstransportSpcthicknessEnum
syn keyword cConstant MasstransportStabilizationEnum
syn keyword cConstant MasstransportVertexPairingEnum
syn keyword cConstant MasstransportNumRequestedOutputsEnum
syn keyword cConstant MasstransportRequestedOutputsEnum
syn keyword cConstant QmuIsdakotaEnum
syn keyword cConstant MassFluxSegmentsEnum
syn keyword cConstant MassFluxSegmentsPresentEnum
syn keyword cConstant QmuMassFluxSegmentsPresentEnum
syn keyword cConstant QmuNumberofpartitionsEnum
syn keyword cConstant QmuNumberofresponsesEnum
syn keyword cConstant QmuPartitionEnum
syn keyword cConstant QmuResponsedescriptorsEnum
syn keyword cConstant QmuVariabledescriptorsEnum
syn keyword cConstant RiftsNumriftsEnum
syn keyword cConstant RiftsRiftstructEnum
syn keyword cConstant SettingsResultsOnNodesEnum
syn keyword cConstant SettingsIoGatherEnum
syn keyword cConstant SettingsLowmemEnum
syn keyword cConstant SettingsOutputFrequencyEnum
syn keyword cConstant SettingsRecordingFrequencyEnum
syn keyword cConstant SettingsWaitonlockEnum
syn keyword cConstant DebugProfilingEnum
syn keyword cConstant ProfilingCurrentMemEnum
syn keyword cConstant ProfilingCurrentFlopsEnum
syn keyword cConstant ProfilingSolutionTimeEnum
syn keyword cConstant SteadystateMaxiterEnum
syn keyword cConstant SteadystateNumRequestedOutputsEnum
syn keyword cConstant SteadystateReltolEnum
syn keyword cConstant SteadystateRequestedOutputsEnum
syn keyword cConstant SurfaceEnum
syn keyword cConstant ThermalIsenthalpyEnum
syn keyword cConstant ThermalIsdynamicbasalspcEnum
syn keyword cConstant ThermalReltolEnum
syn keyword cConstant ThermalMaxiterEnum
syn keyword cConstant ThermalPenaltyFactorEnum
syn keyword cConstant ThermalPenaltyLockEnum
syn keyword cConstant ThermalPenaltyThresholdEnum
syn keyword cConstant ThermalSpctemperatureEnum
syn keyword cConstant ThermalStabilizationEnum
syn keyword cConstant ThermalNumRequestedOutputsEnum
syn keyword cConstant ThermalRequestedOutputsEnum
syn keyword cConstant GiaMantleViscosityEnum
syn keyword cConstant GiaLithosphereThicknessEnum
syn keyword cConstant ThicknessEnum
syn keyword cConstant TimesteppingStartTimeEnum
syn keyword cConstant TimesteppingFinalTimeEnum
syn keyword cConstant TimesteppingCflCoefficientEnum
syn keyword cConstant TimesteppingTimeAdaptEnum
syn keyword cConstant TimesteppingTimeStepEnum
syn keyword cConstant TimesteppingInterpForcingsEnum
syn keyword cConstant TransientIssmbEnum
syn keyword cConstant TransientIscouplerEnum
syn keyword cConstant TransientIsstressbalanceEnum
syn keyword cConstant TransientIsgroundinglineEnum
syn keyword cConstant TransientIsmasstransportEnum
syn keyword cConstant TransientIsthermalEnum
syn keyword cConstant TransientIsgiaEnum
syn keyword cConstant TransientIsdamageevolutionEnum
syn keyword cConstant TransientIshydrologyEnum
syn keyword cConstant TransientIsmovingfrontEnum
syn keyword cConstant TransientIsslrEnum
syn keyword cConstant TransientNumRequestedOutputsEnum
syn keyword cConstant TransientRequestedOutputsEnum
syn keyword cConstant PotentialEnum
syn keyword cConstant BalancethicknessSpcpotentialEnum
syn keyword cConstant BalancethicknessApparentMassbalanceEnum
syn keyword cConstant Balancethickness2MisfitEnum
syn keyword cConstant BalancethicknessDiffusionCoefficientEnum
syn keyword cConstant BalancethicknessCmuEnum
syn keyword cConstant BalancethicknessOmegaEnum
syn keyword cConstant BalancethicknessD0Enum
syn keyword cConstant SmbEnum
syn keyword cConstant SmbAnalysisEnum
syn keyword cConstant SmbSolutionEnum
syn keyword cConstant SmbNumRequestedOutputsEnum
syn keyword cConstant SmbRequestedOutputsEnum
syn keyword cConstant SmbIsInitializedEnum
syn keyword cConstant SMBforcingEnum
syn keyword cConstant SmbMassBalanceEnum
syn keyword cConstant SMBgembEnum
syn keyword cConstant SmbInitDensityScalingEnum
syn keyword cConstant SmbTaEnum
syn keyword cConstant SmbVEnum
syn keyword cConstant SmbDswrfEnum
syn keyword cConstant SmbDlwrfEnum
syn keyword cConstant SmbPEnum
syn keyword cConstant SmbSwfEnum
syn keyword cConstant SmbEAirEnum
syn keyword cConstant SmbPAirEnum
syn keyword cConstant SmbTmeanEnum
syn keyword cConstant SmbCEnum
syn keyword cConstant SmbTzEnum
syn keyword cConstant SmbVzEnum
syn keyword cConstant SmbDtEnum
syn keyword cConstant SmbDzEnum
syn keyword cConstant SmbAIdxEnum
syn keyword cConstant SmbSwIdxEnum
syn keyword cConstant SmbDenIdxEnum
syn keyword cConstant SmbZTopEnum
syn keyword cConstant SmbDzTopEnum
syn keyword cConstant SmbDzMinEnum
syn keyword cConstant SmbZYEnum
syn keyword cConstant SmbZMaxEnum
syn keyword cConstant SmbZMinEnum
syn keyword cConstant SmbOutputFreqEnum
syn keyword cConstant SmbASnowEnum
syn keyword cConstant SmbAIceEnum
syn keyword cConstant SmbCldFracEnum
syn keyword cConstant SmbT0wetEnum
syn keyword cConstant SmbT0dryEnum
syn keyword cConstant SmbKEnum
syn keyword cConstant SmbDEnum
syn keyword cConstant SmbReEnum
syn keyword cConstant SmbGdnEnum
syn keyword cConstant SmbGspEnum
syn keyword cConstant SmbECEnum
syn keyword cConstant SmbCondensationEnum
syn keyword cConstant SmbWEnum
syn keyword cConstant SmbAEnum
syn keyword cConstant SmbTEnum
syn keyword cConstant SmbIsgraingrowthEnum
syn keyword cConstant SmbIsalbedoEnum
syn keyword cConstant SmbIsshortwaveEnum
syn keyword cConstant SmbIsthermalEnum
syn keyword cConstant SmbIsaccumulationEnum
syn keyword cConstant SmbIsmeltEnum
syn keyword cConstant SmbIsdensificationEnum
syn keyword cConstant SmbIsturbulentfluxEnum
syn keyword cConstant SMBpddEnum
syn keyword cConstant SmbDelta18oEnum
syn keyword cConstant SmbDelta18oSurfaceEnum
syn keyword cConstant SmbIsdelta18oEnum
syn keyword cConstant SmbIsmungsmEnum
syn keyword cConstant SmbIsd18opdEnum
syn keyword cConstant SmbPrecipitationsPresentdayEnum
syn keyword cConstant SmbPrecipitationsLgmEnum
syn keyword cConstant SmbTemperaturesPresentdayEnum
syn keyword cConstant SmbTemperaturesLgmEnum
syn keyword cConstant SmbPrecipitationEnum
syn keyword cConstant SmbDesfacEnum
syn keyword cConstant SmbS0pEnum
syn keyword cConstant SmbS0tEnum
syn keyword cConstant SmbRlapsEnum
syn keyword cConstant SmbRlapslgmEnum
syn keyword cConstant SmbPfacEnum
syn keyword cConstant SmbTdiffEnum
syn keyword cConstant SmbSealevEnum
syn keyword cConstant SMBd18opddEnum
syn keyword cConstant SmbDpermilEnum
syn keyword cConstant SMBgradientsEnum
syn keyword cConstant SmbMonthlytemperaturesEnum
syn keyword cConstant SmbHrefEnum
syn keyword cConstant SmbSmbrefEnum
syn keyword cConstant SmbBPosEnum
syn keyword cConstant SmbBNegEnum
syn keyword cConstant SMBhenningEnum
syn keyword cConstant SMBcomponentsEnum
syn keyword cConstant SmbAccumulationEnum
syn keyword cConstant SmbEvaporationEnum
syn keyword cConstant SmbRunoffEnum
syn keyword cConstant SMBmeltcomponentsEnum
syn keyword cConstant SmbMeltEnum
syn keyword cConstant SmbRefreezeEnum
syn keyword cConstant SMBgcmEnum
syn keyword cConstant SmbIspddEnum
syn keyword cConstant SmbIssmbgradientsEnum
syn keyword cConstant SolutionTypeEnum
syn keyword cConstant AnalysisTypeEnum
syn keyword cConstant ConfigurationTypeEnum
syn keyword cConstant AdjointBalancethicknessAnalysisEnum
syn keyword cConstant AdjointBalancethickness2AnalysisEnum
syn keyword cConstant AdjointHorizAnalysisEnum
syn keyword cConstant AnalysisCounterEnum
syn keyword cConstant DefaultAnalysisEnum
syn keyword cConstant BalancethicknessAnalysisEnum
syn keyword cConstant BalancethicknessSolutionEnum
syn keyword cConstant Balancethickness2AnalysisEnum
syn keyword cConstant Balancethickness2SolutionEnum
syn keyword cConstant BalancethicknessSoftAnalysisEnum
syn keyword cConstant BalancethicknessSoftSolutionEnum
syn keyword cConstant BalancevelocityAnalysisEnum
syn keyword cConstant BalancevelocitySolutionEnum
syn keyword cConstant L2ProjectionEPLAnalysisEnum
syn keyword cConstant L2ProjectionBaseAnalysisEnum
syn keyword cConstant BedSlopeSolutionEnum
syn keyword cConstant DamageEvolutionSolutionEnum
syn keyword cConstant DamageEvolutionAnalysisEnum
syn keyword cConstant StressbalanceAnalysisEnum
syn keyword cConstant StressbalanceSIAAnalysisEnum
syn keyword cConstant StressbalanceSolutionEnum
syn keyword cConstant StressbalanceVerticalAnalysisEnum
syn keyword cConstant EnthalpyAnalysisEnum
syn keyword cConstant FlaimAnalysisEnum
syn keyword cConstant FlaimSolutionEnum
syn keyword cConstant HydrologyShreveAnalysisEnum
syn keyword cConstant HydrologyDCInefficientAnalysisEnum
syn keyword cConstant HydrologyDCEfficientAnalysisEnum
syn keyword cConstant HydrologySommersAnalysisEnum
syn keyword cConstant HydrologySolutionEnum
syn keyword cConstant MeltingAnalysisEnum
syn keyword cConstant MasstransportAnalysisEnum
syn keyword cConstant MasstransportSolutionEnum
syn keyword cConstant FreeSurfaceBaseAnalysisEnum
syn keyword cConstant FreeSurfaceTopAnalysisEnum
syn keyword cConstant SurfaceNormalVelocityEnum
syn keyword cConstant ExtrudeFromBaseAnalysisEnum
syn keyword cConstant ExtrudeFromTopAnalysisEnum
syn keyword cConstant DepthAverageAnalysisEnum
syn keyword cConstant SteadystateSolutionEnum
syn keyword cConstant SurfaceSlopeSolutionEnum
syn keyword cConstant SmoothAnalysisEnum
syn keyword cConstant ThermalAnalysisEnum
syn keyword cConstant ThermalSolutionEnum
syn keyword cConstant TransientSolutionEnum
syn keyword cConstant UzawaPressureAnalysisEnum
syn keyword cConstant GiaSolutionEnum
syn keyword cConstant GiaAnalysisEnum
syn keyword cConstant MeshdeformationSolutionEnum
syn keyword cConstant MeshdeformationAnalysisEnum
syn keyword cConstant LevelsetAnalysisEnum
syn keyword cConstant LevelsetStabilizationEnum
syn keyword cConstant ExtrapolationAnalysisEnum
syn keyword cConstant LsfReinitializationAnalysisEnum
syn keyword cConstant ApproximationEnum
syn keyword cConstant NoneApproximationEnum
syn keyword cConstant SIAApproximationEnum
syn keyword cConstant SSAApproximationEnum
syn keyword cConstant SSAHOApproximationEnum
syn keyword cConstant SSAFSApproximationEnum
syn keyword cConstant L1L2ApproximationEnum
syn keyword cConstant HOApproximationEnum
syn keyword cConstant HOFSApproximationEnum
syn keyword cConstant FSApproximationEnum
syn keyword cConstant FSvelocityEnum
syn keyword cConstant FSpressureEnum
syn keyword cConstant DataSetEnum
syn keyword cConstant ConstraintsEnum
syn keyword cConstant LoadsEnum
syn keyword cConstant MaterialsEnum
syn keyword cConstant NodesEnum
syn keyword cConstant ContoursEnum
syn keyword cConstant ParametersEnum
syn keyword cConstant VerticesEnum
syn keyword cConstant ResultsEnum
syn keyword cConstant GenericParamEnum
syn keyword cConstant AdolcParamEnum
syn keyword cConstant BoolInputEnum
syn keyword cConstant BoolParamEnum
syn keyword cConstant ContourEnum
syn keyword cConstant ControlInputEnum
syn keyword cConstant DatasetInputEnum
syn keyword cConstant DoubleInputEnum
syn keyword cConstant DoubleArrayInputEnum
syn keyword cConstant DataSetParamEnum
syn keyword cConstant DoubleMatArrayParamEnum
syn keyword cConstant DoubleMatParamEnum
syn keyword cConstant DoubleParamEnum
syn keyword cConstant DoubleVecParamEnum
syn keyword cConstant ElementEnum
syn keyword cConstant ElementHookEnum
syn keyword cConstant HookEnum
syn keyword cConstant ExternalResultEnum
syn keyword cConstant FileParamEnum
syn keyword cConstant InputEnum
syn keyword cConstant IntInputEnum
syn keyword cConstant InputToExtrudeEnum
syn keyword cConstant InputToL2ProjectEnum
syn keyword cConstant InputToDepthaverageEnum
syn keyword cConstant InputToSmoothEnum
syn keyword cConstant SmoothThicknessMultiplierEnum
syn keyword cConstant IntParamEnum
syn keyword cConstant IntVecParamEnum
syn keyword cConstant TransientParamEnum
syn keyword cConstant MaticeEnum
syn keyword cConstant MatdamageiceEnum
syn keyword cConstant MatestarEnum
syn keyword cConstant MatparEnum
syn keyword cConstant NodeEnum
syn keyword cConstant NumericalfluxEnum
syn keyword cConstant NumericalfluxTypeEnum
syn keyword cConstant NeumannfluxEnum
syn keyword cConstant ParamEnum
syn keyword cConstant MoulinEnum
syn keyword cConstant PengridEnum
syn keyword cConstant PenpairEnum
syn keyword cConstant ProfilerEnum
syn keyword cConstant MatrixParamEnum
syn keyword cConstant MassconEnum
syn keyword cConstant MassconNameEnum
syn keyword cConstant MassconDefinitionenumEnum
syn keyword cConstant MassconLevelsetEnum
syn keyword cConstant MassconaxpbyEnum
syn keyword cConstant MassconaxpbyNameEnum
syn keyword cConstant MassconaxpbyDefinitionenumEnum
syn keyword cConstant MassconaxpbyNamexEnum
syn keyword cConstant MassconaxpbyNameyEnum
syn keyword cConstant MassconaxpbyAlphaEnum
syn keyword cConstant MassconaxpbyBetaEnum
syn keyword cConstant NodeSIdEnum
syn keyword cConstant VectorParamEnum
syn keyword cConstant RiftfrontEnum
syn keyword cConstant RiftfrontTypeEnum
syn keyword cConstant SegmentEnum
syn keyword cConstant SegmentRiftfrontEnum
syn keyword cConstant SpcDynamicEnum
syn keyword cConstant SpcStaticEnum
syn keyword cConstant SpcTransientEnum
syn keyword cConstant StringArrayParamEnum
syn keyword cConstant StringParamEnum
syn keyword cConstant SegEnum
syn keyword cConstant SegInputEnum
syn keyword cConstant TriaEnum
syn keyword cConstant TriaInputEnum
syn keyword cConstant TetraEnum
syn keyword cConstant TetraInputEnum
syn keyword cConstant PentaEnum
syn keyword cConstant PentaInputEnum
syn keyword cConstant VertexEnum
syn keyword cConstant VertexPIdEnum
syn keyword cConstant VertexSIdEnum
syn keyword cConstant AirEnum
syn keyword cConstant IceEnum
syn keyword cConstant MelangeEnum
syn keyword cConstant WaterEnum
syn keyword cConstant ClosedEnum
syn keyword cConstant FreeEnum
syn keyword cConstant OpenEnum
syn keyword cConstant AdjointpEnum
syn keyword cConstant AdjointxEnum
syn keyword cConstant AdjointyEnum
syn keyword cConstant AdjointzEnum
syn keyword cConstant BalancethicknessMisfitEnum
syn keyword cConstant BedSlopeXEnum
syn keyword cConstant BedSlopeYEnum
syn keyword cConstant BoundaryEnum
syn keyword cConstant ConvergedEnum
syn keyword cConstant FillEnum
syn keyword cConstant FractionIncrementEnum
syn keyword cConstant FrictionEnum
syn keyword cConstant InternalEnum
syn keyword cConstant MassFluxEnum
syn keyword cConstant MeltingOffsetEnum
syn keyword cConstant MisfitEnum
syn keyword cConstant PressureEnum
syn keyword cConstant PressurePicardEnum
syn keyword cConstant AndroidFrictionCoefficientEnum
syn keyword cConstant ResetPenaltiesEnum
syn keyword cConstant SegmentOnIceShelfEnum
syn keyword cConstant SurfaceAbsVelMisfitEnum
syn keyword cConstant SurfaceAreaEnum
syn keyword cConstant SurfaceAverageVelMisfitEnum
syn keyword cConstant SurfaceLogVelMisfitEnum
syn keyword cConstant SurfaceLogVxVyMisfitEnum
syn keyword cConstant SurfaceRelVelMisfitEnum
syn keyword cConstant SurfaceSlopeXEnum
syn keyword cConstant SurfaceSlopeYEnum
syn keyword cConstant TemperatureEnum
syn keyword cConstant TemperaturePicardEnum
syn keyword cConstant TemperaturePDDEnum
syn keyword cConstant ThicknessAbsMisfitEnum
syn keyword cConstant SurfaceAbsMisfitEnum
syn keyword cConstant VelEnum
syn keyword cConstant VelocityEnum
syn keyword cConstant VxAverageEnum
syn keyword cConstant VxEnum
syn keyword cConstant VxPicardEnum
syn keyword cConstant VyAverageEnum
syn keyword cConstant VyEnum
syn keyword cConstant VyPicardEnum
syn keyword cConstant VzEnum
syn keyword cConstant VzSSAEnum
syn keyword cConstant VzHOEnum
syn keyword cConstant VzPicardEnum
syn keyword cConstant VzFSEnum
syn keyword cConstant VxMeshEnum
syn keyword cConstant VyMeshEnum
syn keyword cConstant VzMeshEnum
syn keyword cConstant EnthalpyEnum
syn keyword cConstant EnthalpyPicardEnum
syn keyword cConstant ThicknessAbsGradientEnum
syn keyword cConstant ThicknessAlongGradientEnum
syn keyword cConstant ThicknessAcrossGradientEnum
syn keyword cConstant ThicknessPositiveEnum
syn keyword cConstant IntMatParamEnum
syn keyword cConstant RheologyBbarAbsGradientEnum
syn keyword cConstant RheologyBAbsGradientEnum
syn keyword cConstant DragCoefficientAbsGradientEnum
syn keyword cConstant TransientInputEnum
syn keyword cConstant WaterfractionEnum
syn keyword cConstant WatercolumnEnum
syn keyword cConstant BasalFrictionEnum
syn keyword cConstant ViscousHeatingEnum
syn keyword cConstant HydrologyWaterVxEnum
syn keyword cConstant HydrologyWaterVyEnum
syn keyword cConstant DrivingStressXEnum
syn keyword cConstant DrivingStressYEnum
syn keyword cConstant SigmaNNEnum
syn keyword cConstant StressTensorEnum
syn keyword cConstant StressTensorxxEnum
syn keyword cConstant StressTensorxyEnum
syn keyword cConstant StressTensorxzEnum
syn keyword cConstant StressTensoryyEnum
syn keyword cConstant StressTensoryzEnum
syn keyword cConstant StressTensorzzEnum
syn keyword cConstant StressMaxPrincipalEnum
syn keyword cConstant DeviatoricStressEnum
syn keyword cConstant DeviatoricStressxxEnum
syn keyword cConstant DeviatoricStressxyEnum
syn keyword cConstant DeviatoricStressxzEnum
syn keyword cConstant DeviatoricStressyyEnum
syn keyword cConstant DeviatoricStressyzEnum
syn keyword cConstant DeviatoricStresszzEnum
syn keyword cConstant DeviatoricStresseffectiveEnum
syn keyword cConstant LambdaSEnum
syn keyword cConstant StrainRateEnum
syn keyword cConstant StrainRatexxEnum
syn keyword cConstant StrainRatexyEnum
syn keyword cConstant StrainRatexzEnum
syn keyword cConstant StrainRateyyEnum
syn keyword cConstant StrainRateyzEnum
syn keyword cConstant StrainRatezzEnum
syn keyword cConstant DivergenceEnum
syn keyword cConstant MaxDivergenceEnum
syn keyword cConstant GiaCrossSectionShapeEnum
syn keyword cConstant GiadWdtEnum
syn keyword cConstant GiaWEnum
syn keyword cConstant P0Enum
syn keyword cConstant P0ArrayEnum
syn keyword cConstant P1Enum
syn keyword cConstant P1DGEnum
syn keyword cConstant P1bubbleEnum
syn keyword cConstant P1bubblecondensedEnum
syn keyword cConstant P2Enum
syn keyword cConstant P2bubbleEnum
syn keyword cConstant P2bubblecondensedEnum
syn keyword cConstant P2xP1Enum
syn keyword cConstant P1xP2Enum
syn keyword cConstant P1xP3Enum
syn keyword cConstant P2xP4Enum
syn keyword cConstant P1P1Enum
syn keyword cConstant P1P1GLSEnum
syn keyword cConstant MINIEnum
syn keyword cConstant MINIcondensedEnum
syn keyword cConstant TaylorHoodEnum
syn keyword cConstant LATaylorHoodEnum
syn keyword cConstant XTaylorHoodEnum
syn keyword cConstant OneLayerP4zEnum
syn keyword cConstant CrouzeixRaviartEnum
syn keyword cConstant LACrouzeixRaviartEnum
syn keyword cConstant SaveResultsEnum
syn keyword cConstant BoolExternalResultEnum
syn keyword cConstant DoubleExternalResultEnum
syn keyword cConstant DoubleMatExternalResultEnum
syn keyword cConstant IntExternalResultEnum
syn keyword cConstant JEnum
syn keyword cConstant StringExternalResultEnum
syn keyword cConstant StepEnum
syn keyword cConstant TimeEnum
syn keyword cConstant WaterColumnOldEnum
syn keyword cConstant OutputdefinitionEnum
syn keyword cConstant Outputdefinition1Enum
syn keyword cConstant Outputdefinition2Enum
syn keyword cConstant Outputdefinition3Enum
syn keyword cConstant Outputdefinition4Enum
syn keyword cConstant Outputdefinition5Enum
syn keyword cConstant Outputdefinition6Enum
syn keyword cConstant Outputdefinition7Enum
syn keyword cConstant Outputdefinition8Enum
syn keyword cConstant Outputdefinition9Enum
syn keyword cConstant Outputdefinition10Enum
syn keyword cConstant Outputdefinition11Enum
syn keyword cConstant Outputdefinition12Enum
syn keyword cConstant Outputdefinition13Enum
syn keyword cConstant Outputdefinition14Enum
syn keyword cConstant Outputdefinition15Enum
syn keyword cConstant Outputdefinition16Enum
syn keyword cConstant Outputdefinition17Enum
syn keyword cConstant Outputdefinition18Enum
syn keyword cConstant Outputdefinition19Enum
syn keyword cConstant Outputdefinition20Enum
syn keyword cConstant Outputdefinition21Enum
syn keyword cConstant Outputdefinition22Enum
syn keyword cConstant Outputdefinition23Enum
syn keyword cConstant Outputdefinition24Enum
syn keyword cConstant Outputdefinition25Enum
syn keyword cConstant Outputdefinition26Enum
syn keyword cConstant Outputdefinition27Enum
syn keyword cConstant Outputdefinition28Enum
syn keyword cConstant Outputdefinition29Enum
syn keyword cConstant Outputdefinition30Enum
syn keyword cConstant Outputdefinition31Enum
syn keyword cConstant Outputdefinition32Enum
syn keyword cConstant Outputdefinition33Enum
syn keyword cConstant Outputdefinition34Enum
syn keyword cConstant Outputdefinition35Enum
syn keyword cConstant Outputdefinition36Enum
syn keyword cConstant Outputdefinition37Enum
syn keyword cConstant Outputdefinition38Enum
syn keyword cConstant Outputdefinition39Enum
syn keyword cConstant Outputdefinition40Enum
syn keyword cConstant Outputdefinition41Enum
syn keyword cConstant Outputdefinition42Enum
syn keyword cConstant Outputdefinition43Enum
syn keyword cConstant Outputdefinition44Enum
syn keyword cConstant Outputdefinition45Enum
syn keyword cConstant Outputdefinition46Enum
syn keyword cConstant Outputdefinition47Enum
syn keyword cConstant Outputdefinition48Enum
syn keyword cConstant Outputdefinition49Enum
syn keyword cConstant Outputdefinition50Enum
syn keyword cConstant Outputdefinition51Enum
syn keyword cConstant Outputdefinition52Enum
syn keyword cConstant Outputdefinition53Enum
syn keyword cConstant Outputdefinition54Enum
syn keyword cConstant Outputdefinition55Enum
syn keyword cConstant Outputdefinition56Enum
syn keyword cConstant Outputdefinition57Enum
syn keyword cConstant Outputdefinition58Enum
syn keyword cConstant Outputdefinition59Enum
syn keyword cConstant Outputdefinition60Enum
syn keyword cConstant Outputdefinition61Enum
syn keyword cConstant Outputdefinition62Enum
syn keyword cConstant Outputdefinition63Enum
syn keyword cConstant Outputdefinition64Enum
syn keyword cConstant Outputdefinition65Enum
syn keyword cConstant Outputdefinition66Enum
syn keyword cConstant Outputdefinition67Enum
syn keyword cConstant Outputdefinition68Enum
syn keyword cConstant Outputdefinition69Enum
syn keyword cConstant Outputdefinition70Enum
syn keyword cConstant Outputdefinition71Enum
syn keyword cConstant Outputdefinition72Enum
syn keyword cConstant Outputdefinition73Enum
syn keyword cConstant Outputdefinition74Enum
syn keyword cConstant Outputdefinition75Enum
syn keyword cConstant Outputdefinition76Enum
syn keyword cConstant Outputdefinition77Enum
syn keyword cConstant Outputdefinition78Enum
syn keyword cConstant Outputdefinition79Enum
syn keyword cConstant Outputdefinition80Enum
syn keyword cConstant Outputdefinition81Enum
syn keyword cConstant Outputdefinition82Enum
syn keyword cConstant Outputdefinition83Enum
syn keyword cConstant Outputdefinition84Enum
syn keyword cConstant Outputdefinition85Enum
syn keyword cConstant Outputdefinition86Enum
syn keyword cConstant Outputdefinition87Enum
syn keyword cConstant Outputdefinition88Enum
syn keyword cConstant Outputdefinition89Enum
syn keyword cConstant Outputdefinition90Enum
syn keyword cConstant Outputdefinition91Enum
syn keyword cConstant Outputdefinition92Enum
syn keyword cConstant Outputdefinition93Enum
syn keyword cConstant Outputdefinition94Enum
syn keyword cConstant Outputdefinition95Enum
syn keyword cConstant Outputdefinition96Enum
syn keyword cConstant Outputdefinition97Enum
syn keyword cConstant Outputdefinition98Enum
syn keyword cConstant Outputdefinition99Enum
syn keyword cConstant Outputdefinition100Enum
syn keyword cConstant OutputdefinitionListEnum
syn keyword cConstant MassfluxatgateEnum
syn keyword cConstant MassfluxatgateNameEnum
syn keyword cConstant MassfluxatgateDefinitionenumEnum
syn keyword cConstant MassfluxatgateSegmentsEnum
syn keyword cConstant NodalvalueEnum
syn keyword cConstant NodalvalueNameEnum
syn keyword cConstant NodalvalueDefinitionenumEnum
syn keyword cConstant NodalvalueModelEnumEnum
syn keyword cConstant NodalvalueNodeEnum
syn keyword cConstant MisfitNameEnum
syn keyword cConstant MisfitDefinitionenumEnum
syn keyword cConstant MisfitModelEnumEnum
syn keyword cConstant MisfitObservationEnum
syn keyword cConstant MisfitObservationEnumEnum
syn keyword cConstant MisfitLocalEnum
syn keyword cConstant MisfitTimeinterpolationEnum
syn keyword cConstant MisfitWeightsEnum
syn keyword cConstant MisfitWeightsEnumEnum
syn keyword cConstant SurfaceObservationEnum
syn keyword cConstant WeightsSurfaceObservationEnum
syn keyword cConstant VxObsEnum
syn keyword cConstant WeightsVxObsEnum
syn keyword cConstant VyObsEnum
syn keyword cConstant WeightsVyObsEnum
syn keyword cConstant MinVelEnum
syn keyword cConstant MaxVelEnum
syn keyword cConstant MinVxEnum
syn keyword cConstant MaxVxEnum
syn keyword cConstant MaxAbsVxEnum
syn keyword cConstant MinVyEnum
syn keyword cConstant MaxVyEnum
syn keyword cConstant MaxAbsVyEnum
syn keyword cConstant MinVzEnum
syn keyword cConstant MaxVzEnum
syn keyword cConstant MaxAbsVzEnum
syn keyword cConstant FloatingAreaEnum
syn keyword cConstant GroundedAreaEnum
syn keyword cConstant IceMassEnum
syn keyword cConstant IceVolumeEnum
syn keyword cConstant IceVolumeAboveFloatationEnum
syn keyword cConstant TotalFloatingBmbEnum
syn keyword cConstant TotalGroundedBmbEnum
syn keyword cConstant TotalSmbEnum
syn keyword cConstant AbsoluteEnum
syn keyword cConstant IncrementalEnum
syn keyword cConstant AugmentedLagrangianREnum
syn keyword cConstant AugmentedLagrangianRhopEnum
syn keyword cConstant AugmentedLagrangianRlambdaEnum
syn keyword cConstant AugmentedLagrangianRholambdaEnum
syn keyword cConstant AugmentedLagrangianThetaEnum
syn keyword cConstant NoneEnum
syn keyword cConstant AggressiveMigrationEnum
syn keyword cConstant SoftMigrationEnum
syn keyword cConstant SubelementMigrationEnum
syn keyword cConstant SubelementMigration2Enum
syn keyword cConstant ContactEnum
syn keyword cConstant GroundingOnlyEnum
syn keyword cConstant MaskGroundediceLevelsetEnum
syn keyword cConstant GaussSegEnum
syn keyword cConstant GaussTriaEnum
syn keyword cConstant GaussTetraEnum
syn keyword cConstant GaussPentaEnum
syn keyword cConstant FSSolverEnum
syn keyword cConstant AdjointEnum
syn keyword cConstant ColinearEnum
syn keyword cConstant ControlSteadyEnum
syn keyword cConstant FsetEnum
syn keyword cConstant Gradient1Enum
syn keyword cConstant Gradient2Enum
syn keyword cConstant Gradient3Enum
syn keyword cConstant GradientEnum
syn keyword cConstant GroundinglineMigrationEnum
syn keyword cConstant GsetEnum
syn keyword cConstant IndexEnum
syn keyword cConstant IndexedEnum
syn keyword cConstant IntersectEnum
syn keyword cConstant NodalEnum
syn keyword cConstant OldGradientEnum
syn keyword cConstant OutputBufferPointerEnum
syn keyword cConstant OutputBufferSizePointerEnum
syn keyword cConstant OutputFilePointerEnum
syn keyword cConstant ToolkitsFileNameEnum
syn keyword cConstant RootPathEnum
syn keyword cConstant OutputFileNameEnum
syn keyword cConstant InputFileNameEnum
syn keyword cConstant LockFileNameEnum
syn keyword cConstant RestartFileNameEnum
syn keyword cConstant ToolkitsOptionsAnalysesEnum
syn keyword cConstant ToolkitsOptionsStringsEnum
syn keyword cConstant QmuErrNameEnum
syn keyword cConstant QmuInNameEnum
syn keyword cConstant QmuOutNameEnum
syn keyword cConstant RegularEnum
syn keyword cConstant ScaledEnum
syn keyword cConstant SeparateEnum
syn keyword cConstant SsetEnum
syn keyword cConstant VerboseEnum
syn keyword cConstant TriangleInterpEnum
syn keyword cConstant BilinearInterpEnum
syn keyword cConstant NearestInterpEnum
syn keyword cConstant XYEnum
syn keyword cConstant XYZEnum
syn keyword cConstant DenseEnum
syn keyword cConstant MpiDenseEnum
syn keyword cConstant MpiSparseEnum
syn keyword cConstant SeqEnum
syn keyword cConstant MpiEnum
syn keyword cConstant MumpsEnum
syn keyword cConstant GslEnum
syn keyword cConstant OptionEnum
syn keyword cConstant GenericOptionEnum
syn keyword cConstant OptionCellEnum
syn keyword cConstant OptionStructEnum
syn keyword cConstant CuffeyEnum
syn keyword cConstant BuddJackaEnum
syn keyword cConstant CuffeyTemperateEnum
syn keyword cConstant PatersonEnum
syn keyword cConstant ArrheniusEnum
syn keyword cConstant LliboutryDuvalEnum
syn keyword cConstant SpclevelsetEnum
syn keyword cConstant ExtrapolationVariableEnum
syn keyword cConstant IceMaskNodeActivationEnum
syn keyword cConstant LevelsetfunctionSlopeXEnum
syn keyword cConstant LevelsetfunctionSlopeYEnum
syn keyword cConstant LevelsetfunctionPicardEnum
syn keyword cConstant LevelsetReinitFrequencyEnum
syn keyword cConstant SealevelriseSolutionEnum
syn keyword cConstant SealevelriseAnalysisEnum
syn keyword cConstant SealevelEnum
syn keyword cConstant SealevelEustaticEnum
syn keyword cConstant SealevelriseDeltathicknessEnum
syn keyword cConstant SealevelriseMaxiterEnum
syn keyword cConstant SealevelriseReltolEnum
syn keyword cConstant SealevelriseAbstolEnum
syn keyword cConstant SealevelriseLoveHEnum
syn keyword cConstant SealevelriseLoveKEnum
syn keyword cConstant SealevelriseTideLoveHEnum
syn keyword cConstant SealevelriseTideLoveKEnum
syn keyword cConstant SealevelriseRigidEnum
syn keyword cConstant SealevelriseElasticEnum
syn keyword cConstant SealevelriseRotationEnum
syn keyword cConstant SealevelriseGElasticEnum
syn keyword cConstant SealevelriseDegaccEnum
syn keyword cConstant SealevelriseTransitionsEnum
syn keyword cConstant SealevelriseRequestedOutputsEnum
syn keyword cConstant SealevelriseNumRequestedOutputsEnum
syn keyword cConstant ParametersENDEnum
"ISSM's Enums end
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

if !exists("c_no_ansi") || exists("c_ansi_typedefs")
	syn keyword   cType		size_t ssize_t wchar_t ptrdiff_t sig_atomic_t fpos_t
	syn keyword   cType		clock_t time_t va_list jmp_buf FILE DIR div_t ldiv_t
	syn keyword   cType		mbstate_t wctrans_t wint_t wctype_t
	syn keyword   cType		lldiv_t float_t double_t fenv_t fexcept_t
endif
if !exists("c_no_c99") " ISO C99
	syn keyword	cType		bool complex imaginary
	syn keyword	cType		int8_t int16_t int32_t int64_t
	syn keyword	cType		uint8_t uint16_t uint32_t uint64_t
	syn keyword	cType		int_least8_t int_least16_t int_least32_t int_least64_t
	syn keyword	cType		uint_least8_t uint_least16_t uint_least32_t uint_least64_t
	syn keyword	cType		int_fast8_t int_fast16_t int_fast32_t int_fast64_t
	syn keyword	cType		uint_fast8_t uint_fast16_t uint_fast32_t uint_fast64_t
	syn keyword	cType		intptr_t uintptr_t
	syn keyword	cType		intmax_t uintmax_t
endif
if exists("c_gnu")
	syn keyword	cType		__label__ __complex__ __volatile__
endif

syn keyword	cStructure	struct union enum typedef
syn keyword	cStorageClass	static register auto volatile extern const VOL
if exists("c_gnu")
	syn keyword	cStorageClass	inline __attribute__
endif
if !exists("c_no_c99")
	syn keyword	cStorageClass	inline restrict
endif

if !exists("c_no_ansi") || exists("c_ansi_constants") || exists("c_gnu")
	if exists("c_gnu")
		syn keyword cConstant __GNUC__ __FUNCTION__ __PRETTY_FUNCTION__
	endif
	syn keyword cConstant __LINE__ __FILE__ __DATE__ __TIME__ __STDC__ __func__
	syn keyword cConstant __STDC_VERSION__
	syn keyword cConstant CHAR_BIT MB_LEN_MAX MB_CUR_MAX
	syn keyword cConstant UCHAR_MAX UINT_MAX ULONG_MAX USHRT_MAX
	syn keyword cConstant CHAR_MIN INT_MIN LONG_MIN SHRT_MIN
	syn keyword cConstant CHAR_MAX INT_MAX LONG_MAX SHRT_MAX
	syn keyword cConstant SCHAR_MIN SINT_MIN SLONG_MIN SSHRT_MIN
	syn keyword cConstant SCHAR_MAX SINT_MAX SLONG_MAX SSHRT_MAX
	if !exists("c_no_c99")
		syn keyword cConstant LLONG_MIN LLONG_MAX ULLONG_MAX
		syn keyword cConstant INT8_MIN INT16_MIN INT32_MIN INT64_MIN
		syn keyword cConstant INT8_MAX INT16_MAX INT32_MAX INT64_MAX
		syn keyword cConstant UINT8_MAX UINT16_MAX UINT32_MAX UINT64_MAX
		syn keyword cConstant INT_LEAST8_MIN INT_LEAST16_MIN INT_LEAST32_MIN INT_LEAST64_MIN
		syn keyword cConstant INT_LEAST8_MAX INT_LEAST16_MAX INT_LEAST32_MAX INT_LEAST64_MAX
		syn keyword cConstant UINT_LEAST8_MAX UINT_LEAST16_MAX UINT_LEAST32_MAX UINT_LEAST64_MAX
		syn keyword cConstant INT_FAST8_MIN INT_FAST16_MIN INT_FAST32_MIN INT_FAST64_MIN
		syn keyword cConstant INT_FAST8_MAX INT_FAST16_MAX INT_FAST32_MAX INT_FAST64_MAX
		syn keyword cConstant UINT_FAST8_MAX UINT_FAST16_MAX UINT_FAST32_MAX UINT_FAST64_MAX
		syn keyword cConstant INTPTR_MIN INTPTR_MAX UINTPTR_MAX
		syn keyword cConstant INTMAX_MIN INTMAX_MAX UINTMAX_MAX
		syn keyword cConstant PTRDIFF_MIN PTRDIFF_MAX SIG_ATOMIC_MIN SIG_ATOMIC_MAX
		syn keyword cConstant SIZE_MAX WCHAR_MIN WCHAR_MAX WINT_MIN WINT_MAX
	endif
	syn keyword cConstant FLT_RADIX FLT_ROUNDS
	syn keyword cConstant FLT_DIG FLT_MANT_DIG FLT_EPSILON
	syn keyword cConstant DBL_DIG DBL_MANT_DIG DBL_EPSILON
	syn keyword cConstant LDBL_DIG LDBL_MANT_DIG LDBL_EPSILON
	syn keyword cConstant FLT_MIN FLT_MAX FLT_MIN_EXP FLT_MAX_EXP
	syn keyword cConstant FLT_MIN_10_EXP FLT_MAX_10_EXP
	syn keyword cConstant DBL_MIN DBL_MAX DBL_MIN_EXP DBL_MAX_EXP
	syn keyword cConstant DBL_MIN_10_EXP DBL_MAX_10_EXP
	syn keyword cConstant LDBL_MIN LDBL_MAX LDBL_MIN_EXP LDBL_MAX_EXP
	syn keyword cConstant LDBL_MIN_10_EXP LDBL_MAX_10_EXP
	syn keyword cConstant HUGE_VAL CLOCKS_PER_SEC NULL _NO_LEAP_SECONDS _LOCALTIME
	syn keyword cConstant LC_ALL LC_COLLATE LC_CTYPE LC_MONETARY
	syn keyword cConstant LC_NUMERIC LC_TIME
	" syn keyword cConstant SIG_DFL SIG_ERR SIG_IGN
	syn keyword cAnsiFuncPtr SIG_DFL SIG_ERR SIG_IGN
	syn keyword cConstant SIGABRT SIGFPE SIGILL SIGHUP SIGINT SIGSEGV SIGTERM
	syn keyword cConstant INFINITY     FP_SUBNORMAL FP_ILOGB0
	syn keyword cConstant NAN          FP_ZERO      FP_ILOGBNAN
	syn keyword cConstant FP_INFINITE  FP_FAST_FMA
	syn keyword cConstant HUGE_VALF    FP_NAN       FP_FAST_FMAF
	syn keyword cConstant HUGE_VALL    FP_NORMAL    FP_FAST_FMAL
	syn keyword cConstant FE_OVERFLOW      FE_TOWARDZERO
	syn keyword cConstant FE_UNDERFLOW     FE_UPWARD
	syn keyword cConstant FE_DIVBYZERO     FE_ALL_EXCEPT    FE_DFL_ENV
	syn keyword cConstant FE_INEXACT       FE_DOWNWARD
	syn keyword cConstant FE_INVALID       FE_TONEAREST
	syn keyword cConstant	_Complex_I _Imaginary_I 
	" Add POSIX signals as well...
	syn keyword cConstant SIGABRT SIGALRM SIGCHLD SIGCONT SIGFPE SIGHUP
	syn keyword cConstant SIGILL SIGINT SIGKILL SIGPIPE SIGQUIT SIGSEGV
	syn keyword cConstant SIGSTOP SIGTERM SIGTRAP SIGTSTP SIGTTIN SIGTTOU
	syn keyword cConstant SIGUSR1 SIGUSR2
	syn keyword cConstant _IOFBF _IOLBF _IONBF BUFSIZ EOF WEOF
	syn keyword cConstant FOPEN_MAX FILENAME_MAX L_tmpnam
	syn keyword cConstant SEEK_CUR SEEK_END SEEK_SET
	syn keyword cConstant TMP_MAX stderr stdin stdout
	syn keyword cConstant EXIT_FAILURE EXIT_SUCCESS RAND_MAX
	" Add POSIX errors as well
	syn keyword cConstant E2BIG EACCES EAGAIN EBADF EBADMSG EBUSY
	syn keyword cConstant ECANCELED ECHILD EDEADLK EDOM EEXIST EFAULT
	syn keyword cConstant EFBIG EILSEQ EINPROGRESS EINTR EINVAL EIO EISDIR
	syn keyword cConstant EMFILE EMLINK EMSGSIZE ENAMETOOLONG ENFILE ENODEV
	syn keyword cConstant ENOENT ENOEXEC ENOLCK ENOMEM ENOSPC ENOSYS
	syn keyword cConstant ENOTDIR ENOTEMPTY ENOTSUP ENOTTY ENXIO EPERM
	syn keyword cConstant EPIPE ERANGE EROFS ESPIPE ESRCH ETIMEDOUT EXDEV
	" math.h
	syn keyword cConstant M_E M_LOG2E M_LOG10E M_LN2 M_LN10 M_PI M_PI_2 M_PI_4
	syn keyword cConstant M_1_PI M_2_PI M_2_SQRTPI M_SQRT2 M_SQRT1_2
endif
if !exists("c_no_c99") " ISO C99
	syn keyword cConstant true false __bool_true_false_are_defined
endif

" Accept %: for # (C99)
syn region	cPreCondit	start="^\s*\(%:\|#\)\s*\(ifdef\|ifndef\)\>" skip="\\$" end="$" end="//"me=s-1 contains=cComment,cCppString,cCharacter,cCppParen,cParenError,cNumbers,cCommentError,cSpaceError
syn keyword	cDefined	defined contained
"syn match	cDefined	display contained "\<defined(\w\+)" contains=cName,cAnsiName
syn region	cPreConditIf	start="^\s*\(%:\|#\)\s*\(if\|elif\)\>" skip="\\$" end="$" end="//"me=s-1 contains=cDefined,cComment,cCppString,cCharacter,cCppParen,cParenError,cNumbers,cCommentError,cSpaceError
syn match	cPreCondit	display "^\s*\(%:\|#\)\s*\(else\|endif\)\>"
if !exists("c_no_if0")
	syn region	cCppOut		start="^\s*\(%:\|#\)\s*if\s\+0\+\>" end=".\@=\|$" contains=cCppOut2
	syn region	cCppOut2	contained start="0" end="^\s*\(%:\|#\)\s*\(endif\>\|else\>\|elif\>\)" contains=cSpaceError,cCppSkip
		syn region	cCppSkip	contained start="^\s*\(%:\|#\)\s*\(if\>\|ifdef\>\|ifndef\>\)" skip="\\$" end="^\s*\(%:\|#\)\s*endif\>" contains=cSpaceError,cCppSkip
	endif
	syn region	cIncluded	display contained start=+"+ skip=+\\\\\|\\"+ end=+"+
	syn match	cIncluded	display contained "<[^>]*>"
	syn match	cInclude	display "^\s*\(%:\|#\)\s*include\>\s*["<]" contains=cIncluded
	syn match cLineSkip	"\\$"
	syn cluster	cPreProcGroup	contains=cPreConditIf,cPreCondit,cIncluded,cInclude,cDefined,cDefine,cErrInParen,cErrInBracket,cUserLabel,cUserLabel2,cGotoLabel,cSpecial,cOctalZero,cCppOut,cCppOut2,cCppSkip,cFormat,cNumber,cFloat,cOctal,cOctalError,cNumbersCom,cString,cCommentSkip,cCommentString,cComment2String,@cCommentGroup,cCommentStartError,cParen,cBracket,cMulti
	"syn region	cDefine		start="^\s*\(%:\|#\)\s*\(define\|undef\)\>" skip="\\$" end="$" end="//"me=s-1 contains=ALLBUT,@cPreProcGroup,@Spell
	syn region	cDefine		start="^\s*\(%:\|#\)\s*\(define\|undef\)\>" skip="\\$" end="$" end="//"me=s-1 contains=ALLBUT,@cPreProcGroup,cName,cFunction,cAnsiFunction,@Spell
	syn region	cPreProc	start="^\s*\(%:\|#\)\s*\(pragma\>\|line\>\|warning\>\|warn\>\|error\>\)" skip="\\$" end="$" keepend contains=ALLBUT,@cPreProcGroup,@Spell

	" Highlight User Labels
	syn cluster	cMultiGroup	contains=cIncluded,cSpecial,cCommentSkip,cCommentString,cComment2String,@cCommentGroup,cCommentStartError,cUserCont,cUserLabel,cUserLabel2,cGotoLabel,cBitField,cOctalZero,cCppOut,cCppOut2,cCppSkip,cFormat,cNumber,cFloat,cOctal,cOctalError,cNumbersCom,cCppParen,cCppBracket,cCppString
	syn region	cMulti		transparent start='?' skip='::' end=':' contains=ALLBUT,@cMultiGroup,@Spell
	" Avoid matching foo::bar() in C++ by requiring that the next char is not ':'
	syn cluster	cLabelGroup	contains=cUserLabel
	syn match	cUserCont	display "^\s*\I\i*\s*:$" contains=@cLabelGroup
	syn match	cUserCont	display ";\s*\I\i*\s*:$" contains=@cLabelGroup
	syn match	cUserCont	display "^\s*\I\i*\s*:[^:]"me=e-1 contains=@cLabelGroup
	syn match	cUserCont	display ";\s*\I\i*\s*:[^:]"me=e-1 contains=@cLabelGroup

	syn match	cUserLabel	display "\I\i*" contained
	syn match	cUserLabel2	display "\I\i*:;\+"me=e-2
	syn match	cGotoLabel	display "\<goto\s\+\I\i*;"me=e-1,hs=s+5 contains=cGoto
	syn keyword	cGoto		contained goto

	" Avoid recognizing most bitfields as labels
	syn match	cBitField	display "^\s*\I\i*\s*:\s*[1-9]"me=e-1
	syn match	cBitField	display ";\s*\I\i*\s*:\s*[1-9]"me=e-1

	syn match cOperator	"\(<<\|>>\|[-+*/%&^|<>!=]\)="
	syn match cOperator	"<<\|>>\|&&\|||\|++\|--\|->"
	syn match cOperator	"[.!~*&%<>^|=+-]"
	syn match cOperator	"/[^/*=]"me=e-1
	syn match cOperator	"/$"
	syn match cOperator	"[\\]"
	syn match cOperator	"&&\|||"
	syn match cSpecialCharacter	"[,;]"
	syn match cDelimiter    "[][(){}]"
	syn keyword	cType		CHAR byte BYTE uchar ushort uint ulong
	syn keyword	cType		WORD DWORD QWORD INT INT2 INT4 UNS UNS2 UNS4 INT8 UNS8
	syn keyword	cType		CFG_t cfap_t cBYTE cvoid
	syn keyword	cType		_Bool _Complex _Imaginary __int64
	if !exists("c_no_ansi") || exists("c_ansi_typedefs")
		syn keyword   cMC	__near __far FAR __io __direct DIR
	endif

	if exists("c_minlines")
		let b:c_minlines = c_minlines
	else
		if !exists("c_no_if0")
			let b:c_minlines = 50	" #if 0 constructs can be long
		else
			let b:c_minlines = 15	" mostly for () constructs
		endif
	endif
	exec "syn sync ccomment cComment minlines=" . b:c_minlines

	" Define the default highlighting.
	" For version 5.7 and earlier: only when not done already
	" For version 5.8 and later: only when an item doesn't have highlighting yet
	if version >= 508 || !exists("did_c_syn_inits")
		if version < 508
			let did_c_syn_inits = 1
			command -nargs=+ HiLink hi link <args>
		else
			command -nargs=+ HiLink hi def link <args>
		endif

		HiLink cFormat	cSpecial
		HiLink cCppString	cString
		HiLink cCommentL	cComment
		HiLink cCommentStart	cComment
		HiLink cLabel		Label
		" HiLink cUserLabel	Label
		HiLink cUserLabel	UserLabel2
		HiLink cUserLabel2	UserLabel2
		HiLink cGotoLabel	UserLabel2
		HiLink cGoto		Statement
		HiLink cConditional	Conditional
		HiLink cRepeat	Repeat
		HiLink cCharacter	Character
		HiLink cSpecialCharacter cSpecial
		HiLink cNumber	Number
		HiLink cOctal		Number
		HiLink cOctalZero	PreProc		" link this to Error if you want
		HiLink cFloat		Float
		HiLink cOctalError	cError
		HiLink cParenError	cError
		HiLink cErrInParen	cError
		HiLink cErrInBracket	cError
		HiLink cCommentError	cError
		HiLink cCommentStartError	cError
		HiLink cSpaceError	cError
		HiLink cSpecialError	cError
		HiLink cOperator	Operator
		HiLink cOperatorBold	OperatorBold
		HiLink cStructure	Structure
		HiLink cStorageClass	StorageClass
		HiLink cInclude	Include
		HiLink cPreProc	PreProc
		HiLink cDefine	Macro
		HiLink cDefined	PreCondit
		HiLink cIncluded	cString
		HiLink cError		Error
		HiLink cStatement	Statement
		HiLink cPreCondit	PreCondit
		HiLink cPreConditIf	PreCondit
		HiLink cType		Type
		HiLink cConstant	Constant
		HiLink cCommentString cString
		HiLink cComment2String cString
		HiLink cCommentSkip	cComment
		HiLink cString	String
		HiLink cComment	Comment
		HiLink cDelimiter     Delimiter
		HiLink cSpecial	SpecialChar
		HiLink cTodo		Todo
		HiLink cCppSkip	cCppOut
		HiLink cCppOut2	cCppOut
		HiLink cCppOut	Comment
		HiLink cMulti		Operator
		HiLink cMultiMG	Operator
		HiLink cFunction	Function
		HiLink cAnsiFunction	StdFunction
		HiLink cName		Name
		HiLink cBitField	Name
		HiLink cAnsiName	StdName
		"HiLink cBlock	BlockBraces
		HiLink cBraces	BlockBraces
		"HiLink cBraceError	Error
		HiLink cMC		MicroController
		HiLink cAnsiFuncPtr	AnsiFuncPtr

		hi Function		gui=NONE guifg=#e86f00
		"hi StdFunction	gui=bold guifg=#ee0040
		hi StdFunction	gui=bold guifg=#e86f00
		hi Statement		gui=bold guifg=#a06129
		hi UserLabel2		gui=bold guifg=#c96129
		hi Operator		gui=NONE guifg=#000000
		hi OperatorBold	gui=bold guifg=#000000
		hi StdName		gui=bold guifg=#5276e6
		hi Name		gui=NONE guifg=#5276e6
		hi BlockBraces	gui=bold guifg=#000000
		hi Special		gui=NONE guifg=#a000a0
		hi Comment		gui=NONE guifg=grey62
		hi MicroController	gui=bold guifg=#d00000
		hi AnsiFuncPtr	gui=NONE guifg=#ff0000
		" hi PreProc        	gui=NONE guifg=#6a5acd
		hi PreCondit      	gui=NONE guifg=#6a5acd
		" hi Macro          	gui=NONE guifg=#0000ff

		delcommand HiLink
	endif
	hi Normal		gui=NONE guifg=#000000 guibg=Ivory1

	let b:current_syntax = "c"

	" vim: ts=8
