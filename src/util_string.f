module util_string
  use var_type,only:default_string_length
  implicit none
  private
  public::splitAndGetNext,integer_to_string,uppercase,lowercase,is_whitespace,str_trim,str_comp,str_search,remove_word,str_compress,typo_match,glob_match,is_blank_line
  CHARACTER,PARAMETER::whitespace*3=" "//CHAR(9)//CHAR(11),commentChar(3)=(/"#","!","%"/),backslash='\\',star='*',question='?'
contains

  subroutine splitAndGetNext(string,ia,ib)
    CHARACTER(LEN=*), INTENT(IN)             :: string
    INTEGER, INTENT(OUT)                     :: ia, ib

    ia = VERIFY(string,whitespace)
    if (ia.eq.0) then
       ib =0
    else
       ib = SCAN(string(ia+1:),whitespace) + ia - 1
    end if
    if (ib.lt.ia) ib=LEN(string)

  end subroutine splitAndGetNext

! *****************************************************************************
!> \brief   Return a string representing the integer number
! *****************************************************************************
  FUNCTION integer_to_string(number) RESULT(string)
    INTEGER, INTENT(IN)                      :: number
    CHARACTER(LEN=default_string_length)     :: string

    WRITE (UNIT=string,FMT='(I0)') number
  END FUNCTION  integer_to_string

! *****************************************************************************
!> \brief   Convert all lower case characters in a string to upper case.
! *****************************************************************************
  function uppercase(string)
    CHARACTER(LEN=*),INTENT(IN)::string
    character(len=len_trim(string))::uppercase
    INTEGER                                  :: i, iascii
    DO i=1,len(uppercase)
       iascii = IACHAR(string(i:i))
       IF ((iascii >= 97).AND.(iascii <= 122)) THEN
          uppercase(i:i) = ACHAR(iascii - 32)
       ELSE
          uppercase(i:i)=string(i:i)
       END IF
    END DO

  END function uppercase

! *****************************************************************************
  function lowercase(string)
    CHARACTER(LEN=*),INTENT(IN)::string
    character(len=len_trim(string))::lowercase
    INTEGER                                  :: i, iascii
    DO i=1,len(lowercase)
       iascii = ICHAR(string(i:i))
       IF ((iascii >= 65).AND.(iascii <= 90)) THEN
          lowercase(i:i) = CHAR(iascii + 32)
       ELSE
          lowercase(i:i)=string(i:i)
       END IF
    END DO

  END function lowercase

! *****************************************************************************
!> \brief returns .true. if the character passed is a whitespace char.
!> \par History
!>      02.2008 created, AK
! *****************************************************************************
  LOGICAL FUNCTION is_whitespace(testchar)
    CHARACTER(LEN=1), INTENT(IN)             :: testchar

    is_whitespace=.FALSE.
    IF (INDEX(whitespace,testchar).gt.0) is_whitespace = .TRUE.
  END FUNCTION is_whitespace

! *****************************************************************************
!> \brief returns .true. if the line is a comment line or an empty line
! *****************************************************************************
  LOGICAL FUNCTION is_blank_line(line,skipComment)
    CHARACTER(LEN=*), INTENT(IN)             :: line
    logical,intent(in)                       :: skipComment
    INTEGER                                  :: ia,ib

    call str_trim(line,ia,ib)
    if (ia.eq.0) then
       is_blank_line=.true.
    else if (skipComment.and.any(commentChar==line(ia:ia))) then
       is_blank_line=.true.
    else
       is_blank_line=.false.
    end if

  END FUNCTION is_blank_line

! *****************************************************************************
! Return ia, ib, such that string(ia:ib) has no leading or
! trailing spaces
! *****************************************************************************
  SUBROUTINE str_trim(string,ia,ib)
    CHARACTER(LEN=*), INTENT(IN)             :: string
    INTEGER, INTENT(OUT)                     :: ia, ib

    ia = VERIFY(string,whitespace)
    ib = VERIFY(string,whitespace,.true.)

  END SUBROUTINE str_trim

! *****************************************************************************
! Compare two strings ignoring the leading or trailing
! spaces. .true. if str1==str2, .false. otherwise
! *****************************************************************************
  FUNCTION str_comp(str1,str2) RESULT (equal)
    CHARACTER(LEN=*), INTENT(IN)             :: str1, str2
    LOGICAL                                  :: equal

    INTEGER                                  :: i1, i2, j1, j2

    i1 = 0
    i2 = 0
    j1 = 0
    j2 = 0
    CALL str_trim(str1,i1,i2)
    CALL str_trim(str2,j1,j2)
    equal = (str1(i1:i2)==str2(j1:j2))
  END FUNCTION str_comp

! *****************************************************************************
! Return the index, pos, of str2 in str1(1:n), which
! contains an array of strings, using str_comp(str1,str2)
! *****************************************************************************
  FUNCTION str_search(str1,n,str2) RESULT (pos)
    CHARACTER(LEN=*), INTENT(IN)             :: str1(:)
    INTEGER, INTENT(IN)                      :: n
    CHARACTER(LEN=*), INTENT(IN)             :: str2
    INTEGER                                  :: pos

    INTEGER                                  :: i

    pos = 0
    DO i = 1, n
       IF (str_comp(str1(i),str2)) THEN
          pos = i
          EXIT
       END IF
    END DO
  END FUNCTION str_search

! *****************************************************************************
!> \brief   remove a word from a string (words are separated by white spaces)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE remove_word(string)
    CHARACTER(LEN=*), INTENT(INOUT)          :: string

    INTEGER                                  :: i

    i = 1
    ! possibly clean white spaces
    DO WHILE(string(i:i)==" ")
       i = i + 1
    END DO
    ! now remove the word
    DO WHILE(string(i:i)/=" ")
       i = i + 1
    END DO
    string = string(i:)

  END SUBROUTINE remove_word

! *****************************************************************************
!> \brief   Eliminate multiple space characters in a string.
!>          If full is .TRUE., then all spaces are eliminated.  
!> \author  MK
!> \date    23.06.1998
!> \version 1.0
! *****************************************************************************
  SUBROUTINE str_compress(string,full)
    CHARACTER(LEN=*), INTENT(INOUT)          :: string
    LOGICAL, INTENT(IN), OPTIONAL            :: full

    INTEGER                                  :: i, z
    LOGICAL                                  :: remove_all

    IF (PRESENT(full)) THEN
       remove_all = full
    ELSE
       remove_all = .FALSE.
    END IF

    z = 1

    DO i=1,LEN_TRIM(string)
       IF ((z == 1).OR.remove_all) THEN
          IF (string(i:i) /= " ") THEN
             string(z:z) = string(i:i)
             z = z + 1
          END IF
       ELSE
          IF ((string(i:i) /= " ").OR.(string(z-1:z-1) /= " ")) THEN
             string(z:z) = string(i:i)
             z = z + 1
          END IF
       END IF
    END DO

    string(z:) = ""

  END SUBROUTINE str_compress

! *****************************************************************************
!> \brief Convert a sequence of integer numbers (ASCII code) to a string.
!>         Blanks are inserted for invalid ASCII code numbers.  
!> \author  MK
!> \date    19.10.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE ascii_to_string(nascii,string)
    INTEGER, DIMENSION(:), INTENT(IN)        :: nascii
    CHARACTER(LEN=*), INTENT(OUT)            :: string

    INTEGER                                  :: i

    string = ""

    DO i=1,MIN(LEN(string),SIZE(nascii))
       IF ((nascii(i) >= 0).AND.(nascii(i) <= 127)) THEN
          string(i:i) = CHAR(nascii(i))
       ELSE
          string(i:i) = " "
       END IF
    END DO

  END SUBROUTINE ascii_to_string

! *****************************************************************************
!> \brief   Convert a string to sequence of integer numbers.
!> \author  MK
!> \date    19.10.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE string_to_ascii(string,nascii)
    CHARACTER(LEN=*), INTENT(IN)             :: string
    INTEGER, DIMENSION(:), INTENT(OUT)       :: nascii

    INTEGER                                  :: i

    nascii(:) = 0

    DO i=1,MIN(LEN(string),SIZE(nascii))
       nascii(i) = ICHAR(string(i:i))
    END DO

  END SUBROUTINE string_to_ascii

! *****************************************************************************
!> \brief returns a non-zero positive value if typo_string equals string apart from a few typos.
!>     It is case sensitive, apart from typos.
!> \note
!>     could maybe be made a bit smarter
!> \par History
!>      02.2006 created [Joost VandeVondele]
! *****************************************************************************
  FUNCTION typo_match(string,typo_string) RESULT(match)
    CHARACTER(LEN=*), INTENT(IN)             :: string, typo_string
    INTEGER                                  :: match

    CHARACTER(LEN=1)                         :: kind
    CHARACTER(LEN=LEN(string))               :: tmp2
    CHARACTER(LEN=LEN(typo_string))          :: tmp
    INTEGER                                  :: i, j

    match=0
    IF (LEN_TRIM(typo_string).LE.4) THEN
       kind=question 
    ELSE
       kind=star
    ENDIF
    DO i=1,LEN_TRIM(typo_string)
       DO j=i,LEN_TRIM(typo_string)
          tmp=typo_string
          tmp(i:i)=kind
          tmp(j:j)=kind
          IF (i==j .AND. LEN_TRIM(tmp)>2 ) tmp(i:i)=star
          IF (glob_match(string=string,pattern=tmp)) match=match+1
       ENDDO
    ENDDO
    IF (LEN_TRIM(string).LE.4) THEN
       kind=question 
    ELSE
       kind=star
    ENDIF
    DO i=1,LEN_TRIM(string)
       DO j=i,LEN_TRIM(string)
          tmp2=string
          tmp2(i:i)=kind
          tmp2(j:j)=kind
          IF (i==j .AND. LEN_TRIM(tmp2)>2 ) tmp2(i:i)=star
          IF (glob_match(string=typo_string,pattern=tmp2)) match=match+1
       ENDDO
    ENDDO

  END FUNCTION typo_match

! *****************************************************************************
!
!  Imported from Arjen Markus' flibs project
!  http://flibs.sourceforge.net/
!
! globmatch.f90 --
!     Match strings according to (simplified) glob patterns
!
!     The pattern matching is limited to literals, * and ?
!     (character classes are not supported). A backslash escapes
!     any character.
!
! *****************************************************************************
!
! glob_match --
!     Tries to match the given string with the pattern
! Arguments:
!     string     String to be examined
!     pattern    Glob pattern to be used for the matching
! Result:
!     .true. if the entire string matches the pattern, .false.
!     otherwise
! Note:
!     Trailing blanks are ignored
!
! *****************************************************************************
RECURSIVE FUNCTION glob_match( string, pattern ) RESULT(match)
    CHARACTER(len=*), INTENT(in)             :: string, pattern
    LOGICAL                                  :: match

    CHARACTER(len=LEN(pattern))              :: literal
    INTEGER                                  :: k, ll, method, p, ptrim, &
                                                start, strim

    match  = .FALSE.
    method = 0
    ptrim  = LEN_TRIM( pattern )
    strim  = LEN_TRIM( string )
    p      = 1
    ll     = 0
    start  = 1

    !
    ! Split off a piece of the pattern
    !
    DO WHILE ( p <= ptrim )
        SELECT CASE ( pattern(p:p) )
            CASE( star )
                IF ( ll .NE. 0 ) EXIT
                method = 1
            CASE( question )
                IF ( ll .NE. 0 ) EXIT
                method = 2
                start  = start + 1
            CASE( backslash )
                p  = p + 1
                ll = ll + 1
                literal(ll:ll) = pattern(p:p)
            CASE default
                ll = ll + 1
                literal(ll:ll) = pattern(p:p)
        END SELECT

        p = p + 1
    ENDDO

    !
    ! Now look for the literal string (if any!)
    !
    IF ( method == 0 ) THEN
        !
        ! We are at the end of the pattern, and of the string?
        !
        IF ( strim == 0 .AND. ptrim == 0 ) THEN
            match = .TRUE.
        ELSE
            !
            ! The string matches a literal part?
            !
            IF ( ll > 0 ) THEN
                IF ( string(start:MIN(strim,start+ll-1)) == literal(1:ll) ) THEN
                    start = start + ll
                    match = glob_match( string(start:), pattern(p:) )
                ENDIF
            ENDIF
        ENDIF
    ENDIF

    IF ( method == 1 ) THEN
        !
        ! Scan the whole of the remaining string ...
        !
        IF ( ll == 0 ) THEN
            match = .TRUE.
        ELSE
            DO WHILE ( start <= strim )
                k     = INDEX( string(start:), literal(1:ll) )
                IF ( k > 0 ) THEN
                    start = start + k + ll - 1
                    match = glob_match( string(start:), pattern(p:) )
                    IF ( match ) THEN
                        EXIT
                    ENDIF
                ENDIF

                start = start + 1
            ENDDO
        ENDIF
    ENDIF

    IF ( method == 2 .AND. ll > 0 ) THEN
        !
        ! Scan the whole of the remaining string ...
        !
        IF ( string(start:MIN(strim,start+ll-1)) == literal(1:ll) ) THEN
            match = glob_match( string(start+ll:), pattern(p:) )
        ENDIF
    ENDIF
    RETURN
END FUNCTION glob_match

end module util_string
