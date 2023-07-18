MODULE quick_vector_write_m

        IMPLICIT NONE



! *******************************************************************

CONTAINS

SUBROUTINE write_M_by_rows(M, ni, nj, label, unit, file_name)

!       Writes a real matrix, row by row, to formatted file or stdio.
!       Checks optional argument file_name to see if file of thisname is already open,
!       and if so what is unit number.
!       If not, checks to see if unit_number is open.  If open itwrites to this unit.
!       If no file is already open, then:
!       If file_name is present it opens a file with this name
!       If unit_number is present it opens with this unit number,else is uses default_unit
!       If neither unit_number, nor file_name is present it writes to stdio_unit

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: ni, nj
        REAL, INTENT(IN) :: M(ni, nj)
        CHARACTER (len = *) :: label
        INTEGER, INTENT(IN), OPTIONAL :: unit
        CHARACTER (LEN = *), OPTIONAL :: file_name

        INTEGER, PARAMETER :: default_unit = 57, stdio_unit = 6
        LOGICAL :: is_open
        INTEGER :: unit_temp, u_number
        CHARACTER (LEN = 31) :: f_name
        CHARACTER (LEN = 11 ) :: frm

        INTEGER :: i, j

!       Check out status of file

        is_open = .false.
        IF (PRESENT(file_name)) then

                ! Check if file open with this name, get u_number
                INQUIRE (file = TRIM(file_name), opened = is_open, number = unit_temp, &
                                                & formatted = frm)

                IF (is_open) THEN

                        IF (frm == 'UNFORMATTED' ) THEN ! Make sure it's a formatted file.  If not, bail.
                                WRITE (*,*) 'write_M_by_rows: can''t do formatted write to unformatted file -> ', &
                                        & TRIM(file_name)
                                STOP
                        END IF

                        u_number = unit_temp    ! write to open unit number even if unit argument is present

                END IF

        ELSE IF (PRESENT(unit)) THEN

                ! Check if file is open on this unit number
                INQUIRE (unit = unit, opened = is_open, formatted = frm)

                IF (is_open) THEN

                        IF (frm == 'UNFORMATTED' ) THEN ! Make sure it's a formatted file.  If not, bail.
                                WRITE (*,*) 'write_M_by_rows: can''t do formatted write to unformatted file -> ', &
                                        & TRIM(file_name)
                                STOP
                        END IF

                        u_number = unit ! use the unit_number argument
                END IF
        END IF

! If no file open, then open one, or use stdio

        IF (.not. is_open) then

                IF (PRESENT(unit) ) THEN        ! Use unit_number argument

                        u_number = unit

                        IF (PRESENT(file_name)) then    ! use file_name argument if there is one
                                f_name = TRIM(file_name)
                                OPEN (UNIT=u_number, FILE=f_name, FORM='FORMATTED')
                        ELSE
                                OPEN (UNIT=u_number, FORM='FORMATTED') ! use fortran default for file name
                        END IF

                ELSE    ! No unit_number argument, so use default or stdio

                        u_number = default_unit
                        IF (PRESENT(file_name)) then    ! use file_name argument and default unit_number
                                OPEN (UNIT=u_number, FILE=TRIM(file_name), FORM='FORMATTED')
                        ELSE
                                u_number = stdio_unit   ! write to stdio
                        END IF

                END IF

        END IF

!       write data

        WRITE (u_number,*) label

        DO i = 1, ni

                WRITE (u_number,*)
                WRITE (u_number, '("i = ", I6)' ) i
                WRITE (u_number, '(8e15.6)') (M(i,j), j = 1, nj)

        END DO

        IF (.not. is_open) CLOSE(u_number)      ! don't close file if opened outside this routine

END SUBROUTINE write_M_by_rows


SUBROUTINE write_M_by_columns(M, ni, nj, label, unit, file_name)

!       Writes a real matrix, column by column, to formatted file or stdio.
!       Checks optional argument file_name to see if file of this name is already open,
!       and if so what is unit number.
!       If not, checks to see if unit_number is open.  If open it writes to this unit.
!       If no file is already open, then:
!       If file_name is present it opens a file with this name
!       If unit_number is present it opens with this unit number, else is uses default_unit
!       If neither unit_number, nor file_name is present it writes to stdio_unit

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: ni, nj
        REAL, INTENT(IN) :: M(ni, nj)
        CHARACTER (len = *) :: label
        INTEGER, INTENT(IN), OPTIONAL :: unit
        CHARACTER (LEN = *), OPTIONAL :: file_name

        INTEGER, PARAMETER :: default_unit = 57, stdio_unit = 6
        LOGICAL :: is_open
        INTEGER :: unit_temp, u_number
        CHARACTER (LEN = 31) :: f_name
        CHARACTER (LEN = 11 ) :: frm

        INTEGER :: i, j

!       Check out status of file

        is_open = .false.
        IF (PRESENT(file_name)) then

                ! Check if file open with this name, get u_number
                INQUIRE (file = TRIM(file_name), opened = is_open, number = unit_temp, &
                                                & formatted = frm)

                IF (is_open) THEN

                        IF (frm == 'UNFORMATTED' ) THEN ! Make sure it's a formatted file.  If not, bail.
                                WRITE (*,*) &
                                & 'write_M_by_columns: can''t do formatted write to unformatted  file -> ',&
                                & TRIM(file_name)
                                STOP
                        END IF

                        u_number = unit_temp    ! write to open unit number even if unit argument is present

                END IF

        ELSE IF (PRESENT(unit)) THEN

                ! Check if file is open on this unit number
                INQUIRE (unit = unit, opened = is_open, formatted = frm)

                IF (is_open) THEN

                        IF (frm == 'UNFORMATTED' ) THEN ! Make sure it's a formatted file.  If not, bail.
                                WRITE (*,*) 'write_M_by_rows: can''t do formatted write to unformatted file -> ', &
                                        & TRIM(file_name)
                                STOP
                        END IF

                        u_number = unit ! use the unit_number argument
                END IF
        END IF

! If no file open, then open one, or use stdio

        IF (.not. is_open) then

                IF (PRESENT(unit) ) THEN        ! Use unit_number argument

                        u_number = unit

                        IF (PRESENT(file_name)) then    ! use file_name argument if there is one
                                f_name = TRIM(file_name)
                                OPEN (UNIT=u_number, FILE=f_name, FORM='FORMATTED')
                        ELSE
                                OPEN (UNIT=u_number, FORM='FORMATTED') ! use fortran default for file name
                        END IF

                ELSE    ! No unit_number argument, so use default or stdio

                        u_number = default_unit
                        IF (PRESENT(file_name)) then    ! use file_name argument and default unit_number
                                OPEN (UNIT=u_number, FILE=TRIM(file_name), FORM='FORMATTED')
                        ELSE
                                u_number = stdio_unit   ! write to stdio
                        END IF

                END IF

        END IF
!       write data

        WRITE (u_number,*) label

        DO j = 1, nj

                WRITE (u_number,*)
                WRITE (u_number, '("j = ", I6)' ) j
                WRITE (u_number, '(8e15.6)') (M(i,j), i = 1, ni)

        END DO

        IF (.not. is_open) CLOSE(u_number)      ! don't close file if opened outside this routine

END SUBROUTINE write_M_by_columns


SUBROUTINE write_V(V, ni, label, unit, file_name, imin, imax)

!       Writes a real vector to formatted file or stdio.
!       Checks optional argument file_name to see if file of this name is already open,
!       and if so what is unit number.
!       If not, checks to see if unit_number is open.  If open it writes to this unit.
!       If no file is already open, then:
!       If file_name is present it opens a file with this name
!       If unit_number is present it opens with this unit number, else is uses default_unit
!       If neither unit_number, nor file_name is present it writes to stdio_unit

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: ni
        REAL, INTENT(IN) :: V(ni)
        CHARACTER (len = *) :: label
        INTEGER, INTENT(IN), OPTIONAL :: unit
        CHARACTER (LEN = *), OPTIONAL :: file_name
        INTEGER, INTENT(IN), OPTIONAL :: imin, imax

        INTEGER, PARAMETER :: default_unit = 57, stdio_unit = 6
        LOGICAL :: is_open
        INTEGER :: unit_temp, u_number, i_start, i_stop
        CHARACTER (LEN = 31) :: f_name
        CHARACTER (LEN = 11 ) :: frm

        INTEGER :: i

!       Check out status of file

        is_open = .false.
        IF (PRESENT(file_name)) then

                ! Check if file open with this name, get u_number
                INQUIRE (file = TRIM(file_name), opened = is_open, number = unit_temp, &
                                                & formatted = frm)

                IF (is_open) THEN

                        IF (frm == 'UNFORMATTED' ) THEN ! Make sure it's a formatted file.  If not, bail.
                                WRITE (*,*) 'write_V: can''t do formatted write to unformatted file -> ', &
                                        & TRIM(file_name)
                                STOP
                        END IF

                        u_number = unit_temp    ! write to open unit number even if unit argument is present

                END IF

        ELSE IF (PRESENT(unit)) THEN

                ! Check if file is open on this unit number
                INQUIRE (unit = unit, opened = is_open, formatted = frm)

                IF (is_open) THEN

                        IF (frm == 'UNFORMATTED' ) THEN ! Make sure it's a formatted file.  If not, bail.
                                WRITE (*,*) 'write_M_by_rows: can''t do formatted write to unformatted file -> ', &
                                        & TRIM(file_name)
                                STOP
                        END IF

                        u_number = unit ! use the unit_number argument
                END IF
        END IF

! If no file open, then open one, or use stdio

        IF (.not. is_open) then

                IF (PRESENT(unit) ) THEN        ! Use unit_number argument

                        u_number = unit

                        IF (PRESENT(file_name)) then    ! use file_name argument if there is one
                                f_name = TRIM(file_name)
                                OPEN (UNIT=u_number, FILE=f_name, FORM='FORMATTED')
                        ELSE
                                OPEN (UNIT=u_number, FORM='FORMATTED') ! use fortran default for file name
                        END IF

                ELSE    ! No unit_number argument, so use default or stdio

                        u_number = default_unit
                        IF (PRESENT(file_name)) then    ! use file_name argument and default unit_number
                                OPEN (UNIT=u_number, FILE=TRIM(file_name), FORM='FORMATTED')
                        ELSE
                                u_number = stdio_unit   ! write to stdio
                        END IF

                END IF

        END IF

!       write data

        i_start = 1
        i_stop = ni
        WRITE (u_number,*) label

        IF (PRESENT(imin)) then
                i_start = imin
                WRITE (u_number, *) "i Min = ", imin
        END IF

        IF (PRESENT(imax)) then
                i_stop = imax
                WRITE (u_number, *) "i Max = ", imax
        END IF

!       WRITE (u_number,*)
        WRITE (u_number, '(8e15.6)') (V(i), i = i_start, i_stop)

        IF (.not. is_open) CLOSE(u_number)      ! don't close file if opened outside this routine

END SUBROUTINE write_V


! *******************************************************************


END MODULE quick_vector_write_m


