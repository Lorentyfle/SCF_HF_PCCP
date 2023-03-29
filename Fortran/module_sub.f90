module input_HF_PCCP
    implicit none

contains
    subroutine read_data()
        
        use constant, only : sent_coef_Slater,size_matrix_Slater,p_SCF,size_cp_SCF,cp1_SCF,thr_SCF
        implicit none
        ! ***************
        ! * Declaration *
        ! ***************
        character(len=17)                        :: input_fort       = './input/input.txt'
        double precision, allocatable:: sent_coef(:)
        integer                      :: size_matrix
        integer                      :: p
        integer,parameter            :: size_cp          = size_cp_SCF
        integer                      :: cp1(size_cp)
        double precision             :: thr
        character(len=256)      :: sent_prov        = ''
        character(len=9)        :: sent_start       = '--Begin--'
        character(len=1)        :: sent_matrix_size = 'M'
        character(len=7)        :: sent_end         = '--End--'
        character(len=4)        :: sent_cp1         = 'cp1='
        character(len=5)        :: sent_cp1_        = 'cp1 ='
        character(len=2)        :: sent_p1          = 'p='
        character(len=3)        :: sent_p1_         = 'p ='
        character(len=8)        :: sent_thr         = 'thr_SCF='
        character(len=9)        :: sent_thr_        = 'thr_SCF ='
        integer                 :: line_p,line_cp1,line_thr,line_matrix, line_data,line_data_end,have_matrix,line_prov
        integer                 :: max_file         = 100
        integer                 :: i
        integer                 :: nlines
        !
        ! **************
        ! Find the datas
        ! **************
        open(1,file=input_fort,position="rewind")
        !
        nlines = 0
        do
            read(1,*,end=10)
            nlines = nlines + 1
        end do
        10  rewind(1)
        !
        line_p = 0
        do while(sent_prov(1:max_file) /= sent_p1 .and. sent_prov(1:max_file) /= sent_p1_)
            line_p = line_p + 1
            read(1,'(A)') sent_prov
            sent_prov = TRIM(sent_prov)
        end do
        sent_prov = ''
        ! *********************
        rewind(1)
        ! *********************
        line_cp1 = 0
        do while(sent_prov(1:max_file) /= sent_cp1 .and. sent_prov(1:max_file) /= sent_cp1_)
            line_cp1 = line_cp1 + 1
            read(1,'(A)') sent_prov
            sent_prov = TRIM(sent_prov)
        end do
        sent_prov = ''
        ! *********************
        rewind(1)
        ! *********************
        line_thr = 0
        do while(sent_prov(1:max_file) /= sent_thr .and. sent_prov(1:max_file) /= sent_thr_)
            line_thr = line_thr + 1
            read(1,'(A)') sent_prov
            sent_prov = TRIM(sent_prov)
        end do
        sent_prov = ''
        ! *********************
        rewind(1)
        ! *********************


        have_matrix = 0
        line_matrix = 0
        do while(sent_prov(1:max_file) /= sent_matrix_size .and. line_matrix < nlines)
            line_matrix = line_matrix + 1
            read(1,'(A)') sent_prov
            sent_prov = TRIM(sent_prov)
            if ( sent_prov(1:max_file) == sent_matrix_size) then
                have_matrix = 1
            end if
        end do
        sent_prov = ''
        ! *********************
        rewind(1)
        ! *********************
        line_data = 0
        do while(sent_prov(1:max_file) /= sent_start)
            line_data = line_data + 1
            read(1,'(A)') sent_prov
            sent_prov = TRIM(sent_prov)
        end do
        sent_prov = ''
        ! *********************
        rewind(1)
        ! *********************
        if ( have_matrix == 1 ) then
            line_prov = 0
            do i = 1, line_matrix+1
                line_prov = line_prov + 1
                if (line_prov < line_matrix+1) then
                    read(1,'(A)')
                else
                    read(1,'(A)') sent_prov
                    sent_prov = TRIM(sent_prov)
                end if
            end do
            sent_prov = ''
            line_prov = 0
            line_data_end = line_data + line_matrix
            backspace(1)
            read(1,"(I5)") size_matrix
            ! *********************
            rewind(1)
            ! *********************
        else
            line_data_end = 0
            do while(sent_prov(1:max_file) /= sent_end)
                line_data_end = line_data_end + 1
                read(1,'(A)') sent_prov
                sent_prov = TRIM(sent_prov)
            end do
            sent_prov = ''
            ! *********************
            rewind(1)
            ! *********************
            size_matrix = line_data_end - (line_data+1)
        end if
        

        allocate(sent_coef(size_matrix))

        line_prov = 0
        do i = 1, line_data+size_matrix
            if ( i < line_data+1 ) then
                read(1,'(A)')
            else
                read(1,'(1D10.3)') sent_coef(i-(line_data))
                write(*,'(1D10.3)') sent_coef(i-(line_data))
                write(*,'(1F10.3)') sent_coef(i-(line_data))
            end if
        end do
        ! *********************
        rewind(1)
        ! *********************


        line_prov = 0
        do i = 1, line_p+1
            line_prov = line_prov + 1
            if (line_prov < line_p+1) then
                read(1,'(A)')
            else
                read(1,'(A)') sent_prov
                sent_prov = TRIM(sent_prov)
            end if
        end do
        sent_prov = ''
        line_prov = 0
        line_data_end = line_data + line_matrix
        backspace(1)
        read(1,"(I5)") p

        ! *********************
        rewind(1)
        ! *********************

        line_prov = 0
        do i = 1, line_cp1
            line_prov = line_prov + 1
            read(1,'(A)')
        end do
        line_prov = 0
        ! We take the datas back
        do i = 1, size_cp
            read(1,"(I5)") cp1(i)
        end do
        
        ! *********************
        rewind(1)
        ! *********************

        line_prov = 0
        do i = 1, line_thr+1
            line_prov = line_prov + 1
            if (line_prov < line_thr+1) then
                read(1,'(A)')
            else
                read(1,'(A)') sent_prov
                sent_prov = TRIM(sent_prov)
            end if
        end do
        sent_prov = ''
        line_prov = 0
        ! We take the datas back
        backspace(1)
        read(1,"(1E14.10)") thr

        ! Test to see if all the values extracted are correct
        write(*,*) "Slater coefficients:"
        write(*,*) "We have ", size_matrix, " coefficients."
        do i = 1, size(sent_coef)
            write(*,'(1D10.3)') sent_coef(i)
        end do
        write(*,*) "SCF parameters:"
        write(*,*) "SCF threshold = ",thr
        do i = 1, size(cp1)
            if (i <= p) then
                write(*,*) cp1(i), " <=  (p =", p,")"
            else
                write(*,*) cp1(i), " >   (p =", p,")"
            end if
        end do
        !write(*,*) p,cp1,thr,size_matrix,sent_coef
        !size_cp_SCF         = size_cp ! Here the lenght for cp1 is the same.
        

        size_matrix_Slater  = size_matrix
        sent_coef_Slater    = dble(sent_coef)
        

        p_SCF               = p
        cp1_SCF             = cp1
        thr_SCF             = thr


        close(1)
        deallocate(sent_coef)
    end subroutine read_data

    subroutine write_intialisation()
        implicit none
        character(len=19)       :: ofile
        character(len=25)       :: restart_txt
        integer                 :: nlines, end_file,i
        character(len=24)       :: End_text     = "%--End_HF_Calculation--%"
        character(len=25)       :: End_text_    = " %--End_HF_Calculation--%"

        ofile = './output/output.txt'
        open(2,file=ofile,position="rewind")
        !
        nlines = 0
        do
            read(2,*,end=12)
            nlines = nlines + 1
        end do
        12  close(2)
        !
        if ( nlines - 1 > 0 ) then
            end_file = nlines - 1
            open(2,file=ofile)
            do i = 1, end_file
                read(2,*)
            end do
            read(2,"(A)") restart_txt
            close(2)
        else
            end_file = 0
            restart_txt = ""
        end if
        
        restart_txt = trim(restart_txt)
        if ( (restart_txt .eq. End_text) .or. (restart_txt .eq. End_text_) ) then
                open(2,file=ofile,action='write',status='replace')
            else
                open(2,file=ofile,action='write',position='append')
            end if
            write(2,*) "=========================================="
            close(2)
            open(2,file=ofile,action='write',position='append')
            write(2,*) "= Output of the Hartee Fock calculation  ="
            write(2,*) "=========================================="
            close(2)
    end subroutine write_intialisation

    subroutine write_coefficient(Coeff_write, loop)
        use constant, only : size_matrix_Slater
        implicit none
        
        double precision,dimension(size_matrix_Slater), intent(in) :: Coeff_write
        integer,intent(in)      :: loop
        integer                 :: i
        character(len=19)       :: ofile
        !
        ! **************
        ! Size matrix?
        ! **************
        !
        ofile = './output/output.txt'
        !
            open(2,file=ofile,action='write',position='append')
            write(2,*) "=========================================="
            write(2,*) "Loop = ", loop
            do i = 1, size_matrix_Slater
                write(2,*) Coeff_write(i)
            end do
            close(2)
    end subroutine write_coefficient

    subroutine write_energy(Energy_write)
        
        implicit none

        double precision, intent(in)    :: Energy_write
        character(len=500)              :: ofile
        !
        ! **************
        ! Size matrix?
        ! **************
        !
        ofile = './output/output.txt'
        open(2,file=ofile,action='write',position='append')
        write(2,*) "============Converged_energy=============="
        write(2,*) Energy_write
        close(2)

        ! Finish the file
        open(2,file=ofile,action='write',position='append')
        write(2,*) "====================================================="
        write(2,*) "HF program coded by Alex Delhumeau and Timothee Jamin"
        write(2,*) "====================================================="
        write(2,*) "%--End_HF_Calculation--%"
        close(2)
    end subroutine write_energy

end module input_HF_PCCP