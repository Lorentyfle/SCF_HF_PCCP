module input_HF_PCCP
    implicit none

contains
    subroutine read_data()
        
        use constant, only : sent_coef_Slater,size_matrix_Slater,p_SCF,size_cp_SCF,cp1_SCF,thr_SCF,mesh
        implicit none
        ! ***************
        ! * Declaration *
        ! ***************
        character(len=17)            :: input_fort  = './input/input.txt'
        double precision, allocatable:: sent_coef(:)
        integer                      :: size_matrix
        integer                      :: p
        integer,parameter            :: size_cp     = size_cp_SCF
        integer                      :: cp1(size_cp)
        double precision             :: thr
        character(len=256)      :: sent_prov        = ''
        character(len=4)        :: Mesh__           = 'mesh'
        character(len=4)        :: Mesh_            = 'Mesh'
        character(len=9)        :: sent_start       = '--Begin--'
        character(len=1)        :: sent_matrix_size = 'M'
        character(len=7)        :: sent_end         = '--End--'
        character(len=4)        :: sent_cp1         = 'cp1='
        character(len=5)        :: sent_cp1_        = 'cp1 ='
        character(len=2)        :: sent_p1          = 'p='
        character(len=3)        :: sent_p1_         = 'p ='
        character(len=8)        :: sent_thr         = 'thr_SCF='
        character(len=9)        :: sent_thr_        = 'thr_SCF ='
        integer                 :: line_p,line_cp1,line_thr,line_data,line_data_end,line_matrix,line_prov,line_mesh
        integer                 :: have_matrix, have_mesh
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

        have_mesh = 0
        line_mesh = 0
        do while((sent_prov(1:max_file) /= Mesh__ .or. sent_prov(1:max_file) /= Mesh_) .and. line_mesh < nlines)
            line_mesh = line_mesh + 1
            read(1,'(A)') sent_prov
            sent_prov = TRIM(sent_prov)
            if ( sent_prov(1:max_file) .eq. Mesh__ .or. sent_prov(1:max_file) .eq. Mesh_ ) then
                have_mesh = 1
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
        read(1,*) thr


        ! Test to see if all the values extracted are correct
        write(*,*) "Slater coefficients:"
        write(*,*) "We have ", size_matrix, " coefficients."
        do i = 1, size(sent_coef)
            write(*,'(1D10.3)') sent_coef(i)
        end do
        write(*,*) "SCF parameters:"
        do i = 1, size(cp1)
            if (i <= p) then
                write(*,*) cp1(i), " <=  (p =", p,")"
            else
                write(*,*) cp1(i), " >   (p =", p,")"
            end if
        end do


        size_matrix_Slater  = size_matrix
        sent_coef_Slater    = dble(sent_coef)
        

        p_SCF               = p
        cp1_SCF             = cp1
        thr_SCF             = thr
        mesh                = have_mesh


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

    subroutine write_coefficient(Coeff_write, loop,Energy,Density,Eigenvalue,E_diff)
        use constant, only : size_matrix_Slater
        implicit none
        
        double precision,dimension(size_matrix_Slater), intent(in) :: Coeff_write
        double precision,intent(in)     :: Energy, Density, Eigenvalue, E_diff
        integer,intent(in)              :: loop
        integer                         :: i
        character(len=19)               :: ofile
        !
        ! **************
        ! Size matrix?
        ! **************
        !
        ofile = './output/output.txt'
    
        open(2,file=ofile,action='write',position='append')
        write(2,*) "=========================================="
        write(2,*) "Loop = ", loop
        write(2,*) "Energy = ", Energy,"Hartree"
        write(2,*) "Density verification = ", Density
        write(2,*) "Eigenvalue of the ground state = ", Eigenvalue
        write(2,*) "Energy difference = ", E_diff

        do i = 1, size_matrix_Slater
            write(2,*) Coeff_write(i)
        end do
        close(2)
    end subroutine write_coefficient

    subroutine write_energy(Energy_write,Eigenvalue_write)
        
        implicit none

        double precision, intent(in)    :: Energy_write,Eigenvalue_write
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
        write(2,*) "==========Converged_eigenvalue============"
        write(2,*) Eigenvalue_write
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

module Calculation_SCF
    implicit none
    
contains
    subroutine alpha_coef()
        use constant, only: size_matrix_Slater, sent_coef_Slater, alpha_Slater
        implicit none

        double precision, dimension(:), allocatable   :: alpha
        integer         :: M
        integer         :: i

        M = size_matrix_Slater
        allocate(alpha(M))
        do i = 1, M
            alpha(i) = sent_coef_Slater(i)
        end do

        !%%%%%%%%%%%%%%%%
        ! Save the data
        !%%%%%%%%%%%%%%%%

        alpha_Slater = alpha

        deallocate(alpha)
    end subroutine alpha_coef

    subroutine S_matrix(M_size, alpha, S_mat,S_minushalf_sub)
        
        use constant, only: S,S_minushalf
        use diagonalization,only : hshtri, tqli, bubble_sort

        implicit none

        integer, intent(in)                                        :: M_size
        double precision,dimension(:),allocatable,    intent(in)   :: alpha
        double precision,dimension(:,:),allocatable,  intent(out)  :: S_mat,S_minushalf_sub
        
        double precision,dimension(:,:),allocatable                :: Z,T,T_trans,Sbar_sub,Sbar_minushalf_sub,TMP
        double precision,dimension(:),allocatable                  :: e,diag,diag_sorted
        integer,dimension(:),allocatable                           :: diag_id
        integer                                                    :: i,j

        !%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Allocation of the datas
        !%%%%%%%%%%%%%%%%%%%%%%%%%

        allocate(S_mat(M_size,M_size))
        allocate(Z(M_size,M_size),T(M_size,M_size),T_trans(M_size,M_size))
        allocate(Sbar_sub(M_size,M_size),Sbar_minushalf_sub(M_size,M_size))
        allocate(S_minushalf_sub(M_size,M_size))

        allocate(TMP(M_size,M_size))
        allocate(e(M_size),diag(M_size),diag_id(M_size),diag_sorted(M_size))
        
        e = 0.0d0 ; diag = 0.0d0 ; diag_sorted = 0.0d0 ; diag_id = 0
        T = 0.0d0 ; Sbar_sub = 0.0d0 ; Sbar_minushalf_sub = 0.0d0 ; S_minushalf_sub = 0.0d0 ; T_trans = 0.0d0    


        !%%%%%%%%%%%%%
        ! Initilize S
        !%%%%%%%%%%%%%
        
        do i = 1, M_size
            do j = 1, M_size
                if (i .eq. j) then
                    S_mat(i,j) = 1.0d0
                else if (i .ne. j) then
                    S_mat(i,j) = dble((2.0d0*sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**3)
                end if
            end do
        end do

        
        !%%%%%%%%%%%%%%%
        ! Diagonalize S
        !%%%%%%%%%%%%%%%

        
        write(*,*) 'S matrix to diagonalise:'
        do i = 1, M_size
           write(*,'(40F12.8)') (S_mat(i,j), j = 1, M_size)
        end do
        write(*,*)
    
        Z = S_mat
    
        write(*,*) 'Tridiagonalization based on the Householder transformation'
        call hshtri(Z, M_size, diag, e, .false.)
        
        write(*,*) 'Diagonalization based on the QL method'
        call tqli(diag, e, M_size, M_size, Z)
        
        write(*,*) 'Sorting eigenvalues and eigenvectors of S'
        call bubble_sort(diag, M_size, diag_sorted, diag_id)
    
        write(*,*) 'Generating matrix of eigenvectors'
        diag = diag_sorted
    
        do i = 1, M_size
           j = diag_id(i)
           T(:,i) = Z(:,j)
        end do
    
        write(*,*) 'Eigenvectors of S a.k.a. T (column format):'
        do i = 1, size(S_mat,1)
           write(*,'(40f12.8)') (T(i,j), j=1, size(S_mat,2))
        end do
    
        write(*,*)
        write(*,*) 'Generating matrix of associated eigenvalues of S'
        do i = 1, M_size
           Sbar_sub(i,i) = diag(i)
        end do
    
        write(*,*) 'Eigenvalues of S a.k.a. Sbar (column format):'
        do i = 1, size(S_mat,1)
           write(*,'(40f12.8)') (Sbar_sub(i,j), j=1, size(S_mat,2))
        end do
        write(*,*)
    
        !%%%%%%%%%%%%%%%%%%%
        ! Compute Sbar^-1/2
        !%%%%%%%%%%%%%%%%%%%
        
        do i = 1, M_size
           Sbar_minushalf_sub(i,i) = dble(1.0d0/(sqrt(Sbar_sub(i,i))))
        end do
        
        write(*,*) 'Sbar^-1/2 (column format):'
        do i = 1, size(S_mat,1)
           write(*,'(40f12.8)') (Sbar_minushalf_sub(i,j), j=1, size(S_mat,2))
        end do
        write(*,*)
    
        !%%%%%%%%%%%%%%%%
        ! Compute S^-1/2
        !%%%%%%%%%%%%%%%%
        
        T_trans = transpose(T)
        TMP = matmul(Sbar_minushalf_sub, T_trans)
        S_minushalf_sub = matmul(T, TMP)
    
        write(*,*) 'S^-1/2 (column format):'
        do i = 1, size(S_mat,1)
           write(*,'(40f12.8)') (S_minushalf_sub(i,j), j=1, size(S_mat,2))
        end do
        write(*,*)

        !%%%%%%%%%%%%%%%%
        ! Save the data
        !%%%%%%%%%%%%%%%%

        S           = S_mat
        S_minushalf = S_minushalf_sub

        !%%%%%%%%%%%%%%%%%%%
        ! Free the memory
        !%%%%%%%%%%%%%%%%%%%
        !deallocate(S_mat)
        !deallocate(Sbar_sub,Sbar_minushalf_sub)
        !deallocate(S_minushalf_sub)

        !deallocate(Z,T,T_trans)
        !deallocate(TMP)
        !deallocate(e,diag,diag_id,diag_sorted)

    end subroutine S_matrix

    subroutine integrals_He(alpha,hpq,pqrs_sub)
        
        use constant, only: size_matrix_Slater,atomic_charge,h_pq,pqrs

        implicit none

        double precision, dimension(:),intent(in)                    :: alpha
        double precision, dimension(:,:),allocatable, intent(out)    :: hpq
        double precision, dimension(:,:,:,:),allocatable, intent(out):: pqrs_sub
        double precision                                             :: sum_alpha, product_alpha 
        double precision                                             :: pairsum, pairsum_2
        double precision                                             :: f1, f2, f3
        integer         :: M
        integer         :: i,j,k,l

        !%%%%%%%%%%%%%%
        ! Declarations
        !%%%%%%%%%%%%%%

        M = size_matrix_Slater
        allocate(hpq(M,M))
        allocate(pqrs_sub(M,M,M,M))

        !%%%%%%%%%%%%%%%%
        ! h_pq integrals
        !%%%%%%%%%%%%%%%%

        write(*,*) 'Computing numerical values of h_pq integrals:'
        do i = 1, M
            do j = 1, M
                sum_alpha = alpha(i) + alpha(j)
                product_alpha = alpha(i) * alpha(j)
                
                hpq(i,j) = 4.0d0 * dble(((sqrt(product_alpha)/sum_alpha)**3) * (product_alpha - atomic_charge*sum_alpha))
                
                write(*,*) "h_pq with p = ", i, "and q =", j
                write(*,'(40f12.8)') hpq(i,j)
            end do
        end do
        write(*,*)
        
        !%%%%%%%%%%%%%%%%%%%
        ! (pq|rs) integrals
        !%%%%%%%%%%%%%%%%%%%

        write(*,*) 'Computing numerical values of (pq|rs) integrals:'
        do i = 1, M
            do j = 1, M
                do k = 1, M
                    do l = 1, M
                        sum_alpha = alpha(i) + alpha(j) + alpha(k) + alpha(l)
                        pairsum = alpha(i) + alpha(j)
                        pairsum_2 = alpha(k) + alpha(l)
                        product_alpha = alpha(i) * alpha(j) * alpha(k) * alpha(l)

                        f1 = dble(1.0d0/((pairsum)**3 * (pairsum_2)**2))
                        f2 = dble(1.0d0/((pairsum)**3 * (sum_alpha)**2))
                        f3 = dble(1.0d0/((pairsum)**2 * (sum_alpha)**3))
                    
                        pqrs_sub(i,j,k,l) = 32.0d0 * dble((sqrt(product_alpha))**3 * (f1 - f2 - f3))
                        write(*,*) "(pq|rs) with p = ", i, "q =", j, "r = ", k, "s = ", l
                        write(*,'(40f12.8)') pqrs_sub(i,j,k,l)
                    end do
                end do
            end do
        end do
        
        !%%%%%%%%%%%%%%%%
        ! Save the data
        !%%%%%%%%%%%%%%%%

        h_pq = hpq
        pqrs = pqrs_sub

    end subroutine integrals_He

    subroutine Setup_SCF(coefficients_Fock_sub)
        use constant, only: size_matrix_Slater,cp1_SCF,p_SCF,coefficients_Fock
        implicit none

        double precision, dimension(:,:),allocatable,intent(out)   :: coefficients_Fock, coefficients_Fock_sub
        
        integer         :: M
        integer         :: i,j

        M = size_matrix_Slater
        
        allocate(coefficients_Fock_sub(M,M), coefficients_Fock(M,M))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Define initial guesses for Fock coefficients
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do i = 1, M
            do j = 1, M
                if (i .eq. p_SCF) then
                    coefficients_Fock_sub(i,j) = cp1_SCF(1)
                else
                    coefficients_Fock_sub(i,j) = cp1_SCF(2)
                end if
            end do
        end do

        coefficients_Fock = coefficients_Fock_sub
    end subroutine Setup_SCF

    subroutine SCF_HF( E_tot,Rpq_sub,Fock_matrix_sub, I_1_sub, J_11_sub )
        
        use constant, only: size_matrix_Slater,atomic_charge,h_pq,pqrs,thr_SCF,S,S_minushalf
        use constant, only: coefficients_Fock
        use constant, only: Rpq,Fock_matrix,I_1,J_11, E_total, E_verification, eigenvalues
        use input_HF_PCCP, only: write_coefficient
        use Other_function,only: write_matrix

        implicit none

        double precision, intent(out)   ::  E_tot
        double precision, dimension(:,:),allocatable, intent(out) ::  Rpq_sub
        double precision, dimension(:,:),allocatable, intent(out) ::  Fock_matrix_sub
        double precision,intent(out)    :: I_1_sub, J_11_sub

        double precision                                    :: E_tot_old, Verif_Dens, E_check
        double precision, dimension(:,:),allocatable        :: Fock_matrix_prime, Fock_matrix_prime_bar, Density
        double precision, dimension(:,:),allocatable        :: C_Fock
    
        integer                         :: M
        integer                         :: loop
        integer                         :: i,j,k,l

        !%%%%%%%%%%%%%
        ! Declaration
        !%%%%%%%%%%%%%

        M = size_matrix_Slater

        allocate(Fock_matrix_sub(M,M),Fock_matrix_prime(M,M),Fock_matrix_prime_bar(M,M),Rpq_sub(M,M),Density(M,M))
        allocate(C_Fock(M,M))

        E_tot = 0.0d0
        E_tot_old = (E_tot+1.0d0)*5.0d0
        loop = 0
        
        !%%%%%%%%%
        ! Program
        !%%%%%%%%%

        do while(abs(E_tot - E_tot_old) .gt. thr_SCF)
            write(*,*)
            if (loop .eq. 0) then
               write(*,*) "%%%%% Iterative SCF process has begun! %%%%%"
            !end if
            !if (loop .eq. 100) then
            !    write(*,*) "%%%%% Emergency exit, last loop! %%%%%"
            !    go to 10
            else
               write(*,*) "%%%%% Entering a new loop! %%%%%"
            end if
            loop = loop + 1
    
            ! Fock matrix elements
            write(*,*)
            write(*,*) "Loop", loop
            write(*,*) "Coefficients sent:"
            call write_matrix(coefficients_Fock)
    
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the Fock matrix elements and build the F matrix
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            do i = 1, M
                do j = 1, M
                    Fock_matrix_sub(i,j) = h_pq(i,j)
                    do k = 1, M
                        do l = 1, M
                            Fock_matrix_sub(i,j) = Fock_matrix_sub(i,j) + &
                            coefficients_Fock(k,1) * coefficients_Fock(l,1) * pqrs(i,j,k,l)
                        end do
                    end do
                end do
            end do
    
            write(*,*)
            write(*,*) "Fock matrix (column format):"
            call write_matrix(Fock_matrix_sub)
            write(*,*)
    
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Obtain the Fock matrix eigenvectors and eigenvalues
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ! Fock prime matrix
            Fock_matrix_prime = matmul(S_minushalf,matmul(Fock_matrix_sub,S_minushalf))
            
            write(*,*) "F' = S^(-1/2)FS^(-1/2)"
            write(*,*) 'F_prime matrix to diagonalise:'
            call write_matrix(Fock_matrix_prime)
            write(*,*)
            
            !%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Sorting and diagonalize
            !%%%%%%%%%%%%%%%%%%%%%%%%%

            call diag_SCF(Fock_matrix_prime,Fock_matrix_prime_bar,C_Fock)

            !%%%%%%%%%%%%%%%%%%%%%%
            ! Density matrix check
            !%%%%%%%%%%%%%%%%%%%%%%
    
            write(*,*) "Beginning density matrix verification"
            
            ! Formula used:
            ! D = 2R
            ! Rpq = Cp1 * Cq1
            
            do i = 1, M
                do j = 1, M
                    Rpq_sub(i,j) = C_Fock(i,1) * C_Fock(j,1)
                end do
            end do
    
            write(*,*) "R matrix (column format):"
            call write_matrix(Rpq_sub)
            write(*,*)
    
            Density = 2.0d0 * Rpq_sub
            
            write(*,*) "Density matrix (column format):"
            call write_matrix(Density)
            write(*,*)
            
            Verif_Dens = 0.0d0
            do i = 1, M
                do j = 1, M
                    Verif_Dens = Verif_Dens + Density(i,j) * S(i,j)
                end do
            end do
    
            write(*,*) "We have:", Verif_Dens, " =? ", atomic_charge 
    
            !%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Total energy calculation
            !%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            ! Formula for the energy:
            ! E = sum^{M}_{p=1}(sum^{M}_{q=1}(c_p1 * c_q1 *(h_pq + F_pq)))
    
            E_tot_old = E_tot
            E_tot = 0.0d0
            do i = 1, M
                do j = 1, M
    
                    ! Results from the reference papers Roetti and Clementi (see p. 2) and Snow and Bills (see p. 3):
                    ! E_final = -2.8616726                  | E_init = - 2.83308
    
                    E_tot = E_tot + coefficients_Fock(i,1) * coefficients_Fock(j,1) * (h_pq(i,j) + Fock_matrix_sub(i,j))
    
                    ! E_final = -2.8616715937556600 (logic) | E_init = -2.7912499999999998 (logic but don't follow the paper)
                     
                end do
            end do
    
            write(*,*)
            ! write(*,*) "c_11 =", coefficients_Fock(1,1), "c_21 =", coefficients_Fock(2,1)
            ! write(*,*) "F_11 =", Fock_matrix(1,1), "F_12 =", Fock_matrix(1,2), "F_22 =", Fock_matrix(2,2) 
            ! After verification, all this values are correct for the 1st loop. (c11, c21, F11, F12, F22, and epsilon)
            ! Only the total energy is different for the *first* loop.
    
            ! Compare energies
            write(*,*) "The previous energy is: ",  E_tot_old
            write(*,*) "The actual energy is:",     E_tot
            write(*,*) "The difference of energy is:", E_tot_old - E_tot," vs ", thr_SCF
            write(*,*)
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compare Fock coefficients
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%
            write(*,*) "Previous LCAO coefficients for phi_1 a.k.a. the first column of coefficients_Fock (column format):"
            call write_matrix(coefficients_Fock)
            write(*,*)
            write(*,*) "vs"
            write(*,*)
            write(*,*) "New Fock coefficients a.k.a. C_Fock (column format):"
            call write_matrix(C_Fock)
            write(*,*)
            
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Update the Fock coefficients for the next loop
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            write(*,*) "Editing of the starting coefficients"
            do i = 1, M
                coefficients_Fock(i,1) = C_Fock(i,1)
            end do
            write(*,*)
            

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Output actual coefficients to text file
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            call write_coefficient(coefficients_Fock(:,1),loop,E_tot,Verif_Dens,Fock_matrix_prime_bar(1,1),E_tot_old - E_tot)
    
            ! I1 = sum^{M}_{p=1}(sum^{M}_{q=1}(c_pi * c_qi * h_11)
    
            I_1_sub = 0.0d0
            do i = 1, M
               do j = 1, M
                  I_1_sub = I_1_sub + coefficients_Fock(i,1) * coefficients_Fock(j,1) * h_pq(i,j)
               end do
            end do
    
            write(*,*) "I_1 = ", I_1_sub
            write(*,*)
    
            ! J11 = sum=i_1^M sum=i_1^M sum=i_1^M sum=i_1^M (c_p1 * c_q1 * c_r1 * c_s1 * (pq|rs))
    
            J_11_sub = 0.0d0
            do i = 1, M
               do j = 1, M
                  do k = 1, M
                     do l = 1, M
                        J_11_sub = J_11_sub + (coefficients_Fock(i,1) * coefficients_Fock(j,1) &
                             * coefficients_Fock(k,1) * coefficients_Fock(l,1) * pqrs(i,j,k,l))
                     end do
                  end do
               end do
            end do
    
            write(*,*) "J_11 = ", J_11_sub
            write(*,*)
    
            ! Total energy with E = 2 * I_{1} + J_{11}
    
            E_check = 2.0d0 * I_1_sub + J_11_sub
            write(*,*) "Check: E_check = ", E_check    
        end do

        !%%%%%%%%%%%%%%%%
        ! Save the data
        !%%%%%%%%%%%%%%%%
        
        10 I_1      = I_1_sub
        J_11        = J_11_sub
        Rpq         = Rpq_sub
        Fock_matrix = Fock_matrix_sub
        E_total     = E_tot
        E_verification= E_check
        eigenvalues = Fock_matrix_prime_bar

    end subroutine SCF_HF

    subroutine diag_SCF(Fock_matrix_prime,Fock_matrix_prime_bar,C_Fock)
        use constant, only: size_matrix_Slater,S_minushalf
        use diagonalization,only : hshtri, tqli, bubble_sort
        use Other_function,only: write_matrix

        implicit none

        double precision, dimension(:,:),allocatable,intent(out)   :: Fock_matrix_prime_bar
        double precision, dimension(:,:),allocatable,intent(out)   :: C_Fock
        
        double precision, dimension(:,:),allocatable,intent(in)    :: Fock_matrix_prime

        double precision,dimension(:,:),allocatable                :: Z,T,T_trans,TMP
        double precision,dimension(:),allocatable                  :: e,diag,diag_sorted
        integer,dimension(:),allocatable                           :: diag_id

        integer             :: M
        integer             :: i,k

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Initialisation
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        M = size_matrix_Slater
        allocate(Z(M,M),T(M,M),T_trans(M,M))
        allocate(TMP(M,M))
        allocate(e(M),diag(M),diag_id(M),diag_sorted(M))
        allocate(C_Fock(M,M),Fock_matrix_prime_bar(M,M))
        C_Fock = 0.0d0 ; Z = Fock_matrix_prime

        diag = 0.0d0 ; e = 0.0d0 ; diag_sorted = 0.0d0
        T = 0.0d0  ; diag_id = 0

        !%%%%%%%%%%%%%%%%%
        ! Diagonalisation
        !%%%%%%%%%%%%%%%%%

        write(*,*) 'Tridiagonalization based on the Householder transformation'
        call hshtri(Z, M, diag, e, .false.)
        
        write(*,*) 'Diagonalization based on the QL method'
        call tqli(diag, e, M, M, Z)
        
        write(*,*) 'Sorting eigenvalues and eigenvectors of F_prime'
        call bubble_sort(diag, M, diag_sorted, diag_id)
    
        write(*,*) 'Generating matrix of eigenvectors'
        diag = diag_sorted
    
        do i = 1, M
           k = diag_id(i)
           T(:,i) = Z(:,k)
        end do
        write(*,*)
    
        write(*,*) 'Eigenvectors of F_prime (column format):'
        call write_matrix(Fock_matrix_prime)
        
        ! C = S^(-1/2)C'
        
        C_Fock = matmul(S_minushalf,T)
        write(*,*)

        write(*,*) 'Eigenvectors of C (column format):'
        call write_matrix(C_Fock)

        write(*,*)
        write(*,*) 'Generating matrix of associated eigenvalues of F_prime'
        do i = 1, M
           Fock_matrix_prime_bar(i,i) = diag(i)
        end do
    
        write(*,*) 'Eigenvalues of F_prime a.k.a. F_prime_bar (column format):'
        call write_matrix(Fock_matrix_prime_bar)
        write(*,*)

        ! epsilon_1 = -0.918164 (paper) |  epsilon_1 = -0.91816350674952829 (logic)
        write(*,*) "Lowest energy eigenvalue a.k.a. epsilon_1 = ", Fock_matrix_prime_bar(1,1)
        write(*,*)

        deallocate(Z,T,T_trans)
        deallocate(TMP)
        deallocate(e,diag,diag_id,diag_sorted)

    end subroutine diag_SCF

    subroutine after_SCF()
        use input_HF_PCCP, only: write_energy
        use constant, only: size_matrix_Slater,thr_SCF,S
        use constant, only: Rpq,Fock_matrix,J_11,E_total, E_verification, eigenvalues
        use Other_function, only: write_matrix
        implicit none

        double precision                                    :: E_consist
        double precision, dimension(:,:),allocatable        :: indempotency, commutation1, commutation2

        integer         :: M

        M = size_matrix_Slater



        write(*,*) "%%%%% Iterative SCF process has terminated! %%%%%"
        write(*,*)
        write(*,*) "The final energy is: ", E_total
    
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Check consistency between two formulae for the energy
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ! E = 2* epsilon_{1} - J_{11}
    
        E_consist = 2.0d0 * eigenvalues(1,1)
        E_consist = E_consist - J_11
    
        write(*,*) "Check: E_consist = ", E_consist
        write(*,*)
    
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Idempotency of density matrix
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        write(*,*) "Idempotency property:"
        write(*,*)
        
        ! Compute density matrix
        write(*,*) "R matrix (column format):"
        call write_matrix(Rpq)
         
        ! Idempotency property
    
        allocate(indempotency(M,M))
        indempotency = matmul(Rpq,matmul(S,Rpq))
        write(*,*)
        write(*,*) "RSR (column format):"
        call write_matrix(indempotency)
        write(*,*)
        write(*,*) "If equal we have idempotency."
        write(*,*)
    
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Density and Fock matrix commutation
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        write(*,*) "Commutation property:"
        write(*,*)
        allocate(commutation1(M,M),commutation2(M,M))
        commutation1 = matmul(Fock_matrix,matmul(Rpq,S))
        commutation2 = matmul(S,matmul(Rpq,Fock_matrix))
        write(*,*) "FRS ="
        call write_matrix(commutation1)
        write(*,*)
        write(*,*) "SRF ="
        call write_matrix(commutation2)
        write(*,*)
        write(*,*) "If equal we have commutation."
        write(*,*)
        deallocate(commutation2,commutation1,indempotency)
    
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Addition of the Koopman theorem
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        ! EI = - h_{N} - sum^{N}_{j=1} (2J_{Nj} - K_{Nj}) = - epsilon_{HOMO}
        ! Approximative value of He+ ~ -1.99916 Ha
        write(*,"(A54)") "%%%% Energy following the Koompan theorem for He+ %%%%"
        write(*,"(A17,(1D30.20),A3)") "Value computed = ", - eigenvalues(1,1) ," Ha"
        write(*,"(A54,(1D30.20),A3)") "Total energy - energy of the hydrogenoide calculated =", -1.99916 - E_total ," Ha"
        write(*,"(A29)") "Experimental value = 0.904 Ha" 
        ! We compare our eigenvalue with the experimental one: 0.904 (a.u.)
    
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Write final energy output to text file
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        call write_energy(E_total,eigenvalues(1,1))
        
        ! Energy checks
        write(*,*) "%%%%% Energy checks (temporary for debugging) %%%%%"
        write(*,*)
        write(*,"(A85,(1D30.20))") "Total energy calculated from sum^{M}_{p=1}(sum^{M}_{q=1}(c_p1 * c_q1 *(h_pq + F_pq)):", E_total
        write(*,"(A44,(1D30.20))") "Total energy calculated from 2 * I_1 + J_11:", E_verification
        write(*,"(A50,(1D30.20))") "Total energy calculated from 2 * epsilon_1 - J_11:", E_consist
        write(*,*) "Threshold", thr_SCF
		
		write(*,*)
		write(*,*) "%%%%% Datas asked for in Part C %%%%%"
		write(*,*) "Total energy:", E_total
		write(*,*) "Occupied orbital energy", eigenvalues(1,1)
        
    end subroutine after_SCF
end module Calculation_SCF

module Other_function
    implicit none
    
contains
    subroutine write_matrix(matrix)
        ! Only 2D matrix can be printed with this subroutine
        implicit none
        double precision,dimension(:,:),allocatable,intent(in) :: matrix
        integer         :: i,j

        do i = 1, size(matrix,1)
            write(*,'(40f12.8)') (matrix(i,j), j=1, size(matrix,2))
         end do

    end subroutine write_matrix
end module Other_function