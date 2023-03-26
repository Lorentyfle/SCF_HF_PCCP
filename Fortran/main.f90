program main
  
    !%%%%%%%%%%%%%
    ! Declaration
    !%%%%%%%%%%%%%

    use constant,       only : sent_coef_Slater,size_matrix_Slater,p_SCF,cp1_SCF,thr_SCF
    use input_HF_PCCP,  only : read_data,write_coefficient,write_energy,write_intialisation
    use diagonalization,only : hshtri, tqli, bubble_sort
    
    implicit none
    
    integer                                             :: M, i, j, k, l, atomic_charge,loop
    double precision, dimension(:,:), allocatable       :: S, Z, T, Sbar, Sbar_minushalf, S_minushalf, T_trans, TMP, h
    double precision, dimension(:), allocatable         :: alpha, e, diag, diag_sorted
    integer, dimension(:), allocatable                  :: diag_id

    double precision                                    :: sum_alpha, product_alpha, pairsum, pairsum_2, f1, f2, f3
    double precision, dimension(:,:,:,:), allocatable   :: pqrs
    double precision, dimension(:,:),allocatable        :: coefficients_Fock
    
    double precision                                    :: E_tot,E_tot_old,Verif_Dens,E_consist
    double precision, dimension(:,:),allocatable        :: Fock_matrix, Fock_matrix_prime,Fock_matrix_prime_bar,Rpq,Density
    double precision, dimension(:,:),allocatable        :: indempotency, commutation1,commutation2
    double precision, dimension(:,:),allocatable        :: C_Fock
    ! We take the datas from the files
    
    call read_data()
    write(*,*) "All the constant extracted from the .txt file, are:"
    write(*,*) "sent_coef_Slater", sent_coef_Slater
    write(*,*) "size_matrix_Slater", size_matrix_Slater
    write(*,*) "p_SCF", p_SCF
    write(*,*) "cp1_SCF", cp1_SCF
    write(*,*) "thr_SCF", thr_SCF
    call write_intialisation()

    !%%%%%%%%%%%%%%%%
    ! Initialization
    !%%%%%%%%%%%%%%%%

    write(*,*)
    write(*,*) 'Initializing variables'

    atomic_charge = 2
    M = size_matrix_Slater

    allocate(S(M,M), Z(M,M), T(M,M), Sbar(M,M), Sbar_minushalf(M,M), S_minushalf(M,M), T_trans(M,M), TMP(M,M), h(M,M))
    allocate(alpha(M), e(M), diag(M), diag_sorted(M), diag_id(M))
    allocate(pqrs(M,M,M,M))

    e = 0.0d0 ; diag = 0.0d0 ; diag_sorted = 0.0d0 ; diag_id = 0
    T = 0.0d0 ; Sbar = 0.0d0 ; Sbar_minushalf = 0.0d0 ; S_minushalf = 0.0d0 ; T_trans = 0.0d0

    ! initialize alpha
    do i = 1, M
       alpha(i) = sent_coef_Slater(i)
    end do

    ! initialize S
    do i = 1, M
       do j = 1, M
            if (i .eq. j) then
                S(i,j) = 1.0d0
            else if (i .ne. j) then
                S(i,j) = dble((2*sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**3)
            end if
        end do
    end do

    !%%%%%%%%%%%%%%%
    ! Diagonalize S
    !%%%%%%%%%%%%%%%

    write(*,*) 'S matrix to diagonalise:'
    do i = 1, M
       write(*,'(100F10.4)') (S(i,j), j = 1, M)
    end do
    write(*,*)

    Z = S

    write(*,*) 'Tridiagonalization based on the Householder transformation'
    call hshtri(Z, M, diag, e, .false.)
    
    write(*,*) 'Diagonalization based on the QL method'
    call tqli(diag, e, M, M, Z)
    
    write(*,*) 'Sorting eigenvalues and eigenvectors of S'
    call bubble_sort(diag, M, diag_sorted, diag_id)

    write(*,*) 'Generating matrix of eigenvectors'
    diag = diag_sorted

    do i = 1, M
       j = diag_id(i)
       T(:,i) = Z(:,j)
    end do

    write(*,*) 'Eigenvectors of S a.k.a. T (column format):'
    do i = 1, size(S,1)
       write(*,'(40f12.8)') (T(i,j), j=1, size(S,2))
    end do

    write(*,*)
    write(*,*) 'Generating matrix of associated eigenvalues of S'
    do i = 1, M
       Sbar(i,i) = diag(i)
    end do

    write(*,*) 'Eigenvalues of S a.k.a. Sbar (column format):'
    do i = 1, size(S,1)
       write(*,'(40f12.8)') (Sbar(i,j), j=1, size(S,2))
    end do
    write(*,*)

    !%%%%%%%%%%%%%%%%%%%
    ! Compute Sbar^-1/2
    !%%%%%%%%%%%%%%%%%%%
    
    do i = 1, M
       Sbar_minushalf(i,i) = dble(1/(sqrt(Sbar(i,i))))
    end do
	
    write(*,*) 'Sbar^-1/2 (column format):'
    do i = 1, size(S,1)
       write(*,'(40f12.8)') (Sbar_minushalf(i,j), j=1, size(S,2))
    end do
    write(*,*)

    !%%%%%%%%%%%%%%%%
    ! Compute S^-1/2
    !%%%%%%%%%%%%%%%%
    
    T_trans = transpose(T)
    TMP = matmul(Sbar_minushalf, T_trans)
    S_minushalf = matmul(T, TMP)

    write(*,*) 'S^-1/2 (column format):'
    do i = 1, size(S,1)
       write(*,'(40f12.8)') (S_minushalf(i,j), j=1, size(S,2))
    end do
    write(*,*)

    !%%%%%%%%%%%%%%%%
    ! h_pq integrals
    !%%%%%%%%%%%%%%%%

    write(*,*) 'Computing numerical values of h_pq integrals:'
    do i = 1, M
        do j = 1, M
            sum_alpha = alpha(i) + alpha(j)
            product_alpha = alpha(i) * alpha(j)
            
            h(i,j) = 4 * dble(((sqrt(product_alpha)/sum_alpha)**3) * (product_alpha - atomic_charge*sum_alpha))
            
            write(*,*) "h_pq with p = ", i, "and q =", j
            write(*,'(40f12.8)') h(i,j)
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

                    f1 = dble(1/((pairsum)**3 * (pairsum_2)**2))
                    f2 = dble(1/((pairsum)**3 * (sum_alpha)**2))
                    f3 = dble(1/((pairsum)**2 * (sum_alpha)**3))
                
                    pqrs(i,j,k,l) = 32 * ((sqrt(product_alpha))**3 * (f1 - f2 - f3))
                    write(*,*) "(pq|rs) with p = ", i, "q =", j, "r = ", k, "s = ", l
                    write(*,'(40f12.8)') pqrs(i,j,k,l)
                end do
            end do
        end do
    end do

    !%%%%%%%%%%%%%%%
    ! Iterative SCF
    !%%%%%%%%%%%%%%%

    write(*,*)
    write(*,*) "Beginning iterative SCF"
    
    ! Setup
    allocate(coefficients_Fock(M,1))
    write(*,*) "p_SCF", p_SCF
    write(*,*) "cp1_SCF", cp1_SCF
    write(*,*) "Matrix size", M
    write(*,*) "thr_SCF", thr_SCF

    ! Define initial guesses for Fock coefficients
    do i = 1, M
        do j = 1, M
            if (i .eq. p_SCF) then
                coefficients_Fock(i,j) = cp1_SCF(1)
            else
                coefficients_Fock(i,j) = cp1_SCF(2)
            end if
        end do
    end do
    
    ! Write initial guesses for Fock coefficients
    write(*,*)
    write(*,*) 'Initial values of the coefficients'
    do i = 1, size(coefficients_Fock,1)
        write(*,'(40f12.8)') (coefficients_Fock(i,j), j=1, size(coefficients_Fock,2))
    end do
    
    allocate(Fock_matrix(M,M),Fock_matrix_prime(M,M),Fock_matrix_prime_bar(M,M),Rpq(M,M),Density(M,M))
    allocate(C_Fock(M,M))

    !%%%%%%%%%%%%%%%%%%%%%%%
    ! Iterative SCF Process
    !%%%%%%%%%%%%%%%%%%%%%%%
    
    E_tot = 0
    E_tot_old = (E_tot+1)*5
    loop = 0
    
    do while(abs(E_tot - E_tot_old) .gt. thr_SCF)
        loop = loop + 1

        ! Fock matrix elements
        write(*,*)
        write(*,*) "loop", loop
        write(*,*) "Coefficients sent:"
        do i = 1, size(coefficients_Fock,1)
            write(*,'(40f12.8)') (coefficients_Fock(i,j), j=1, 1)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the Fock matrix elements and build the F matrix
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        do i = 1, M
            do j = 1, M
                Fock_matrix(i,j) = h(i,j)
                !if ( loop == 1 ) then
                !    write(*,*) "h(i,j) =", Fock_matrix(i,j), h(i,j)
                !end if
                do k = 1, M
                    do l = 1, M
                        Fock_matrix(i,j) = Fock_matrix(i,j) + coefficients_Fock(k,1) * coefficients_Fock(l,1) * pqrs(i,j,k,l)
                        !if ( loop == 1 ) then
                        !    write(*,*) "k =",k, coefficients_Fock(k,1),"| l =",l,coefficients_Fock(l,1), &
                        !    "| pqrs =" ,pqrs(i,j,k,l) ,"| F =", Fock_matrix(i,j)
                        !end if
                    end do
                end do
            end do
        end do

        write(*,*)
        write(*,*) "Fock matrix (column format):"
        do i = 1, size(Fock_matrix,1)
            write(*,'(40f12.8)') (Fock_matrix(i,j), j=1, size(Fock_matrix,2))
        end do
        write(*,*)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Obtain the Fock matrix eigenvectors and eigenvalues
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ! Fock prime matrix
        Fock_matrix_prime = matmul(S_minushalf,matmul(Fock_matrix,S_minushalf))
        
        write(*,*) "F' = S^(-1/2)FS^(-1/2)"
        write(*,*) 'F_prime matrix to diagonalise:'
        do i = 1, M
           write(*,'(100F10.4)') (Fock_matrix_prime(i,j), j = 1, M)
        end do
        write(*,*)
        
        ! Re-Init variables for sorting
        Z = Fock_matrix_prime
        diag = 0.0d0
        e = 0.0d0
        diag_sorted = 0.0d0
        T = 0.0d0
        C_Fock = 0.0d0
        diag_id = 0

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
        do i = 1, size(Fock_matrix_prime,1)
           write(*,'(40f12.8)') (T(i,j), j=1, size(Fock_matrix_prime,2))
        end do
        
        ! C = S^(-1/2)C'
        
        !write(*,*) S_minushalf
        !write(*,*) T
        C_Fock = matmul(S_minushalf,T)

        write(*,*) 'Eigenvectors of C:'
        do i = 1, size(C_Fock,1)
           write(*,'(40f12.8)') (C_Fock(i,j), j=1, size(C_Fock,2))
        end do

        write(*,*)
        write(*,*) 'Generating matrix of associated eigenvalues of F_prime'
        do i = 1, M
           Fock_matrix_prime_bar(i,i) = diag(i)
        end do
    
        write(*,*) 'Eigenvalues of F_prime a.k.a. F_prime_bar (column format):'
        do i = 1, size(Fock_matrix_prime_bar,1)
           write(*,'(40f12.8)') (Fock_matrix_prime_bar(i,j), j=1, size(Fock_matrix_prime_bar,2))
        end do
        write(*,*)
        
        write(*,*) "Lowest energy eigenvalue a.k.a. epsilon_1 = ", Fock_matrix_prime_bar(1,1)
        write(*,*)

        !%%%%%%%%%%%%%%%%%%%%%%
        ! Density matrix check
        !%%%%%%%%%%%%%%%%%%%%%%

        write(*,*) "Beginning density matrix verification"
        
        ! D = 2R
        ! Rpq = Cp1 * Cq1
        do i = 1, M
            do j = 1, M
                Rpq(i,j) = C_Fock(i,1) * C_Fock(j,1)
            end do
        end do

        write(*,*) "R matrix (column format):"
        do i = 1, size(Rpq,1)
            write(*,'(40f12.8)') (Rpq(i,j), j=1, size(Rpq,2))
        end do
        write(*,*)

        Density = 2 * Rpq
        
        write(*,*) "Density matrix (column format):"
        do i = 1, size(Density,1)
            write(*,'(40f12.8)') (Density(i,j), j=1, size(Density,2))
        end do
        write(*,*)
        
        Verif_Dens = 0.0d0
        do i = 1, M
            do j = 1, M
                Verif_Dens = Verif_Dens + Density(i,j) * S(i,j)
            end do
        end do

        write(*,*) "We have:", Verif_Dens, " =? ", atomic_charge 

        !%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Total energy from eq. 9
        !%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Formula for the energy
        ! E = sum^{M}_{p=1}(sum^{M}_{q=1}(C_p1 * C_q1 *(h_pq + F_pq)))

        E_tot_old = E_tot
        E_tot = 0.0d0
        do i = 1, M
            do j = 1, M
                ! Calculation needs to be properly done for logical results

                ! not-logic attempt #1
               
                ! E_tot = E_tot + Rpq(i,1) * Rpq(j,1) * (h(i,j) + Fock_matrix(i,j))
                ! E_final =~ -2.023296 (not logic)      | E_init = - 1.92495

                ! not-logic attempt #2
               
                ! E_tot = E_tot + Rpq(i,j) * (h(i,j) + Fock_matrix(i,j))
                ! E_final =~ -2.861672                  | E_init = - 2.93938 (not logic)

                ! logic attempt #3
               
                E_tot = E_tot + coefficients_Fock(i,1) * coefficients_Fock(j,1) * (h(i,j) + Fock_matrix(i,j))
                ! E_final = -2.8616715937556600 (logic) | E_init = -2.7912499999999998 (logic but don't follow the paper)

                ! Results from the reference papers Snow and Bills (see p. 3):
                ! E_final = -2.86167                    | E_init = - 2.83308 
            end do
        end do

        write(*,*)
        ! write(*,*) "c11 =", coefficients_Fock(1,1), "c21 =", coefficients_Fock(2,1)
        ! write(*,*) "F11 =", Fock_matrix(1,1),"F12 =", Fock_matrix(1,2), "F22 =", Fock_matrix(2,2) 
        ! After verification, all this values are correct for the 1st loop. (c11,c21,F11,F12,F22 and epsilon)
        ! Only the total energy is different for the first loop.

        ! Compare energies
        write(*,*) "The previous energy is: ",  E_tot_old
        write(*,*) "The actual energy is:",     E_tot
        write(*,*)

        ! Compare Fock coefficients
        write(*,*) "Previous Fock coefficients a.k.a. coefficients_Fock (column format):"
        do i = 1, size(coefficients_Fock,1)
            write(*,'(40f12.8)') (coefficients_Fock(i,j), j=1, size(coefficients_Fock,2))
        end do
        write(*,*)
        write(*,*) "vs"
        write(*,*)
        write(*,*) "New Fock coefficients a.k.a. C_Fock (column format):"
        do i = 1, size(C_Fock,1)
            write(*,'(40f12.8)') (C_Fock(i,j), j=1, size(C_Fock,2))
        end do
        write(*,*)

        ! Update the Fock coefficients for the next loop
        write(*,*) "Editing of the starting coefficients"
        do i = 1, M
            coefficients_Fock(i,1) = C_Fock(i,1)
        end do
        
        !write(*,*) "New Fock coefficients (column format):"
        !do i = 1, size(coefficients_Fock,1)
        !    write(*,'(40f12.8)') (coefficients_Fock(i,j), j=1, size(coefficients_Fock,2))
        !end do

        ! Output actual coefficients to text file
        call write_coefficient(coefficients_Fock(:,1),loop)
    end do
    write(*,*)

    !%%%%%%%%%%%%%%%%%%%%%%%
    ! After SCF convergence
    !%%%%%%%%%%%%%%%%%%%%%%%

    write(*,*) "Iterative SCF process has terminated."
    write(*,*) "The final energy is: ", E_tot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Check consistency between two formulae for the energy
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! E = 2* epsilon_{1} - J_{11}

    E_consist = 2 * Fock_matrix_prime_bar(1,1)
    do i = 1, M
        do j = 1, M
            do k = 1, M
                do l = 1, M
                    ! E not right, E_consist found is too high,                E_consist = -2.3488957671479835
                    ! E_consist = E_consist - (Rpq(i,1) * Rpq(j,1) * Rpq(k,1) * Rpq(l,1) * pqrs(i,j,k,l))

                    ! This also doesn't work... but it's closer to working,    E_consist = -2.8236751900322368
                    E_consist = E_consist - (coefficients_Fock(i,1) * coefficients_Fock(j,1) &
                                * coefficients_Fock(k,2) * coefficients_Fock(l,2) * pqrs(i,j,k,l))
                end do                
            end do
        end do
    end do
    !write(*,*) 2*Fock_matrix_prime_bar(1,1) !- (Rpq(1,1)*Rpq(2,1)* pqrs(1,1,1,1))
    write(*,*) "Check: E consist = ", E_consist

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Idempotency of density matrix
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! Compute density matrix
    write(*,*) "R matrix (column format):"
    do i = 1, size(Rpq,1)
        write(*,'(40f12.8)') (Rpq(i,j), j=1, size(Rpq,2))
    end do
     
    ! Idempotency property

    allocate(indempotency(M,M))
    indempotency = matmul(Rpq,matmul(S,Rpq))
    write(*,*)
    write(*,*) "RSR (column format):"
    do i = 1, size(indempotency,1)
        write(*,'(40f12.8)') (indempotency(i,j), j=1, size(indempotency,2))
    end do
    write(*,*)
    write(*,*) "If equal we have idempotenty."
    write(*,*)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Density and Fock matrix commutation
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    write(*,*) "Commutation property, FRS = SRF"
    allocate(commutation1(M,M),commutation2(M,M))
    commutation1 = matmul(Fock_matrix,matmul(Rpq,S))
    commutation2 = matmul(S,matmul(Rpq,Fock_matrix))
    write(*,*) "FRS ="
    do i = 1, size(commutation1,1)
        write(*,'(40f12.8)') (commutation1(i,j), j=1, size(commutation1,2))
    end do
    write(*,*)
    write(*,*) "SRF ="
    do i = 1, size(commutation2,1)
        write(*,'(40f12.8)') (commutation2(i,j), j=1, size(commutation2,2))
    end do
    write(*,*)
    write(*,*) "If equal we have commutation."
    write(*,*)
    deallocate(commutation2,commutation1,indempotency)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Write final energy output to text file
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call write_energy(E_tot)
    
    

    
end program main

