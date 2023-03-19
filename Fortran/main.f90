program main
  
    !%%%%%%%%%%%%%
    ! Declaration
    !%%%%%%%%%%%%%

    use constant, only : sent_coef_Slater,size_matrix_Slater,p_SCF,cp1_SCF,thr_SCF
    use input_HF_PCCP, only : read_data
	use diagonalization, only: hshtri, tqli, bubble_sort
    
    implicit none

    integer :: M, i, j, k, l, atomic_charge
    double precision, dimension(:,:), allocatable :: S, Z, T, Sbar, Sbar_minushalf, S_minushalf, T_trans, TMP, h
    double precision, dimension(:), allocatable :: alpha, e, diag, diag_sorted
    integer, dimension(:), allocatable :: diag_id

    double precision :: sum_alpha, product_alpha, pairsum, pairsum_2, f1, f2, f3
    double precision, dimension(:,:,:,:), allocatable :: pqrs

    ! We take the datas from the files
    
    call read_data()
    write(*,*) "All the constant extracted from the .txt file, are:"
    write(*,*) "sent_coef_Slater", sent_coef_Slater
    write(*,*) "size_matrix_Slater", size_matrix_Slater
    write(*,*) "p_SCF", p_SCF
    write(*,*) "cp1_SCF", cp1_SCF
    write(*,*) "thr_SCF", thr_SCF

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

    write(*,*) 'eigenvectors of S a.k.a. T (column format):'
    do i = 1, size(S,1)
       write(*,'(40f12.8)') (T(i,j), j=1, size(S,2))
    end do

    write(*,*)
    write(*,*) 'generating matrix of associated eigenvalues of S'
    do i = 1, M
       Sbar(i,i) = diag(i)
    end do

    write(*,*) 'eigenvalues of S a.k.a. Sbar (column format):'
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
                f3 = dble(1/(pairsum)**2 * (sum_alpha)**3)
                
				pqrs(i,j,k,l) = 32 * dble((sqrt(product_alpha))**3 * (f1 - f2 - f3))
				
				write(*,*) "(pq|rs) with p = ", i, "q =", j, "r = ", k, "s = ", l
				write(*,'(40f12.8)') pqrs(i,j,k,l)
             end do
          end do
       end do
    end do
    
end program main

