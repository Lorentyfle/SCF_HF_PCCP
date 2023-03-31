module diagonalization
implicit none

contains
	subroutine hshtri( A_, N, d_diag, d_off , verbose)
  
	  implicit none

	  ! passed variables
	  integer, intent(in) :: N
	  double precision, dimension(N,N), intent(inout) :: A_ 
	  double precision, dimension(N),   intent(out)   :: d_off 
	  double precision, dimension(N),   intent(out)   :: d_diag

	  logical :: verbose

	  ! local variables  
	  double precision, dimension(N,N) :: A_tri, identity, w, U_, P, Q
	  double precision, dimension(N)   :: vector
	  double precision, dimension(N)   :: a, e, u

	  double precision :: norm_a, sign_test, error, threshold, old_max

	  integer :: i, j, k, iter

	  double precision :: cpu_init
	  double precision :: cpu_fini

	  ! setup identity matrix
	  identity = 0.0d0
	  do i = 1, N
		 identity(i,i) = 1.0d0
	  end do
	  !
	  if (verbose) then

		 write(*,*)
		 write(*,*) '##################################################'
		 write(*,*) ' TRIDIAG start'
		 write(*,*) '##################################################'

	  end if
	  !
	  call cpu_time(cpu_init)
	  !
	  ! Initialisation:
	  A_tri = A_
	  Q     = identity  
	  !
	  A_matrix_tridiag: do k = 1, N - 2
		 !
		 ! build the first unit vector
		 e      = 0.0d0
		 e(k+1) = sign( 1.0d0, a(k+2) )
		 !
		 ! get the vector(a) :=
		 ! the lasts (N-k) elements of the first
		 ! column of A_tri
		 a        = 0.0d0
		 a(k+1:N) = A_tri(k+1:N,k)
		 norm_a   = sqrt( dot_product(a,a) )
		 !norm_a   = norm2( a )

		 !
		 ! compute the vector(u)
		 u = a - norm_a * e
		 !
		 w = 0.0d0
		 ! compute the matrix(w)
		 do i = 1, N 
			do j = 1, N
			   w(i,j) = u(i) * u(j)
			end do
		 end do
		 !
		 ! compute the projector matrix(P)
		 w = w / dot_product(u,u)
		 P = identity - 2.0d0 * w
		 !
		 if (verbose) then
			write(*,*)
			write(*,*) '#### iteration:', k, '####'
			!
			write(*,*) 'projector matrix P'
			do i = 1, N
			   write(*,'(100f16.8)') (P(i,j), j = 1, N)
			end do
		 end if
		 !
		 ! check is P is unitary
		 U_ = matmul( transpose(P), P )
		 !
		 if (verbose) then        
			write(*,*)
			write(*,*) 'P is unitary?'
			!
			do i = 1, N
			   write(*,'(100f16.8)') (U_(i,j), j = 1, N)
			end do
		 end if
		 !
		 ! Backup P in Q such that: Q = P_1 P_2 P_3 ... P_N-2 
		 Q = matmul(Q,P)     
		 !
		 ! applied P to the left of A
		 A_tri = matmul( transpose(P), matmul(A_tri,P) )
		 !
		 if (verbose) then     
			!if ( k /= (N-2) ) then
			write(*,*)
			write(*,*) 'temporary matrix A_tri'
			!
			do i = 1, N
			   write(*,'(100f16.8)') (A_tri(i,j), j = 1, N)
			end do
			!
			!end if
			!
		 end if
		 !
	  end do A_matrix_tridiag

	  call cpu_time(cpu_fini)

	  if (verbose) then

		 write(*,*)
		 write(*,*) '##################################################'
		 write(*,*) ' TRIDIAG stop'
		 write(*,*) '##################################################'
		 !
		 write(*,*)             
		 write(*,*) 'final tridiagonal matrix'
		 !
		 do i = 1, N
			write(*,'(100f16.8)') (A_tri(i,j), j = 1, N)
		 end do
		 !
		 write(*,*)
		 write(*,*) 'final Q matrix'
		 !
		 do i = 1, N
			write(*,'(100f16.8)') (Q(i,j), j = 1, N)
		 end do
		 !       
		 ! check is Q is unitary
		 U_ = matmul( transpose(Q), Q )
		 !
		 write(*,*)
		 write(*,*) 'Q is unitary?'
		 !
		 do i = 1, N
			write(*,'(100f16.8)') (U_(i,j), j = 1, N)
		 end do
		 !
		 ! check is Q A_tri Qt = A?
		 U_ = matmul( Q, matmul(A_tri, transpose(Q) ) )
		 !
		 write(*,*)
		 write(*,*) 'Qt A_tri Q  = A?'
		 !
		 do i = 1, N
			write(*,'(100f16.8)') (U_(i,j), j = 1, N)
		 end do
		 !
	  end if
	  !
	  write(*,*)
	  write (*,'(" CPU time HSHTRI = ",f12.6," seconds")') cpu_fini - cpu_init
	  write(*,*)
	  !
	  ! backup d_diag, offdiag and Q
	  d_diag = 0.0d0
	  do i = 1, N
		 d_diag(i) = A_tri(i,i)
	  end do
	  !
	  d_off = 0.0d0
	  do i = 1, N-1
		 d_off(i+1) = A_tri(i,i+1)
	  end do
	  !
	  A_ = Q

	  return
	end subroutine hshtri
	
	SUBROUTINE TQLI(D,E,N,NP,Z)

	  IMPLICIT NONE

	  INTEGER                             :: N, NP
	  DOUBLE PRECISION, DIMENSION(NP, NP) :: Z
	  DOUBLE PRECISION, DIMENSION(NP)     :: D, E

	  INTEGER          :: I, ITER, K, L, M
	  DOUBLE PRECISION :: B, C, DD, F, G, P, R, S

	  IF ( N .GT. 1 ) THEN
		 DO 11 I = 2, N
			E(I-1) = E(I)
	11   END DO
		 E(N) = 0.0D0
		 DO 15 L = 1, N
			ITER = 0
	1       DO 12 M = L, N-1
			   DD = ABS(D(M))+ABS(D(M+1))
			   IF ( ABS(E(M))+DD .EQ. DD ) GO TO 2
	12      END DO
			M=N
	2       IF ( M .NE. L) THEN
			   IF ( ITER .EQ. 30 ) EXIT
			   !IF(ITER.EQ.30) EXIT 'too many iterations'
			   ITER = ITER + 1
			   G = (D(L+1)-D(L))/(2.0D0*E(L))
			   R = SQRT(G**2+1.0D0)
			   G = D(M)-D(L)+E(L)/(G+SIGN(R,G))
			   S = 1.0D0
			   C = 1.0D0
			   P = 0.0D0
			   DO 14 I = M-1, L, -1
				  F = S*E(I)
				  B = C*E(I)
				  IF ( ABS(F) .GE. ABS(G) ) THEN
					 C = G/F
					 R = SQRT(C**2 + 1.0D0)
					 E(I+1) = F*R
					 S = 1.0D0/R
					 C = C*S
				  ELSE
					 S = F/G
					 R = SQRT(S**2 + 1.0D0)
					 E(I+1) = G*R
					 C = 1.0D0/R  
					 S = S*C
				  ENDIF
				  G     = D(I+1) - P
				  R     = (D(I)-G)*S + 2.0D0*C*B
				  P     = S*R
				  D(I+1)= G + P
				  G     = C*R - B
				  DO 13 K = 1, N
					 F        = Z(K,I+1)
					 Z(K,I+1) = S*Z(K,I) + C*F
					 Z(K,I  ) = C*Z(K,I) - S*F
	13            END DO
	14         END DO
			   D(L)=D(L)-P
			   E(L)=G
			   E(M)=0.0D0
			   GO TO 1
			ENDIF
	15   END DO
	  ENDIF

	  RETURN
	END SUBROUTINE TQLI
	
	subroutine bubble_sort(u_inp, n, u_out, id)
	  !
	  implicit none
	  !
	  integer, intent(in)   :: n  ! input variable
	  double precision, dimension(n), intent(in)  :: u_inp  ! input variable to sort
	  double precision, dimension(n), intent(out) :: u_out  ! output variable sorted
	  integer,          dimension(n), intent(out) :: id     ! output new indexation
	  double precision, dimension(n) :: u ! working variable
	  !
	  double precision :: tmp
	  integer :: i, tmp_id
	  logical :: swap

	  ! initialisation
	  u    = u_inp
	  swap = .true.
	  do i = 1, n
		 id(i) = i
	  end do
	  !  
	  do while ( swap .eqv. .true. )
		 !
		 swap = .false.
		 !
		 loop_do: do i = 1, n - 1
			!
			if ( u(i) > u(i+1) ) then
			   !
			   ! swap!
			   tmp    = u(i  ) ; tmp_id  = id(i)
			   u(i  ) = u(i+1) ; id(i)   = id(i+1)
			   u(i+1) = tmp    ; id(i+1) = tmp_id
			   swap   = .true.
			end if
			!
		 end do loop_do

	  end do

	  u_out = u

	end subroutine bubble_sort

end module diagonalization
