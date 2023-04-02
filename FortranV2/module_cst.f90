module constant
    implicit none
    
    save

    ! Variables
    double precision, allocatable:: sent_coef_Slater(:)
    integer                      :: size_matrix_Slater
    integer                      :: p_SCF
    integer,parameter            :: size_cp_SCF = 2
    integer                      :: cp1_SCF(size_cp_SCF)
    double precision             :: thr_SCF
    integer                      :: mesh
    integer                      :: atomic_charge
    double precision, dimension(:,:), allocatable :: S, S_minushalf
    double precision, dimension(:,:), allocatable :: h_pq
    double precision, dimension(:,:,:,:), allocatable :: pqrs
    double precision, dimension(:), allocatable   :: alpha_Slater
    double precision, dimension(:,:),allocatable  :: coefficients_Fock
    double precision, dimension(:,:),allocatable  :: eigenvalues
    double precision, dimension(:,:),allocatable  ::  Rpq
    double precision, dimension(:,:),allocatable  ::  Fock_matrix
    double precision            :: E_total,E_verification,I_1, J_11

    
end module constant