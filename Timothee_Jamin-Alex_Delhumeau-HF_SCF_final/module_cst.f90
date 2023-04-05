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
    
end module constant