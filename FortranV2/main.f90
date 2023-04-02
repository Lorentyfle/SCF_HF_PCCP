program main
  
    !%%%%%%%%%%%%%
    ! Declaration
    !%%%%%%%%%%%%%

    use constant,       only : sent_coef_Slater,size_matrix_Slater,thr_SCF,atomic_charge
    use constant,       only : S,S_minushalf,h_pq,pqrs,alpha_Slater
    use constant,       only : coefficients_Fock,Rpq,Fock_matrix,E_total,I_1,J_11
    use input_HF_PCCP,  only : read_data,write_coefficient,write_energy,write_intialisation
    use diagonalization,only : hshtri, tqli, bubble_sort
    use Other_function, only : write_matrix
    use Calculation_SCF, only: alpha_coef,S_matrix,integrals_He,Setup_SCF,SCF_HF,diag_SCF,after_SCF

    implicit none

    ! We take the datas from the files
    write(*,*) "Extracting input parameters from input.txt"
    write(*,*)
    call read_data()
    write(*,*) "Coefficients extracted: ", sent_coef_Slater
    write(*,*) "Size of the Slater matrix =", size_matrix_Slater
    write(*,*) "SCF threshold =", thr_SCF
    call write_intialisation()

    !%%%%%%%%%%%%%%%%
    ! Initialization
    !%%%%%%%%%%%%%%%%

    write(*,*)
    write(*,*) 'Initializing variables'
    write(*,*)
    ! Because we have an He atom with two electrons, we have:
    atomic_charge = 2

    !%%%%%%%%%%%%%%%%%%
    ! initialize alpha
    !%%%%%%%%%%%%%%%%%%

    call alpha_coef()

    !%%%%%%%%%%%%%%%%%%%%%%
    ! initialize S, S^-1/2
    !%%%%%%%%%%%%%%%%%%%%%%
    
    call S_matrix(size_matrix_Slater,alpha_Slater,S,S_minushalf)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! h_pq and (pq|rs) integrals
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call integrals_He(alpha_Slater,h_pq,pqrs)

    !%%%%%%%%%%%%%%%
    ! Iterative SCF
    !%%%%%%%%%%%%%%%

    call Setup_SCF(coefficients_Fock)

    !%%%%%%%%%%%%%%%%%%%%%%%
    ! Iterative SCF Process
    !%%%%%%%%%%%%%%%%%%%%%%%
    write(*,*) "Hey"    
    call SCF_HF(E_total,Rpq,Fock_matrix,I_1,J_11)

    !%%%%%%%%%%%%%%%%%%%%%%%
    ! After SCF convergence
    !%%%%%%%%%%%%%%%%%%%%%%%

    call after_SCF()
    
end program main

