program main
    ! ***************
    ! * Declaration *
    ! ***************
    use constant, only : sent_coef_Slater,size_matrix_Slater,p_SCF,cp1_SCF,thr_SCF
    use input_HF_PCCP, only : read_data
    implicit none


    ! We take the datas from the files
    call read_data()
    write(*,*) "All the constant extracted from the .txt file, are:"
    write(*,*) sent_coef_Slater,size_matrix_Slater,p_SCF,cp1_SCF,thr_SCF

end program main


