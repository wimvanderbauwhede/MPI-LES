module singleton_module_press_map_87

contains

subroutine press_map_87_scal(rhs,rhsav)
    ! local vars: i,j,k
    ! parallelfortran: synthesised loop variable decls
! READ
    real(kind=4), intent(in) :: rhsav
! WRITTEN
! READ & WRITTEN
    real, dimension(0:(ip+1),0:(jp+1),0:(kp+1)), intent(inout) :: rhs
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
                rhs(i,j,k) = rhs(i,j,k) - rhsav
end subroutine press_map_87

end module singleton_module_press_map_87

