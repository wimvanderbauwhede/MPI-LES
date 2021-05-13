module singleton_module_press_map_138

contains

subroutine press_map_138_scal(p0,pav)
    ! local vars: i,j,k
    ! parallelfortran: synthesised loop variable decls
! READ
    real(kind=4), intent(in) :: pav
! WRITTEN
! READ & WRITTEN
    real, dimension(0:(ip+2),0:(jp+2),0:(kp+1)), intent(inout) :: p0
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
                p0(i,j,k) = p0(i,j,k) - pav
end subroutine press_map_138

end module singleton_module_press_map_138

