module singleton_module_press_reduce_129

contains

subroutine press_reduce_129_scal(p0,dx1,dy1,dzn,pav,pco)
    ! local vars: i,j,k
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension(0:(ip+2),0:(jp+2),0:(kp+1)), intent(in) :: p0
    real, dimension((-1):(ip+1)), intent(in) :: dx1
    real, dimension(0:(jp+1)), intent(in) :: dy1
    real, dimension((-1):(kp+2)), intent(in) :: dzn
! WRITTEN
! READ & WRITTEN
    real(kind=4), intent(inout) :: pav
    real(kind=4), intent(InOut) :: pco
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
                pav = pav + p0(i,j,k) * dx1(i) * dy1(j) * dzn(k)
                pco = pco + dx1(i) * dy1(j) * dzn(k)
end subroutine press_reduce_129

end module singleton_module_press_reduce_129

