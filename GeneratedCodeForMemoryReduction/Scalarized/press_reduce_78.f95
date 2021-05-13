module singleton_module_press_reduce_78

contains

subroutine press_reduce_78_scal(dx1,dy1,dzn,rhs,rhsav,area)
    ! local vars: i,j,k
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension((-1):(ip+1)), intent(in) :: dx1
    real, dimension(0:(jp+1)), intent(in) :: dy1
    real, dimension((-1):(kp+2)), intent(in) :: dzn
    real, dimension(0:(ip+1),0:(jp+1),0:(kp+1)), intent(in) :: rhs
! WRITTEN
! READ & WRITTEN
    real(kind=4), intent(inout) :: rhsav
    real(kind=4), intent(InOut) :: area
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
                rhsav = rhsav + dx1(i) * dy1(j) * dzn(k) * rhs(i,j,k)
                area = area + dx1(i) * dy1(j) * dzn(k)
end subroutine press_reduce_78

end module singleton_module_press_reduce_78

