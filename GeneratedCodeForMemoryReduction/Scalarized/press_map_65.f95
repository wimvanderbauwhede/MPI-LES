module singleton_module_press_map_65

contains

subroutine press_map_65_scal(u,dx1,v,dy1,w,dzn,f,g,h,rhs,dt)
    ! local vars: i,j,k
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: u
    real, dimension((-1):(ip+1)), intent(in) :: dx1
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: v
    real, dimension(0:(jp+1)), intent(in) :: dy1
    real, dimension(0:(ip+1),(-1):(jp+1),(-1):(kp+1)), intent(in) :: w
    real, dimension((-1):(kp+2)), intent(in) :: dzn
    real, dimension(0:ip,0:jp,0:kp), intent(in) :: f
    real, dimension(0:ip,0:jp,0:kp), intent(in) :: g
    real, dimension(0:ip,0:jp,0:kp), intent(in) :: h
    real(kind=4), intent(in) :: dt
! WRITTEN
! READ & WRITTEN
    real, dimension(0:(ip+1),0:(jp+1),0:(kp+1)), intent(inout) :: rhs
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
                rhs(i,j,k) = (-u(i - 1,j,k) + u(i,j,k)) / dx1(i) + (-v(i,j - 1,k) + v(i,j, &
      k)) / dy1(j) + (-w(i,j,k - 1) + w(i,j,k)) / dzn(k)
                rhs(i,j,k) = (f(i,j,k) - f(i - 1,j,k)) / dx1(i) + (g(i,j,k) - g(i,j - 1, &
      k)) / dy1(j) + (h(i,j,k) - h(i,j,k - 1)) / dzn(k) + rhs(i,j,k) / dt
end subroutine press_map_65

end module singleton_module_press_map_65

