module singleton_module_velnw_map_36

contains

subroutine velnw_map_36_scal(p0,ro,dxs,u,dt,f)
    ! local vars: pz,i,j,k
    real(kind=4) :: pz
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension(0:(ip+2),0:(jp+2),0:(kp+1)), intent(in) :: p0
    real, dimension(0:ip), intent(in) :: dxs
    real, dimension(0:ip,0:jp,0:kp), intent(in) :: f
    real(kind=4), intent(in) :: ro
    real(kind=4), intent(in) :: dt
! WRITTEN
! READ & WRITTEN
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: u
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
        pz = (-p0(i,j,k) + p0(i + 1,j,k)) / ro / dxs(i)
        u(i,j,k) = u(i,j,k) + dt * (f(i,j,k) - pz)
end subroutine velnw_map_36

end module singleton_module_velnw_map_36

