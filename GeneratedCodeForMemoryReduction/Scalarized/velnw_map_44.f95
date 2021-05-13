module singleton_module_velnw_map_44

contains

subroutine velnw_map_44_scal(p0,ro,dys,v,dt,g)
    ! local vars: pz,i,j,k
    real(kind=4) :: pz
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension(0:(ip+2),0:(jp+2),0:(kp+1)), intent(in) :: p0
    real, dimension(0:jp), intent(in) :: dys
    real, dimension(0:ip,0:jp,0:kp), intent(in) :: g
    real(kind=4), intent(in) :: ro
    real(kind=4), intent(in) :: dt
! WRITTEN
! READ & WRITTEN
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: v
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
        pz = (-p0(i,j,k) + p0(i,j + 1,k)) / ro / dys(j)
        v(i,j,k) = v(i,j,k) + dt * (g(i,j,k) - pz)
end subroutine velnw_map_44

end module singleton_module_velnw_map_44

