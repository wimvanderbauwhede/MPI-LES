module singleton_module_velnw_map_52

contains

subroutine velnw_map_52_scal(p0,ro,dzs,w,dt,h)
    ! local vars: pz,i,j,k
    real(kind=4) :: pz
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension(0:(ip+2),0:(jp+2),0:(kp+1)), intent(in) :: p0
    real, dimension((-1):(kp+2)), intent(in) :: dzs
    real, dimension(0:ip,0:jp,0:kp), intent(in) :: h
    real(kind=4), intent(in) :: ro
    real(kind=4), intent(in) :: dt
! WRITTEN
! READ & WRITTEN
    real, dimension(0:(ip+1),(-1):(jp+1),(-1):(kp+1)), intent(in) :: w
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
        pz = (-p0(i,j,k) + p0(i,j,k + 1)) / ro / dzs(k)
        w(i,j,k) = w(i,j,k) + dt * (h(i,j,k) - pz)
end subroutine velnw_map_52

end module singleton_module_velnw_map_52

