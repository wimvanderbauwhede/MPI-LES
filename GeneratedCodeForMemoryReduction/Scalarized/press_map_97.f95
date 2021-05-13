module singleton_module_press_map_97

contains

subroutine press_map_97_scal(dzs,dys,dxs,nrd,p0,rhs,p1)
        parameter(omega=1.)
    ! local vars: i,j,k,dz1,dz2,cn4s,cn4l,cn3s,cn3l,cn2s,cn2l,cn1,reltmp
    real(kind=4) :: dz1
    real(kind=4) :: dz2
    real(kind=4) :: cn4s
    real(kind=4) :: cn4l
    real(kind=4) :: cn3s
    real(kind=4) :: cn3l
    real(kind=4) :: cn2s
    real(kind=4) :: cn2l
    real(kind=4) :: cn1
    real(kind=4) :: reltmp
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension((-1):(kp+2)), intent(in) :: dzs
    real, dimension(0:jp), intent(in) :: dys
    real, dimension(0:ip), intent(in) :: dxs
    real, dimension(0:(ip+1),0:(jp+1),0:(kp+1)), intent(in) :: rhs
    integer, intent(in) :: nrd
! WRITTEN
! READ & WRITTEN
    real, dimension(0:(ip+2),0:(jp+2),0:(kp+1)), intent(inout) :: p1
    real, dimension(0:(ip+2),0:(jp+2),0:(kp+1)), intent(inout) :: p0
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
                dz1 = dzs(k - 1)
                dz2 = dzs(k)
                cn4s = 2. / (dz1 * (dz1 + dz2))
                cn4l = 2. / (dz2 * (dz1 + dz2))
                    cn3s = 2. / (dys(j - 1) * (dys(j - 1) + dys(j)))
                    cn3l = 2. / (dys(j) * (dys(j - 1) + dys(j)))
                        cn2s = 2. / (dxs(i - 1) * (dxs(i - 1) + dxs(i)))
                        cn2l = 2. / (dxs(i) * (dxs(i - 1) + dxs(i)))
                        cn1 = 1. / (2. / (dxs(i - 1) * dxs(i)) + 2. / (dys(j - 1) * dys(j)) + 2. /  &
      (dz1 * dz2))
                      if (nrd == 0) then
                        reltmp = omega * (cn1 * (cn2l * p0(i + 1,j,k) + cn2s * p0(i - 1,j, &
      k) + cn3l * p0(i,j + 1,k) + cn3s * p0(i,j - 1,k) + cn4l * p0(i,j,k + 1) + cn4s * p0(i,j, &
      k - 1) - rhs(i,j,k)) - p0(i,j,k))
                        p1(i,j,k) = p0(i,j,k) + reltmp
                      else
                        reltmp = omega * (cn1 * (cn2l * p1(i + 1,j,k) + cn2s * p1(i - 1,j, &
      k) + cn3l * p1(i,j + 1,k) + cn3s * p1(i,j - 1,k) + cn4l * p1(i,j,k + 1) + cn4s * p1(i,j, &
      k - 1) - rhs(i,j,k)) - p1(i,j,k))
                        p0(i,j,k) = p1(i,j,k) + reltmp
                      end if
end subroutine press_map_97

end module singleton_module_press_map_97

