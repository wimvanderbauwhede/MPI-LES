module singleton_module_velfg_map_175

contains

subroutine velfg_map_175_scal(dzn,u,w,dx1,v,dy1,h)
        parameter(ip=300)
        parameter(jp=300)
            parameter(u0=0)
    ! local vars: i,j,k,covx1,covy1,covz1,covc,nou7_,diu7_,cov7_i,nou7_ip1,diu7_ip1,cov7_ip1,nou8_,diu8_,cov8_j,nou8_jp1,diu8_jp1,cov8_jp1,nou9_,diu9_,cov9_k,nou9_kp1,diu9_kp1,cov9_kp1
    real(kind=4) :: covx1
    real(kind=4) :: covy1
    real(kind=4) :: covz1
    real(kind=4) :: covc
    real(kind=4) :: nou7_
    real(kind=4) :: diu7_
    real(kind=4) :: cov7_i
    real(kind=4) :: cov7_ip1
    real(kind=4) :: nou8_
    real(kind=4) :: diu8_
    real(kind=4) :: cov8_j
    real(kind=4) :: nou8_jp1
    real(kind=4) :: diu8_jp1
    real(kind=4) :: cov8_jp1
    real(kind=4) :: nou9_
    real(kind=4) :: diu9_
    real(kind=4) :: cov9_k
    real(kind=4) :: nou9_kp1
    real(kind=4) :: diu9_kp1
    real(kind=4) :: cov9_kp1
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension((-1):(kp+2)), intent(in) :: dzn
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: u
    real, dimension(0:(ip+1),(-1):(jp+1),(-1):(kp+1)), intent(in) :: w
    real, dimension((-1):(ip+1)), intent(in) :: dx1
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: v
    real, dimension(0:(jp+1)), intent(in) :: dy1
! WRITTEN
    real, dimension(0:ip,0:jp,0:kp), intent(out) :: h
! READ & WRITTEN
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
      nou7_ = (dzn(k + 1) * u(i - 1,j,k) + dzn(k) * u(i - 1,j,k + 1)) / (dzn(k) + dzn(k + 1))
      diu7_ = 2. * (-w(i - 1,j,k) + w(i,j,k)) / (dx1(i - 1) + dx1(i))
      cov7_i = (nou7_ - u0) * diu7_
      cov7_ip1 = (nou7_ - u0) * diu7_
      if (i == ip)        cov7_ip1 = cov7_i
      nou8_ = (dzn(k + 1) * v(i,j - 1,k) + dzn(k) * v(i,j - 1,k + 1)) / (dzn(k) + dzn(k + 1))
      diu8_ = 2. * (-w(i,j - 1,k) + w(i,j,k)) / (dy1(j - 1) + dy1(j))
      cov8_j = nou8_ * diu8_
      nou8_jp1 = (dzn(k + 1) * v(i,j,k) + dzn(k) * v(i,j,k + 1)) / (dzn(k) + dzn(k + 1))
      diu8_jp1 = 2. * (-w(i,j,k) + w(i,j + 1,k)) / (dy1(j) + dy1(j + 1))
      cov8_jp1 = nou8_jp1 * diu8_jp1
      if (j == jp) then
        nou8_ = (dzn(k + 1) * v(i,0,k) + dzn(k) * v(i,0,k + 1)) / (dzn(k) + dzn(k + 1))
        diu8_ = 2. * (-w(i,0,k) + w(i,1,k)) / (dy1(0) + dy1(1))
        cov8_jp1 = nou8_ * diu8_
      end if
      nou9_ = (w(i,j,k - 1) + w(i,j,k)) / 2.
      diu9_ = (-w(i,j,k - 1) + w(i,j,k)) / dzn(k)
      cov9_k = nou9_ * diu9_
      nou9_kp1 = (w(i,j,k) + w(i,j,k + 1)) / 2.
      diu9_kp1 = (-w(i,j,k) + w(i,j,k + 1)) / dzn(k + 1)
      cov9_kp1 = nou9_kp1 * diu9_kp1
       covx1 = (cov7_i + cov7_ip1) / 2.
       covy1 = (cov8_j + cov8_jp1) / 2.
       covz1 = (dzn(k + 1) * cov9_k + dzn(k) * cov9_kp1) / (dzn(k) + dzn(k + 1))
       covc = covx1 + covy1 + covz1
        h(i,j,k) = (-covc)
end subroutine velfg_map_175

end module singleton_module_velfg_map_175

