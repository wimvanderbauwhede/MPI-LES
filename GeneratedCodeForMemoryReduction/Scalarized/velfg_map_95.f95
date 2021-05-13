module singleton_module_velfg_map_95

contains

subroutine velfg_map_95_scal(u,dx1,v,dy1,w,dzs,dzn,f)
        parameter(ip=300)
        parameter(jp=300)
    ! local vars: nou1_,i,j,k,diu1_,cov1_i,nou1_ip1,diu1_ip1,cov1_ip1,nou2_,diu2_,cov2_j,nou2_jp1,diu2_jp1,cov2_jp1,nou3_,diu3_,cov3_k,nou3_kp1,diu3_kp1,cov3_kp1,covx1,covy1,covz1,covc
    real(kind=4) :: nou1_
    real(kind=4) :: diu1_
    real(kind=4) :: cov1_i
    real(kind=4) :: nou1_ip1
    real(kind=4) :: diu1_ip1
    real(kind=4) :: cov1_ip1
    real(kind=4) :: nou2_
    real(kind=4) :: diu2_
    real(kind=4) :: cov2_j
    real(kind=4) :: nou2_jp1
    real(kind=4) :: diu2_jp1
    real(kind=4) :: cov2_jp1
    real(kind=4) :: nou3_
    real(kind=4) :: diu3_
    real(kind=4) :: cov3_k
    real(kind=4) :: nou3_kp1
    real(kind=4) :: diu3_kp1
    real(kind=4) :: cov3_kp1
    real(kind=4) :: covx1
    real(kind=4) :: covy1
    real(kind=4) :: covz1
    real(kind=4) :: covc
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: u
    real, dimension((-1):(ip+1)), intent(in) :: dx1
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: v
    real, dimension(0:(jp+1)), intent(in) :: dy1
    real, dimension(0:(ip+1),(-1):(jp+1),(-1):(kp+1)), intent(in) :: w
    real, dimension((-1):(kp+2)), intent(in) :: dzs
    real, dimension((-1):(kp+2)), intent(in) :: dzn
! WRITTEN
    real, dimension(0:ip,0:jp,0:kp), intent(out) :: f
! READ & WRITTEN
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
      nou1_ = (u(i - 1,j,k) + u(i,j,k)) / 2.
      diu1_ = (-u(i - 1,j,k) + u(i,j,k)) / dx1(i)
      cov1_i = nou1_ * diu1_
      nou1_ip1 = (u(i,j,k) + u(i + 1,j,k)) / 2.
      diu1_ip1 = (-u(i,j,k) + u(i + 1,j,k)) / dx1(i + 1)
      cov1_ip1 = nou1_ip1 * diu1_ip1
      if (i == ip)        cov1_ip1 = cov1_i
      nou2_ = (dx1(i + 1) * v(i,j - 1,k) + dx1(i) * v(i + 1,j - 1,k)) / (dx1(i) + dx1(i + 1))
      diu2_ = 2. * (-u(i,j - 1,k) + u(i,j,k)) / (dy1(j - 1) + dy1(j))
      cov2_j = nou2_ * diu2_
      nou2_jp1 = (dx1(i + 1) * v(i,j,k) + dx1(i) * v(i + 1,j,k)) / (dx1(i) + dx1(i + 1))
      diu2_jp1 = 2. * (-u(i,j,k) + u(i,j + 1,k)) / (dy1(j) + dy1(j + 1))
      cov2_jp1 = nou2_jp1 * diu2_jp1
      if (j == jp) then
          nou2_ = (dx1(i + 1) * v(i,0,k) + dx1(i) * v(i + 1,0,k)) / (dx1(i) + dx1(i + 1))
          diu2_ = 2. * (-u(i,0,k) + u(i,1,k)) / (dy1(0) + dy1(1))
          cov2_jp1 = nou2_ * diu2_
      end if
        nou3_ = (dx1(i + 1) * w(i,j,k - 1) + dx1(i) * w(i + 1,j,k - 1)) / (dx1(i) + dx1(i + 1))
        diu3_ = (-u(i,j,k - 1) + u(i,j,k)) / dzs(k - 1)
        cov3_k = nou3_ * diu3_
        nou3_kp1 = (dx1(i + 1) * w(i,j,k) + dx1(i) * w(i + 1,j,k)) / (dx1(i) + dx1(i + 1))
        diu3_kp1 = (-u(i,j,k) + u(i,j,k + 1)) / dzs(k)
        cov3_kp1 = nou3_kp1 * diu3_kp1
        if (k == 1) then
            nou3_ = 0.5 * (dx1(i + 1) * w(i,j,1) + dx1(i) * w(i + 1,j,1)) / (dx1(i) + dx1(i + 1))
            diu3_ = 0.4 * u(i,j,1) / (alog(5.0 * dzn(1)) * 0.2 * dzn(1))
            cov3_k = nou3_ * diu3_
        end if
        covx1 = (dx1(i + 1) * cov1_i + dx1(i) * cov1_ip1) / (dx1(i) + dx1(i + 1))
        covy1 = (cov2_j + cov2_jp1) / 2.
        covz1 = (cov3_k + cov3_kp1) / 2.
        covc = covx1 + covy1 + covz1
        f(i,j,k) = (-covc)
end subroutine velfg_map_95

end module singleton_module_velfg_map_95

