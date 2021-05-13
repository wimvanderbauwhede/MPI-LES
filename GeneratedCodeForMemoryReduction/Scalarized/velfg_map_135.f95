module singleton_module_velfg_map_135

contains

subroutine velfg_map_135_scal(dy1,u,v,dx1,w,dzs,dzn,g)
        parameter(ip=300)
        parameter(jp=300)
            parameter(u0=0)
    ! local vars: i,j,k,covx1,covy1,covz1,covc,nou4_,diu4_,cov4_i,nou4_ip1,diu4_ip1,cov4_ip1,nou5_,diu5_,cov5_j,nou5_jp1,diu5_jp1,cov5_jp1,nou6_,diu6_,cov6_k,nou6_kp1,diu6_kp1,cov6_kp1
    real(kind=4) :: covx1
    real(kind=4) :: covy1
    real(kind=4) :: covz1
    real(kind=4) :: covc
    real(kind=4) :: nou4_
    real(kind=4) :: diu4_
    real(kind=4) :: cov4_i
    real(kind=4) :: nou4_ip1
    real(kind=4) :: diu4_ip1
    real(kind=4) :: cov4_ip1
    real(kind=4) :: nou5_
    real(kind=4) :: diu5_
    real(kind=4) :: cov5_j
    real(kind=4) :: nou5_jp1
    real(kind=4) :: diu5_jp1
    real(kind=4) :: cov5_jp1
    real(kind=4) :: nou6_
    real(kind=4) :: diu6_
    real(kind=4) :: cov6_k
    real(kind=4) :: nou6_kp1
    real(kind=4) :: diu6_kp1
    real(kind=4) :: cov6_kp1
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension(0:(jp+1)), intent(in) :: dy1
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: u
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: v
    real, dimension((-1):(ip+1)), intent(in) :: dx1
    real, dimension(0:(ip+1),(-1):(jp+1),(-1):(kp+1)), intent(in) :: w
    real, dimension((-1):(kp+2)), intent(in) :: dzs
    real, dimension((-1):(kp+2)), intent(in) :: dzn
! WRITTEN
    real, dimension(0:ip,0:jp,0:kp), intent(out) :: g
! READ & WRITTEN
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
      nou4_ = (dy1(j + 1) * u(i - 1,j,k) + dy1(j) * u(i - 1,j + 1,k)) / (dy1(j) + dy1(j + 1))
      diu4_ = 2. * (-v(i - 1,j,k) + v(i,j,k)) / (dx1(i - 1) + dx1(i))
      cov4_i = (nou4_ - u0) * diu4_
      nou4_ip1 = (dy1(j + 1) * u(i,j,k) + dy1(j) * u(i,j + 1,k)) / (dy1(j) + dy1(j + 1))
      diu4_ip1 = 2. * (-v(i,j,k) + v(i + 1,j,k)) / (dx1(i) + dx1(i + 1))
      cov4_ip1 = (nou4_ip1 - u0) * diu4_ip1
      if (i == ip)        cov4_ip1 = cov4_i
      nou5_ = (v(i,j - 1,k) + v(i,j,k)) / 2.
      diu5_ = (-v(i,j - 1,k) + v(i,j,k)) / dy1(j)
      cov5_j = nou5_ * diu5_
      nou5_jp1 = (v(i,j,k) + v(i,j + 1,k)) / 2.
      diu5_jp1 = (-v(i,j,k) + v(i,j + 1,k)) / dy1(j + 1)
      cov5_jp1 = nou5_jp1 * diu5_jp1
      if (j == jp) then
          nou5_ = (v(i,1,k) + v(i,2,k)) / 2.
          diu5_ = (-v(i,1,k) + v(i,2,k)) / dy1(j)
          cov5_jp1 = nou5_ * diu5_
      end if  
        nou6_ = (dy1(j + 1) * w(i,j,k - 1) + dy1(j) * w(i,j + 1,k - 1)) / (dy1(j) + dy1(j + 1))
        diu6_ = (-v(i,j,k - 1) + v(i,j,k)) / dzs(k - 1)
        cov6_k = nou6_ * diu6_
        nou6_kp1 = (dy1(j + 1) * w(i,j,k) + dy1(j) * w(i,j + 1,k)) / (dy1(j) + dy1(j + 1))
        diu6_kp1 = (-v(i,j,k) + v(i,j,k + 1)) / dzs(k)
        cov6_kp1 = nou6_kp1 * diu6_kp1
        if (k == 1) then
            nou6_ = 0.5 * (dy1(j + 1) * w(i,j,1) + dy1(j) * w(i,j + 1,1)) / (dy1(j) + dy1(j + 1))
            diu6_ = 0.4 * v(i,j,1) / (alog(5.0 * dzn(1)) * 0.2 * dzn(1))
            cov6_k = nou6_ * diu6_
        end if
        covx1 = (cov4_i + cov4_ip1) / 2.
        covy1 = (dy1(j + 1) * cov5_j + dy1(j) * cov5_jp1) / (dy1(j) + dy1(j + 1))
        covz1 = (cov6_k + cov6_kp1) / 2.
        covc = covx1 + covy1 + covz1
        g(i,j,k) = (-covc)
end subroutine velfg_map_135

end module singleton_module_velfg_map_135

