module singleton_module_les_map_164

contains

subroutine les_map_164_scal(sm,dy1,dx1,dzn,u,v,dzs,w,dxs,f)
    ! local vars: diu1_i_j_k,i,j,k,diu2_i_j_k,diu2_i_jp1_k,diu3_i_j_k,diu3_i_j_kp1,diu4_ip1_j_k,diu4_ip1_jm1_k,diu7_ip1_j_k,diu7_ip1_j_km1,evsx2,evsx1,evsy2,evsy1,evsz2,evsz1,diu1_ip1_j_k,visux2,visux1,visuy2,visuy1,visuz2,visuz1,vfu
    real(kind=4) :: diu1_i_j_k
    real(kind=4) :: diu2_i_j_k
    real(kind=4) :: diu2_i_jp1_k
    real(kind=4) :: diu3_i_j_k
    real(kind=4) :: diu3_i_j_kp1
    real(kind=4) :: diu4_ip1_j_k
    real(kind=4) :: diu4_ip1_jm1_k
    real(kind=4) :: diu7_ip1_j_k
    real(kind=4) :: diu7_ip1_j_km1
    real(kind=4) :: evsx2
    real(kind=4) :: evsx1
    real(kind=4) :: evsy2
    real(kind=4) :: evsy1
    real(kind=4) :: evsz2
    real(kind=4) :: evsz1
    real(kind=4) :: diu1_ip1_j_k
    real(kind=4) :: visux2
    real(kind=4) :: visux1
    real(kind=4) :: visuy2
    real(kind=4) :: visuy1
    real(kind=4) :: visuz2
    real(kind=4) :: visuz1
    real(kind=4) :: vfu
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension((-1):(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: sm
    real, dimension(0:(jp+1)), intent(in) :: dy1
    real, dimension((-1):(ip+1)), intent(in) :: dx1
    real, dimension((-1):(kp+2)), intent(in) :: dzn
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: u
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: v
    real, dimension((-1):(kp+2)), intent(in) :: dzs
    real, dimension(0:(ip+1),(-1):(jp+1),(-1):(kp+1)), intent(in) :: w
    real, dimension(0:ip), intent(in) :: dxs
! WRITTEN
! READ & WRITTEN
    real, dimension(0:ip,0:jp,0:kp), intent(inout) :: f
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
      evsx2 = sm(i + 1,j,k)
      evsx1 = sm(i,j,k)
      evsy2 = (dy1(j + 1) * ((dx1(i + 1) * sm(i,j,k) + dx1(i) * sm(i + 1,j, &
      k)) / (dx1(i) + dx1(i + 1))) + dy1(j) * ((dx1(i + 1) * sm(i,j + 1,k) + dx1(i) * sm(i + 1, &
      j + 1,k)) / (dx1(i) + dx1(i + 1)))) / (dy1(j) + dy1(j + 1))
      evsy1 = (dy1(j + 1) * ((dx1(i + 1) * sm(i,j - 1,k) + dx1(i) * sm(i + 1,j - 1, &
      k)) / (dx1(i) + dx1(i + 1))) + dy1(j) * ((dx1(i + 1) * sm(i,j,k) + dx1(i) * sm(i + 1,j, &
      k)) / (dx1(i) + dx1(i + 1)))) / (dy1(j) + dy1(j + 1))
      evsz2 = (dzn(k + 1) * ((dx1(i + 1) * sm(i,j,k) + dx1(i) * sm(i + 1,j, &
      k)) / (dx1(i) + dx1(i + 1))) + dzn(k) * ((dx1(i + 1) * sm(i,j,k + 1) + dx1(i) * sm(i + 1,j, &
      k + 1)) / (dx1(i) + dx1(i + 1)))) / (dzn(k) + dzn(k + 1))
      evsz1 = (dzn(k) * ((dx1(i + 1) * sm(i,j,k - 1) + dx1(i) * sm(i + 1,j, &
      k - 1)) / (dx1(i) + dx1(i + 1))) + dzn(k - 1) * ((dx1(i + 1) * sm(i,j, &
      k) + dx1(i) * sm(i + 1,j,k)) / (dx1(i) + dx1(i + 1)))) / (dzn(k - 1) + dzn(k))
      diu1_i_j_k = (-u(i - 1,j,k) + u(i,j,k)) / dx1(i)
      diu1_ip1_j_k = (-u(i,j,k) + u(i + 1,j,k)) / dx1(i + 1)
      visux2 = (evsx2) * 2. * diu1_ip1_j_k
      visux1 = (evsx1) * 2. * diu1_i_j_k
      diu2_i_jp1_k = 2. * (-u(i,j,k) + u(i,j + 1,k)) / (dy1(j) + dy1(j + 1))
      diu4_ip1_j_k = 2. * (-v(i,j,k) + v(i + 1,j,k)) / (dx1(i) + dx1(i + 1))
      visuy2 = (evsy2) * (diu2_i_jp1_k + diu4_ip1_j_k)
      diu2_i_j_k = 2. * (-u(i,j - 1,k) + u(i,j,k)) / (dy1(j - 1) + dy1(j))
      diu4_ip1_jm1_k = 2. * (-v(i,j - 1,k) + v(i + 1,j - 1,k)) / (dx1(i) + dx1(i + 1))
      visuy1 = (evsy1) * (diu2_i_j_k + diu4_ip1_jm1_k)
      diu3_i_j_kp1 = (-u(i,j,k) + u(i,j,k + 1)) / dzs(k)
      diu7_ip1_j_k = 2. * (-w(i,j,k) + w(i + 1,j,k)) / (dx1(i) + dx1(i + 1))
      visuz2 = (evsz2) * (diu3_i_j_kp1 + diu7_ip1_j_k)
      diu3_i_j_k = (-u(i,j,k - 1) + u(i,j,k)) / dzs(k - 1)
      diu7_ip1_j_km1 = 2. * (-w(1,j,k - 1) + w(i + 1,j,k - 1)) / (dx1(i) + dx1(i + 1))
      visuz1 = (evsz1) * (diu3_i_j_k + diu7_ip1_j_km1)
      vfu = (visux2 - visux1) / dxs(i) + (visuy2 - visuy1) / dy1(j) + (visuz2 - visuz1) / dzn(k)
      f(i,j,k) = (f(i,j,k) + vfu)
end subroutine les_map_164

end module singleton_module_les_map_164

