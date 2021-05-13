module singleton_module_les_map_202

contains

subroutine les_map_202_scal(sm,dy1,dx1,dzn,u,v,dzs,w,dys,g)
    ! local vars: i,j,k,diu2_im1_jp1_k,diu2_i_jp1_k,diu4_i_j_k,diu4_ip1_j_k,diu5_i_j_k,diu6_i_j_k,diu6_i_j_kp1,diu8_i_jp1_k,diu8_i_jp1_km1,evsx2,evsx1,evsy2,evsy1,evsz2,evsz1,visvx2,visvx1,diu5_i_jp1_k,visvy2,visvy1,visvz2,visvz1,vfv
    real(kind=4) :: diu2_im1_jp1_k
    real(kind=4) :: diu2_i_jp1_k
    real(kind=4) :: diu4_i_j_k
    real(kind=4) :: diu4_ip1_j_k
    real(kind=4) :: diu5_i_j_k
    real(kind=4) :: diu6_i_j_k
    real(kind=4) :: diu6_i_j_kp1
    real(kind=4) :: diu8_i_jp1_k
    real(kind=4) :: diu8_i_jp1_km1
    real(kind=4) :: evsx2
    real(kind=4) :: evsx1
    real(kind=4) :: evsy2
    real(kind=4) :: evsy1
    real(kind=4) :: evsz2
    real(kind=4) :: evsz1
    real(kind=4) :: visvx2
    real(kind=4) :: visvx1
    real(kind=4) :: diu5_i_jp1_k
    real(kind=4) :: visvy2
    real(kind=4) :: visvy1
    real(kind=4) :: visvz2
    real(kind=4) :: visvz1
    real(kind=4) :: vfv
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
    real, dimension(0:jp), intent(in) :: dys
! WRITTEN
! READ & WRITTEN
    real, dimension(0:ip,0:jp,0:kp), intent(inout) :: g
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
      evsy2 = sm(i,j + 1,k)
      evsy1 = sm(i,j,k)
      evsx2 = (dy1(j + 1) * ((dx1(i + 1) * sm(i,j,k) + dx1(i) * sm(i + 1,j, &
      k)) / (dx1(i) + dx1(i + 1))) + dy1(j) * ((dx1(i + 1) * sm(i,j + 1,k) + dx1(i) * sm(i + 1, &
      j + 1,k)) / (dx1(i) + dx1(i + 1)))) / (dy1(j) + dy1(j + 1))
      evsx1 = (dy1(j + 1) * ((dx1(i) * sm(i - 1,j,k) + dx1(i - 1) * sm(i,j, &
      k)) / (dx1(i - 1) + dx1(i))) + dy1(j) * ((dx1(i) * sm(i - 1,j + 1,k) + dx1(i - 1) * sm(i, &
      j + 1,k)) / (dx1(i - 1) + dx1(i)))) / (dy1(j) + dy1(j + 1))
      evsz2 = (dzn(k + 1) * ((dx1(i + 1) * sm(i,j,k) + dx1(i) * sm(i + 1,j, &
      k)) / (dx1(i) + dx1(i + 1))) + dzn(k) * ((dx1(i + 1) * sm(i,j,k + 1) + dx1(i) * sm(i + 1,j, &
      k + 1)) / (dx1(i) + dx1(i + 1)))) / (dzn(k) + dzn(k + 1))
      evsz1 = (dzn(k) * ((dx1(i + 1) * sm(i,j,k - 1) + dx1(i) * sm(i + 1,j, &
      k - 1)) / (dx1(i) + dx1(i + 1))) + dzn(k - 1) * ((dx1(i + 1) * sm(i,j, &
      k) + dx1(i) * sm(i + 1,j,k)) / (dx1(i) + dx1(i + 1)))) / (dzn(k - 1) + dzn(k))
      diu2_i_jp1_k = 2. * (-u(i,j,k) + u(i,j + 1,k)) / (dy1(j) + dy1(j + 1))
      diu4_ip1_j_k = 2. * (-v(i,j,k) + v(i + 1,j,k)) / (dx1(i) + dx1(i + 1))
      visvx2 = (evsx2) * (diu2_i_jp1_k + diu4_ip1_j_k)
      diu2_im1_jp1_k = 2. * (-u(i - 1,j,k) + u(i - 1,j + 1,k)) / (dy1(j) + dy1(j + 1))
      diu4_i_j_k = 2. * (-v(i - 1,j,k) + v(i,j,k)) / (dx1(i - 1) + dx1(i))
      visvx1 = (evsx1) * (diu2_im1_jp1_k + diu4_i_j_k)
      diu5_i_jp1_k = (-v(i,j,k) + v(i,j + 1,k)) / dy1(j + 1)
      visvy2 = (evsy2) * 2. * diu5_i_jp1_k
      diu5_i_j_k = (-v(i,j - 1,k) + v(i,j,k)) / dy1(j)
      visvy1 = (evsy1) * 2. * diu5_i_j_k
      diu6_i_j_kp1 = (-v(i,j,k) + v(i,j,k + 1)) / dzs(k)
      diu8_i_jp1_k = 2. * (-w(i,j,k) + w(i,j + 1,k)) / (dy1(j) + dy1(j + 1))
      visvz2 = (evsz2) * (diu6_i_j_kp1 + diu8_i_jp1_k)
      diu6_i_j_k = (-v(i,j,k - 1) + v(i,j,k)) / dzs(k - 1)
      diu8_i_jp1_km1 = 2. * (-w(i,j,k - 1) + w(i,j + 1,k - 1)) / (dy1(j) + dy1(j + 1))
      visvz1 = (evsz1) * (diu6_i_j_k + diu8_i_jp1_km1)
      vfv = (visvx2 - visvx1) / dx1(i) + (visvy2 - visvy1) / dys(j) + (visvz2 - visvz1) / dzn(k)
      g(i,j,k) = (g(i,j,k) + vfv)
end subroutine les_map_202

end module singleton_module_les_map_202

