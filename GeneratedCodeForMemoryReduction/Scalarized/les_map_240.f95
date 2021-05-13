module singleton_module_les_map_240

contains

subroutine les_map_240_scal(sm,dzn,dx1,dy1,u,dzs,w,v,h)
    ! local vars: i,j,k,diu3_im1_j_kp1,diu3_i_j_kp1,diu6_i_jm1_kp1,diu6_i_j_kp1,diu7_i_j_k,diu7_ip1_j_k,diu8_i_j_k,diu8_i_jp1_k,diu9_i_j_k,evsx2,evsx1,evsy2,evsy1,evsz2,evsz1,viswx2,viswx1,viswy2,viswy1,diu9_i_j_kp1,viswz2,viswz1,vfw
    real(kind=4) :: diu3_im1_j_kp1
    real(kind=4) :: diu3_i_j_kp1
    real(kind=4) :: diu6_i_jm1_kp1
    real(kind=4) :: diu6_i_j_kp1
    real(kind=4) :: diu7_i_j_k
    real(kind=4) :: diu7_ip1_j_k
    real(kind=4) :: diu8_i_j_k
    real(kind=4) :: diu8_i_jp1_k
    real(kind=4) :: diu9_i_j_k
    real(kind=4) :: evsx2
    real(kind=4) :: evsx1
    real(kind=4) :: evsy2
    real(kind=4) :: evsy1
    real(kind=4) :: evsz2
    real(kind=4) :: evsz1
    real(kind=4) :: viswx2
    real(kind=4) :: viswx1
    real(kind=4) :: viswy2
    real(kind=4) :: viswy1
    real(kind=4) :: diu9_i_j_kp1
    real(kind=4) :: viswz2
    real(kind=4) :: viswz1
    real(kind=4) :: vfw
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension((-1):(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: sm
    real, dimension((-1):(kp+2)), intent(in) :: dzn
    real, dimension((-1):(ip+1)), intent(in) :: dx1
    real, dimension(0:(jp+1)), intent(in) :: dy1
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: u
    real, dimension((-1):(kp+2)), intent(in) :: dzs
    real, dimension(0:(ip+1),(-1):(jp+1),(-1):(kp+1)), intent(in) :: w
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: v
! WRITTEN
! READ & WRITTEN
    real, dimension(0:ip,0:jp,0:kp), intent(inout) :: h
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
      evsz2 = sm(i,j,k + 1)
      evsz1 = sm(i,j,k)
      evsx2 = (dzn(k + 1) * ((dx1(i + 1) * sm(i,j,k) + dx1(i) * sm(i + 1,j, &
      k)) / (dx1(i) + dx1(i + 1))) + dzn(k) * ((dx1(i + 1) * sm(i,j,k + 1) + dx1(i) * sm(i + 1,j, &
      k + 1)) / (dx1(i) + dx1(i + 1)))) / (dzn(k) + dzn(k + 1))
      evsx1 = (dzn(k + 1) * ((dx1(i) * sm(i - 1,j,k) + dx1(i - 1) * sm(i,j, &
      k)) / (dx1(i - 1) + dx1(i))) + dzn(k) * ((dx1(i) * sm(i - 1,j,k + 1) + dx1(i - 1) * sm(i,j, &
      k + 1)) / (dx1(i - 1) + dx1(i)))) / (dzn(k) + dzn(k + 1))
      evsy2 = (dzn(k + 1) * ((dy1(j + 1) * sm(i,j,k) + dy1(j) * sm(i,j + 1, &
      k)) / (dy1(j) + dy1(j + 1))) + dzn(k) * ((dy1(j + 1) * sm(i,j,k + 1) + dy1(j) * sm(i,j + 1, &
      k + 1)) / (dy1(j) + dy1(j + 1)))) / (dzn(k) + dzn(k + 1))
      evsy1 = (dzn(k + 1) * ((dy1(j) * sm(i,j - 1,k) + dy1(j - 1) * sm(i,j, &
      k)) / (dy1(j - 1) + dy1(j))) + dzn(k) * ((dy1(j) * sm(i,j - 1,k + 1) + dy1(j - 1) * sm(i,j, &
      k + 1)) / (dy1(j - 1) + dy1(j)))) / (dzn(k) + dzn(k + 1))
      diu3_i_j_kp1 = (-u(i,j,k) + u(i,j,k + 1)) / dzs(k)
      diu7_ip1_j_k = 2. * (-w(i,j,k) + w(i + 1,j,k)) / (dx1(i) + dx1(i + 1))
      viswx2 = (evsx2) * (diu3_i_j_kp1 + diu7_ip1_j_k)
      diu3_im1_j_kp1 = (-u(i - 1,j,k) + u(i - 1,j,k + 1)) / dzs(k)
      diu7_i_j_k = 2. * (-w(i - 1,j,k) + w(i,j,k)) / (dx1(i - 1) + dx1(i))
      viswx1 = (evsx1) * (diu3_im1_j_kp1 + diu7_i_j_k)
      diu6_i_j_kp1 = (-v(i,j,k) + v(i,j,k + 1)) / dzs(k)
      diu8_i_jp1_k = 2. * (-w(i,j,k) + w(i,j + 1,k)) / (dy1(j) + dy1(j + 1))
      viswy2 = (evsy2) * (diu6_i_j_kp1 + diu8_i_jp1_k)
      diu6_i_jm1_kp1 = (-v(i,j - 1,k) + v(i,j - 1,k + 1)) / dzs(k)
      diu8_i_j_k = 2. * (-w(i,j - 1,k) + w(i,j,k)) / (dy1(j - 1) + dy1(j))
      viswy1 = (evsy1) * (diu6_i_jm1_kp1 + diu8_i_j_k)
      diu9_i_j_kp1 = (-w(i,j,k) + w(i,j,k + 1)) / dzn(k + 1)
      viswz2 = (evsz2) * 2. * diu9_i_j_kp1
      diu9_i_j_k = (-w(i,j,k - 1) + w(i,j,k)) / dzn(k)
      viswz1 = (evsz1) * 2. * diu9_i_j_k
      vfw = (viswx2 - viswx1) / dx1(i) + (viswy2 - viswy1) / dy1(j) + (viswz2 - viswz1) / dzn(k)
      h(i,j,k) = (h(i,j,k) + vfw)
end subroutine les_map_240

end module singleton_module_les_map_240

