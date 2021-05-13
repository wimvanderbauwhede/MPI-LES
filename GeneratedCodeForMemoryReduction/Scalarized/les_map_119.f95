module singleton_module_les_map_119

contains

subroutine les_map_119_scal(u,dx1,dy1,dzs,v,w,dzn,delx1,sm)
        parameter(cs0=0.14)
    ! local vars: diu1_i_j_k,i,j,k,dudxx1,diu2_i_j_k,diu2_im1_j_k,diu2_im1_jp1_k,diu2_i_jp1_k,dudyx1,diu3_i_j_k,diu3_im1_j_k,diu3_im1_j_kp1,diu3_i_j_kp1,dudzx1,diu4_i_j_k,diu4_i_jm1_k,diu4_ip1_j_k,diu4_ip1_jm1_k,dvdxx1,diu5_i_j_k,dvdyx1,diu6_i_j_k,diu6_i_jm1_k,diu6_i_jm1_kp1,diu6_i_j_kp1,dvdzx1,diu7_i_j_k,diu7_i_j_km1,diu7_ip1_j_k,diu7_ip1_j_km1,dwdxx1,diu8_i_j_k,diu8_i_j_km1,diu8_i_jp1_k,diu8_i_jp1_km1,dwdyx1,diu9_i_j_k,dwdzx1,csx1
    real(kind=4) :: diu1_i_j_k
    real(kind=4) :: dudxx1
    real(kind=4) :: diu2_i_j_k
    real(kind=4) :: diu2_im1_j_k
    real(kind=4) :: diu2_im1_jp1_k
    real(kind=4) :: diu2_i_jp1_k
    real(kind=4) :: dudyx1
    real(kind=4) :: diu3_i_j_k
    real(kind=4) :: diu3_im1_j_k
    real(kind=4) :: diu3_im1_j_kp1
    real(kind=4) :: diu3_i_j_kp1
    real(kind=4) :: dudzx1
    real(kind=4) :: diu4_i_j_k
    real(kind=4) :: diu4_i_jm1_k
    real(kind=4) :: diu4_ip1_j_k
    real(kind=4) :: diu4_ip1_jm1_k
    real(kind=4) :: dvdxx1
    real(kind=4) :: diu5_i_j_k
    real(kind=4) :: dvdyx1
    real(kind=4) :: diu6_i_j_k
    real(kind=4) :: diu6_i_jm1_k
    real(kind=4) :: diu6_i_jm1_kp1
    real(kind=4) :: diu6_i_j_kp1
    real(kind=4) :: dvdzx1
    real(kind=4) :: diu7_i_j_k
    real(kind=4) :: diu7_i_j_km1
    real(kind=4) :: diu7_ip1_j_k
    real(kind=4) :: diu7_ip1_j_km1
    real(kind=4) :: dwdxx1
    real(kind=4) :: diu8_i_j_k
    real(kind=4) :: diu8_i_j_km1
    real(kind=4) :: diu8_i_jp1_k
    real(kind=4) :: diu8_i_jp1_km1
    real(kind=4) :: dwdyx1
    real(kind=4) :: diu9_i_j_k
    real(kind=4) :: dwdzx1
    real(kind=4) :: csx1
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: u
    real, dimension((-1):(ip+1)), intent(in) :: dx1
    real, dimension(0:(jp+1)), intent(in) :: dy1
    real, dimension((-1):(kp+2)), intent(in) :: dzs
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: v
    real, dimension(0:(ip+1),(-1):(jp+1),(-1):(kp+1)), intent(in) :: w
    real, dimension((-1):(kp+2)), intent(in) :: dzn
    real, dimension(1:kp), intent(in) :: delx1
! WRITTEN
    real, dimension((-1):(ip+1),(-1):(jp+1),0:(kp+1)), intent(out) :: sm
! READ & WRITTEN
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
        diu1_i_j_k = (-u(i - 1,j,k) + u(i,j,k)) / dx1(i)
        dudxx1 = diu1_i_j_k
        diu2_i_j_k = 2. * (-u(i,j - 1,k) + u(i,j,k)) / (dy1(j - 1) + dy1(j))
        diu2_im1_j_k = 2. * (-u(i - 1,j - 1,k) + u(i - 1,j,k)) / (dy1(j - 1) + dy1(j))
        diu2_im1_jp1_k = 2. * (-u(i - 1,j,k) + u(i - 1,j + 1,k)) / (dy1(j) + dy1(j + 1))
        diu2_i_jp1_k = 2. * (-u(i,j,k) + u(i,j + 1,k)) / (dy1(j) + dy1(j + 1))
        dudyx1 = (diu2_im1_j_k + diu2_im1_jp1_k + diu2_i_j_k + diu2_i_jp1_k) * .25
        diu3_i_j_k = (-u(i,j,k - 1) + u(i,j,k)) / dzs(k - 1)
        diu3_im1_j_k = (-u(i - 1,j,k - 1) + u(i - 1,j,k)) / dzs(k - 1)
        diu3_im1_j_kp1 = (-u(i - 1,j,k) + u(i - 1,j,k + 1)) / dzs(k)
        diu3_i_j_kp1 = (-u(i,j,k) + u(i,j,k + 1)) / dzs(k)
        dudzx1 = (diu3_im1_j_k + diu3_im1_j_kp1 + diu3_i_j_k + diu3_i_j_kp1) * .25
        diu4_i_j_k = 2. * (-v(i - 1,j,k) + v(i,j,k)) / (dx1(i - 1) + dx1(i))
        diu4_i_jm1_k = 2. * (-v(i - 1,j - 1,k) + v(i,j - 1,k)) / (dx1(i - 1) + dx1(i))
        diu4_ip1_j_k = 2. * (-v(i,j,k) + v(i + 1,j,k)) / (dx1(i) + dx1(i + 1))
        diu4_ip1_jm1_k = 2. * (-v(i,j - 1,k) + v(i + 1,j - 1,k)) / (dx1(i) + dx1(i + 1))
        dvdxx1 = (diu4_i_j_k + diu4_i_jm1_k + diu4_ip1_j_k + diu4_ip1_jm1_k) * .25
        diu5_i_j_k = (-v(i,j - 1,k) + v(i,j,k)) / dy1(j)
        dvdyx1 = diu5_i_j_k
        diu6_i_j_k = (-v(i,j,k - 1) + v(i,j,k)) / dzs(k - 1)
        diu6_i_jm1_k = (-v(i,j - 1,k - 1) + v(i,j - 1,k)) / dzs(k - 1)
        diu6_i_jm1_kp1 = (-v(i,j - 1,k) + v(i,j - 1,k + 1)) / dzs(k)
        diu6_i_j_kp1 = (-v(i,j,k) + v(i,j,k + 1)) / dzs(k)
        dvdzx1 = (diu6_i_jm1_k + diu6_i_jm1_kp1 + diu6_i_j_k + diu6_i_j_kp1) * .25
        diu7_i_j_k = 2. * (-w(i - 1,j,k) + w(i,j,k)) / (dx1(i - 1) + dx1(i))
        diu7_i_j_km1 = 2. * (-w(i - 1,j,k - 1) + w(i,j,k - 1)) / (dx1(i - 1) + dx1(i))
        diu7_ip1_j_k = 2. * (-w(i,j,k) + w(i + 1,j,k)) / (dx1(i) + dx1(i + 1))
        diu7_ip1_j_km1 = 2. * (-w(1,j,k - 1) + w(i + 1,j,k - 1)) / (dx1(i) + dx1(i + 1))
        dwdxx1 = (diu7_i_j_k + diu7_i_j_km1 + diu7_ip1_j_k + diu7_ip1_j_km1) * .25
        diu8_i_j_k = 2. * (-w(i,j - 1,k) + w(i,j,k)) / (dy1(j - 1) + dy1(j))
        diu8_i_j_km1 = 2. * (-w(i,j - 1,k - 1) + w(i,j,k - 1)) / (dy1(j - 1) + dy1(j))
        diu8_i_jp1_k = 2. * (-w(i,j,k) + w(i,j + 1,k)) / (dy1(j) + dy1(j + 1))
        diu8_i_jp1_km1 = 2. * (-w(i,j,k - 1) + w(i,j + 1,k - 1)) / (dy1(j) + dy1(j + 1))
        dwdyx1 = (diu8_i_j_k + diu8_i_j_km1 + diu8_i_jp1_k + diu8_i_jp1_km1) * .25
        diu9_i_j_k = (-w(i,j,k - 1) + w(i,j,k)) / dzn(k)
        dwdzx1 = diu9_i_j_k
      csx1 = cs0
      sm(i,j, &
      k) = (csx1 * delx1(k)) ** 2 * sqrt(2. * (dudxx1 ** 2 + dvdyx1 ** 2 + dwdzx1 ** 2) + (dudyx1  &
      + dvdxx1) ** 2 + (dwdyx1 + dvdzx1) ** 2 + (dudzx1 + dwdxx1) ** 2)
end subroutine les_map_119

end module singleton_module_les_map_119

