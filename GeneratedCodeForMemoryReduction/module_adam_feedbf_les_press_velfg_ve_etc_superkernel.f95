module module_adam_feedbf_les_press_velfg_ve_etc_superkernel


    contains


subroutine adam_map_36(f,g,h,fold,gold,hold)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
    ! Local vars: fd,i,j,k,gd,hd
    real(kind=4) :: fd
    integer :: i
    integer :: j
    integer :: k
    real(kind=4) :: gd
    real(kind=4) :: hd
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: f
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: g
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: h
    real(kind=4), dimension(ip,jp,kp), intent(InOut) :: fold
    real(kind=4), dimension(ip,jp,kp), intent(InOut) :: gold
    real(kind=4), dimension(ip,jp,kp), intent(InOut) :: hold
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
                fd = f(i,j,k)
                gd = g(i,j,k)
                hd = h(i,j,k)
                f(i,j,k) = 1.5*f(i,j,k)-0.5*fold(i,j,k)
                g(i,j,k) = 1.5*g(i,j,k)-0.5*gold(i,j,k)
                h(i,j,k) = 1.5*h(i,j,k)-0.5*hold(i,j,k)
                fold(i,j,k) = fd
                gold(i,j,k) = gd
                hold(i,j,k) = hd

end subroutine adam_map_36


subroutine feedbf_map_49(usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,fx,fy,fz)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
    ! Local vars: i,j,k,f1x,f1y,f1z,f2x,f2y,f2z
    integer :: i
    integer :: j
    integer :: k
    real(kind=4) :: f1x
    real(kind=4) :: f1y
    real(kind=4) :: f1z
    real(kind=4) :: f2x
    real(kind=4) :: f2y
    real(kind=4) :: f2z
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: u
    real(kind=4), dimension((-1):(ip + 1),0:(jp + 1),0:(kp + 1)), intent(In) :: bmask1
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: v
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: cmask1
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(In) :: w
    real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)), intent(In) :: dmask1
    real(kind=4) :: alpha
    real(kind=4) :: dt
    real(kind=4) :: beta
! WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(Out) :: fx
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(Out) :: fy
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(Out) :: fz
! READ & WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: usum
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: vsum
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: wsum
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
                usum(i,j,k) = (usum(i,j,k)+u(i,j,k))*bmask1(i,j,k)
                vsum(i,j,k) = (vsum(i,j,k)+v(i,j,k))*cmask1(i,j,k)
                wsum(i,j,k) = (wsum(i,j,k)+w(i,j,k))*dmask1(i,j,k)
                f1x = alpha*usum(i,j,k)*dt
                f1y = alpha*vsum(i,j,k)*dt
                f1z = alpha*wsum(i,j,k)*dt
                f2x = beta*u(i,j,k)*bmask1(i,j,k)
                f2y = beta*v(i,j,k)*cmask1(i,j,k)
                f2z = beta*w(i,j,k)*dmask1(i,j,k)
                fx(i,j,k) = f1x+f2x
                fy(i,j,k) = f1y+f2y
                fz(i,j,k) = f1z+f2z

end subroutine feedbf_map_49


subroutine feedbf_map_67(f,fx,g,fy,h,fz)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
    ! Local vars: i,j,k
    integer :: i
    integer :: j
    integer :: k
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In) :: fx
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In) :: fy
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In) :: fz
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: f
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: g
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: h
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
                f(i,j,k) = f(i,j,k)+fx(i,j,k)
                g(i,j,k) = g(i,j,k)+fy(i,j,k)
                h(i,j,k) = h(i,j,k)+fz(i,j,k)

end subroutine feedbf_map_67


subroutine les_map_119(u,dx1,dy1,dzs,v,w,dzn,delx1,sm)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
    ! Local vars: diu1_i_j_k,i,j,k,dudxx1,diu2_i_j_k,diu2_im1_j_k,diu2_im1_jp1_k,diu2_i_jp1_k,dudyx1,diu3_i_j_k,diu3_im1_j_k,diu3_im1_j_kp1,diu3_i_j_kp1,dudzx1,diu4_i_j_k,diu4_i_jm1_k,diu4_ip1_j_k,diu4_ip1_jm1_k,dvdxx1,diu5_i_j_k,dvdyx1,diu6_i_j_k,diu6_i_jm1_k,diu6_i_jm1_kp1,diu6_i_j_kp1,dvdzx1,diu7_i_j_k,diu7_i_j_km1,diu7_ip1_j_k,diu7_ip1_j_km1,dwdxx1,diu8_i_j_k,diu8_i_j_km1,diu8_i_jp1_k,diu8_i_jp1_km1,dwdyx1,diu9_i_j_k,dwdzx1,csx1
    real(kind=4) :: diu1_i_j_k
    integer :: i
    integer :: j
    integer :: k
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
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: u
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzs
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: v
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(In) :: w
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
    real(kind=4), dimension(kp), intent(In) :: delx1
! WRITTEN
    real(kind=4), dimension((-1):(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(Out) :: sm
! READ & WRITTEN
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
        diu1_i_j_k = (-u(i-1,j,k)+u(i,j,k))/dx1(i) 
        dudxx1 = diu1_i_j_k
        diu2_i_j_k = 2.*(-u(i,j-1,k)+u(i,j,k))/(dy1(j-1)+dy1(j))
        diu2_im1_j_k = 2.*(-u(i-1,j-1,k)+u(i-1,j,k))/(dy1(j-1)+dy1(j))
        diu2_im1_jp1_k = 2.*(-u(i-1,j,k)+u(i-1,j+1,k))/(dy1(j)+dy1(j+1))
        diu2_i_jp1_k = 2.*(-u(i,j,k)+u(i,j+1,k))/(dy1(j)+dy1(j+1))
        dudyx1 =  (diu2_im1_j_k+diu2_im1_jp1_k +diu2_i_j_k+diu2_i_jp1_k) *.25
        diu3_i_j_k = (-u(i,j,k-1)+u(i,j,k))/dzs(k-1)
        diu3_im1_j_k = (-u(i-1,j,k-1)+u(i-1,j,k))/dzs(k-1)
        diu3_im1_j_kp1 = (-u(i-1,j,k)+u(i-1,j,k+1))/dzs(k)
        diu3_i_j_kp1 = (-u(i,j,k)+u(i,j,k+1))/dzs(k)
        dudzx1 =  (diu3_im1_j_k+diu3_im1_j_kp1 +diu3_i_j_k+diu3_i_j_kp1) *.25
        diu4_i_j_k = 2.*(-v(i-1,j,k)+v(i,j,k))/(dx1(i-1)+dx1(i))
        diu4_i_jm1_k = 2.*(-v(i-1,j-1,k)+v(i,j-1,k))/(dx1(i-1)+dx1(i))
        diu4_ip1_j_k = 2.*(-v(i,j,k)+v(i+1,j,k))/(dx1(i)+dx1(i+1))
        diu4_ip1_jm1_k = 2.*(-v(i,j-1,k)+v(i+1,j-1,k))/(dx1(i)+dx1(i+1))
        dvdxx1 =  (diu4_i_j_k+diu4_i_jm1_k +diu4_ip1_j_k+diu4_ip1_jm1_k ) *.25
        diu5_i_j_k = (-v(i,j-1,k)+v(i,j,k))/dy1(j)
        dvdyx1 =  diu5_i_j_k
        diu6_i_j_k = (-v(i,j,k-1)+v(i,j,k))/dzs(k-1)
        diu6_i_jm1_k = (-v(i,j-1,k-1)+v(i,j-1,k))/dzs(k-1)
        diu6_i_jm1_kp1 = (-v(i,j-1,k)+v(i,j-1,k+1))/dzs(k)
        diu6_i_j_kp1 = (-v(i,j,k)+v(i,j,k+1))/dzs(k)
        dvdzx1 =  (diu6_i_jm1_k+diu6_i_jm1_kp1 +diu6_i_j_k+diu6_i_j_kp1 )  *.25
        diu7_i_j_k = 2.*(-w(i-1,j,k)+w(i,j,k))/(dx1(i-1)+dx1(i))
        diu7_i_j_km1 = 2.*(-w(i-1,j,k-1)+w(i,j,k-1))/(dx1(i-1)+dx1(i))
        diu7_ip1_j_k = 2.*(-w(i,j,k)+w(i+1,j,k))/(dx1(i)+dx1(i+1))
        diu7_ip1_j_km1 = 2.*(-w(1,j,k-1)+w(i+1,j,k-1))/(dx1(i)+dx1(i+1))
        dwdxx1 =  (diu7_i_j_k+diu7_i_j_km1 +diu7_ip1_j_k+diu7_ip1_j_km1 ) *.25
        diu8_i_j_k = 2.*(-w(i,j-1,k)+w(i,j,k))/(dy1(j-1)+dy1(j))
        diu8_i_j_km1 = 2.*(-w(i,j-1,k-1)+w(i,j,k-1))/(dy1(j-1)+dy1(j))
        diu8_i_jp1_k = 2.*(-w(i,j,k)+w(i,j+1,k))/(dy1(j)+dy1(j+1))
        diu8_i_jp1_km1 = 2.*(-w(i,j,k-1)+w(i,j+1,k-1))/(dy1(j)+dy1(j+1))
        dwdyx1 =  (diu8_i_j_k+diu8_i_j_km1 +diu8_i_jp1_k+diu8_i_jp1_km1 ) *.25
        diu9_i_j_k = (-w(i,j,k-1)+w(i,j,k))/dzn(k)
        dwdzx1 =  diu9_i_j_k
      csx1 = cs0
      sm(i,j,k) = ( csx1*delx1(k) )**2  * sqrt( 2.*( dudxx1**2+dvdyx1**2+dwdzx1**2 ) +( dudyx1+dvdxx1 )**2  &
      +( dwdyx1+dvdzx1 )**2 +( dudzx1+dwdxx1 )**2 )

end subroutine les_map_119


subroutine les_map_164(sm,dy1,dx1,dzn,u,v,dzs,w,dxs,f)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
    ! Local vars: diu1_i_j_k,i,j,k,diu2_i_j_k,diu2_i_jp1_k,diu3_i_j_k,diu3_i_j_kp1,diu4_ip1_j_k,diu4_ip1_jm1_k,diu7_ip1_j_k,diu7_ip1_j_km1,evsx2,evsx1,evsy2,evsy1,evsz2,evsz1,diu1_ip1_j_k,visux2,visux1,visuy2,visuy1,visuz2,visuz1,vfu
    real(kind=4) :: diu1_i_j_k
    integer :: i
    integer :: j
    integer :: k
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
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension((-1):(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: sm
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: u
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: v
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzs
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(In) :: w
    real(kind=4), dimension(0:ip), intent(In) :: dxs
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: f
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 2) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 2)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
      evsx2 = sm(i+1,j,k)
      evsx1 = sm(i,j,k)
      evsy2 = (dy1(j+1)*((dx1(i+1)*sm(i,j,k)+dx1(i)*sm(i+1,j, &
      k)) /(dx1(i)+dx1(i+1))) +dy1(j)*((dx1(i+1)*sm(i,j+1,k)+dx1(i)*sm(i+1,j+1, &
      k)) /(dx1(i)+dx1(i+1)))) /(dy1(j)+dy1(j+1))
      evsy1 = (dy1(j+1)*((dx1(i+1)*sm(i,j-1,k)+dx1(i)*sm(i+1,j-1, &
      k)) /(dx1(i)+dx1(i+1))) +dy1(j)*((dx1(i+1)*sm(i,j,k)+dx1(i)*sm(i+1,j, &
      k)) /(dx1(i)+dx1(i+1)))) /(dy1(j)+dy1(j+1))
      evsz2 = (dzn(k+1)*((dx1(i+1)*sm(i,j,k)+dx1(i)*sm(i+1,j, &
      k)) /(dx1(i)+dx1(i+1))) +dzn(k)*((dx1(i+1)*sm(i,j,k+1)+dx1(i)*sm(i+1,j, &
      k+1)) /(dx1(i)+dx1(i+1)))) /(dzn(k)+dzn(k+1))
      evsz1 = (dzn(k)*((dx1(i+1)*sm(i,j,k-1)+dx1(i)*sm(i+1,j, &
      k-1)) /(dx1(i)+dx1(i+1))) +dzn(k-1)*((dx1(i+1)*sm(i,j,k)+dx1(i)*sm(i+1,j, &
      k)) /(dx1(i)+dx1(i+1)))) /(dzn(k-1)+dzn(k))
      diu1_i_j_k = (-u(i-1,j,k)+u(i,j,k))/dx1(i) 
      diu1_ip1_j_k = (-u(i,j,k)+u(i+1,j,k))/dx1(i+1) 
      visux2 = (evsx2)*2.*diu1_ip1_j_k
      visux1 = (evsx1)*2.*diu1_i_j_k
      diu2_i_jp1_k = 2.*(-u(i,j,k)+u(i,j+1,k))/(dy1(j)+dy1(j+1))
      diu4_ip1_j_k = 2.*(-v(i,j,k)+v(i+1,j,k))/(dx1(i)+dx1(i+1))
      visuy2 = (evsy2)* ( diu2_i_jp1_k + diu4_ip1_j_k  )
      diu2_i_j_k = 2.*(-u(i,j-1,k)+u(i,j,k))/(dy1(j-1)+dy1(j))
      diu4_ip1_jm1_k = 2.*(-v(i,j-1,k)+v(i+1,j-1,k))/(dx1(i)+dx1(i+1))
      visuy1 = (evsy1)* ( diu2_i_j_k  + diu4_ip1_jm1_k  )
      diu3_i_j_kp1 = (-u(i,j,k)+u(i,j,k+1))/dzs(k)
      diu7_ip1_j_k = 2.*(-w(i,j,k)+w(i+1,j,k))/(dx1(i)+dx1(i+1))
      visuz2 = (evsz2)* ( diu3_i_j_kp1 + diu7_ip1_j_k  )
      diu3_i_j_k = (-u(i,j,k-1)+u(i,j,k))/dzs(k-1)
      diu7_ip1_j_km1 = 2.*(-w(1,j,k-1)+w(i+1,j,k-1))/(dx1(i)+dx1(i+1))
      visuz1 = (evsz1)* ( diu3_i_j_k + diu7_ip1_j_km1)
      vfu = (visux2-visux1)/dxs(i) +(visuy2-visuy1)/dy1(j) +(visuz2-visuz1)/dzn(k)
      f(i,j,k) = (f(i,j,k)+vfu)

end subroutine les_map_164


subroutine les_map_202(sm,dy1,dx1,dzn,u,v,dzs,w,dys,g)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
    ! Local vars: i,j,k,diu2_im1_jp1_k,diu2_i_jp1_k,diu4_i_j_k,diu4_ip1_j_k,diu5_i_j_k,diu6_i_j_k,diu6_i_j_kp1,diu8_i_jp1_k,diu8_i_jp1_km1,evsx2,evsx1,evsy2,evsy1,evsz2,evsz1,visvx2,visvx1,diu5_i_jp1_k,visvy2,visvy1,visvz2,visvz1,vfv
    integer :: i
    integer :: j
    integer :: k
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
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension((-1):(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: sm
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: u
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: v
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzs
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(In) :: w
    real(kind=4), dimension(0:jp), intent(In) :: dys
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: g
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 2) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 2)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
      evsy2 = sm(i,j+1,k)
      evsy1 = sm(i,j,k)
      evsx2 = (dy1(j+1)*((dx1(i+1)*sm(i,j,k)+dx1(i)*sm(i+1,j, &
      k)) /(dx1(i)+dx1(i+1))) +dy1(j)*((dx1(i+1)*sm(i,j+1,k)+dx1(i)*sm(i+1,j+1, &
      k)) /(dx1(i)+dx1(i+1)))) /(dy1(j)+dy1(j+1))
      evsx1 = (dy1(j+1)*((dx1(i)*sm(i-1,j,k)+dx1(i-1)*sm(i,j, &
      k)) /(dx1(i-1)+dx1(i))) +dy1(j)*((dx1(i)*sm(i-1,j+1,k)+dx1(i-1)*sm(i,j+1, &
      k)) /(dx1(i-1)+dx1(i)))) /(dy1(j)+dy1(j+1))
      evsz2 = (dzn(k+1)*((dx1(i+1)*sm(i,j,k)+dx1(i)*sm(i+1,j, &
      k)) /(dx1(i)+dx1(i+1))) +dzn(k)*((dx1(i+1)*sm(i,j,k+1)+dx1(i)*sm(i+1,j, &
      k+1)) /(dx1(i)+dx1(i+1)))) /(dzn(k)+dzn(k+1))
      evsz1 = (dzn(k)*((dx1(i+1)*sm(i,j,k-1)+dx1(i)*sm(i+1,j, &
      k-1)) /(dx1(i)+dx1(i+1))) +dzn(k-1)*((dx1(i+1)*sm(i,j,k)+dx1(i)*sm(i+1,j, &
      k)) /(dx1(i)+dx1(i+1)))) /(dzn(k-1)+dzn(k))
      diu2_i_jp1_k = 2.*(-u(i,j,k)+u(i,j+1,k))/(dy1(j)+dy1(j+1))
      diu4_ip1_j_k = 2.*(-v(i,j,k)+v(i+1,j,k))/(dx1(i)+dx1(i+1))
      visvx2 = (evsx2)* ( diu2_i_jp1_k+diu4_ip1_j_k )
      diu2_im1_jp1_k = 2.*(-u(i-1,j,k)+u(i-1,j+1,k))/(dy1(j)+dy1(j+1))
      diu4_i_j_k = 2.*(-v(i-1,j,k)+v(i,j,k))/(dx1(i-1)+dx1(i))
      visvx1 = (evsx1)* ( diu2_im1_jp1_k+diu4_i_j_k )
      diu5_i_jp1_k = (-v(i,j,k)+v(i,j+1,k))/dy1(j+1)
      visvy2 = (evsy2)*2.*diu5_i_jp1_k
      diu5_i_j_k = (-v(i,j-1,k)+v(i,j,k))/dy1(j)
      visvy1 = (evsy1)*2.*diu5_i_j_k
      diu6_i_j_kp1 = (-v(i,j,k)+v(i,j,k+1))/dzs(k)
      diu8_i_jp1_k = 2.*(-w(i,j,k)+w(i,j+1,k))/(dy1(j)+dy1(j+1))
      visvz2 = (evsz2)* ( diu6_i_j_kp1+diu8_i_jp1_k )
      diu6_i_j_k = (-v(i,j,k-1)+v(i,j,k))/dzs(k-1)
      diu8_i_jp1_km1 = 2.*(-w(i,j,k-1)+w(i,j+1,k-1))/(dy1(j)+dy1(j+1))
      visvz1 = (evsz1)* ( diu6_i_j_k+diu8_i_jp1_km1 )
      vfv = (visvx2-visvx1)/dx1(i) +(visvy2-visvy1)/dys(j) +(visvz2-visvz1)/dzn(k)
      g(i,j,k) = (g(i,j,k)+vfv)

end subroutine les_map_202


subroutine les_map_240(sm,dzn,dx1,dy1,u,dzs,w,v,h)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
    ! Local vars: i,j,k,diu3_im1_j_kp1,diu3_i_j_kp1,diu6_i_jm1_kp1,diu6_i_j_kp1,diu7_i_j_k,diu7_ip1_j_k,diu8_i_j_k,diu8_i_jp1_k,diu9_i_j_k,evsx2,evsx1,evsy2,evsy1,evsz2,evsz1,viswx2,viswx1,viswy2,viswy1,diu9_i_j_kp1,viswz2,viswz1,vfw
    integer :: i
    integer :: j
    integer :: k
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
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension((-1):(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: sm
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: u
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzs
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(In) :: w
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: v
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: h
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
      evsz2 = sm(i,j,k+1)
      evsz1 = sm(i,j,k)
      evsx2 = (dzn(k+1)*((dx1(i+1)*sm(i,j,k)+dx1(i)*sm(i+1,j, &
      k)) /(dx1(i)+dx1(i+1))) +dzn(k)*((dx1(i+1)*sm(i,j,k+1)+dx1(i)*sm(i+1,j, &
      k+1)) /(dx1(i)+dx1(i+1)))) /(dzn(k)+dzn(k+1))
      evsx1 = (dzn(k+1)*((dx1(i)*sm(i-1,j,k)+dx1(i-1)*sm(i,j, &
      k)) /(dx1(i-1)+dx1(i))) +dzn(k)*((dx1(i)*sm(i-1,j,k+1)+dx1(i-1)*sm(i,j, &
      k+1)) /(dx1(i-1)+dx1(i)))) /(dzn(k)+dzn(k+1))
      evsy2 = (dzn(k+1)*((dy1(j+1)*sm(i,j,k)+dy1(j)*sm(i,j+1, &
      k)) /(dy1(j)+dy1(j+1))) +dzn(k)*((dy1(j+1)*sm(i,j,k+1)+dy1(j)*sm(i,j+1, &
      k+1)) /(dy1(j)+dy1(j+1)))) /(dzn(k)+dzn(k+1))
      evsy1 = (dzn(k+1)*((dy1(j)*sm(i,j-1,k)+dy1(j-1)*sm(i,j, &
      k)) /(dy1(j-1)+dy1(j))) +dzn(k)*((dy1(j)*sm(i,j-1,k+1)+dy1(j-1)*sm(i,j, &
      k+1)) /(dy1(j-1)+dy1(j)))) /(dzn(k)+dzn(k+1))
      diu3_i_j_kp1 = (-u(i,j,k)+u(i,j,k+1))/dzs(k)
      diu7_ip1_j_k = 2.*(-w(i,j,k)+w(i+1,j,k))/(dx1(i)+dx1(i+1))
      viswx2 = (evsx2)* ( diu3_i_j_kp1 +diu7_ip1_j_k )
      diu3_im1_j_kp1 = (-u(i-1,j,k)+u(i-1,j,k+1))/dzs(k)
      diu7_i_j_k = 2.*(-w(i-1,j,k)+w(i,j,k))/(dx1(i-1)+dx1(i))
      viswx1 = (evsx1)* ( diu3_im1_j_kp1 +diu7_i_j_k )
      diu6_i_j_kp1 = (-v(i,j,k)+v(i,j,k+1))/dzs(k)
      diu8_i_jp1_k = 2.*(-w(i,j,k)+w(i,j+1,k))/(dy1(j)+dy1(j+1))
      viswy2 = (evsy2)* ( diu6_i_j_kp1 +diu8_i_jp1_k )
      diu6_i_jm1_kp1 = (-v(i,j-1,k)+v(i,j-1,k+1))/dzs(k)
      diu8_i_j_k = 2.*(-w(i,j-1,k)+w(i,j,k))/(dy1(j-1)+dy1(j))
      viswy1 = (evsy1)* ( diu6_i_jm1_kp1 +diu8_i_j_k )
      diu9_i_j_kp1 = (-w(i,j,k)+w(i,j,k+1))/dzn(k+1)
      viswz2 = (evsz2)*2.*diu9_i_j_kp1
      diu9_i_j_k = (-w(i,j,k-1)+w(i,j,k))/dzn(k)
      viswz1 = (evsz1)*2.*diu9_i_j_k
      vfw = (viswx2-viswx1)/dx1(i) +(viswy2-viswy1)/dy1(j) +(viswz2-viswz1)/dzn(k)
      h(i,j,k) = (h(i,j,k)+vfw)

end subroutine les_map_240


subroutine press_map_65(u,dx1,v,dy1,w,dzn,f,g,h,rhs,dt)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
        real, parameter  :: pjuge = 0.0001
        integer, parameter  :: nmaxp = 50 
        real, parameter  :: omega = 1.
    ! Local vars: i,j,k
    integer :: i
    integer :: j
    integer :: k
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: u
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: v
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(In) :: w
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In) :: f
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In) :: g
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In) :: h
    real(kind=4) :: dt
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)), intent(InOut) :: rhs
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
                rhs(i,j,k) = (-u(i-1,j,k)+u(i,j,k))/dx1(i) +(-v(i,j-1,k)+ &
                             v(i,j,k))/dy1(j) +(-w(i,j,k-1)+w(i,j,k))/dzn(k)
                rhs(i,j,k) = (f(i,j,k)-f(i-1,j,k))/dx1(i) +(g(i,j,k)- &
                              g(i,j-1,k))/dy1(j) +(h(i,j,k)-h(i,j,k-1))/dzn(k) &
                              +rhs(i,j,k)/dt

end subroutine press_map_65


subroutine press_reduce_78(dx1,dy1,dzn,rhs,rhsav,area)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
        real, parameter  :: pjuge = 0.0001
        integer, parameter  :: nmaxp = 50 
        real, parameter  :: omega = 1.
    ! Local vars: i,j,k
    integer :: i
    integer :: j
    integer :: k
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
    real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)), intent(In) :: rhs
! WRITTEN
! READ & WRITTEN
    real(kind=4) :: rhsav
    real(kind=4) :: area
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
                rhsav = rhsav+dx1(i)*dy1(j)*dzn(k)*rhs(i,j,k)
                area = area +dx1(i)*dy1(j)*dzn(k)

end subroutine press_reduce_78


subroutine press_map_87(rhs,rhsav)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
        real, parameter  :: pjuge = 0.0001
        integer, parameter  :: nmaxp = 50 
        real, parameter  :: omega = 1.
    ! Local vars: i,j,k
    integer :: i
    integer :: j
    integer :: k
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4) :: rhsav
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)), intent(InOut) :: rhs
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
                rhs(i,j,k) = rhs(i,j,k)-rhsav

end subroutine press_map_87


subroutine press_map_97(dzs,dys,dxs,nrd,p0,rhs,p1)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
        real, parameter  :: pjuge = 0.0001
        integer, parameter  :: nmaxp = 50 
        real, parameter  :: omega = 1.
    ! Local vars: i,j,k,dz1,dz2,cn4s,cn4l,cn3s,cn3l,cn2s,cn2l,cn1,reltmp
    integer :: i
    integer :: j
    integer :: k
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
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzs
    real(kind=4), dimension(0:jp), intent(In) :: dys
    real(kind=4), dimension(0:ip), intent(In) :: dxs
    real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)), intent(In) :: rhs
    integer :: nrd
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p1
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p0
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
                dz1 = dzs(k-1)
                dz2 = dzs(k)
                cn4s = 2./(dz1*(dz1+dz2))
                cn4l = 2./(dz2*(dz1+dz2))
                    cn3s = 2./(dys(j-1)*(dys(j-1)+dys(j)))
                    cn3l = 2./(dys(j)*(dys(j-1)+dys(j)))
                        cn2s = 2./(dxs(i-1)*(dxs(i-1)+dxs(i)))
                        cn2l = 2./(dxs(i)*(dxs(i-1)+dxs(i)))
                        cn1 = 1./ (2./(dxs(i-1)*dxs(i))  + 2./(dys(j-1)*dys(j)) + 2./(dz1*dz2))
                      if (nrd==0) then
                        reltmp = omega*(cn1 *(cn2l*p0(i+1,j,k) + &
                                 cn2s*p0(i-1,j,k) +cn3l*p0(i,j+1,k) + &
                                 cn3s*p0(i,j-1,k) +cn4l*p0(i,j,k+1) + &
                                 cn4s*p0(i,j,k-1) -rhs(i,j,k))-p0(i,j,k))
                        p1(i,j,k) = p0(i,j,k) +reltmp
                      else
                        reltmp = omega*(cn1 *(cn2l*p1(i+1,j,k) + &
                                 cn2s*p1(i-1,j,k) +cn3l*p1(i,j+1,k) + &
                                 cn3s*p1(i,j-1,k) +cn4l*p1(i,j,k+1) + &
                                 cn4s*p1(i,j,k-1) -rhs(i,j,k))-p1(i,j,k))
                        p0(i,j,k) = p1(i,j,k) +reltmp
                      end if

end subroutine press_map_97


subroutine press_reduce_129(p0,dx1,dy1,dzn,pav,pco)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
        real, parameter  :: pjuge = 0.0001
        integer, parameter  :: nmaxp = 50 
        real, parameter  :: omega = 1.
    ! Local vars: i,j,k
    integer :: i
    integer :: j
    integer :: k
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(In) :: p0
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
! WRITTEN
! READ & WRITTEN
    real(kind=4) :: pav
    real(kind=4) :: pco
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
                pav = pav+p0(i,j,k)*dx1(i)*dy1(j)*dzn(k)
                pco = pco+dx1(i)*dy1(j)*dzn(k)

end subroutine press_reduce_129


subroutine press_map_138(p0,pav)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
        real, parameter  :: pjuge = 0.0001
        integer, parameter  :: nmaxp = 50 
        real, parameter  :: omega = 1.
    ! Local vars: i,j,k
    integer :: i
    integer :: j
    integer :: k
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4) :: pav
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p0
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
                p0(i,j,k) = p0(i,j,k)-pav

end subroutine press_map_138


subroutine velFG_map_95(u,dx1,v,dy1,w,dzs,dzn,f)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
            integer, parameter  :: u0 = 0
    ! Local vars: nou1_,i,j,k,diu1_,cov1_i,nou1_ip1,diu1_ip1,cov1_ip1,nou2_,diu2_,cov2_j,nou2_jp1,diu2_jp1,cov2_jp1,nou3_,diu3_,cov3_k,nou3_kp1,diu3_kp1,cov3_kp1,covx1,covy1,covz1,covc
    real(kind=4) :: nou1_
    integer :: i
    integer :: j
    integer :: k
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
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: u
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: v
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(In) :: w
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzs
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
! WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(Out) :: f
! READ & WRITTEN
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
      nou1_ = ( u(i-1,j,k)+u(i,j,k))/2. 
      diu1_ = (-u(i-1,j,k)+u(i,j,k))/dx1(i) 
      cov1_i = nou1_*diu1_ 
      nou1_ip1 = ( u(i,j,k)+u(i+1,j,k))/2. 
      diu1_ip1 = (-u(i,j,k)+u(i+1,j,k))/dx1(i+1) 
      cov1_ip1 = nou1_ip1*diu1_ip1 
      if (i==ip) cov1_ip1 = cov1_i
      nou2_ = (dx1(i+1)*v(i,j-1,k)+dx1(i)*v(i+1,j-1,k)) /(dx1(i)+dx1(i+1))
      diu2_ = 2.*(-u(i,j-1,k)+u(i,j,k))/(dy1(j-1)+dy1(j))
      cov2_j = nou2_*diu2_
      nou2_jp1 = (dx1(i+1)*v(i,j,k)+dx1(i)*v(i+1,j,k)) /(dx1(i)+dx1(i+1))
      diu2_jp1 = 2.*(-u(i,j,k)+u(i,j+1,k))/(dy1(j)+dy1(j+1))
      cov2_jp1 = nou2_jp1*diu2_jp1
      if (j==jp) then
          nou2_ = (dx1(i+1)*v(i,0,k)+dx1(i)*v(i+1,0,k)) /(dx1(i)+dx1(i+1))
          diu2_ = 2.*(-u(i,0,k)+u(i,1,k))/(dy1(0)+dy1(1))
          cov2_jp1 = nou2_*diu2_
      end if
        nou3_ = (dx1(i+1)*w(i,j,k-1)+dx1(i)*w(i+1,j,k-1)) /(dx1(i)+dx1(i+1))
        diu3_ = (-u(i,j,k-1)+u(i,j,k))/dzs(k-1)
        cov3_k = nou3_*diu3_
        nou3_kp1 = (dx1(i+1)*w(i,j,k)+dx1(i)*w(i+1,j,k)) /(dx1(i)+dx1(i+1))
        diu3_kp1 = (-u(i,j,k)+u(i,j,k+1))/dzs(k)
        cov3_kp1 = nou3_kp1*diu3_kp1
        if (k==1) then
            nou3_ = 0.5*(dx1(i+1)*w(i,j,1)+dx1(i)*w(i+1,j,1))/(dx1(i)+dx1(i+1))
            diu3_ = 0.4*u(i,j,1) / ( alog( 5.0 * dzn(1) ) * 0.2 * dzn(1) )
            cov3_k = nou3_*diu3_
        end if
        covx1 = (dx1(i+1)*cov1_i+dx1(i)*cov1_ip1) /(dx1(i)+dx1(i+1))
        covy1 = (cov2_j+cov2_jp1)/2.
        covz1 = (cov3_k+cov3_kp1)/2.
        covc = covx1+covy1+covz1
        f(i,j,k) = (-covc)

end subroutine velFG_map_95


subroutine velFG_map_135(dy1,u,v,dx1,w,dzs,dzn,g)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
            integer, parameter  :: u0 = 0
    ! Local vars: i,j,k,covx1,covy1,covz1,covc,nou4_,diu4_,cov4_i,nou4_ip1,diu4_ip1,cov4_ip1,nou5_,diu5_,cov5_j,nou5_jp1,diu5_jp1,cov5_jp1,nou6_,diu6_,cov6_k,nou6_kp1,diu6_kp1,cov6_kp1
    integer :: i
    integer :: j
    integer :: k
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
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: u
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: v
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(In) :: w
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzs
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
! WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(Out) :: g
! READ & WRITTEN
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
      nou4_ = (dy1(j+1)*u(i-1,j,k)+dy1(j)*u(i-1,j+1,k)) /(dy1(j)+dy1(j+1))
      diu4_ = 2.*(-v(i-1,j,k)+v(i,j,k))/(dx1(i-1)+dx1(i))
      cov4_i = (nou4_-u0)*diu4_
      nou4_ip1 = (dy1(j+1)*u(i,j,k)+dy1(j)*u(i,j+1,k)) /(dy1(j)+dy1(j+1))
      diu4_ip1 = 2.*(-v(i,j,k)+v(i+1,j,k))/(dx1(i)+dx1(i+1))
      cov4_ip1 = (nou4_ip1-u0)*diu4_ip1
      if (i==ip) cov4_ip1 = cov4_i
      nou5_ = ( v(i,j-1,k)+v(i,j,k))/2.
      diu5_ = (-v(i,j-1,k)+v(i,j,k))/dy1(j)
      cov5_j = nou5_*diu5_
      nou5_jp1 = ( v(i,j,k)+v(i,j+1,k))/2.
      diu5_jp1 = (-v(i,j,k)+v(i,j+1,k))/dy1(j+1)
      cov5_jp1 = nou5_jp1*diu5_jp1
      if (j==jp) then
          nou5_ = ( v(i,1,k)+v(i,2,k))/2.
          diu5_ = (-v(i,1,k)+v(i,2,k))/dy1(j)
          cov5_jp1 = nou5_*diu5_
      end if  
        nou6_ = (dy1(j+1)*w(i,j,k-1)+dy1(j)*w(i,j+1,k-1)) /(dy1(j)+dy1(j+1))
        diu6_ = (-v(i,j,k-1)+v(i,j,k))/dzs(k-1)
        cov6_k = nou6_*diu6_
        nou6_kp1 = (dy1(j+1)*w(i,j,k)+dy1(j)*w(i,j+1,k)) /(dy1(j)+dy1(j+1))
        diu6_kp1 = (-v(i,j,k)+v(i,j,k+1))/dzs(k)
        cov6_kp1 = nou6_kp1*diu6_kp1
        if (k==1) then
            nou6_ = 0.5*(dy1(j+1)*w(i,j,1)+dy1(j)*w(i,j+1,1))/(dy1(j)+dy1(j+1))
            diu6_= 0.4*v(i,j,1) / (alog(5.0*dzn(1)) * 0.2 * dzn(1))
            cov6_k = nou6_*diu6_
        end if
        covx1 = (cov4_i+cov4_ip1)/2.
        covy1 = (dy1(j+1)*cov5_j+dy1(j)*cov5_jp1) /(dy1(j)+dy1(j+1))
        covz1 = (cov6_k+cov6_kp1)/2.
        covc = covx1+covy1+covz1
        g(i,j,k) = (-covc)

end subroutine velFG_map_135


subroutine velFG_map_175(dzn,u,w,dx1,v,dy1,h)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
            integer, parameter  :: u0 = 0
    ! Local vars: i,j,k,covx1,covy1,covz1,covc,nou7_,diu7_,cov7_i,nou7_ip1,diu7_ip1,cov7_ip1,nou8_,diu8_,cov8_j,nou8_jp1,diu8_jp1,cov8_jp1,nou9_,diu9_,cov9_k,nou9_kp1,diu9_kp1,cov9_kp1
    integer :: i
    integer :: j
    integer :: k
    real(kind=4) :: covx1
    real(kind=4) :: covy1
    real(kind=4) :: covz1
    real(kind=4) :: covc
    real(kind=4) :: nou7_
    real(kind=4) :: diu7_
    real(kind=4) :: cov7_i
    real(kind=4) :: nou7_ip1
    real(kind=4) :: diu7_ip1
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
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: u
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(In) :: w
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: v
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
! WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(Out) :: h
! READ & WRITTEN
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = (((80 - 1) - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
      nou7_ = (dzn(k+1)*u(i-1,j,k)+dzn(k)*u(i-1,j,k+1)) /(dzn(k)+dzn(k+1))
      diu7_ = 2.*(-w(i-1,j,k)+w(i,j,k))/(dx1(i-1)+dx1(i))
      cov7_i = (nou7_-u0)*diu7_
      nou7_ip1 = (dzn(k+1)*u(i,j,k)+dzn(k)*u(i,j,k+1)) /(dzn(k)+dzn(k+1))
      diu7_ip1 = 2.*(-w(i,j,k)+w(i+1,j,k))/(dx1(i)+dx1(i+1))
      cov7_ip1 = (nou7_-u0)*diu7_
      if (i==ip) cov7_ip1 = cov7_i
      nou8_ = (dzn(k+1)*v(i,j-1,k)+dzn(k)*v(i,j-1,k+1)) /(dzn(k)+dzn(k+1))
      diu8_ = 2.*(-w(i,j-1,k)+w(i,j,k))/(dy1(j-1)+dy1(j))
      cov8_j = nou8_*diu8_
      nou8_jp1 = (dzn(k+1)*v(i,j,k)+dzn(k)*v(i,j,k+1)) /(dzn(k)+dzn(k+1))
      diu8_jp1 = 2.*(-w(i,j,k)+w(i,j+1,k))/(dy1(j)+dy1(j+1))
      cov8_jp1 = nou8_jp1*diu8_jp1
      if (j==jp) then
        nou8_ = (dzn(k+1)*v(i,0,k)+dzn(k)*v(i,0,k+1)) /(dzn(k)+dzn(k+1))
        diu8_ = 2.*(-w(i,0,k)+w(i,1,k))/(dy1(0)+dy1(1))
        cov8_jp1 = nou8_*diu8_
      end if
      nou9_ = ( w(i,j,k-1)+w(i,j,k))/2.
      diu9_ = (-w(i,j,k-1)+w(i,j,k))/dzn(k)
      cov9_k = nou9_*diu9_
      nou9_kp1 = ( w(i,j,k)+w(i,j,k+1))/2.
      diu9_kp1 = (-w(i,j,k)+w(i,j,k+1))/dzn(k+1)
      cov9_kp1 = nou9_kp1*diu9_kp1
       covx1 = (cov7_i+cov7_ip1)/2.
       covy1 = (cov8_j+cov8_jp1)/2.
       covz1 = (dzn(k+1)*cov9_k+dzn(k)*cov9_kp1) /(dzn(k)+dzn(k+1))
       covc = covx1+covy1+covz1
        h(i,j,k) = (-covc)      

end subroutine velFG_map_175


subroutine velnw_map_36(p0,ro,dxs,u,dt,f)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
    ! Local vars: pz,i,j,k
    real(kind=4) :: pz
    integer :: i
    integer :: j
    integer :: k
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(In) :: p0
    real(kind=4), dimension(0:ip), intent(In) :: dxs
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In) :: f
    real(kind=4) :: ro
    real(kind=4) :: dt
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: u
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
        pz = (-p0(i,j,k)+p0(i+1,j,k))/ro/dxs(i)
        u(i,j,k) = u(i,j,k)+dt*(f(i,j,k)-pz)

end subroutine velnw_map_36


subroutine velnw_map_44(p0,ro,dys,v,dt,g)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
    ! Local vars: pz,i,j,k
    real(kind=4) :: pz
    integer :: i
    integer :: j
    integer :: k
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(In) :: p0
    real(kind=4), dimension(0:jp), intent(In) :: dys
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In) :: g
    real(kind=4) :: ro
    real(kind=4) :: dt
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: v
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
        pz = (-p0(i,j,k)+p0(i,j+1,k))/ro/dys(j)
        v(i,j,k) = v(i,j,k)+dt*(g(i,j,k)-pz)

end subroutine velnw_map_44


subroutine velnw_map_52(p0,ro,dzs,w,dt,h)

    integer, parameter :: kp=80
        integer, parameter :: ip = 300
        integer, parameter :: jp = 300
        integer, parameter :: ipmax = ip
        integer, parameter :: jpmax = jp
        real, parameter :: dxgrid = 4.
        real, parameter :: dygrid = 4.
        real, parameter :: cs0 = 0.14
        integer, parameter :: i_anime=1
        integer, parameter :: avetime=2
        integer, parameter :: km_sl=80
        integer, parameter :: i_aveflow=0
        integer, parameter :: i_ifdata_out=0
        real, parameter :: dt_orig = 0.05 
    ! Local vars: pz,i,j,k
    real(kind=4) :: pz
    integer :: i
    integer :: j
    integer :: k
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(In) :: p0
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzs
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In) :: h
    real(kind=4) :: ro
    real(kind=4) :: dt
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(InOut) :: w
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = (((80 - 1) - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
        pz = (-p0(i,j,k)+p0(i,j,k+1))/ro/dzs(k)
        w(i,j,k) = w(i,j,k)+dt*(h(i,j,k)-pz)

end subroutine velnw_map_52


subroutine adam_feedbf_les_press_velfg_ve_etc_superkernel(f,g,h,fold,gold,hold,usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,fx,fy,fz,dx1,dy1,dzs,dzn,delx1,sm,dxs,dys,rhs,global_rhsav_array,global_area_array,rhsav,nrd,p0,p1,global_pav_array,global_pco_array,pav,ro,state_ptr)
use module_adam_feedbf_les_press_velfg_ve_etc_superkernel_init
  real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: u
  real(kind=4), dimension((-1):(ip + 1),0:(jp + 1),0:(kp + 1)), intent(In) :: bmask1
  real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: v
  real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: cmask1
  real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(InOut) :: w
  real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)), intent(In) :: dmask1
  real(kind=4), intent(In), dimension(1) :: alpha
  real(kind=4), intent(In), dimension(1) :: dt
  real(kind=4), intent(In), dimension(1) :: beta
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: fx
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: fy
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: fz
  real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
  real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
  real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzs
  real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
  real(kind=4), dimension(kp), intent(In) :: delx1
  real(kind=4), dimension((-1):(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: sm
  real(kind=4), dimension(0:ip), intent(In) :: dxs
  real(kind=4), dimension(0:jp), intent(In) :: dys
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: f
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: g
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: h
  real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)), intent(InOut) :: rhs
  real(kind=4), intent(In), dimension(1) :: rhsav
  integer, intent(In), dimension(1) :: nrd
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p0
  real(kind=4), intent(In), dimension(1) :: pav
  integer, intent(In), dimension(1) :: k
  integer, intent(In), dimension(1) :: j
  real(kind=4), intent(In), dimension(1) :: ro
  real(kind=4), dimension(ip,jp,kp), intent(InOut) :: fold
  real(kind=4), dimension(ip,jp,kp), intent(InOut) :: gold
  real(kind=4), dimension(ip,jp,kp), intent(InOut) :: hold
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: usum
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: vsum
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: wsum
  real(kind=4), dimension(1:NUNITS), intent(Out) :: global_rhsav_array
  real(kind=4), dimension(1:NUNITS), intent(Out) :: global_area_array
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p1
  real(kind=4), dimension(1:NUNITS), intent(Out) :: global_pav_array
  real(kind=4), dimension(1:NUNITS), intent(Out) :: global_pco_array

  integer :: state
  integer, dimension(1) :: state_ptr
integer, parameter :: ST_ADAM_MAP_36 = 0 !  adam_map_36
integer, parameter :: ST_FEEDBF_MAP_49 = 1 !  feedbf_map_49
integer, parameter :: ST_FEEDBF_MAP_67 = 2 !  feedbf_map_67
integer, parameter :: ST_LES_MAP_119 = 3 !  les_map_119
integer, parameter :: ST_LES_MAP_164 = 4 !  les_map_164
integer, parameter :: ST_LES_MAP_202 = 5 !  les_map_202
integer, parameter :: ST_LES_MAP_240 = 6 !  les_map_240
integer, parameter :: ST_PRESS_MAP_65 = 7 !  press_map_65
integer, parameter :: ST_PRESS_REDUCE_78 = 8 !  press_reduce_78
integer, parameter :: ST_PRESS_MAP_87 = 9 !  press_map_87
integer, parameter :: ST_PRESS_MAP_97 = 10 !  press_map_97
integer, parameter :: ST_PRESS_REDUCE_129 = 11 !  press_reduce_129
integer, parameter :: ST_PRESS_MAP_138 = 12 !  press_map_138
integer, parameter :: ST_VELFG_MAP_95 = 13 !  velFG_map_95
integer, parameter :: ST_VELFG_MAP_135 = 14 !  velFG_map_135
integer, parameter :: ST_VELFG_MAP_175 = 15 !  velFG_map_175
integer, parameter :: ST_VELNW_MAP_36 = 16 !  velnw_map_36
integer, parameter :: ST_VELNW_MAP_44 = 17 !  velnw_map_44
integer, parameter :: ST_VELNW_MAP_52 = 18 !  velnw_map_52
  state = state_ptr(1) ! state 
! SUPERKERNEL BODY
  select case(state)
    case (ST_ADAM_MAP_36)
      call adam_map_36(f,g,h,fold,gold,hold)
    case (ST_FEEDBF_MAP_49)
      call feedbf_map_49(usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,fx,fy,fz)
    case (ST_FEEDBF_MAP_67)
      call feedbf_map_67(f,fx,g,fy,h,fz)
    case (ST_LES_MAP_119)
      call les_map_119(u,dx1,dy1,dzs,v,w,dzn,delx1,sm)
    case (ST_LES_MAP_164)
      call les_map_164(sm,dy1,dx1,dzn,u,v,dzs,w,dxs,f)
    case (ST_LES_MAP_202)
      call les_map_202(sm,dy1,dx1,dzn,u,v,dzs,w,dys,g)
    case (ST_LES_MAP_240)
      call les_map_240(sm,dzn,dx1,dy1,u,dzs,w,v,h)
    case (ST_PRESS_MAP_65)
      call press_map_65(u,dx1,v,dy1,w,dzn,f,g,h,rhs,dt)
    case (ST_PRESS_REDUCE_78)
      call press_reduce_78(dx1,dy1,dzn,rhs,global_rhsav_array,global_area_array)
    case (ST_PRESS_MAP_87)
      call press_map_87(rhs,rhsav)
    case (ST_PRESS_MAP_97)
      call press_map_97(dzs,dys,dxs,nrd,p0,rhs,p1)
    case (ST_PRESS_REDUCE_129)
      call press_reduce_129(p0,dx1,dy1,dzn,global_pav_array,global_pco_array)
    case (ST_PRESS_MAP_138)
      call press_map_138(p0,pav)
    case (ST_VELFG_MAP_95)
      call velFG_map_95(u,dx1,v,dy1,w,dzs,dzn,f)
    case (ST_VELFG_MAP_135)
      call velFG_map_135(dy1,u,v,dx1,w,dzs,dzn,g)
    case (ST_VELFG_MAP_175)
      call velFG_map_175(dzn,u,w,dx1,v,dy1,h)
    case (ST_VELNW_MAP_36)
      call velnw_map_36(p0,ro,dxs,u,dt,f)
    case (ST_VELNW_MAP_44)
      call velnw_map_44(p0,ro,dys,v,dt,g)
    case (ST_VELNW_MAP_52)
      call velnw_map_52(p0,ro,dzs,w,dt,h)
  end select
end subroutine adam_feedbf_les_press_velfg_ve_etc_superkernel
end module module_adam_feedbf_les_press_velfg_ve_etc_superkernel