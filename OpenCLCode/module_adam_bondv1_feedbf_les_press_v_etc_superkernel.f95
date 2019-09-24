module module_adam_bondv1_feedbf_les_press_v_etc_superkernel


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


subroutine bondv1_map_38(z2,dzn,u,v,w)

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
    ! Local vars: u_val,k,i,j
    real(kind=4) :: u_val
    integer :: k
    integer :: i
    integer :: j
    ! ParallelFortran: Synthesised loop variable decls
    integer :: i_range
    integer :: k_range
    integer :: j_range
    integer :: i_rel
    integer :: k_rel
    integer :: j_rel
! READ
    real(kind=4), dimension(0:(kp + 2)), intent(In) :: z2
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
! WRITTEN
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(Out) :: u
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(Out) :: v
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(Out) :: w
! READ & WRITTEN
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    i_range = ((1 - 0) + 1)
    k_range = ((78 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_rel = (global_id / (k_range * j_range))
    i = (i_rel + 0)
    k_rel = ((global_id - (i_rel * (k_range * j_range))) / j_range)
    k = (k_rel + 1)
    j_rel = ((global_id - (i_rel * (k_range * j_range))) - (k_rel * j_range))
    j = (j_rel + 1)


    ! ParallelFortran: Original code
                    u_val = 5.*((z2(k)+0.5*dzn(k))/600.)**0.2
                    u(i,j,k) = u_val
                    v(i,j,k) = 0.0
                    w(i,j,k) = 0.0

end subroutine bondv1_map_38


subroutine bondv1_map_48(u,v,w)

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
    ! Local vars: k,i,j
    integer :: k
    integer :: i
    integer :: j
    ! ParallelFortran: Synthesised loop variable decls
    integer :: i_range
    integer :: k_range
    integer :: j_range
    integer :: i_rel
    integer :: k_rel
    integer :: j_rel
! READ
! WRITTEN
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(Out) :: v
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(Out) :: w
! READ & WRITTEN
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: u
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    i_range = ((1 - 0) + 1)
    k_range = ((80 - 79) + 1)
    j_range = ((300 - 1) + 1)
    i_rel = (global_id / (k_range * j_range))
    i = (i_rel + 0)
    k_rel = ((global_id - (i_rel * (k_range * j_range))) / j_range)
    k = (k_rel + 79)
    j_rel = ((global_id - (i_rel * (k_range * j_range))) - (k_rel * j_range))
    j = (j_rel + 1)


    ! ParallelFortran: Original code
                    u(i,j,k) = u(i,j,77)
                    v(i,j,k) = 0.0
                    w(i,j,k) = 0.0

end subroutine bondv1_map_48


subroutine bondv1_map_58(u,v,w)

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
    ! Local vars: k,i,j
    integer :: k
    integer :: i
    integer :: j
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
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: u
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: v
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(InOut) :: w
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 2) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 2)


    ! ParallelFortran: Original code
                    u(i,j,k) = u(1,j,k)
                    v(i,j,k) = v(1,j,k)
                    w(i,j,k) = w(1,j,k)

end subroutine bondv1_map_58


subroutine bondv1_reduce_69(u,global_aaa_array)

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
    ! Missing args: 
    ! Local vars: k,j

    integer :: k
    integer :: j

    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: k_rel
    integer :: j_rel

    integer :: local_id
    integer :: local_id_fortran
    integer :: group_id
    integer :: group_id_fortran
    integer :: global_id
    integer :: r_iter
    integer :: local_chunk_size
    integer :: start_position

    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: u
    real(kind=4), dimension(1:NUNITS), intent(Out) :: global_aaa_array

    real(kind=4) :: local_aaa

    local_chunk_size = (80 * 300)
    start_position = 0
    local_aaa = 0

    do r_iter=start_position, ((start_position + local_chunk_size) - 1)
        k_range = ((80 - 1) + 1)
        j_range = ((300 - 1) + 1)
        k_rel = (r_iter / j_range)
        k = (k_rel + 1)
        j_rel = (r_iter - (k_rel * j_range))
        j = (j_rel + 1)
        local_aaa = amax1(local_aaa,u(300,j,k))
    end do

    global_aaa_array(1) = local_aaa

end subroutine bondv1_reduce_69


subroutine bondv1_reduce_76(u,global_bbb_array)

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
    ! Missing args: 
    ! Local vars: k,j

    integer :: k
    integer :: j

    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: k_rel
    integer :: j_rel

    integer :: local_id
    integer :: local_id_fortran
    integer :: group_id
    integer :: group_id_fortran
    integer :: global_id
    integer :: r_iter
    integer :: local_chunk_size
    integer :: start_position

    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: u
    real(kind=4), dimension(1:NUNITS), intent(Out) :: global_bbb_array

    real(kind=4) :: local_bbb

    local_chunk_size = (80 * 300)
    start_position = 0
    local_bbb = 1e38

    do r_iter=start_position, ((start_position + local_chunk_size) - 1)
        k_range = ((80 - 1) + 1)
        j_range = ((300 - 1) + 1)
        k_rel = (r_iter / j_range)
        k = (k_rel + 1)
        j_rel = (r_iter - (k_rel * j_range))
        j = (j_rel + 1)
        local_bbb = amin1(local_bbb,u(300,j,k))
    end do

    global_bbb_array(1) = local_bbb

end subroutine bondv1_reduce_76


subroutine bondv1_map_83(u,dt,uout,dxs,v,w)

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
    ! Local vars: k,j
    integer :: k
    integer :: j
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: k_rel
    integer :: j_rel
! READ
    real(kind=4), dimension(0:ip), intent(In) :: dxs
    real(kind=4) :: dt
    real(kind=4) :: uout
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: u
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: v
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(InOut) :: w
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    k_rel = (global_id / j_range)
    k = (k_rel + 1)
    j_rel = (global_id - (k_rel * j_range))
    j = (j_rel + 1)


    ! ParallelFortran: Original code
              u(ip,j,k) = u(ip,j,k)-dt*uout *(u(ip,j,k)-u(ip-1,j,k))/dxs(ip)
              v(ip+1,j,k) = v(ip+1,j,k)-dt*uout *(v(ip+1,j,k)-v(ip,j,k))/dxs(ip)
              w(ip+1,j,k) = w(ip+1,j,k)-dt*uout *(w(ip+1,j,k)-w(ip,j,k))/dxs(ip)

end subroutine bondv1_map_83


subroutine bondv1_map_98(u,v,w)

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
    ! Local vars: k,i
    integer :: k
    integer :: i
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: i_range
    integer :: k_rel
    integer :: i_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: u
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: v
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(InOut) :: w
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = (((80 + 1) - 0) + 1)
    i_range = (((300 + 1) - 0) + 1)
    k_rel = (global_id / i_range)
    k = (k_rel + 0)
    i_rel = (global_id - (k_rel * i_range))
    i = (i_rel + 0)


    ! ParallelFortran: Original code
            u(i,   0,k) = u(i,jp  ,k)
            u(i,jp+1,k) = u(i,   1,k)
            v(i,   0,k) = v(i,jp  ,k)
            v(i,jp+1,k) = v(i,   1,k)
    if ((k < 80)) then
            w(i,   0,k) = w(i,jp  ,k)
            w(i,jp+1,k) = w(i,   1,k)
    end if

end subroutine bondv1_map_98


subroutine bondv1_map_116(u,v)

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
    ! Local vars: i,j
    integer :: i
    integer :: j
    ! ParallelFortran: Synthesised loop variable decls
    integer :: j_range
    integer :: i_range
    integer :: j_rel
    integer :: i_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: u
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: v
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    j_range = (((300 + 1) - 0) + 1)
    i_range = (((300 + 1) - 0) + 1)
    j_rel = (global_id / i_range)
    j = (j_rel + 0)
    i_rel = (global_id - (j_rel * i_range))
    i = (i_rel + 0)


    ! ParallelFortran: Original code
            u(i,j,   0) = -u(i,j, 1)
            u(i,j,kp+1) = u(i,j,kp)
            v(i,j,   0) = -v(i,j, 1)
            v(i,j,kp+1) = v(i,j,kp)

end subroutine bondv1_map_116


subroutine bondv1_map_128(w)

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
    ! Local vars: i,j
    integer :: i
    integer :: j
    ! ParallelFortran: Synthesised loop variable decls
    integer :: j_range
    integer :: i_range
    integer :: j_rel
    integer :: i_rel
! READ
! WRITTEN
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(Out) :: w
! READ & WRITTEN
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    j_range = (((300 + 1) - (-1)) + 1)
    i_range = (((300 + 1) - 0) + 1)
    j_rel = (global_id / i_range)
    j = (j_rel + (-1))
    i_rel = (global_id - (j_rel * i_range))
    i = (i_rel + 0)


    ! ParallelFortran: Original code
            w(i,j, 0) = 0.0
            w(i,j,kp) = 0.0

end subroutine bondv1_map_128


subroutine feedbf_map_35(kp,usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,jp,ip,f,fx,g,fy,h,fz)
    use common_sn

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
    real :: kp
    real(kind=4) :: alpha
    real(kind=4) :: dt
    real(kind=4) :: beta
    real :: jp
    real :: ip
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: usum
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: vsum
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: wsum
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: fx
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: fy
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: fz
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: f
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: g
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: h
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((kp - 1) + 1)
    j_range = ((jp - 1) + 1)
    i_range = ((ip - 1) + 1)
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
                f(i,j,k) = f(i,j,k)+fx(i,j,k)
                g(i,j,k) = g(i,j,k)+fy(i,j,k)
                h(i,j,k) = h(i,j,k)+fz(i,j,k)

end subroutine feedbf_map_35


subroutine les_map_88(dx1,dy1,dzn,delx1)

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
    ! Local vars: k
    integer :: k
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: k_rel
! READ
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
! WRITTEN
    real(kind=4), dimension(kp), intent(Out) :: delx1
! READ & WRITTEN
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_rel = global_id
    k = (k_rel + 1)


    ! ParallelFortran: Original code
          delx1(k) = (dx1(0)*dy1(0)*dzn(k))**(1./3.)

end subroutine les_map_88


subroutine les_map_91(diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,delx1,sm)

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
    ! Local vars: k,dudxx1,i,j,dudyx1,dudzx1,dvdxx1,dvdyx1,dvdzx1,dwdxx1,dwdyx1,dwdzx1,csx1
    integer :: k
    real(kind=4) :: dudxx1
    integer :: i
    integer :: j
    real(kind=4) :: dudyx1
    real(kind=4) :: dudzx1
    real(kind=4) :: dvdxx1
    real(kind=4) :: dvdyx1
    real(kind=4) :: dvdzx1
    real(kind=4) :: dwdxx1
    real(kind=4) :: dwdyx1
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
    real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu1
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu2
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu3
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu4
    real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu5
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu6
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu7
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu8
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu9
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
      dudxx1 =  diu1(i,j,k)
      dudyx1 =  (diu2(i-1,j,k)+diu2(i-1,j+1,k) +diu2(i  ,j,k)+diu2(i  ,j+1,k) ) *.25
      dudzx1 =  (diu3(i-1,j,k)+diu3(i-1,j,k+1) +diu3(i  ,j,k)+diu3(i  ,j,k+1) ) *.25
      dvdxx1 =  (diu4(i  ,j,k)+diu4(i  ,j-1,k) +diu4(i+1,j,k)+diu4(i+1,j-1,k) ) *.25
      dvdyx1 =  diu5(i,j,k)
      dvdzx1 =  (diu6(i,j-1,k)+diu6(i,j-1,k+1) +diu6(i,j  ,k)+diu6(i,j  ,k+1) )  *.25
      dwdxx1 =  (diu7(i  ,j,k)+diu7(i  ,j,k-1) +diu7(i+1,j,k)+diu7(i+1,j,k-1) ) *.25
      dwdyx1 =  (diu8(i,j  ,k)+diu8(i,j  ,k-1) +diu8(i,j+1,k)+diu8(i,j+1,k-1) ) *.25
      dwdzx1 =  diu9(i,j,k)
      csx1 = cs0
      sm(i,j,k) = ( csx1*delx1(k) )**2  * sqrt( 2.*( dudxx1**2+dvdyx1**2+dwdzx1**2 ) +( dudyx1+dvdxx1 )**2  &
      +( dwdyx1+dvdzx1 )**2 +( dudzx1+dwdxx1 )**2 )

end subroutine les_map_91


subroutine les_map_109(sm)

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
    ! Local vars: k,j
    integer :: k
    integer :: j
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: k_rel
    integer :: j_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension((-1):(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: sm
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = (((80 + 1) - 0) + 1)
    j_range = (((300 + 1) - (-1)) + 1)
    k_rel = (global_id / j_range)
    k = (k_rel + 0)
    j_rel = (global_id - (k_rel * j_range))
    j = (j_rel + (-1))


    ! ParallelFortran: Original code
                    sm(   0,j,k) = sm(1 ,j,k) 
                    sm(ip+1,j,k) = sm(ip,j,k)

end subroutine les_map_109


subroutine les_map_115(sm)

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
    ! Local vars: k,i
    integer :: k
    integer :: i
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: i_range
    integer :: k_rel
    integer :: i_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension((-1):(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: sm
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = (((80 + 1) - 0) + 1)
    i_range = (((300 + 1) - 0) + 1)
    k_rel = (global_id / i_range)
    k = (k_rel + 0)
    i_rel = (global_id - (k_rel * i_range))
    i = (i_rel + 0)


    ! ParallelFortran: Original code
                    sm(i,jp+1,k) = sm(i,jp  ,k)
                    sm(i,0,k) = sm(i,1   ,k) 

end subroutine les_map_115


subroutine les_map_121(sm)

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
    ! Local vars: i,j
    integer :: i
    integer :: j
    ! ParallelFortran: Synthesised loop variable decls
    integer :: j_range
    integer :: i_range
    integer :: j_rel
    integer :: i_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension((-1):(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: sm
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    j_range = (((300 + 1) - (-1)) + 1)
    i_range = (((300 + 1) - 0) + 1)
    j_rel = (global_id / i_range)
    j = (j_rel + (-1))
    i_rel = (global_id - (j_rel * i_range))
    i = (i_rel + 0)


    ! ParallelFortran: Original code
            sm(i,j,   0) = -sm(i,j, 1)
            sm(i,j,kp+1) = sm(i,j,kp)

end subroutine les_map_121


subroutine les_map_127(sm,dy1,dx1,dzn,diu1,diu2,diu4,diu3,diu7,dxs,f)

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
    ! Local vars: k,i,j,evsx2,evsx1,evsy2,evsy1,evsz2,evsz1,visux2,visux1,visuy2,visuy1,visuz2,visuz1,vfu
    integer :: k
    integer :: i
    integer :: j
    real(kind=4) :: evsx2
    real(kind=4) :: evsx1
    real(kind=4) :: evsy2
    real(kind=4) :: evsy1
    real(kind=4) :: evsz2
    real(kind=4) :: evsz1
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
    real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu1
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu2
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu4
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu3
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu7
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
      visux2 = (evsx2)*2.*diu1(i+1,j  ,k  )
      visux1 = (evsx1)*2.*diu1(i  ,j,  k  )
      visuy2 = (evsy2)* ( diu2(i  ,j+1,k  )+diu4(i+1,j  ,k  ) )
      visuy1 = (evsy1)* ( diu2(i  ,j  ,k  )+diu4(i+1,j-1,k  ) )
      visuz2 = (evsz2)* ( diu3(i  ,j  ,k+1)+diu7(i+1,j  ,k  ) )
      visuz1 = (evsz1)* ( diu3(i  ,j  ,k  )+diu7(i+1,j  ,k-1) )
      vfu = (visux2-visux1)/dxs(i) +(visuy2-visuy1)/dy1(j) +(visuz2-visuz1)/dzn(k)
      f(i,j,k) = (f(i,j,k)+vfu)

end subroutine les_map_127


subroutine les_map_155(sm,dy1,dx1,dzn,diu1,diu2,diu4,diu3,diu7,uspd,u,dxs,f)

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
    ! Local vars: i,j,evsx2,evsx1,evsy2,evsy1,evsz2,visux2,visux1,visuy2,visuy1,visuz2,visuz1,vfu
    integer :: i
    integer :: j
    real(kind=4) :: evsx2
    real(kind=4) :: evsx1
    real(kind=4) :: evsy2
    real(kind=4) :: evsy1
    real(kind=4) :: evsz2
    real(kind=4) :: visux2
    real(kind=4) :: visux1
    real(kind=4) :: visuy2
    real(kind=4) :: visuy1
    real(kind=4) :: visuz2
    real(kind=4) :: visuz1
    real(kind=4) :: vfu
    ! ParallelFortran: Synthesised loop variable decls
    integer :: j_range
    integer :: i_range
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension((-1):(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: sm
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
    real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu1
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu2
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu4
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu3
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu7
    real(kind=4), dimension(0:(ip + 1),0:(jp + 1)), intent(In) :: uspd
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: u
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
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    j_rel = (global_id / i_range)
    j = (j_rel + 1)
    i_rel = (global_id - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
      evsx2=sm(i+1,j,1)
      evsx1=sm(i,j,1)
      evsy2=(dy1(j+1)*((dx1(i+1)*sm(i,j,1)+dx1(i)*sm(i+1,j,1))/(dx1(i)+dx1(i+1)))&
      +dy1(j)*((dx1(i+1)*sm(i,j+1,1)+dx1(i)*sm(i+1,j+1,1))/(dx1(i)+dx1(i+1))))/(dy1(j)+dy1(j+1))
      evsy1=(dy1(j+1)*((dx1(i+1)*sm(i,j-1,1)+dx1(i)*sm(i+1,j-1,1))/(dx1(i)+dx1(i+1)))&
      +dy1(j)*((dx1(i+1)*sm(i,j,1)+dx1(i)*sm(i+1,j,1))/(dx1(i)+dx1(i+1))))/(dy1(j)+dy1(j+1))
      evsz2=(dzn(2)*((dx1(i+1)*sm(i,j,1)+dx1(i)*sm(i+1,j,1))/(dx1(i)+dx1(i+1)))&
      +dzn(1)*((dx1(i+1)*sm(i,j,2)+dx1(i)*sm(i+1,j,2))/(dx1(i)+dx1(i+1))))/(dzn(1)+dzn(2))
      visux2=(evsx2)*2.*diu1(i+1,j  ,1  )
      visux1=(evsx1)*2.*diu1(i  ,j,  1  )
      visuy2=(evsy2)* ( diu2(i  ,j+1,1  )+diu4(i+1,j  ,1 ) )
      visuy1=(evsy1)* ( diu2(i  ,j  ,1  )+diu4(i+1,j-1,1 ) )
      visuz2=(evsz2)* ( diu3(i  ,j  ,2  )+diu7(i+1,j  ,1 ) )
      visuz1=(0.4*uspd(i,j)/alog(0.5*dzn(1)/0.1))**2*(u(i,j,1)/uspd(i,j))
      vfu= (visux2-visux1)/dxs(i)+(visuy2-visuy1)/dy1(j)+(visuz2-visuz1)/dzn(1)
      f(i,j,1)=(f(i,j,1)+vfu)

end subroutine les_map_155


subroutine les_map_175(sm,dy1,dx1,dzn,diu2,diu4,diu5,diu6,diu8,dys,g)

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
    ! Local vars: k,i,j,evsx2,evsx1,evsy2,evsy1,evsz2,evsz1,visvx2,visvx1,visvy2,visvy1,visvz2,visvz1,vfv
    integer :: k
    integer :: i
    integer :: j
    real(kind=4) :: evsx2
    real(kind=4) :: evsx1
    real(kind=4) :: evsy2
    real(kind=4) :: evsy1
    real(kind=4) :: evsz2
    real(kind=4) :: evsz1
    real(kind=4) :: visvx2
    real(kind=4) :: visvx1
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
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu2
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu4
    real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu5
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu6
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu8
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
      visvx2 = (evsx2)* ( diu2(i  ,j+1,k  )+diu4(i+1,j  ,k  ) )
      visvx1 = (evsx1)* ( diu2(i-1,j+1,k  )+diu4(i  ,j  ,k  ) )
      visvy2 = (evsy2)*2.*diu5(i  ,j+1,k  )
      visvy1 = (evsy1)*2.*diu5(i  ,j  ,k  )
      visvz2 = (evsz2)* ( diu6(i  ,j  ,k+1)+diu8(i  ,j+1,k  ) )
      visvz1 = (evsz1)* ( diu6(i  ,j  ,k  )+diu8(i  ,j+1,k-1) )
      vfv = (visvx2-visvx1)/dx1(i) +(visvy2-visvy1)/dys(j) +(visvz2-visvz1)/dzn(k)
      g(i,j,k) = (g(i,j,k)+vfv)

end subroutine les_map_175


subroutine les_map_203(sm,dy1,dx1,dzn,diu2,diu4,diu5,diu6,diu8,vspd,v,dys,g)

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
    ! Local vars: i,j,evsx2,evsx1,evsy2,evsy1,evsz2,visvx2,visvx1,visvy2,visvy1,visvz2,visvz1,vfv
    integer :: i
    integer :: j
    real(kind=4) :: evsx2
    real(kind=4) :: evsx1
    real(kind=4) :: evsy2
    real(kind=4) :: evsy1
    real(kind=4) :: evsz2
    real(kind=4) :: visvx2
    real(kind=4) :: visvx1
    real(kind=4) :: visvy2
    real(kind=4) :: visvy1
    real(kind=4) :: visvz2
    real(kind=4) :: visvz1
    real(kind=4) :: vfv
    ! ParallelFortran: Synthesised loop variable decls
    integer :: j_range
    integer :: i_range
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension((-1):(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: sm
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu2
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu4
    real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu5
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu6
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu8
    real(kind=4), dimension(0:(ip + 1),0:(jp + 1)), intent(In) :: vspd
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: v
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
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    j_rel = (global_id / i_range)
    j = (j_rel + 1)
    i_rel = (global_id - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
      evsy2=sm(i,j+1,1)
      evsy1=sm(i,j,1)
      evsx2=(dy1(j+1)*((dx1(i+1)*sm(i,j,1)+dx1(i)*sm(i+1,j,1))/(dx1(i)+dx1(i+1)))&
      +dy1(j)*((dx1(i+1)*sm(i,j+1,1)+dx1(i)*sm(i+1,j+1,1))/(dx1(i)+dx1(i+1))))/(dy1(j)+dy1(j+1))
      evsx1=(dy1(j+1)*((dx1(i)*sm(i-1,j,1)+dx1(i-1)*sm(i,j,1))/(dx1(i-1)+dx1(i)))&
      +dy1(j)*((dx1(i)*sm(i-1,j+1,1)+dx1(i-1)*sm(i,j+1,1))/(dx1(i-1)+dx1(i))))/(dy1(j)+dy1(j+1))
      evsz2=(dzn(2)*((dx1(i+1)*sm(i,j,1)+dx1(i)*sm(i+1,j,1))/(dx1(i)+dx1(i+1)))&
      +dzn(1)*((dx1(i+1)*sm(i,j,2)+dx1(i)*sm(i+1,j,2))/(dx1(i)+dx1(i+1))))/(dzn(1)+dzn(2))
      visvx2=(evsx2)* ( diu2(i  ,j+1,1  )+diu4(i+1,j  ,1  ) )
      visvx1=(evsx1)* ( diu2(i-1,j+1,1  )+diu4(i  ,j  ,1  ) )
      visvy2=(evsy2)*2.*diu5(i  ,j+1,1  )
      visvy1=(evsy1)*2.*diu5(i  ,j  ,1  )
      visvz2=(evsz2)* ( diu6(i  ,j  ,2  )+diu8(i  ,j+1,1  ) )
      visvz1=(0.4*vspd(i,j)/alog(0.5*dzn(1)/0.1))**2*(v(i,j,1)/vspd(i,j))
      vfv=(visvx2-visvx1)/dx1(i)+(visvy2-visvy1)/dys(j)+(visvz2-visvz1)/dzn(1)
      g(i,j,1)=(g(i,j,1)+vfv)

end subroutine les_map_203


subroutine les_map_223(sm,dzn,dx1,dy1,diu3,diu7,diu6,diu8,diu9,h)

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
    ! Local vars: k,i,j,evsx2,evsx1,evsy2,evsy1,evsz2,evsz1,viswx2,viswx1,viswy2,viswy1,viswz2,viswz1,vfw
    integer :: k
    integer :: i
    integer :: j
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
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu3
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu7
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu6
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu8
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu9
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
      viswx2 = (evsx2)* ( diu3(i  ,j  ,k+1)+diu7(i+1,j  ,k  ) )
      viswx1 = (evsx1)* ( diu3(i-1,j  ,k+1)+diu7(i  ,j  ,k  ) )
      viswy2 = (evsy2)* ( diu6(i  ,j  ,k+1)+diu8(i  ,j+1,k  ) )
      viswy1 = (evsy1)* ( diu6(i  ,j-1,k+1)+diu8(i  ,j  ,k  ) )
      viswz2 = (evsz2)*2.*diu9(i  ,j  ,k+1)
      viswz1 = (evsz1)*2.*diu9(i  ,j  ,k  )
      vfw = (viswx2-viswx1)/dx1(i) +(viswy2-viswy1)/dy1(j) +(viswz2-viswz1)/dzn(k)
      h(i,j,k) = (h(i,j,k)+vfw)

end subroutine les_map_223


subroutine press_map_64(f)

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
    ! Local vars: j,k
    integer :: j
    integer :: k
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: k_rel
    integer :: j_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: f
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    j_range = ((300 - 1) + 1)
    k_rel = (global_id / j_range)
    k = (k_rel + 1)
    j_rel = (global_id - (k_rel * j_range))
    j = (j_rel + 1)


    ! ParallelFortran: Original code
        f( 0,j,k)=f(1  ,j,k)

end subroutine press_map_64


subroutine press_map_69(g)

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
    ! Local vars: k,i
    integer :: k
    integer :: i
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: i_range
    integer :: k_rel
    integer :: i_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: g
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = ((80 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / i_range)
    k = (k_rel + 1)
    i_rel = (global_id - (k_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
        g(i, 0,k)=g(i,jp  ,k)

end subroutine press_map_69


subroutine press_map_74(h)

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
    ! Local vars: j,i
    integer :: j
    integer :: i
    ! ParallelFortran: Synthesised loop variable decls
    integer :: j_range
    integer :: i_range
    integer :: j_rel
    integer :: i_rel
! READ
! WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(Out) :: h
! READ & WRITTEN
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    j_rel = (global_id / i_range)
    j = (j_rel + 1)
    i_rel = (global_id - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
        h(i,j, 0)=0.0
        h(i,j,kp)=0.0

end subroutine press_map_74


subroutine press_map_80(u,dx1,v,dy1,w,dzn,f,g,h,rhs,dt)

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
    ! Local vars: j,k,i
    integer :: j
    integer :: k
    integer :: i
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

end subroutine press_map_80


subroutine press_reduce_93(dx1,dy1,dzn,rhs,global_rhsav_array,global_area_array)

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
    ! Missing args: 
    ! Local vars: j,k,i

    integer :: j
    integer :: k
    integer :: i

    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel

    integer :: local_id
    integer :: local_id_fortran
    integer :: group_id
    integer :: group_id_fortran
    integer :: global_id
    integer :: r_iter
    integer :: local_chunk_size
    integer :: start_position

    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
    real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)), intent(In) :: rhs
    real(kind=4), dimension(1:NUNITS), intent(Out) :: global_rhsav_array
    real(kind=4), dimension(1:NUNITS), intent(Out) :: global_area_array

    real(kind=4) :: local_rhsav
    real(kind=4) :: local_area

    local_chunk_size = (80 * (300 * 300))
    start_position = 0
    local_rhsav = 0
    local_area = 0

    do r_iter=start_position, ((start_position + local_chunk_size) - 1)
        k_range = ((80 - 1) + 1)
        j_range = ((300 - 1) + 1)
        i_range = ((300 - 1) + 1)
        k_rel = (r_iter / (j_range * i_range))
        k = (k_rel + 1)
        j_rel = ((r_iter - (k_rel * (j_range * i_range))) / i_range)
        j = (j_rel + 1)
        i_rel = ((r_iter - (k_rel * (j_range * i_range))) - (j_rel * i_range))
        i = (i_rel + 1)
        local_rhsav = (local_rhsav + (((dx1(i) * dy1(j)) * dzn(k)) * rhs(i,j,k)))
        local_area = (local_area + ((dx1(i) * dy1(j)) * dzn(k)))
    end do

    global_rhsav_array(1) = local_rhsav
    global_area_array(1) = local_area

end subroutine press_reduce_93


subroutine press_map_102(rhs,rhsav)

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
    ! Local vars: j,k,i
    integer :: j
    integer :: k
    integer :: i
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

end subroutine press_map_102


subroutine press_map_112(dzs,dys,dxs,nrd,p,rhs)

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
    ! Local vars: j,k,i,dz1,dz2,cn4s,cn4l,cn3s,cn3l,cn2s,cn2l,cn1,reltmp
    integer :: j
    integer :: k
    integer :: i
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
    real(kind=4), dimension(0:1,0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p
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
                        reltmp = omega*(cn1 *(cn2l*p(0,i+1,j,k) + &
                                 cn2s*p(0,i-1,j,k) +cn3l*p(0,i,j+1,k) + &
                                 cn3s*p(0,i,j-1,k) +cn4l*p(0,i,j,k+1) + &
                                 cn4s*p(0,i,j,k-1) -rhs(i,j,k))-p(0,i,j,k))
                        p(1,i,j,k) = p(0,i,j,k) +reltmp
                      else
                        reltmp = omega*(cn1 *(cn2l*p(1,i+1,j,k) + &
                                 cn2s*p(1,i-1,j,k) +cn3l*p(1,i,j+1,k) + &
                                 cn3s*p(1,i,j-1,k) +cn4l*p(1,i,j,k+1) + &
                                 cn4s*p(1,i,j,k-1) -rhs(i,j,k))-p(1,i,j,k))
                        p(0,i,j,k) = p(1,i,j,k) +reltmp
                      end if

end subroutine press_map_112


subroutine press_map_140(p)

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
    ! Local vars: j,k
    integer :: j
    integer :: k
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: k_rel
    integer :: j_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:1,0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = (((80 + 1) - 0) + 1)
    j_range = (((300 + 1) - 0) + 1)
    k_rel = (global_id / j_range)
    k = (k_rel + 0)
    j_rel = (global_id - (k_rel * j_range))
    j = (j_rel + 0)


    ! ParallelFortran: Original code
            p(0,   0,j,k) = p(0,1 ,j,k)
            p(0,ip+1,j,k) = p(0,ip,j,k)

end subroutine press_map_140


subroutine press_map_146(p)

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
    ! Local vars: k,i
    integer :: k
    integer :: i
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: i_range
    integer :: k_rel
    integer :: i_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:1,0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = (((80 + 1) - 0) + 1)
    i_range = (((300 + 1) - 0) + 1)
    k_rel = (global_id / i_range)
    k = (k_rel + 0)
    i_rel = (global_id - (k_rel * i_range))
    i = (i_rel + 0)


    ! ParallelFortran: Original code
            p(0,i,   0,k) = p(0,i,jp,k)
            p(0,i,jp+1,k) = p(0,i, 1,k)

end subroutine press_map_146


subroutine press_map_153(p)

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
    ! Local vars: j,i
    integer :: j
    integer :: i
    ! ParallelFortran: Synthesised loop variable decls
    integer :: j_range
    integer :: i_range
    integer :: j_rel
    integer :: i_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:1,0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    j_range = (((300 + 1) - 0) + 1)
    i_range = (((300 + 1) - 0) + 1)
    j_rel = (global_id / i_range)
    j = (j_rel + 0)
    i_rel = (global_id - (j_rel * i_range))
    i = (i_rel + 0)


    ! ParallelFortran: Original code
        p(0,i,j,   0) = p(0,i,j,1)
        p(0,i,j,kp+1) = p(0,i,j,kp)

end subroutine press_map_153


subroutine press_reduce_162(p,dx1,dy1,dzn,global_pav_array,global_pco_array)

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
    ! Missing args: 
    ! Local vars: j,k,i

    integer :: j
    integer :: k
    integer :: i

    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: i_range
    integer :: k_rel
    integer :: j_rel
    integer :: i_rel

    integer :: local_id
    integer :: local_id_fortran
    integer :: group_id
    integer :: group_id_fortran
    integer :: global_id
    integer :: r_iter
    integer :: local_chunk_size
    integer :: start_position

    real(kind=4), dimension(0:1,0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(In) :: p
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
    real(kind=4), dimension(1:NUNITS), intent(Out) :: global_pav_array
    real(kind=4), dimension(1:NUNITS), intent(Out) :: global_pco_array

    real(kind=4) :: local_pav
    real(kind=4) :: local_pco

    local_chunk_size = (80 * (300 * 300))
    start_position = 0
    local_pav = 0
    local_pco = 0

    do r_iter=start_position, ((start_position + local_chunk_size) - 1)
        k_range = ((80 - 1) + 1)
        j_range = ((300 - 1) + 1)
        i_range = ((300 - 1) + 1)
        k_rel = (r_iter / (j_range * i_range))
        k = (k_rel + 1)
        j_rel = ((r_iter - (k_rel * (j_range * i_range))) / i_range)
        j = (j_rel + 1)
        i_rel = ((r_iter - (k_rel * (j_range * i_range))) - (j_rel * i_range))
        i = (i_rel + 1)
        local_pav = (local_pav + (((p(0,i,j,k) * dx1(i)) * dy1(j)) * dzn(k)))
        local_pco = (local_pco + ((dx1(i) * dy1(j)) * dzn(k)))
    end do

    global_pav_array(1) = local_pav
    global_pco_array(1) = local_pco

end subroutine press_reduce_162


subroutine press_map_171(p,pav)

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
    ! Local vars: j,k,i
    integer :: j
    integer :: k
    integer :: i
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
    real(kind=4), dimension(0:1,0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p
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
                p(0,i,j,k) = p(0,i,j,k)-pav

end subroutine press_map_171


subroutine press_map_178(p)

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
    ! Local vars: j,k
    integer :: j
    integer :: k
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: j_range
    integer :: k_rel
    integer :: j_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:1,0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = (((80 + 1) - 0) + 1)
    j_range = (((300 + 1) - 0) + 1)
    k_rel = (global_id / j_range)
    k = (k_rel + 0)
    j_rel = (global_id - (k_rel * j_range))
    j = (j_rel + 0)


    ! ParallelFortran: Original code
            p(0,   0,j,k) = p(0,1 ,j,k)
            p(0,ip+1,j,k) = p(0,ip,j,k)

end subroutine press_map_178


subroutine press_map_184(p)

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
    ! Local vars: k,i
    integer :: k
    integer :: i
    ! ParallelFortran: Synthesised loop variable decls
    integer :: k_range
    integer :: i_range
    integer :: k_rel
    integer :: i_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:1,0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    k_range = (((80 + 1) - 0) + 1)
    i_range = (((300 + 1) - 0) + 1)
    k_rel = (global_id / i_range)
    k = (k_rel + 0)
    i_rel = (global_id - (k_rel * i_range))
    i = (i_rel + 0)


    ! ParallelFortran: Original code
            p(0,i,   0,k) = p(0,i,jp,k)
            p(0,i,jp+1,k) = p(0,i, 1,k)

end subroutine press_map_184


subroutine press_map_190(p)

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
    ! Local vars: j,i
    integer :: j
    integer :: i
    ! ParallelFortran: Synthesised loop variable decls
    integer :: j_range
    integer :: i_range
    integer :: j_rel
    integer :: i_rel
! READ
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:1,0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    j_range = (((300 + 1) - 0) + 1)
    i_range = (((300 + 1) - 0) + 1)
    j_rel = (global_id / i_range)
    j = (j_rel + 0)
    i_rel = (global_id - (j_rel * i_range))
    i = (i_rel + 0)


    ! ParallelFortran: Original code
        p(0,i,j,   0) = p(0,i,j,1)
        p(0,i,j,kp+1) = p(0,i,j,kp)

end subroutine press_map_190


subroutine velFG_map_135(u,v,dx1,dy1,uspd,vspd)

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
    ! Local vars: i,j
    integer :: i
    integer :: j
    ! ParallelFortran: Synthesised loop variable decls
    integer :: j_range
    integer :: i_range
    integer :: j_rel
    integer :: i_rel
! READ
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: u
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: v
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
! WRITTEN
    real(kind=4), dimension(0:(ip + 1),0:(jp + 1)), intent(Out) :: uspd
    real(kind=4), dimension(0:(ip + 1),0:(jp + 1)), intent(Out) :: vspd
! READ & WRITTEN
! globalIdDeclaration
    integer :: global_id
! globalIdInitialisation
    call get_global_id(global_id,0)
! ptrAssignments_fseq

    ! ParallelFortran: Synthesised loop variables
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    j_rel = (global_id / i_range)
    j = (j_rel + 1)
    i_rel = (global_id - (j_rel * i_range))
    i = (i_rel + 1)


    ! ParallelFortran: Original code
         uspd(i,j)=(u(i,j,1)**2+((0.5*(v(i,j-1,1)+v(i,j,1))*dx1(i+1)&
     +0.5*(v(i+1,j-1,1)+v(i+1,j,1))*dx1(i))/(dx1(i)+dx1(i+1)))**2)**0.5
         vspd(i,j)=(v(i,j,1)**2+((0.5*(u(i-1,j,1)+u(i,j,1))*dy1(j+1)&
     +0.5*(u(i-1,j+1,1)+u(i,j+1,1))*dy1(j))/(dy1(j)+dy1(j+1)))**2)**0.5

end subroutine velFG_map_135


subroutine velFG_map_148(dx1,cov1,cov2,cov3,cov4,dy1,cov5,cov6,cov7,cov8,dzn,cov9,f,g,h)

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
    ! Local vars: i,j,covx1,k,covy1,covz1,covc
    integer :: i
    integer :: j
    real(kind=4) :: covx1
    integer :: k
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
    real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
    real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov1
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov2
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov3
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov4
    real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
    real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov5
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov6
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov7
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov8
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
    real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov9
! WRITTEN
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(Out) :: f
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(Out) :: g
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(Out) :: h
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
        covx1 = (dx1(i+1)*cov1(i,j,k)+dx1(i)*cov1(i+1,j,k)) /(dx1(i)+dx1(i+1))
        covy1 = (cov2(i,j,k)+cov2(i,j+1,k))/2.
        covz1 = (cov3(i,j,k)+cov3(i,j,k+1))/2.
        covc = covx1+covy1+covz1
        f(i,j,k) = (-covc)
        covx1 = (cov4(i,j,k)+cov4(i+1,j,k))/2.
        covy1 = (dy1(j+1)*cov5(i,j,k)+dy1(j)*cov5(i,j+1,k)) /(dy1(j)+dy1(j+1))
        covz1 = (cov6(i,j,k)+cov6(i,j,k+1))/2.
        covc = covx1+covy1+covz1
        g(i,j,k) = (-covc)
    if ((k < (80 - 1))) then
       covx1 = (cov7(i,j,k)+cov7(i+1,j,k))/2.
       covy1 = (cov8(i,j,k)+cov8(i,j+1,k))/2.
       covz1 = (dzn(k+1)*cov9(i,j,k)+dzn(k)*cov9(i,j,k+1)) /(dzn(k)+dzn(k+1))
       covc = covx1+covy1+covz1
        h(i,j,k) = (-covc)      
    end if

end subroutine velFG_map_148


subroutine velnw_map_36(p,ro,dxs,u,dt,f,dys,v,g,dzs,w,h)

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
    real(kind=4), dimension(0:1,0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(In) :: p
    real(kind=4), dimension(0:ip), intent(In) :: dxs
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In) :: f
    real(kind=4), dimension(0:jp), intent(In) :: dys
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In) :: g
    real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzs
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In) :: h
    real(kind=4) :: ro
    real(kind=4) :: dt
! WRITTEN
! READ & WRITTEN
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: u
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: v
    real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(InOut) :: w
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
        pz = (-p(0,i,j,k)+p(0,i+1,j,k))/ro/dxs(i)
        u(i,j,k) = u(i,j,k)+dt*(f(i,j,k)-pz)
        pz = (-p(0,i,j,k)+p(0,i,j+1,k))/ro/dys(j)
        v(i,j,k) = v(i,j,k)+dt*(g(i,j,k)-pz)
    if ((k < (80 - 1))) then
        pz = (-p(0,i,j,k)+p(0,i,j,k+1))/ro/dzs(k)
        w(i,j,k) = w(i,j,k)+dt*(h(i,j,k)-pz)
    end if

end subroutine velnw_map_36


subroutine adam_bondv1_feedbf_les_press_v_etc_superkernel(f,g,h,fold,gold,hold,z2,dzn,u,v,w,global_aaa_array,global_bbb_array,dt,uout,dxs,kp,usum,bmask1,vsum,cmask1,wsum,dmask1,alpha,beta,jp,ip,fx,fy,fz,dx1,dy1,delx1,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,sm,uspd,dys,vspd,rhs,global_rhsav_array,global_area_array,rhsav,dzs,nrd,p,global_pav_array,global_pco_array,pav,cov1,cov2,cov3,cov4,cov5,cov6,cov7,cov8,cov9,ro,state_ptr)
use module_adam_bondv1_feedbf_les_press_v_etc_superkernel_init
use common_sn
  real(kind=4), dimension(0:(kp + 2)), intent(In) :: z2
  real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzn
  real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: u
  real(kind=4), intent(In), dimension(1) :: dt
  real(kind=4), intent(In), dimension(1) :: uout
  real(kind=4), dimension(0:ip), intent(In) :: dxs
  real, dimension(1) :: kp
  real(kind=4), dimension((-1):(ip + 1),0:(jp + 1),0:(kp + 1)), intent(In) :: bmask1
  real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: v
  real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(In) :: cmask1
  real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)), intent(InOut) :: w
  real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)), intent(In) :: dmask1
  real(kind=4), intent(In), dimension(1) :: alpha
  real(kind=4), intent(In), dimension(1) :: beta
  real, dimension(1) :: jp
  real, dimension(1) :: ip
  real(kind=4), dimension((-1):(ip + 1)), intent(In) :: dx1
  real(kind=4), dimension(0:(jp + 1)), intent(In) :: dy1
  real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu1
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu2
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu3
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu4
  real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu5
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu6
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu7
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu8
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: diu9
  real(kind=4), dimension(kp), intent(InOut) :: delx1
  real(kind=4), dimension((-1):(ip + 1),(-1):(jp + 1),0:(kp + 1)), intent(InOut) :: sm
  real(kind=4), dimension(0:(ip + 1),0:(jp + 1)), intent(InOut) :: uspd
  real(kind=4), dimension(0:jp), intent(In) :: dys
  real(kind=4), dimension(0:(ip + 1),0:(jp + 1)), intent(InOut) :: vspd
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: f
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: g
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: h
  real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)), intent(InOut) :: rhs
  real(kind=4), intent(In), dimension(1) :: rhsav
  real(kind=4), dimension((-1):(kp + 2)), intent(In) :: dzs
  integer, intent(In), dimension(1) :: nrd
  real(kind=4), dimension(0:1,0:(ip + 2),0:(jp + 2),0:(kp + 1)), intent(InOut) :: p
  real(kind=4), intent(In), dimension(1) :: pav
  real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov1
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov2
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov3
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov4
  real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov5
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov6
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov7
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov8
  real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)), intent(In) :: cov9
  real(kind=4), intent(In), dimension(1) :: ro
  real(kind=4), dimension(ip,jp,kp), intent(InOut) :: fold
  real(kind=4), dimension(ip,jp,kp), intent(InOut) :: gold
  real(kind=4), dimension(ip,jp,kp), intent(InOut) :: hold
  real(kind=4), dimension(1:NUNITS), intent(Out) :: global_aaa_array
  real(kind=4), dimension(1:NUNITS), intent(Out) :: global_bbb_array
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: usum
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: vsum
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: wsum
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: fx
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: fy
  real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut) :: fz
  real(kind=4), dimension(1:NUNITS), intent(Out) :: global_rhsav_array
  real(kind=4), dimension(1:NUNITS), intent(Out) :: global_area_array
  real(kind=4), dimension(1:NUNITS), intent(Out) :: global_pav_array
  real(kind=4), dimension(1:NUNITS), intent(Out) :: global_pco_array

  integer :: state
  integer, dimension(1) :: state_ptr
integer, parameter :: ST_ADAM_MAP_36 = 0 !  adam_map_36
integer, parameter :: ST_BONDV1_MAP_38 = 1 !  bondv1_map_38
integer, parameter :: ST_BONDV1_MAP_48 = 2 !  bondv1_map_48
integer, parameter :: ST_BONDV1_MAP_58 = 3 !  bondv1_map_58
integer, parameter :: ST_BONDV1_REDUCE_69 = 4 !  bondv1_reduce_69
integer, parameter :: ST_BONDV1_REDUCE_76 = 5 !  bondv1_reduce_76
integer, parameter :: ST_BONDV1_MAP_83 = 6 !  bondv1_map_83
integer, parameter :: ST_BONDV1_MAP_98 = 7 !  bondv1_map_98
integer, parameter :: ST_BONDV1_MAP_116 = 8 !  bondv1_map_116
integer, parameter :: ST_BONDV1_MAP_128 = 9 !  bondv1_map_128
integer, parameter :: ST_FEEDBF_MAP_35 = 10 !  feedbf_map_35
integer, parameter :: ST_LES_MAP_88 = 11 !  les_map_88
integer, parameter :: ST_LES_MAP_91 = 12 !  les_map_91
integer, parameter :: ST_LES_MAP_109 = 13 !  les_map_109
integer, parameter :: ST_LES_MAP_115 = 14 !  les_map_115
integer, parameter :: ST_LES_MAP_121 = 15 !  les_map_121
integer, parameter :: ST_LES_MAP_127 = 16 !  les_map_127
integer, parameter :: ST_LES_MAP_155 = 17 !  les_map_155
integer, parameter :: ST_LES_MAP_175 = 18 !  les_map_175
integer, parameter :: ST_LES_MAP_203 = 19 !  les_map_203
integer, parameter :: ST_LES_MAP_223 = 20 !  les_map_223
integer, parameter :: ST_PRESS_MAP_64 = 21 !  press_map_64
integer, parameter :: ST_PRESS_MAP_69 = 22 !  press_map_69
integer, parameter :: ST_PRESS_MAP_74 = 23 !  press_map_74
integer, parameter :: ST_PRESS_MAP_80 = 24 !  press_map_80
integer, parameter :: ST_PRESS_REDUCE_93 = 25 !  press_reduce_93
integer, parameter :: ST_PRESS_MAP_102 = 26 !  press_map_102
integer, parameter :: ST_PRESS_MAP_112 = 27 !  press_map_112
integer, parameter :: ST_PRESS_MAP_140 = 28 !  press_map_140
integer, parameter :: ST_PRESS_MAP_146 = 29 !  press_map_146
integer, parameter :: ST_PRESS_MAP_153 = 30 !  press_map_153
integer, parameter :: ST_PRESS_REDUCE_162 = 31 !  press_reduce_162
integer, parameter :: ST_PRESS_MAP_171 = 32 !  press_map_171
integer, parameter :: ST_PRESS_MAP_178 = 33 !  press_map_178
integer, parameter :: ST_PRESS_MAP_184 = 34 !  press_map_184
integer, parameter :: ST_PRESS_MAP_190 = 35 !  press_map_190
integer, parameter :: ST_VELFG_MAP_135 = 36 !  velFG_map_135
integer, parameter :: ST_VELFG_MAP_148 = 37 !  velFG_map_148
integer, parameter :: ST_VELNW_MAP_36 = 38 !  velnw_map_36
  state = state_ptr(1) ! state 
! SUPERKERNEL BODY
  select case(state)
    case (ST_ADAM_MAP_36)
      call adam_map_36(f,g,h,fold,gold,hold)
    case (ST_BONDV1_MAP_38)
      call bondv1_map_38(z2,dzn,u,v,w)
    case (ST_BONDV1_MAP_48)
      call bondv1_map_48(u,v,w)
    case (ST_BONDV1_MAP_58)
      call bondv1_map_58(u,v,w)
    case (ST_BONDV1_REDUCE_69)
      call bondv1_reduce_69(u,global_aaa_array)
    case (ST_BONDV1_REDUCE_76)
      call bondv1_reduce_76(u,global_bbb_array)
    case (ST_BONDV1_MAP_83)
      call bondv1_map_83(u,dt,uout,dxs,v,w)
    case (ST_BONDV1_MAP_98)
      call bondv1_map_98(u,v,w)
    case (ST_BONDV1_MAP_116)
      call bondv1_map_116(u,v)
    case (ST_BONDV1_MAP_128)
      call bondv1_map_128(w)
    case (ST_FEEDBF_MAP_35)
      call feedbf_map_35(kp,usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,jp,ip,f,fx,g,fy,h,fz)
    case (ST_LES_MAP_88)
      call les_map_88(dx1,dy1,dzn,delx1)
    case (ST_LES_MAP_91)
      call les_map_91(diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,delx1,sm)
    case (ST_LES_MAP_109)
      call les_map_109(sm)
    case (ST_LES_MAP_115)
      call les_map_115(sm)
    case (ST_LES_MAP_121)
      call les_map_121(sm)
    case (ST_LES_MAP_127)
      call les_map_127(sm,dy1,dx1,dzn,diu1,diu2,diu4,diu3,diu7,dxs,f)
    case (ST_LES_MAP_155)
      call les_map_155(sm,dy1,dx1,dzn,diu1,diu2,diu4,diu3,diu7,uspd,u,dxs,f)
    case (ST_LES_MAP_175)
      call les_map_175(sm,dy1,dx1,dzn,diu2,diu4,diu5,diu6,diu8,dys,g)
    case (ST_LES_MAP_203)
      call les_map_203(sm,dy1,dx1,dzn,diu2,diu4,diu5,diu6,diu8,vspd,v,dys,g)
    case (ST_LES_MAP_223)
      call les_map_223(sm,dzn,dx1,dy1,diu3,diu7,diu6,diu8,diu9,h)
    case (ST_PRESS_MAP_64)
      call press_map_64(f)
    case (ST_PRESS_MAP_69)
      call press_map_69(g)
    case (ST_PRESS_MAP_74)
      call press_map_74(h)
    case (ST_PRESS_MAP_80)
      call press_map_80(u,dx1,v,dy1,w,dzn,f,g,h,rhs,dt)
    case (ST_PRESS_REDUCE_93)
      call press_reduce_93(dx1,dy1,dzn,rhs,global_rhsav_array,global_area_array)
    case (ST_PRESS_MAP_102)
      call press_map_102(rhs,rhsav)
    case (ST_PRESS_MAP_112)
      call press_map_112(dzs,dys,dxs,nrd,p,rhs)
    case (ST_PRESS_MAP_140)
      call press_map_140(p)
    case (ST_PRESS_MAP_146)
      call press_map_146(p)
    case (ST_PRESS_MAP_153)
      call press_map_153(p)
    case (ST_PRESS_REDUCE_162)
      call press_reduce_162(p,dx1,dy1,dzn,global_pav_array,global_pco_array)
    case (ST_PRESS_MAP_171)
      call press_map_171(p,pav)
    case (ST_PRESS_MAP_178)
      call press_map_178(p)
    case (ST_PRESS_MAP_184)
      call press_map_184(p)
    case (ST_PRESS_MAP_190)
      call press_map_190(p)
    case (ST_VELFG_MAP_135)
      call velFG_map_135(u,v,dx1,dy1,uspd,vspd)
    case (ST_VELFG_MAP_148)
      call velFG_map_148(dx1,cov1,cov2,cov3,cov4,dy1,cov5,cov6,cov7,cov8,dzn,cov9,f,g,h)
    case (ST_VELNW_MAP_36)
      call velnw_map_36(p,ro,dxs,u,dt,f,dys,v,g,dzs,w,h)
  end select
end subroutine adam_bondv1_feedbf_les_press_v_etc_superkernel
end module module_adam_bondv1_feedbf_les_press_v_etc_superkernel