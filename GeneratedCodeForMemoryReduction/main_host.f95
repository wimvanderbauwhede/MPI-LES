program main
    use module_init
    use module_grid
    use module_set
    use module_aveflow
    use module_ifdata
    use module_anime
    use module_velnw
    use module_velFG
    use module_les
    use module_press
    use module_adam

    use oclWrapper
    use module_adam_feedbf_les_press_velfg_ve_etc_superkernel_init
      use params_common_sn
      implicit none
    
! Original declarations
      real(4) :: alpha !!
      integer :: ical !!
      integer :: im !!
      integer :: jm !!
      integer :: km !!
      integer :: n !!
      integer :: n0 !!
      integer :: n1 !!
      integer :: nif !!
      integer :: nmax !!
      real(4) :: beta !!
      character*70 :: data10 !!
      character*70 :: data11 !!
      character*70 :: data12 !!
      character*70 :: data13 !!
      character*70 :: data14 !!
      character*70 :: data15 !!
      character*70 :: data20 !!
      character*70 :: data21 !!
      character*70 :: data22 !!
      character*70 :: data23 !!
      character*70 :: data24 !!
      character*70 :: data25 !!
      character*70 :: data26 !!
      character*70 :: data27 !!
      character*70 :: data30 !!
      character*70 :: data31 !!
      real(4) :: dt !!
      real(4) :: ro !!
      real(4) :: time !!
      real(4) :: vn !!
      real(4), dimension(0:ip+1,0:jp+1,0:kp+1) :: amask1 !!
      real(4), dimension(-1:ip+1,0:jp+1,0:kp+1) :: bmask1 !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1) :: cmask1 !!
      real(4), dimension(0:ip+1,0:jp+1,0:kp+1) :: dmask1 !!
      real(4), dimension(ip,jp,kp) :: cn1 !!
      real(4), dimension(ip) :: cn2l !!
      real(4), dimension(ip) :: cn2s !!
      real(4), dimension(jp) :: cn3l !!
      real(4), dimension(jp) :: cn3s !!
      real(4), dimension(kp) :: cn4l !!
      real(4), dimension(kp) :: cn4s !!
      real(4), dimension(kp) :: delx1 !!
      real(4), dimension(-1:ip+1) :: dx1 !!
      real(4), dimension(0:ip) :: dxl !!
      real(4), dimension(0:ip) :: dxs !!
      real(4), dimension(0:jp+1) :: dy1 !!
      real(4), dimension(0:jp) :: dyl !!
      real(4), dimension(0:jp) :: dys !!
      real(4), dimension(-1:kp+2) :: dzn !!
      real(4), dimension(-1:kp+2) :: dzs !!
      real(4), dimension(0:ip,0:jp,0:kp) :: f !!
      real(4), dimension(ip,jp,kp) :: fold !!
      real(4), dimension(0:ip,0:jp,0:kp) :: fx !!
      real(4), dimension(0:ip,0:jp,0:kp) :: fy !!
      real(4), dimension(0:ip,0:jp,0:kp) :: fz !!
      real(4), dimension(0:ip,0:jp,0:kp) :: g !!
      real(4), dimension(ip,jp,kp) :: gold !!
      real(4), dimension(0:ip,0:jp,0:kp) :: h !!
      real(4), dimension(ip,jp,kp) :: hold !!
      real(4), dimension(0:1,0:ip+2,0:jp+2,0:kp+1) :: p !!
      real(4), dimension(0:ip+1,0:jp+1,0:kp+1) :: rhs !!
      real(4), dimension(-1:ip+1,-1:jp+1,0:kp+1) :: sm !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1) :: u !!
      real(4), dimension(0:ip,0:jp,0:kp) :: usum !!
      real(4), dimension(ip,jp,kp) :: uwfx !!
      real(4), dimension(ip,kp) :: uwfxs !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1) :: v !!
      real(4), dimension(0:ip,0:jp,0:kp) :: vsum !!
      real(4), dimension(0:ip+1,-1:jp+1,-1:kp+1) :: w !!
      real(4), dimension(0:ip,0:jp,0:kp) :: wsum !!
      real(4), dimension(0:kp+2) :: z2 !!
      real(4), dimension(-1:ipmax+1,-1:jpmax+1) :: zbm !!
      integer :: clock_rate !!
      integer(4), dimension(0:9) :: timestamp !!
      integer(4) :: i !!



!

! otherStatements

! remainingDecls

    ! Extra declarations
!NOTHING!
    ! Buffer declarations
      integer(8) :: state_ptr_buf
      integer(8) :: ro_buf
      integer(8) :: dys_buf
      integer(8) :: dxs_buf
      integer(8) :: dt_buf
      integer(8) :: f_buf
      integer(8) :: g_buf
      integer(8) :: dzs_buf
      integer(8) :: h_buf
      integer(8) :: dx1_buf
      integer(8) :: dy1_buf
      integer(8) :: dzn_buf
      integer(8) :: delx1_buf
      integer(8) :: fold_buf
      integer(8) :: gold_buf
      integer(8) :: hold_buf
      integer(8) :: usum_buf
      integer(8) :: vsum_buf
      integer(8) :: wsum_buf
      integer(8) :: fx_buf
      integer(8) :: fy_buf
      integer(8) :: fz_buf
      integer(8) :: sm_buf

    integer, dimension(1) :: state_ptr

    ! Size declarations
      integer, dimension(1) :: state_ptr_sz
      integer, dimension(1) :: ro_ptr_sz
      integer, dimension(1) :: dys_sz
      integer, dimension(1) :: dxs_sz
      integer, dimension(1) :: dt_ptr_sz
      integer, dimension(3) :: f_sz
      integer, dimension(3) :: g_sz
      integer, dimension(1) :: dzs_sz
      integer, dimension(3) :: h_sz
      integer, dimension(1) :: dx1_sz
      integer, dimension(1) :: dy1_sz
      integer, dimension(1) :: dzn_sz
      integer, dimension(1) :: delx1_sz
      integer, dimension(3) :: fold_sz
      integer, dimension(3) :: gold_sz
      integer, dimension(3) :: hold_sz
      integer, dimension(3) :: usum_sz
      integer, dimension(3) :: vsum_sz
      integer, dimension(3) :: wsum_sz
      integer, dimension(3) :: fx_sz
      integer, dimension(3) :: fy_sz
      integer, dimension(3) :: fz_sz
      integer, dimension(3) :: sm_sz
      real(kind=4), dimension(1) :: ro_ptr
      real(kind=4), dimension(1) :: dt_ptr

    call adam_feedbf_les_press_velfg_ve_etc_superkernel_init()
    
! Size assignments
      state_ptr_sz = shape(state_ptr)
      ro_ptr_sz = shape(ro_ptr)
      dys_sz = shape(dys)
      dxs_sz = shape(dxs)
      dt_ptr_sz = shape(dt_ptr)
      f_sz = shape(f)
      g_sz = shape(g)
      dzs_sz = shape(dzs)
      h_sz = shape(h)
      dx1_sz = shape(dx1)
      dy1_sz = shape(dy1)
      dzn_sz = shape(dzn)
      delx1_sz = shape(delx1)
      fold_sz = shape(fold)
      gold_sz = shape(gold)
      hold_sz = shape(hold)
      usum_sz = shape(usum)
      vsum_sz = shape(vsum)
      wsum_sz = shape(wsum)
      fx_sz = shape(fx)
      fy_sz = shape(fy)
      fz_sz = shape(fz)
      sm_sz = shape(sm)

    ! Buffer loads
      call oclLoadBuffer(STATE_PTR_BUF_IDX, state_ptr_buf)
      call oclLoadBuffer(RO_BUF_IDX, ro_buf)
      call oclLoadBuffer(DYS_BUF_IDX, dys_buf)
      call oclLoadBuffer(DXS_BUF_IDX, dxs_buf)
      call oclLoadBuffer(DT_BUF_IDX, dt_buf)
      call oclLoadBuffer(F_BUF_IDX, f_buf)
      call oclLoadBuffer(G_BUF_IDX, g_buf)
      call oclLoadBuffer(DZS_BUF_IDX, dzs_buf)
      call oclLoadBuffer(H_BUF_IDX, h_buf)
      call oclLoadBuffer(DX1_BUF_IDX, dx1_buf)
      call oclLoadBuffer(DY1_BUF_IDX, dy1_buf)
      call oclLoadBuffer(DZN_BUF_IDX, dzn_buf)
      call oclLoadBuffer(DELX1_BUF_IDX, delx1_buf)
      call oclLoadBuffer(FOLD_BUF_IDX, fold_buf)
      call oclLoadBuffer(GOLD_BUF_IDX, gold_buf)
      call oclLoadBuffer(HOLD_BUF_IDX, hold_buf)
      call oclLoadBuffer(USUM_BUF_IDX, usum_buf)
      call oclLoadBuffer(VSUM_BUF_IDX, vsum_buf)
      call oclLoadBuffer(WSUM_BUF_IDX, wsum_buf)
      call oclLoadBuffer(FX_BUF_IDX, fx_buf)
      call oclLoadBuffer(FY_BUF_IDX, fy_buf)
      call oclLoadBuffer(FZ_BUF_IDX, fz_buf)
      call oclLoadBuffer(SM_BUF_IDX, sm_buf)

    ! Original code with buffer writes and reads
    call set(data10,data11,data20,data21,data22,data23,data24,data25,data26,&
             data27,data30,data31,ical,nif,n0,n1,nmax,dt,ro,&
             vn,alpha,beta,data12,data13,data14,data15)
    ro_ptr(1) = ro
    call oclWrite1DFloatArrayBuffer(ro_buf,ro_ptr_sz,ro_ptr)! Automatic conversion to array
    call grid(dx1,dxl,dy1,dyl,z2,dzn,dzs,dxs,dys)
    call init(u,v,w,p,cn2s,dxs,cn2l,cn3s,dys,cn3l,dzs,cn4s,cn4l,cn1,&
              amask1,bmask1,cmask1,dmask1, &
              zbm,z2,dzn)
    call oclWrite1DFloatArrayBuffer(dys_buf,dys_sz,dys)
    call ifdata( &
                fold,gold,hold,&
                 time, &
                n,u,v,w,p,usum,vsum,wsum,delx1,dx1,dy1,dzn,&
                sm,f,g,h,z2,dt,dxs,vn, &
                dzs,&
                amask1,bmask1,cmask1,dmask1,&
                alpha,beta,fx,fy,fz,zbm,ical,nif)
    call oclWrite1DFloatArrayBuffer(dxs_buf,dxs_sz,dxs)
    dt_ptr(1) = dt
    call oclWrite1DFloatArrayBuffer(dt_buf,dt_ptr_sz,dt_ptr)! Automatic conversion to array
    call oclWrite3DFloatArrayBuffer(f_buf,f_sz,f)
    call oclWrite3DFloatArrayBuffer(g_buf,g_sz,g)
    call oclWrite1DFloatArrayBuffer(dzs_buf,dzs_sz,dzs)
    call oclWrite3DFloatArrayBuffer(h_buf,h_sz,h)
    call oclWrite1DFloatArrayBuffer(dx1_buf,dx1_sz,dx1)
    call oclWrite1DFloatArrayBuffer(dy1_buf,dy1_sz,dy1)
    call oclWrite1DFloatArrayBuffer(dzn_buf,dzn_sz,dzn)
    call oclWrite1DFloatArrayBuffer(delx1_buf,delx1_sz,delx1)
    call oclWrite3DFloatArrayBuffer(fold_buf,fold_sz,fold)
    call oclWrite3DFloatArrayBuffer(gold_buf,gold_sz,gold)
    call oclWrite3DFloatArrayBuffer(hold_buf,hold_sz,hold)
    call oclWrite3DFloatArrayBuffer(usum_buf,usum_sz,usum)
    call oclWrite3DFloatArrayBuffer(vsum_buf,vsum_sz,vsum)
    call oclWrite3DFloatArrayBuffer(wsum_buf,wsum_sz,wsum)
    call oclWrite3DFloatArrayBuffer(fx_buf,fx_sz,fx)
    call oclWrite3DFloatArrayBuffer(fy_buf,fy_sz,fy)
    call oclWrite3DFloatArrayBuffer(fz_buf,fz_sz,fz)
    call oclWrite3DFloatArrayBuffer(sm_buf,sm_sz,sm)
#ifdef TIMINGS
!    nmax=201
    call system_clock(timestamp(8), clock_rate)
#endif

    do n = n0,nmax
        time = float(n-n0)*dt
#ifdef TIMINGS
        print *, 'run_LES_reference: time step = ',n
        call system_clock(timestamp(0), clock_rate)
#endif

        call velnw(p,ro,dxs,u,dt,f,dys,v,g,dzs,w,h) 
#ifdef TIMINGS
        call system_clock(timestamp(1), clock_rate)
#endif

#ifdef TIMINGS
        call system_clock(timestamp(2), clock_rate)
#endif

      call velfg(dx1,dy1,dzn,f, &
      g,dzs,h,u,v,w)
#ifdef TIMINGS
        call system_clock(timestamp(3), clock_rate)
#endif

#ifdef TIMINGS
        call system_clock(timestamp(4), clock_rate)
#endif

        call les(u,v,w,f,g,h,sm,dx1,dy1,dzn,dzs,dxs,dys)
#ifdef TIMINGS
        call system_clock(timestamp(5), clock_rate)
#endif

        call adam(n,nmax,data21,fold,gold,hold,&
        f,g,h) 
#ifdef TIMINGS
        call system_clock(timestamp(6), clock_rate)
#endif

        call press(u,v,w,p,rhs,f,g,h,dx1,dy1,dzn,dxs,dys,dzs,dt,n,nmax&
        )
#ifdef TIMINGS
        call system_clock(timestamp(7), clock_rate)
        do i=1, 7
            print '("Time for state ",i2," = ",f6.3," s")',i, &
                  (timestamp(i)-timestamp(i-1))/ real(clock_rate)
        end do
#endif

        call anime(n,n0,n1,&
                    u,w,v,&
                    p(0,:,:,:) &
                    ) 
     end do
#ifdef TIMINGS
    call system_clock(timestamp(9))
    print *,"Total time:" ,(timestamp(9)-timestamp(8))/real(clock_rate), &
          "s for ",nmax-n0,"iterations"
#endif

    end program main
