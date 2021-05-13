module module_adam_feedbf_les_press_velfg_ve_etc_superkernel_init

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
        integer, parameter ::ALPHA_BUF_IDX = 16
        integer, parameter ::BETA_BUF_IDX = 18
        integer, parameter ::BMASK1_BUF_IDX = 9
        integer, parameter ::CMASK1_BUF_IDX = 12
        integer, parameter ::DELX1_BUF_IDX = 26
        integer, parameter ::DMASK1_BUF_IDX = 15
        integer, parameter ::DT_BUF_IDX = 17
        integer, parameter ::DX1_BUF_IDX = 22
        integer, parameter ::DXS_BUF_IDX = 28
        integer, parameter ::DY1_BUF_IDX = 23
        integer, parameter ::DYS_BUF_IDX = 29
        integer, parameter ::DZN_BUF_IDX = 25
        integer, parameter ::DZS_BUF_IDX = 24
        integer, parameter ::F_BUF_IDX = 1
        integer, parameter ::FOLD_BUF_IDX = 4
        integer, parameter ::FX_BUF_IDX = 19
        integer, parameter ::FY_BUF_IDX = 20
        integer, parameter ::FZ_BUF_IDX = 21
        integer, parameter ::G_BUF_IDX = 2
        integer, parameter ::GLOBAL_AREA_ARRAY_BUF_IDX = 32
        integer, parameter ::GLOBAL_PAV_ARRAY_BUF_IDX = 37
        integer, parameter ::GLOBAL_PCO_ARRAY_BUF_IDX = 38
        integer, parameter ::GLOBAL_RHSAV_ARRAY_BUF_IDX = 31
        integer, parameter ::GOLD_BUF_IDX = 5
        integer, parameter ::H_BUF_IDX = 3
        integer, parameter ::HOLD_BUF_IDX = 6
        integer, parameter ::NRD_BUF_IDX = 34
        integer, parameter ::P0_BUF_IDX = 35
        integer, parameter ::P1_BUF_IDX = 36
        integer, parameter ::PAV_BUF_IDX = 39
        integer, parameter ::RHS_BUF_IDX = 30
        integer, parameter ::RHSAV_BUF_IDX = 33
        integer, parameter ::RO_BUF_IDX = 40
        integer, parameter ::SM_BUF_IDX = 27
        integer, parameter ::STATE_PTR_BUF_IDX = 41
        integer, parameter ::U_BUF_IDX = 8
        integer, parameter ::USUM_BUF_IDX = 7
        integer, parameter ::V_BUF_IDX = 11
        integer, parameter ::VSUM_BUF_IDX = 10
        integer, parameter ::W_BUF_IDX = 14
        integer, parameter ::WSUM_BUF_IDX = 13

contains

    subroutine adam_feedbf_les_press_velfg_ve_etc_superkernel_init()

        use oclWrapper
        character(len=*), parameter :: srcstr = "module_adam_feedbf_les_press_velfg_ve_etc_superkernel.cl"
        character(len=*), parameter :: kstr   = "adam_feedbf_les_press_velfg_ve_etc_superkernel"
! parameters
              integer, parameter :: kp = 80 
              integer, parameter :: ip = 300 
              integer, parameter :: jp = 300 
              integer, parameter :: ipmax = ip 
              integer, parameter :: jpmax = jp 
              real, parameter :: dxgrid = 4. 
              real, parameter :: dygrid = 4. 
              real, parameter :: cs0 = 0.14 
              integer, parameter :: i_anime = 1 
              integer, parameter :: avetime = 2 
              integer, parameter :: km_sl = 80 
              integer, parameter :: i_aveflow = 0 
              integer, parameter :: i_ifdata_out = 0 
              real, parameter :: dt_orig = 0.05 
              real, parameter :: pjuge = 0.0001 
              integer, parameter :: nmaxp = 50 
              real, parameter :: omega = 1. 
              integer, parameter :: u0 = 0 
! declarations
        real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)) :: u
        real(kind=4), dimension((-1):(ip + 1),0:(jp + 1),0:(kp + 1)) :: bmask1
        real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)) :: v
        real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)) :: cmask1
        real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)) :: w
        real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)) :: dmask1
        real(kind=4) :: alpha
        real(kind=4) :: dt
        real(kind=4) :: beta
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: fx
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: fy
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: fz
        real(kind=4), dimension((-1):(ip + 1)) :: dx1
        real(kind=4), dimension(0:(jp + 1)) :: dy1
        real(kind=4), dimension((-1):(kp + 2)) :: dzs
        real(kind=4), dimension((-1):(kp + 2)) :: dzn
        real(kind=4), dimension(kp) :: delx1
        real(kind=4), dimension((-1):(ip + 1),(-1):(jp + 1),0:(kp + 1)) :: sm
        real(kind=4), dimension(0:ip) :: dxs
        real(kind=4), dimension(0:jp) :: dys
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: f
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: g
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: h
        real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)) :: rhs
        real(kind=4) :: rhsav
        integer :: nrd
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 1)) :: p0
        real(kind=4) :: pav
        real(kind=4) :: ro
        real(kind=4), dimension(ip,jp,kp) :: fold
        real(kind=4), dimension(ip,jp,kp) :: gold
        real(kind=4), dimension(ip,jp,kp) :: hold
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: usum
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: vsum
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: wsum
        real(kind=4), dimension(1:NUNITS) :: global_rhsav_array
        real(kind=4), dimension(1:NUNITS) :: global_area_array
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 1)) :: p1
        real(kind=4), dimension(1:NUNITS) :: global_pav_array
        real(kind=4), dimension(1:NUNITS) :: global_pco_array
        integer, dimension(1) :: state_ptr
! buffer declarations
        integer(8) :: u_buf
        integer(8) :: bmask1_buf
        integer(8) :: v_buf
        integer(8) :: cmask1_buf
        integer(8) :: w_buf
        integer(8) :: dmask1_buf
        integer(8) :: alpha_buf
        integer(8) :: dt_buf
        integer(8) :: beta_buf
        integer(8) :: fx_buf
        integer(8) :: fy_buf
        integer(8) :: fz_buf
        integer(8) :: dx1_buf
        integer(8) :: dy1_buf
        integer(8) :: dzs_buf
        integer(8) :: dzn_buf
        integer(8) :: delx1_buf
        integer(8) :: sm_buf
        integer(8) :: dxs_buf
        integer(8) :: dys_buf
        integer(8) :: f_buf
        integer(8) :: g_buf
        integer(8) :: h_buf
        integer(8) :: rhs_buf
        integer(8) :: rhsav_buf
        integer(8) :: nrd_buf
        integer(8) :: p0_buf
        integer(8) :: pav_buf
        integer(8) :: ro_buf
        integer(8) :: fold_buf
        integer(8) :: gold_buf
        integer(8) :: hold_buf
        integer(8) :: usum_buf
        integer(8) :: vsum_buf
        integer(8) :: wsum_buf
        integer(8) :: global_rhsav_array_buf
        integer(8) :: global_area_array_buf
        integer(8) :: p1_buf
        integer(8) :: global_pav_array_buf
        integer(8) :: global_pco_array_buf
        integer(8) :: state_ptr_buf
        integer, dimension(3) :: u_sz
        integer, dimension(3) :: bmask1_sz
        integer, dimension(3) :: v_sz
        integer, dimension(3) :: cmask1_sz
        integer, dimension(3) :: w_sz
        integer, dimension(3) :: dmask1_sz
        integer, dimension(1) :: alpha_ptr_sz
        integer, dimension(1) :: dt_ptr_sz
        integer, dimension(1) :: beta_ptr_sz
        integer, dimension(3) :: fx_sz
        integer, dimension(3) :: fy_sz
        integer, dimension(3) :: fz_sz
        integer, dimension(1) :: dx1_sz
        integer, dimension(1) :: dy1_sz
        integer, dimension(1) :: dzs_sz
        integer, dimension(1) :: dzn_sz
        integer, dimension(1) :: delx1_sz
        integer, dimension(3) :: sm_sz
        integer, dimension(1) :: dxs_sz
        integer, dimension(1) :: dys_sz
        integer, dimension(3) :: f_sz
        integer, dimension(3) :: g_sz
        integer, dimension(3) :: h_sz
        integer, dimension(3) :: rhs_sz
        integer, dimension(1) :: rhsav_ptr_sz
        integer, dimension(1) :: nrd_ptr_sz
        integer, dimension(3) :: p0_sz
        integer, dimension(1) :: pav_ptr_sz
        integer, dimension(1) :: ro_ptr_sz
        integer, dimension(3) :: fold_sz
        integer, dimension(3) :: gold_sz
        integer, dimension(3) :: hold_sz
        integer, dimension(3) :: usum_sz
        integer, dimension(3) :: vsum_sz
        integer, dimension(3) :: wsum_sz
        integer, dimension(1) :: global_rhsav_array_sz
        integer, dimension(1) :: global_area_array_sz
        integer, dimension(3) :: p1_sz
        integer, dimension(1) :: global_pav_array_sz
        integer, dimension(1) :: global_pco_array_sz
        integer, dimension(1) :: state_ptr_sz
        real(kind=4), dimension(1) :: alpha_ptr
        real(kind=4), dimension(1) :: dt_ptr
        real(kind=4), dimension(1) :: beta_ptr
        real(kind=4), dimension(1) :: rhsav_ptr
        integer, dimension(1) :: nrd_ptr
        real(kind=4), dimension(1) :: pav_ptr
        real(kind=4), dimension(1) :: ro_ptr

        call oclInit(srcstr,kstr)

        u_sz = shape(u)
        bmask1_sz = shape(bmask1)
        v_sz = shape(v)
        cmask1_sz = shape(cmask1)
        w_sz = shape(w)
        dmask1_sz = shape(dmask1)
        alpha_ptr_sz = shape(alpha_ptr)
        dt_ptr_sz = shape(dt_ptr)
        beta_ptr_sz = shape(beta_ptr)
        fx_sz = shape(fx)
        fy_sz = shape(fy)
        fz_sz = shape(fz)
        dx1_sz = shape(dx1)
        dy1_sz = shape(dy1)
        dzs_sz = shape(dzs)
        dzn_sz = shape(dzn)
        delx1_sz = shape(delx1)
        sm_sz = shape(sm)
        dxs_sz = shape(dxs)
        dys_sz = shape(dys)
        f_sz = shape(f)
        g_sz = shape(g)
        h_sz = shape(h)
        rhs_sz = shape(rhs)
        rhsav_ptr_sz = shape(rhsav_ptr)
        nrd_ptr_sz = shape(nrd_ptr)
        p0_sz = shape(p0)
        pav_ptr_sz = shape(pav_ptr)
        ro_ptr_sz = shape(ro_ptr)
        fold_sz = shape(fold)
        gold_sz = shape(gold)
        hold_sz = shape(hold)
        usum_sz = shape(usum)
        vsum_sz = shape(vsum)
        wsum_sz = shape(wsum)
        global_rhsav_array_sz = shape(global_rhsav_array)
        global_area_array_sz = shape(global_area_array)
        p1_sz = shape(p1)
        global_pav_array_sz = shape(global_pav_array)
        global_pco_array_sz = shape(global_pco_array)
        state_ptr_sz = shape(state_ptr)

        call oclMake3DFloatArrayReadWriteBuffer(u_buf,u_sz,u)
        call oclMake3DFloatArrayReadWriteBuffer(bmask1_buf,bmask1_sz,bmask1)
        call oclMake3DFloatArrayReadWriteBuffer(v_buf,v_sz,v)
        call oclMake3DFloatArrayReadWriteBuffer(cmask1_buf,cmask1_sz,cmask1)
        call oclMake3DFloatArrayReadWriteBuffer(w_buf,w_sz,w)
        call oclMake3DFloatArrayReadWriteBuffer(dmask1_buf,dmask1_sz,dmask1)
        call oclMake1DFloatArrayReadWriteBuffer(alpha_buf,alpha_ptr_sz,alpha_ptr)! Automatic conversion to array
        call oclMake1DFloatArrayReadWriteBuffer(dt_buf,dt_ptr_sz,dt_ptr)! Automatic conversion to array
        call oclMake1DFloatArrayReadWriteBuffer(beta_buf,beta_ptr_sz,beta_ptr)! Automatic conversion to array
        call oclMake3DFloatArrayReadWriteBuffer(fx_buf,fx_sz,fx)
        call oclMake3DFloatArrayReadWriteBuffer(fy_buf,fy_sz,fy)
        call oclMake3DFloatArrayReadWriteBuffer(fz_buf,fz_sz,fz)
        call oclMake1DFloatArrayReadWriteBuffer(dx1_buf,dx1_sz,dx1)
        call oclMake1DFloatArrayReadWriteBuffer(dy1_buf,dy1_sz,dy1)
        call oclMake1DFloatArrayReadWriteBuffer(dzs_buf,dzs_sz,dzs)
        call oclMake1DFloatArrayReadWriteBuffer(dzn_buf,dzn_sz,dzn)
        call oclMake1DFloatArrayReadWriteBuffer(delx1_buf,delx1_sz,delx1)
        call oclMake3DFloatArrayReadWriteBuffer(sm_buf,sm_sz,sm)
        call oclMake1DFloatArrayReadWriteBuffer(dxs_buf,dxs_sz,dxs)
        call oclMake1DFloatArrayReadWriteBuffer(dys_buf,dys_sz,dys)
        call oclMake3DFloatArrayReadWriteBuffer(f_buf,f_sz,f)
        call oclMake3DFloatArrayReadWriteBuffer(g_buf,g_sz,g)
        call oclMake3DFloatArrayReadWriteBuffer(h_buf,h_sz,h)
        call oclMake3DFloatArrayReadWriteBuffer(rhs_buf,rhs_sz,rhs)
        call oclMake1DFloatArrayReadWriteBuffer(rhsav_buf,rhsav_ptr_sz,rhsav_ptr)! Automatic conversion to array
        call oclMake1DIntArrayReadWriteBuffer(nrd_buf,nrd_ptr_sz,nrd_ptr)! Automatic conversion to array
        call oclMake3DFloatArrayReadWriteBuffer(p0_buf,p0_sz,p0)
        call oclMake1DFloatArrayReadWriteBuffer(pav_buf,pav_ptr_sz,pav_ptr)! Automatic conversion to array
        call oclMake1DFloatArrayReadWriteBuffer(ro_buf,ro_ptr_sz,ro_ptr)! Automatic conversion to array
        call oclMake3DFloatArrayReadWriteBuffer(fold_buf,fold_sz,fold)
        call oclMake3DFloatArrayReadWriteBuffer(gold_buf,gold_sz,gold)
        call oclMake3DFloatArrayReadWriteBuffer(hold_buf,hold_sz,hold)
        call oclMake3DFloatArrayReadWriteBuffer(usum_buf,usum_sz,usum)
        call oclMake3DFloatArrayReadWriteBuffer(vsum_buf,vsum_sz,vsum)
        call oclMake3DFloatArrayReadWriteBuffer(wsum_buf,wsum_sz,wsum)
        call oclMake1DFloatArrayReadWriteBuffer(global_rhsav_array_buf,global_rhsav_array_sz,global_rhsav_array)
        call oclMake1DFloatArrayReadWriteBuffer(global_area_array_buf,global_area_array_sz,global_area_array)
        call oclMake3DFloatArrayReadWriteBuffer(p1_buf,p1_sz,p1)
        call oclMake1DFloatArrayReadWriteBuffer(global_pav_array_buf,global_pav_array_sz,global_pav_array)
        call oclMake1DFloatArrayReadWriteBuffer(global_pco_array_buf,global_pco_array_sz,global_pco_array)
        call oclMake1DIntArrayReadWriteBuffer(state_ptr_buf,state_ptr_sz,state_ptr)

        call oclSetFloatArrayArg(7, u_buf)
        call oclSetFloatArrayArg(8, bmask1_buf)
        call oclSetFloatArrayArg(10, v_buf)
        call oclSetFloatArrayArg(11, cmask1_buf)
        call oclSetFloatArrayArg(13, w_buf)
        call oclSetFloatArrayArg(14, dmask1_buf)
        call oclSetFloatArrayArg(15, alpha_buf)
        call oclSetFloatArrayArg(16, dt_buf)
        call oclSetFloatArrayArg(17, beta_buf)
        call oclSetFloatArrayArg(18, fx_buf)
        call oclSetFloatArrayArg(19, fy_buf)
        call oclSetFloatArrayArg(20, fz_buf)
        call oclSetFloatArrayArg(21, dx1_buf)
        call oclSetFloatArrayArg(22, dy1_buf)
        call oclSetFloatArrayArg(23, dzs_buf)
        call oclSetFloatArrayArg(24, dzn_buf)
        call oclSetFloatArrayArg(25, delx1_buf)
        call oclSetFloatArrayArg(26, sm_buf)
        call oclSetFloatArrayArg(27, dxs_buf)
        call oclSetFloatArrayArg(28, dys_buf)
        call oclSetFloatArrayArg(0, f_buf)
        call oclSetFloatArrayArg(1, g_buf)
        call oclSetFloatArrayArg(2, h_buf)
        call oclSetFloatArrayArg(29, rhs_buf)
        call oclSetFloatArrayArg(32, rhsav_buf)
        call oclSetIntArrayArg(33, nrd_buf)
        call oclSetFloatArrayArg(34, p0_buf)
        call oclSetFloatArrayArg(38, pav_buf)
        call oclSetFloatArrayArg(39, ro_buf)
        call oclSetFloatArrayArg(3, fold_buf)
        call oclSetFloatArrayArg(4, gold_buf)
        call oclSetFloatArrayArg(5, hold_buf)
        call oclSetFloatArrayArg(6, usum_buf)
        call oclSetFloatArrayArg(9, vsum_buf)
        call oclSetFloatArrayArg(12, wsum_buf)
        call oclSetFloatArrayArg(30, global_rhsav_array_buf)
        call oclSetFloatArrayArg(31, global_area_array_buf)
        call oclSetFloatArrayArg(35, p1_buf)
        call oclSetFloatArrayArg(36, global_pav_array_buf)
        call oclSetFloatArrayArg(37, global_pco_array_buf)
        call oclSetIntArrayArg(40, state_ptr_buf)

        call oclStoreBuffer(ALPHA_BUF_IDX, alpha_buf)
        call oclStoreBuffer(BETA_BUF_IDX, beta_buf)
        call oclStoreBuffer(BMASK1_BUF_IDX, bmask1_buf)
        call oclStoreBuffer(CMASK1_BUF_IDX, cmask1_buf)
        call oclStoreBuffer(DELX1_BUF_IDX, delx1_buf)
        call oclStoreBuffer(DMASK1_BUF_IDX, dmask1_buf)
        call oclStoreBuffer(DT_BUF_IDX, dt_buf)
        call oclStoreBuffer(DX1_BUF_IDX, dx1_buf)
        call oclStoreBuffer(DXS_BUF_IDX, dxs_buf)
        call oclStoreBuffer(DY1_BUF_IDX, dy1_buf)
        call oclStoreBuffer(DYS_BUF_IDX, dys_buf)
        call oclStoreBuffer(DZN_BUF_IDX, dzn_buf)
        call oclStoreBuffer(DZS_BUF_IDX, dzs_buf)
        call oclStoreBuffer(F_BUF_IDX, f_buf)
        call oclStoreBuffer(FOLD_BUF_IDX, fold_buf)
        call oclStoreBuffer(FX_BUF_IDX, fx_buf)
        call oclStoreBuffer(FY_BUF_IDX, fy_buf)
        call oclStoreBuffer(FZ_BUF_IDX, fz_buf)
        call oclStoreBuffer(G_BUF_IDX, g_buf)
        call oclStoreBuffer(GLOBAL_AREA_ARRAY_BUF_IDX, global_area_array_buf)
        call oclStoreBuffer(GLOBAL_PAV_ARRAY_BUF_IDX, global_pav_array_buf)
        call oclStoreBuffer(GLOBAL_PCO_ARRAY_BUF_IDX, global_pco_array_buf)
        call oclStoreBuffer(GLOBAL_RHSAV_ARRAY_BUF_IDX, global_rhsav_array_buf)
        call oclStoreBuffer(GOLD_BUF_IDX, gold_buf)
        call oclStoreBuffer(H_BUF_IDX, h_buf)
        call oclStoreBuffer(HOLD_BUF_IDX, hold_buf)
        call oclStoreBuffer(NRD_BUF_IDX, nrd_buf)
        call oclStoreBuffer(P0_BUF_IDX, p0_buf)
        call oclStoreBuffer(P1_BUF_IDX, p1_buf)
        call oclStoreBuffer(PAV_BUF_IDX, pav_buf)
        call oclStoreBuffer(RHS_BUF_IDX, rhs_buf)
        call oclStoreBuffer(RHSAV_BUF_IDX, rhsav_buf)
        call oclStoreBuffer(RO_BUF_IDX, ro_buf)
        call oclStoreBuffer(SM_BUF_IDX, sm_buf)
        call oclStoreBuffer(STATE_PTR_BUF_IDX, state_ptr_buf)
        call oclStoreBuffer(U_BUF_IDX, u_buf)
        call oclStoreBuffer(USUM_BUF_IDX, usum_buf)
        call oclStoreBuffer(V_BUF_IDX, v_buf)
        call oclStoreBuffer(VSUM_BUF_IDX, vsum_buf)
        call oclStoreBuffer(W_BUF_IDX, w_buf)
        call oclStoreBuffer(WSUM_BUF_IDX, wsum_buf)


    end subroutine adam_feedbf_les_press_velfg_ve_etc_superkernel_init
end module module_adam_feedbf_les_press_velfg_ve_etc_superkernel_init