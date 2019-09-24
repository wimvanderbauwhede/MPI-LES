module module_adam_bondv1_feedbf_les_press_v_etc_superkernel_init

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
        integer, parameter ::ALPHA_BUF_IDX = 24
        integer, parameter ::BETA_BUF_IDX = 25
        integer, parameter ::BMASK1_BUF_IDX = 19
        integer, parameter ::CMASK1_BUF_IDX = 21
        integer, parameter ::COV1_BUF_IDX = 57
        integer, parameter ::COV2_BUF_IDX = 58
        integer, parameter ::COV3_BUF_IDX = 59
        integer, parameter ::COV4_BUF_IDX = 60
        integer, parameter ::COV5_BUF_IDX = 61
        integer, parameter ::COV6_BUF_IDX = 62
        integer, parameter ::COV7_BUF_IDX = 63
        integer, parameter ::COV8_BUF_IDX = 64
        integer, parameter ::COV9_BUF_IDX = 65
        integer, parameter ::DELX1_BUF_IDX = 33
        integer, parameter ::DIU1_BUF_IDX = 34
        integer, parameter ::DIU2_BUF_IDX = 35
        integer, parameter ::DIU3_BUF_IDX = 36
        integer, parameter ::DIU4_BUF_IDX = 37
        integer, parameter ::DIU5_BUF_IDX = 38
        integer, parameter ::DIU6_BUF_IDX = 39
        integer, parameter ::DIU7_BUF_IDX = 40
        integer, parameter ::DIU8_BUF_IDX = 41
        integer, parameter ::DIU9_BUF_IDX = 42
        integer, parameter ::DMASK1_BUF_IDX = 23
        integer, parameter ::DT_BUF_IDX = 14
        integer, parameter ::DX1_BUF_IDX = 31
        integer, parameter ::DXS_BUF_IDX = 16
        integer, parameter ::DY1_BUF_IDX = 32
        integer, parameter ::DYS_BUF_IDX = 45
        integer, parameter ::DZN_BUF_IDX = 8
        integer, parameter ::DZS_BUF_IDX = 51
        integer, parameter ::F_BUF_IDX = 1
        integer, parameter ::FOLD_BUF_IDX = 4
        integer, parameter ::FX_BUF_IDX = 28
        integer, parameter ::FY_BUF_IDX = 29
        integer, parameter ::FZ_BUF_IDX = 30
        integer, parameter ::G_BUF_IDX = 2
        integer, parameter ::GLOBAL_AAA_ARRAY_BUF_IDX = 12
        integer, parameter ::GLOBAL_AREA_ARRAY_BUF_IDX = 49
        integer, parameter ::GLOBAL_BBB_ARRAY_BUF_IDX = 13
        integer, parameter ::GLOBAL_PAV_ARRAY_BUF_IDX = 54
        integer, parameter ::GLOBAL_PCO_ARRAY_BUF_IDX = 55
        integer, parameter ::GLOBAL_RHSAV_ARRAY_BUF_IDX = 48
        integer, parameter ::GOLD_BUF_IDX = 5
        integer, parameter ::H_BUF_IDX = 3
        integer, parameter ::HOLD_BUF_IDX = 6
        integer, parameter ::IP_BUF_IDX = 27
        integer, parameter ::JP_BUF_IDX = 26
        integer, parameter ::KP_BUF_IDX = 17
        integer, parameter ::NRD_BUF_IDX = 52
        integer, parameter ::P_BUF_IDX = 53
        integer, parameter ::PAV_BUF_IDX = 56
        integer, parameter ::RHS_BUF_IDX = 47
        integer, parameter ::RHSAV_BUF_IDX = 50
        integer, parameter ::RO_BUF_IDX = 66
        integer, parameter ::SM_BUF_IDX = 43
        integer, parameter ::STATE_PTR_BUF_IDX = 67
        integer, parameter ::U_BUF_IDX = 9
        integer, parameter ::UOUT_BUF_IDX = 15
        integer, parameter ::USPD_BUF_IDX = 44
        integer, parameter ::USUM_BUF_IDX = 18
        integer, parameter ::V_BUF_IDX = 10
        integer, parameter ::VSPD_BUF_IDX = 46
        integer, parameter ::VSUM_BUF_IDX = 20
        integer, parameter ::W_BUF_IDX = 11
        integer, parameter ::WSUM_BUF_IDX = 22
        integer, parameter ::Z2_BUF_IDX = 7

contains

    subroutine adam_bondv1_feedbf_les_press_v_etc_superkernel_init()

        use oclWrapper
        use common_sn
        character(len=*), parameter :: srcstr = "module_adam_bondv1_feedbf_les_press_v_etc_superkernel.cl"
        character(len=*), parameter :: kstr   = "adam_bondv1_feedbf_les_press_v_etc_superkernel"
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
! declarations
        real(kind=4), dimension(0:(kp + 2)) :: z2
        real(kind=4), dimension((-1):(kp + 2)) :: dzn
        real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)) :: u
        real(kind=4) :: dt
        real(kind=4) :: uout
        real(kind=4), dimension(0:ip) :: dxs
        real :: kp
        real(kind=4), dimension((-1):(ip + 1),0:(jp + 1),0:(kp + 1)) :: bmask1
        real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)) :: v
        real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),0:(kp + 1)) :: cmask1
        real(kind=4), dimension(0:(ip + 1),(-1):(jp + 1),(-1):(kp + 1)) :: w
        real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)) :: dmask1
        real(kind=4) :: alpha
        real(kind=4) :: beta
        real :: jp
        real :: ip
        real(kind=4), dimension((-1):(ip + 1)) :: dx1
        real(kind=4), dimension(0:(jp + 1)) :: dy1
        real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)) :: diu1
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: diu2
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: diu3
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: diu4
        real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)) :: diu5
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: diu6
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: diu7
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: diu8
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: diu9
        real(kind=4), dimension(kp) :: delx1
        real(kind=4), dimension((-1):(ip + 1),(-1):(jp + 1),0:(kp + 1)) :: sm
        real(kind=4), dimension(0:(ip + 1),0:(jp + 1)) :: uspd
        real(kind=4), dimension(0:jp) :: dys
        real(kind=4), dimension(0:(ip + 1),0:(jp + 1)) :: vspd
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: f
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: g
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: h
        real(kind=4), dimension(0:(ip + 1),0:(jp + 1),0:(kp + 1)) :: rhs
        real(kind=4) :: rhsav
        real(kind=4), dimension((-1):(kp + 2)) :: dzs
        integer :: nrd
        real(kind=4), dimension(0:1,0:(ip + 2),0:(jp + 2),0:(kp + 1)) :: p
        real(kind=4) :: pav
        real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)) :: cov1
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: cov2
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: cov3
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: cov4
        real(kind=4), dimension((-1):(ip + 2),0:(jp + 2),0:(kp + 2)) :: cov5
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: cov6
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: cov7
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: cov8
        real(kind=4), dimension(0:(ip + 2),0:(jp + 2),0:(kp + 2)) :: cov9
        real(kind=4) :: ro
        real(kind=4), dimension(ip,jp,kp) :: fold
        real(kind=4), dimension(ip,jp,kp) :: gold
        real(kind=4), dimension(ip,jp,kp) :: hold
        real(kind=4), dimension(1:NUNITS) :: global_aaa_array
        real(kind=4), dimension(1:NUNITS) :: global_bbb_array
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: usum
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: vsum
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: wsum
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: fx
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: fy
        real(kind=4), dimension(0:ip,0:jp,0:kp) :: fz
        real(kind=4), dimension(1:NUNITS) :: global_rhsav_array
        real(kind=4), dimension(1:NUNITS) :: global_area_array
        real(kind=4), dimension(1:NUNITS) :: global_pav_array
        real(kind=4), dimension(1:NUNITS) :: global_pco_array
        integer, dimension(1) :: state_ptr
! buffer declarations
        integer(8) :: z2_buf
        integer(8) :: dzn_buf
        integer(8) :: u_buf
        integer(8) :: dt_buf
        integer(8) :: uout_buf
        integer(8) :: dxs_buf
        integer(8) :: kp_buf
        integer(8) :: bmask1_buf
        integer(8) :: v_buf
        integer(8) :: cmask1_buf
        integer(8) :: w_buf
        integer(8) :: dmask1_buf
        integer(8) :: alpha_buf
        integer(8) :: beta_buf
        integer(8) :: jp_buf
        integer(8) :: ip_buf
        integer(8) :: dx1_buf
        integer(8) :: dy1_buf
        integer(8) :: diu1_buf
        integer(8) :: diu2_buf
        integer(8) :: diu3_buf
        integer(8) :: diu4_buf
        integer(8) :: diu5_buf
        integer(8) :: diu6_buf
        integer(8) :: diu7_buf
        integer(8) :: diu8_buf
        integer(8) :: diu9_buf
        integer(8) :: delx1_buf
        integer(8) :: sm_buf
        integer(8) :: uspd_buf
        integer(8) :: dys_buf
        integer(8) :: vspd_buf
        integer(8) :: f_buf
        integer(8) :: g_buf
        integer(8) :: h_buf
        integer(8) :: rhs_buf
        integer(8) :: rhsav_buf
        integer(8) :: dzs_buf
        integer(8) :: nrd_buf
        integer(8) :: p_buf
        integer(8) :: pav_buf
        integer(8) :: cov1_buf
        integer(8) :: cov2_buf
        integer(8) :: cov3_buf
        integer(8) :: cov4_buf
        integer(8) :: cov5_buf
        integer(8) :: cov6_buf
        integer(8) :: cov7_buf
        integer(8) :: cov8_buf
        integer(8) :: cov9_buf
        integer(8) :: ro_buf
        integer(8) :: fold_buf
        integer(8) :: gold_buf
        integer(8) :: hold_buf
        integer(8) :: global_aaa_array_buf
        integer(8) :: global_bbb_array_buf
        integer(8) :: usum_buf
        integer(8) :: vsum_buf
        integer(8) :: wsum_buf
        integer(8) :: fx_buf
        integer(8) :: fy_buf
        integer(8) :: fz_buf
        integer(8) :: global_rhsav_array_buf
        integer(8) :: global_area_array_buf
        integer(8) :: global_pav_array_buf
        integer(8) :: global_pco_array_buf
        integer(8) :: state_ptr_buf
        integer, dimension(1) :: z2_sz
        integer, dimension(1) :: dzn_sz
        integer, dimension(3) :: u_sz
        integer, dimension(1) :: dt_ptr_sz
        integer, dimension(1) :: uout_ptr_sz
        integer, dimension(1) :: dxs_sz
        integer, dimension(1) :: kp_ptr_sz
        integer, dimension(3) :: bmask1_sz
        integer, dimension(3) :: v_sz
        integer, dimension(3) :: cmask1_sz
        integer, dimension(3) :: w_sz
        integer, dimension(3) :: dmask1_sz
        integer, dimension(1) :: alpha_ptr_sz
        integer, dimension(1) :: beta_ptr_sz
        integer, dimension(1) :: jp_ptr_sz
        integer, dimension(1) :: ip_ptr_sz
        integer, dimension(1) :: dx1_sz
        integer, dimension(1) :: dy1_sz
        integer, dimension(3) :: diu1_sz
        integer, dimension(3) :: diu2_sz
        integer, dimension(3) :: diu3_sz
        integer, dimension(3) :: diu4_sz
        integer, dimension(3) :: diu5_sz
        integer, dimension(3) :: diu6_sz
        integer, dimension(3) :: diu7_sz
        integer, dimension(3) :: diu8_sz
        integer, dimension(3) :: diu9_sz
        integer, dimension(1) :: delx1_sz
        integer, dimension(3) :: sm_sz
        integer, dimension(2) :: uspd_sz
        integer, dimension(1) :: dys_sz
        integer, dimension(2) :: vspd_sz
        integer, dimension(3) :: f_sz
        integer, dimension(3) :: g_sz
        integer, dimension(3) :: h_sz
        integer, dimension(3) :: rhs_sz
        integer, dimension(1) :: rhsav_ptr_sz
        integer, dimension(1) :: dzs_sz
        integer, dimension(1) :: nrd_ptr_sz
        integer, dimension(4) :: p_sz
        integer, dimension(1) :: pav_ptr_sz
        integer, dimension(3) :: cov1_sz
        integer, dimension(3) :: cov2_sz
        integer, dimension(3) :: cov3_sz
        integer, dimension(3) :: cov4_sz
        integer, dimension(3) :: cov5_sz
        integer, dimension(3) :: cov6_sz
        integer, dimension(3) :: cov7_sz
        integer, dimension(3) :: cov8_sz
        integer, dimension(3) :: cov9_sz
        integer, dimension(1) :: ro_ptr_sz
        integer, dimension(3) :: fold_sz
        integer, dimension(3) :: gold_sz
        integer, dimension(3) :: hold_sz
        integer, dimension(1) :: global_aaa_array_sz
        integer, dimension(1) :: global_bbb_array_sz
        integer, dimension(3) :: usum_sz
        integer, dimension(3) :: vsum_sz
        integer, dimension(3) :: wsum_sz
        integer, dimension(3) :: fx_sz
        integer, dimension(3) :: fy_sz
        integer, dimension(3) :: fz_sz
        integer, dimension(1) :: global_rhsav_array_sz
        integer, dimension(1) :: global_area_array_sz
        integer, dimension(1) :: global_pav_array_sz
        integer, dimension(1) :: global_pco_array_sz
        integer, dimension(1) :: state_ptr_sz
        real(kind=4), dimension(1) :: dt_ptr
        real(kind=4), dimension(1) :: uout_ptr
        real, dimension(1) :: kp_ptr
        real(kind=4), dimension(1) :: alpha_ptr
        real(kind=4), dimension(1) :: beta_ptr
        real, dimension(1) :: jp_ptr
        real, dimension(1) :: ip_ptr
        real(kind=4), dimension(1) :: rhsav_ptr
        integer, dimension(1) :: nrd_ptr
        real(kind=4), dimension(1) :: pav_ptr
        real(kind=4), dimension(1) :: ro_ptr

        call oclInit(srcstr,kstr)

        z2_sz = shape(z2)
        dzn_sz = shape(dzn)
        u_sz = shape(u)
        dt_ptr_sz = shape(dt_ptr)
        uout_ptr_sz = shape(uout_ptr)
        dxs_sz = shape(dxs)
        kp_ptr_sz = shape(kp_ptr)
        bmask1_sz = shape(bmask1)
        v_sz = shape(v)
        cmask1_sz = shape(cmask1)
        w_sz = shape(w)
        dmask1_sz = shape(dmask1)
        alpha_ptr_sz = shape(alpha_ptr)
        beta_ptr_sz = shape(beta_ptr)
        jp_ptr_sz = shape(jp_ptr)
        ip_ptr_sz = shape(ip_ptr)
        dx1_sz = shape(dx1)
        dy1_sz = shape(dy1)
        diu1_sz = shape(diu1)
        diu2_sz = shape(diu2)
        diu3_sz = shape(diu3)
        diu4_sz = shape(diu4)
        diu5_sz = shape(diu5)
        diu6_sz = shape(diu6)
        diu7_sz = shape(diu7)
        diu8_sz = shape(diu8)
        diu9_sz = shape(diu9)
        delx1_sz = shape(delx1)
        sm_sz = shape(sm)
        uspd_sz = shape(uspd)
        dys_sz = shape(dys)
        vspd_sz = shape(vspd)
        f_sz = shape(f)
        g_sz = shape(g)
        h_sz = shape(h)
        rhs_sz = shape(rhs)
        rhsav_ptr_sz = shape(rhsav_ptr)
        dzs_sz = shape(dzs)
        nrd_ptr_sz = shape(nrd_ptr)
        p_sz = shape(p)
        pav_ptr_sz = shape(pav_ptr)
        cov1_sz = shape(cov1)
        cov2_sz = shape(cov2)
        cov3_sz = shape(cov3)
        cov4_sz = shape(cov4)
        cov5_sz = shape(cov5)
        cov6_sz = shape(cov6)
        cov7_sz = shape(cov7)
        cov8_sz = shape(cov8)
        cov9_sz = shape(cov9)
        ro_ptr_sz = shape(ro_ptr)
        fold_sz = shape(fold)
        gold_sz = shape(gold)
        hold_sz = shape(hold)
        global_aaa_array_sz = shape(global_aaa_array)
        global_bbb_array_sz = shape(global_bbb_array)
        usum_sz = shape(usum)
        vsum_sz = shape(vsum)
        wsum_sz = shape(wsum)
        fx_sz = shape(fx)
        fy_sz = shape(fy)
        fz_sz = shape(fz)
        global_rhsav_array_sz = shape(global_rhsav_array)
        global_area_array_sz = shape(global_area_array)
        global_pav_array_sz = shape(global_pav_array)
        global_pco_array_sz = shape(global_pco_array)
        state_ptr_sz = shape(state_ptr)

        call oclMake1DFloatArrayReadWriteBuffer(z2_buf,z2_sz,z2)
        call oclMake1DFloatArrayReadWriteBuffer(dzn_buf,dzn_sz,dzn)
        call oclMake3DFloatArrayReadWriteBuffer(u_buf,u_sz,u)
        call oclMake1DFloatArrayReadWriteBuffer(dt_buf,dt_ptr_sz,dt_ptr)! Automatic conversion to array
        call oclMake1DFloatArrayReadWriteBuffer(uout_buf,uout_ptr_sz,uout_ptr)! Automatic conversion to array
        call oclMake1DFloatArrayReadWriteBuffer(dxs_buf,dxs_sz,dxs)
        call oclMake1DFloatArrayReadWriteBuffer(kp_buf,kp_ptr_sz,kp_ptr)! Automatic conversion to array
        call oclMake3DFloatArrayReadWriteBuffer(bmask1_buf,bmask1_sz,bmask1)
        call oclMake3DFloatArrayReadWriteBuffer(v_buf,v_sz,v)
        call oclMake3DFloatArrayReadWriteBuffer(cmask1_buf,cmask1_sz,cmask1)
        call oclMake3DFloatArrayReadWriteBuffer(w_buf,w_sz,w)
        call oclMake3DFloatArrayReadWriteBuffer(dmask1_buf,dmask1_sz,dmask1)
        call oclMake1DFloatArrayReadWriteBuffer(alpha_buf,alpha_ptr_sz,alpha_ptr)! Automatic conversion to array
        call oclMake1DFloatArrayReadWriteBuffer(beta_buf,beta_ptr_sz,beta_ptr)! Automatic conversion to array
        call oclMake1DFloatArrayReadWriteBuffer(jp_buf,jp_ptr_sz,jp_ptr)! Automatic conversion to array
        call oclMake1DFloatArrayReadWriteBuffer(ip_buf,ip_ptr_sz,ip_ptr)! Automatic conversion to array
        call oclMake1DFloatArrayReadWriteBuffer(dx1_buf,dx1_sz,dx1)
        call oclMake1DFloatArrayReadWriteBuffer(dy1_buf,dy1_sz,dy1)
        call oclMake3DFloatArrayReadWriteBuffer(diu1_buf,diu1_sz,diu1)
        call oclMake3DFloatArrayReadWriteBuffer(diu2_buf,diu2_sz,diu2)
        call oclMake3DFloatArrayReadWriteBuffer(diu3_buf,diu3_sz,diu3)
        call oclMake3DFloatArrayReadWriteBuffer(diu4_buf,diu4_sz,diu4)
        call oclMake3DFloatArrayReadWriteBuffer(diu5_buf,diu5_sz,diu5)
        call oclMake3DFloatArrayReadWriteBuffer(diu6_buf,diu6_sz,diu6)
        call oclMake3DFloatArrayReadWriteBuffer(diu7_buf,diu7_sz,diu7)
        call oclMake3DFloatArrayReadWriteBuffer(diu8_buf,diu8_sz,diu8)
        call oclMake3DFloatArrayReadWriteBuffer(diu9_buf,diu9_sz,diu9)
        call oclMake1DFloatArrayReadWriteBuffer(delx1_buf,delx1_sz,delx1)
        call oclMake3DFloatArrayReadWriteBuffer(sm_buf,sm_sz,sm)
        call oclMake2DFloatArrayReadWriteBuffer(uspd_buf,uspd_sz,uspd)
        call oclMake1DFloatArrayReadWriteBuffer(dys_buf,dys_sz,dys)
        call oclMake2DFloatArrayReadWriteBuffer(vspd_buf,vspd_sz,vspd)
        call oclMake3DFloatArrayReadWriteBuffer(f_buf,f_sz,f)
        call oclMake3DFloatArrayReadWriteBuffer(g_buf,g_sz,g)
        call oclMake3DFloatArrayReadWriteBuffer(h_buf,h_sz,h)
        call oclMake3DFloatArrayReadWriteBuffer(rhs_buf,rhs_sz,rhs)
        call oclMake1DFloatArrayReadWriteBuffer(rhsav_buf,rhsav_ptr_sz,rhsav_ptr)! Automatic conversion to array
        call oclMake1DFloatArrayReadWriteBuffer(dzs_buf,dzs_sz,dzs)
        call oclMake1DIntArrayReadWriteBuffer(nrd_buf,nrd_ptr_sz,nrd_ptr)! Automatic conversion to array
        call oclMake4DFloatArrayReadWriteBuffer(p_buf,p_sz,p)
        call oclMake1DFloatArrayReadWriteBuffer(pav_buf,pav_ptr_sz,pav_ptr)! Automatic conversion to array
        call oclMake3DFloatArrayReadWriteBuffer(cov1_buf,cov1_sz,cov1)
        call oclMake3DFloatArrayReadWriteBuffer(cov2_buf,cov2_sz,cov2)
        call oclMake3DFloatArrayReadWriteBuffer(cov3_buf,cov3_sz,cov3)
        call oclMake3DFloatArrayReadWriteBuffer(cov4_buf,cov4_sz,cov4)
        call oclMake3DFloatArrayReadWriteBuffer(cov5_buf,cov5_sz,cov5)
        call oclMake3DFloatArrayReadWriteBuffer(cov6_buf,cov6_sz,cov6)
        call oclMake3DFloatArrayReadWriteBuffer(cov7_buf,cov7_sz,cov7)
        call oclMake3DFloatArrayReadWriteBuffer(cov8_buf,cov8_sz,cov8)
        call oclMake3DFloatArrayReadWriteBuffer(cov9_buf,cov9_sz,cov9)
        call oclMake1DFloatArrayReadWriteBuffer(ro_buf,ro_ptr_sz,ro_ptr)! Automatic conversion to array
        call oclMake3DFloatArrayReadWriteBuffer(fold_buf,fold_sz,fold)
        call oclMake3DFloatArrayReadWriteBuffer(gold_buf,gold_sz,gold)
        call oclMake3DFloatArrayReadWriteBuffer(hold_buf,hold_sz,hold)
        call oclMake1DFloatArrayReadWriteBuffer(global_aaa_array_buf,global_aaa_array_sz,global_aaa_array)
        call oclMake1DFloatArrayReadWriteBuffer(global_bbb_array_buf,global_bbb_array_sz,global_bbb_array)
        call oclMake3DFloatArrayReadWriteBuffer(usum_buf,usum_sz,usum)
        call oclMake3DFloatArrayReadWriteBuffer(vsum_buf,vsum_sz,vsum)
        call oclMake3DFloatArrayReadWriteBuffer(wsum_buf,wsum_sz,wsum)
        call oclMake3DFloatArrayReadWriteBuffer(fx_buf,fx_sz,fx)
        call oclMake3DFloatArrayReadWriteBuffer(fy_buf,fy_sz,fy)
        call oclMake3DFloatArrayReadWriteBuffer(fz_buf,fz_sz,fz)
        call oclMake1DFloatArrayReadWriteBuffer(global_rhsav_array_buf,global_rhsav_array_sz,global_rhsav_array)
        call oclMake1DFloatArrayReadWriteBuffer(global_area_array_buf,global_area_array_sz,global_area_array)
        call oclMake1DFloatArrayReadWriteBuffer(global_pav_array_buf,global_pav_array_sz,global_pav_array)
        call oclMake1DFloatArrayReadWriteBuffer(global_pco_array_buf,global_pco_array_sz,global_pco_array)
        call oclMake1DIntArrayReadWriteBuffer(state_ptr_buf,state_ptr_sz,state_ptr)

        call oclSetFloatArrayArg(6, z2_buf)
        call oclSetFloatArrayArg(7, dzn_buf)
        call oclSetFloatArrayArg(8, u_buf)
        call oclSetFloatArrayArg(13, dt_buf)
        call oclSetFloatArrayArg(14, uout_buf)
        call oclSetFloatArrayArg(15, dxs_buf)
        call oclSetFloatArrayArg(16, kp_buf)
        call oclSetFloatArrayArg(18, bmask1_buf)
        call oclSetFloatArrayArg(9, v_buf)
        call oclSetFloatArrayArg(20, cmask1_buf)
        call oclSetFloatArrayArg(10, w_buf)
        call oclSetFloatArrayArg(22, dmask1_buf)
        call oclSetFloatArrayArg(23, alpha_buf)
        call oclSetFloatArrayArg(24, beta_buf)
        call oclSetFloatArrayArg(25, jp_buf)
        call oclSetFloatArrayArg(26, ip_buf)
        call oclSetFloatArrayArg(30, dx1_buf)
        call oclSetFloatArrayArg(31, dy1_buf)
        call oclSetFloatArrayArg(33, diu1_buf)
        call oclSetFloatArrayArg(34, diu2_buf)
        call oclSetFloatArrayArg(35, diu3_buf)
        call oclSetFloatArrayArg(36, diu4_buf)
        call oclSetFloatArrayArg(37, diu5_buf)
        call oclSetFloatArrayArg(38, diu6_buf)
        call oclSetFloatArrayArg(39, diu7_buf)
        call oclSetFloatArrayArg(40, diu8_buf)
        call oclSetFloatArrayArg(41, diu9_buf)
        call oclSetFloatArrayArg(32, delx1_buf)
        call oclSetFloatArrayArg(42, sm_buf)
        call oclSetFloatArrayArg(43, uspd_buf)
        call oclSetFloatArrayArg(44, dys_buf)
        call oclSetFloatArrayArg(45, vspd_buf)
        call oclSetFloatArrayArg(0, f_buf)
        call oclSetFloatArrayArg(1, g_buf)
        call oclSetFloatArrayArg(2, h_buf)
        call oclSetFloatArrayArg(46, rhs_buf)
        call oclSetFloatArrayArg(49, rhsav_buf)
        call oclSetFloatArrayArg(50, dzs_buf)
        call oclSetIntArrayArg(51, nrd_buf)
        call oclSetFloatArrayArg(52, p_buf)
        call oclSetFloatArrayArg(55, pav_buf)
        call oclSetFloatArrayArg(56, cov1_buf)
        call oclSetFloatArrayArg(57, cov2_buf)
        call oclSetFloatArrayArg(58, cov3_buf)
        call oclSetFloatArrayArg(59, cov4_buf)
        call oclSetFloatArrayArg(60, cov5_buf)
        call oclSetFloatArrayArg(61, cov6_buf)
        call oclSetFloatArrayArg(62, cov7_buf)
        call oclSetFloatArrayArg(63, cov8_buf)
        call oclSetFloatArrayArg(64, cov9_buf)
        call oclSetFloatArrayArg(65, ro_buf)
        call oclSetFloatArrayArg(3, fold_buf)
        call oclSetFloatArrayArg(4, gold_buf)
        call oclSetFloatArrayArg(5, hold_buf)
        call oclSetFloatArrayArg(11, global_aaa_array_buf)
        call oclSetFloatArrayArg(12, global_bbb_array_buf)
        call oclSetFloatArrayArg(17, usum_buf)
        call oclSetFloatArrayArg(19, vsum_buf)
        call oclSetFloatArrayArg(21, wsum_buf)
        call oclSetFloatArrayArg(27, fx_buf)
        call oclSetFloatArrayArg(28, fy_buf)
        call oclSetFloatArrayArg(29, fz_buf)
        call oclSetFloatArrayArg(47, global_rhsav_array_buf)
        call oclSetFloatArrayArg(48, global_area_array_buf)
        call oclSetFloatArrayArg(53, global_pav_array_buf)
        call oclSetFloatArrayArg(54, global_pco_array_buf)
        call oclSetIntArrayArg(66, state_ptr_buf)

        call oclStoreBuffer(ALPHA_BUF_IDX, alpha_buf)
        call oclStoreBuffer(BETA_BUF_IDX, beta_buf)
        call oclStoreBuffer(BMASK1_BUF_IDX, bmask1_buf)
        call oclStoreBuffer(CMASK1_BUF_IDX, cmask1_buf)
        call oclStoreBuffer(COV1_BUF_IDX, cov1_buf)
        call oclStoreBuffer(COV2_BUF_IDX, cov2_buf)
        call oclStoreBuffer(COV3_BUF_IDX, cov3_buf)
        call oclStoreBuffer(COV4_BUF_IDX, cov4_buf)
        call oclStoreBuffer(COV5_BUF_IDX, cov5_buf)
        call oclStoreBuffer(COV6_BUF_IDX, cov6_buf)
        call oclStoreBuffer(COV7_BUF_IDX, cov7_buf)
        call oclStoreBuffer(COV8_BUF_IDX, cov8_buf)
        call oclStoreBuffer(COV9_BUF_IDX, cov9_buf)
        call oclStoreBuffer(DELX1_BUF_IDX, delx1_buf)
        call oclStoreBuffer(DIU1_BUF_IDX, diu1_buf)
        call oclStoreBuffer(DIU2_BUF_IDX, diu2_buf)
        call oclStoreBuffer(DIU3_BUF_IDX, diu3_buf)
        call oclStoreBuffer(DIU4_BUF_IDX, diu4_buf)
        call oclStoreBuffer(DIU5_BUF_IDX, diu5_buf)
        call oclStoreBuffer(DIU6_BUF_IDX, diu6_buf)
        call oclStoreBuffer(DIU7_BUF_IDX, diu7_buf)
        call oclStoreBuffer(DIU8_BUF_IDX, diu8_buf)
        call oclStoreBuffer(DIU9_BUF_IDX, diu9_buf)
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
        call oclStoreBuffer(GLOBAL_AAA_ARRAY_BUF_IDX, global_aaa_array_buf)
        call oclStoreBuffer(GLOBAL_AREA_ARRAY_BUF_IDX, global_area_array_buf)
        call oclStoreBuffer(GLOBAL_BBB_ARRAY_BUF_IDX, global_bbb_array_buf)
        call oclStoreBuffer(GLOBAL_PAV_ARRAY_BUF_IDX, global_pav_array_buf)
        call oclStoreBuffer(GLOBAL_PCO_ARRAY_BUF_IDX, global_pco_array_buf)
        call oclStoreBuffer(GLOBAL_RHSAV_ARRAY_BUF_IDX, global_rhsav_array_buf)
        call oclStoreBuffer(GOLD_BUF_IDX, gold_buf)
        call oclStoreBuffer(H_BUF_IDX, h_buf)
        call oclStoreBuffer(HOLD_BUF_IDX, hold_buf)
        call oclStoreBuffer(IP_BUF_IDX, ip_buf)
        call oclStoreBuffer(JP_BUF_IDX, jp_buf)
        call oclStoreBuffer(KP_BUF_IDX, kp_buf)
        call oclStoreBuffer(NRD_BUF_IDX, nrd_buf)
        call oclStoreBuffer(P_BUF_IDX, p_buf)
        call oclStoreBuffer(PAV_BUF_IDX, pav_buf)
        call oclStoreBuffer(RHS_BUF_IDX, rhs_buf)
        call oclStoreBuffer(RHSAV_BUF_IDX, rhsav_buf)
        call oclStoreBuffer(RO_BUF_IDX, ro_buf)
        call oclStoreBuffer(SM_BUF_IDX, sm_buf)
        call oclStoreBuffer(STATE_PTR_BUF_IDX, state_ptr_buf)
        call oclStoreBuffer(U_BUF_IDX, u_buf)
        call oclStoreBuffer(UOUT_BUF_IDX, uout_buf)
        call oclStoreBuffer(USPD_BUF_IDX, uspd_buf)
        call oclStoreBuffer(USUM_BUF_IDX, usum_buf)
        call oclStoreBuffer(V_BUF_IDX, v_buf)
        call oclStoreBuffer(VSPD_BUF_IDX, vspd_buf)
        call oclStoreBuffer(VSUM_BUF_IDX, vsum_buf)
        call oclStoreBuffer(W_BUF_IDX, w_buf)
        call oclStoreBuffer(WSUM_BUF_IDX, wsum_buf)
        call oclStoreBuffer(Z2_BUF_IDX, z2_buf)


    end subroutine adam_bondv1_feedbf_les_press_v_etc_superkernel_init
end module module_adam_bondv1_feedbf_les_press_v_etc_superkernel_init