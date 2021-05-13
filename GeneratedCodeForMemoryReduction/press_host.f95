module module_press
    use module_bondFG 
    use module_boundp 
 contains 
subroutine press(u,v,w,p0,p1,rhs,f,g,h,dx1,dy1,dzn,dxs,dys,dzs,dt,n,nmax &
)

        use module_adam_feedbf_les_press_velfg_ve_etc_superkernel_init
        use oclWrapper
!      implicit none
      use params_common_sn

      implicit none

! otherStatements

! remainingDecls
      real(4), dimension(0:ip), intent(in) :: dxs !!
      real(4), dimension(0:jp), intent(in) :: dys !!
      real(4), dimension(-1:kp+2), intent(in) :: dzs !!
      real(4) :: cn1 !!
      real(4) :: cn2l !!
      real(4) :: cn2s !!
      real(4) :: cn3l !!
      real(4) :: cn3s !!
      real(4) :: cn4l !!
      real(4) :: cn4s !!
      real(4) :: dz1 !!
      real(4) :: dz2 !!
      real(4), intent(in) :: dt !!
      real(4), dimension(-1:ip+1), intent(in) :: dx1 !!
      real(4), dimension(0:jp+1), intent(in) :: dy1 !!
      real(4), dimension(-1:kp+2), intent(in) :: dzn !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: f !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: g !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: h !!
      integer, intent(in) :: n !!
      integer, intent(in) :: nmax !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+1) :: p0 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+1) :: p1 !!
      real(4), dimension(0:ip+1,0:jp+1,0:kp+1), intent(out) :: rhs !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(in) :: u !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(in) :: v !!
      real(4), dimension(0:ip+1,-1:jp+1,-1:kp+1), intent(in) :: w !!
      integer :: nn !!
      integer :: i !!
      integer :: j !!
      integer :: k !!
      integer :: l !!
      integer :: nrd !!
      real(4) :: rhsav !!
      real(4) :: pav !!
      real(4) :: area !!
      real(4) :: pco !!
      real(4) :: sor !!
      real(4) :: reltmp !!
      real, parameter :: pjuge = 0.0001  !!
      integer, parameter :: nmaxp = 50  !!
      real, parameter :: omega = 1.  !!

        ! Extra declarations
        real (kind=4) :: exectime
        real(kind=4), dimension(1:NUNITS) :: global_pav_array
        real(kind=4), dimension(1:NUNITS) :: global_pco_array
        real(kind=4), dimension(1:NUNITS) :: global_rhsav_array
        real(kind=4), dimension(1:NUNITS) :: global_area_array

        ! Buffer declarations
        integer(8) :: state_ptr_buf
        integer(8) :: u_buf
        integer(8) :: dx1_buf
        integer(8) :: v_buf
        integer(8) :: dy1_buf
        integer(8) :: w_buf
        integer(8) :: dzn_buf
        integer(8) :: f_buf
        integer(8) :: g_buf
        integer(8) :: h_buf
        integer(8) :: rhs_buf
        integer(8) :: dt_buf
        integer(8) :: global_rhsav_array_buf
        integer(8) :: global_area_array_buf
        integer(8) :: rhsav_buf
        integer(8) :: dzs_buf
        integer(8) :: dys_buf
        integer(8) :: dxs_buf
        integer(8) :: nrd_buf
        integer(8) :: p0_buf
        integer(8) :: p1_buf
        integer(8) :: global_pav_array_buf
        integer(8) :: global_pco_array_buf
        integer(8) :: pav_buf

        integer, dimension(1) :: state_ptr

        ! Size declarations
        integer, dimension(1) :: state_ptr_sz
        integer, dimension(3) :: u_sz
        integer, dimension(1) :: dx1_sz
        integer, dimension(3) :: v_sz
        integer, dimension(1) :: dy1_sz
        integer, dimension(3) :: w_sz
        integer, dimension(1) :: dzn_sz
        integer, dimension(3) :: f_sz
        integer, dimension(3) :: g_sz
        integer, dimension(3) :: h_sz
        integer, dimension(3) :: rhs_sz
        integer, dimension(1) :: dt_ptr_sz
        integer, dimension(1) :: rhsav_ptr_sz
        integer, dimension(1) :: dzs_sz
        integer, dimension(1) :: dys_sz
        integer, dimension(1) :: dxs_sz
        integer, dimension(1) :: nrd_ptr_sz
        integer, dimension(3) :: p0_sz
        integer, dimension(3) :: p1_sz
        integer, dimension(1) :: pav_ptr_sz
        integer, dimension(1) :: global_rhsav_array_sz
        integer, dimension(1) :: global_area_array_sz
        integer, dimension(1) :: global_pav_array_sz
        integer, dimension(1) :: global_pco_array_sz
        real(kind=4), dimension(1) :: dt_ptr
        real(kind=4), dimension(1) :: rhsav_ptr
        integer, dimension(1) :: nrd_ptr
        real(kind=4), dimension(1) :: pav_ptr

        Integer :: r_iter
        
! Size assignments
        state_ptr_sz = shape(state_ptr)
        u_sz = shape(u)
        dx1_sz = shape(dx1)
        v_sz = shape(v)
        dy1_sz = shape(dy1)
        w_sz = shape(w)
        dzn_sz = shape(dzn)
        f_sz = shape(f)
        g_sz = shape(g)
        h_sz = shape(h)
        rhs_sz = shape(rhs)
        dt_ptr_sz = shape(dt_ptr)
        rhsav_ptr_sz = shape(rhsav_ptr)
        dzs_sz = shape(dzs)
        dys_sz = shape(dys)
        dxs_sz = shape(dxs)
        nrd_ptr_sz = shape(nrd_ptr)
        p0_sz = shape(p0)
        p1_sz = shape(p1)
        pav_ptr_sz = shape(pav_ptr)
        global_rhsav_array_sz = shape(global_rhsav_array)
        global_area_array_sz = shape(global_area_array)
        global_pav_array_sz = shape(global_pav_array)
        global_pco_array_sz = shape(global_pco_array)

        ! Buffer loads
        call oclLoadBuffer(STATE_PTR_BUF_IDX, state_ptr_buf)
        call oclLoadBuffer(U_BUF_IDX, u_buf)
        call oclLoadBuffer(DX1_BUF_IDX, dx1_buf)
        call oclLoadBuffer(V_BUF_IDX, v_buf)
        call oclLoadBuffer(DY1_BUF_IDX, dy1_buf)
        call oclLoadBuffer(W_BUF_IDX, w_buf)
        call oclLoadBuffer(DZN_BUF_IDX, dzn_buf)
        call oclLoadBuffer(F_BUF_IDX, f_buf)
        call oclLoadBuffer(G_BUF_IDX, g_buf)
        call oclLoadBuffer(H_BUF_IDX, h_buf)
        call oclLoadBuffer(RHS_BUF_IDX, rhs_buf)
        call oclLoadBuffer(DT_BUF_IDX, dt_buf)
        call oclLoadBuffer(GLOBAL_RHSAV_ARRAY_BUF_IDX, global_rhsav_array_buf)
        call oclLoadBuffer(GLOBAL_AREA_ARRAY_BUF_IDX, global_area_array_buf)
        call oclLoadBuffer(RHSAV_BUF_IDX, rhsav_buf)
        call oclLoadBuffer(DZS_BUF_IDX, dzs_buf)
        call oclLoadBuffer(DYS_BUF_IDX, dys_buf)
        call oclLoadBuffer(DXS_BUF_IDX, dxs_buf)
        call oclLoadBuffer(NRD_BUF_IDX, nrd_buf)
        call oclLoadBuffer(P0_BUF_IDX, p0_buf)
        call oclLoadBuffer(P1_BUF_IDX, p1_buf)
        call oclLoadBuffer(GLOBAL_PAV_ARRAY_BUF_IDX, global_pav_array_buf)
        call oclLoadBuffer(GLOBAL_PCO_ARRAY_BUF_IDX, global_pco_array_buf)
        call oclLoadBuffer(PAV_BUF_IDX, pav_buf)

        ! Original code with buffer writes and reads
! ---- BEGIN press_map_65 -----------------------------------------------------------------------------------------------------
        oclGlobalRange = (80 * (300 * 300))
        oclLocalRange = 0
        state_ptr(1) = ST_PRESS_MAP_65
        
        call oclWrite3DFloatArrayBuffer(u_buf,u_sz,u)
        call oclWrite1DFloatArrayBuffer(dx1_buf,dx1_sz,dx1)
        call oclWrite3DFloatArrayBuffer(v_buf,v_sz,v)
        call oclWrite1DFloatArrayBuffer(dy1_buf,dy1_sz,dy1)
        call oclWrite3DFloatArrayBuffer(w_buf,w_sz,w)
        call oclWrite1DFloatArrayBuffer(dzn_buf,dzn_sz,dzn)
        call oclWrite3DFloatArrayBuffer(f_buf,f_sz,f)
        call oclWrite3DFloatArrayBuffer(g_buf,g_sz,g)
        call oclWrite3DFloatArrayBuffer(h_buf,h_sz,h)
        dt_ptr(1) = dt
        call oclWrite1DFloatArrayBuffer(dt_buf,dt_ptr_sz,dt_ptr)! Automatic conversion to array
        call oclWrite3DFloatArrayBuffer(rhs_buf,rhs_sz,rhs)
        call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
        call runOcl(oclGlobalRange,oclLocalRange,exectime)
        ! call press_map_65

        call oclRead3DFloatArrayBuffer(rhs_buf,rhs_sz,rhs)
! ---- END --------------------------------------------------------------------------------------------------------------------
    rhsav = 0.0
    area = 0.0
! ---- BEGIN press_reduce_78 --------------------------------------------------------------------------------------------------
        oclGlobalRange = 1
        oclLocalRange = 1
        state_ptr(1) = ST_PRESS_REDUCE_78
        
        call oclWrite1DFloatArrayBuffer(dx1_buf,dx1_sz,dx1)
        call oclWrite1DFloatArrayBuffer(dy1_buf,dy1_sz,dy1)
        call oclWrite1DFloatArrayBuffer(dzn_buf,dzn_sz,dzn)
        call oclWrite3DFloatArrayBuffer(rhs_buf,rhs_sz,rhs)
        call oclWrite1DFloatArrayBuffer(global_rhsav_array_buf,global_rhsav_array_sz,global_rhsav_array)
        call oclWrite1DFloatArrayBuffer(global_area_array_buf,global_area_array_sz,global_area_array)
        call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)

        call runOcl(oclGlobalRange,oclLocalRange,exectime)
        ! call press_reduce_78
        call oclRead1DFloatArrayBuffer(global_rhsav_array_buf,global_rhsav_array_sz,global_rhsav_array)
        call oclRead1DFloatArrayBuffer(global_area_array_buf,global_area_array_sz,global_area_array)

        rhsav = (rhsav + global_rhsav_array(1))
        area = (area + global_area_array(1))

! ---- END --------------------------------------------------------------------------------------------------------------------
    rhsav = rhsav/area
! ---- BEGIN press_map_87 -----------------------------------------------------------------------------------------------------
        oclGlobalRange = (80 * (300 * 300))
        oclLocalRange = 0
        state_ptr(1) = ST_PRESS_MAP_87
        
        rhsav_ptr(1) = rhsav
        call oclWrite1DFloatArrayBuffer(rhsav_buf,rhsav_ptr_sz,rhsav_ptr)! Automatic conversion to array
        call oclWrite3DFloatArrayBuffer(rhs_buf,rhs_sz,rhs)
        call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
        call runOcl(oclGlobalRange,oclLocalRange,exectime)
        ! call press_map_87

        call oclRead3DFloatArrayBuffer(rhs_buf,rhs_sz,rhs)
! ---- END --------------------------------------------------------------------------------------------------------------------
        do l=1, 50
        sor = 0.0
            do nrd=0, 1
! ---- BEGIN press_map_97 -----------------------------------------------------------------------------------------------------
                oclGlobalRange = (80 * (300 * 300))
                oclLocalRange = 0
                state_ptr(1) = ST_PRESS_MAP_97
                
                call oclWrite1DFloatArrayBuffer(dzs_buf,dzs_sz,dzs)
                call oclWrite1DFloatArrayBuffer(dys_buf,dys_sz,dys)
                call oclWrite1DFloatArrayBuffer(dxs_buf,dxs_sz,dxs)
                nrd_ptr(1) = nrd
                call oclWrite1DIntArrayBuffer(nrd_buf,nrd_ptr_sz,nrd_ptr)! Automatic conversion to array
                call oclWrite3DFloatArrayBuffer(rhs_buf,rhs_sz,rhs)
                call oclWrite3DFloatArrayBuffer(p1_buf,p1_sz,p1)
                call oclWrite3DFloatArrayBuffer(p0_buf,p0_sz,p0)
                call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
                call runOcl(oclGlobalRange,oclLocalRange,exectime)
                ! call press_map_97

                call oclRead3DFloatArrayBuffer(p1_buf,p1_sz,p1)
                call oclRead3DFloatArrayBuffer(p0_buf,p0_sz,p0)
! ---- END --------------------------------------------------------------------------------------------------------------------
            end do
        end do
    pav = 0.0
    pco = 0.0
! ---- BEGIN press_reduce_129 -------------------------------------------------------------------------------------------------
        oclGlobalRange = 1
        oclLocalRange = 1
        state_ptr(1) = ST_PRESS_REDUCE_129
        
        call oclWrite3DFloatArrayBuffer(p0_buf,p0_sz,p0)
        call oclWrite1DFloatArrayBuffer(dx1_buf,dx1_sz,dx1)
        call oclWrite1DFloatArrayBuffer(dy1_buf,dy1_sz,dy1)
        call oclWrite1DFloatArrayBuffer(dzn_buf,dzn_sz,dzn)
        call oclWrite1DFloatArrayBuffer(global_pav_array_buf,global_pav_array_sz,global_pav_array)
        call oclWrite1DFloatArrayBuffer(global_pco_array_buf,global_pco_array_sz,global_pco_array)
        call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)

        call runOcl(oclGlobalRange,oclLocalRange,exectime)
        ! call press_reduce_129
        call oclRead1DFloatArrayBuffer(global_pav_array_buf,global_pav_array_sz,global_pav_array)
        call oclRead1DFloatArrayBuffer(global_pco_array_buf,global_pco_array_sz,global_pco_array)

        pav = (pav + global_pav_array(1))
        pco = (pco + global_pco_array(1))

! ---- END --------------------------------------------------------------------------------------------------------------------
    pav = pav/pco
! ---- BEGIN press_map_138 ----------------------------------------------------------------------------------------------------
        oclGlobalRange = (80 * (300 * 300))
        oclLocalRange = 0
        state_ptr(1) = ST_PRESS_MAP_138
        
        pav_ptr(1) = pav
        call oclWrite1DFloatArrayBuffer(pav_buf,pav_ptr_sz,pav_ptr)! Automatic conversion to array
        call oclWrite3DFloatArrayBuffer(p0_buf,p0_sz,p0)
        call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
        call runOcl(oclGlobalRange,oclLocalRange,exectime)
        ! call press_map_138

        call oclRead3DFloatArrayBuffer(p0_buf,p0_sz,p0)
! ---- END --------------------------------------------------------------------------------------------------------------------
end subroutine press
! Footer (produceCode_progUnit c)
end module module_press
