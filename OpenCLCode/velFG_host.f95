module module_velFG
    use module_vel2 
 contains 
        subroutine velfg(dx1,dy1,dzn,f,g,h,u,v,w, &
        dfu1,dfv1,dfw1,vn,dzs, &
        diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9, &
        cov1,cov2,cov3,cov4,cov5,cov6,cov7,cov8,cov9, &
        nou1,nou2,nou3,nou4,nou5,nou6,nou7,nou8,nou9, &
        uspd,vspd) 

            use module_adam_bondv1_feedbf_les_press_v_etc_superkernel_init
            use oclWrapper
!      implicit none
      use params_common_sn

      implicit none

! otherStatements

! remainingDecls
      real(4), dimension(-1:ip+2,0:jp+2,0:kp+2), intent(out) :: cov1 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: cov2 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: cov3 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: cov4 !!
      real(4), dimension(-1:ip+2,0:jp+2,0:kp+2), intent(out) :: cov5 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: cov6 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: cov7 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: cov8 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: cov9 !!
      real(4), dimension(0:ip,jp,kp), intent(out) :: dfu1 !!
      real(4), dimension(ip,0:jp,kp), intent(out) :: dfv1 !!
      real(4), dimension(ip,jp,kp), intent(out) :: dfw1 !!
      real(4), intent(in) :: vn !!
      real(4), dimension(-1:ip+2,0:jp+2,0:kp+2), intent(out) :: diu1 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: diu2 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: diu3 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: diu4 !!
      real(4), dimension(-1:ip+2,0:jp+2,0:kp+2), intent(out) :: diu5 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: diu6 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: diu7 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: diu8 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: diu9 !!
      real(4), dimension(-1:ip+1), intent(in) :: dx1 !!
      real(4), dimension(0:jp+1), intent(in) :: dy1 !!
      real(4), dimension(-1:kp+2), intent(in) :: dzn !!
      real(4), dimension(-1:kp+2), intent(in) :: dzs !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(out) :: f !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(out) :: g !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(out) :: h !!
      real(4), dimension(-1:ip+2,0:jp+2,0:kp+2), intent(out) :: nou1 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: nou2 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: nou3 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: nou4 !!
      real(4), dimension(-1:ip+2,0:jp+2,0:kp+2), intent(out) :: nou5 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: nou6 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: nou7 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: nou8 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(out) :: nou9 !!
      real(4) :: nou1_ !!
      real(4) :: nou2_ !!
      real(4) :: nou3_ !!
      real(4) :: nou4_ !!
      real(4) :: nou5_ !!
      real(4) :: nou6_ !!
      real(4) :: nou7_ !!
      real(4) :: nou8_ !!
      real(4) :: nou9_ !!
      real(4) :: diu1_ !!
      real(4) :: diu2_ !!
      real(4) :: diu3_ !!
      real(4) :: diu4_ !!
      real(4) :: diu5_ !!
      real(4) :: diu6_ !!
      real(4) :: diu7_ !!
      real(4) :: diu8_ !!
      real(4) :: diu9_ !!
      real(4) :: cov1_i !!
      real(4) :: cov2_j !!
      real(4) :: cov3_k !!
      real(4) :: cov4_i !!
      real(4) :: cov5_j !!
      real(4) :: cov6_k !!
      real(4) :: cov7_i !!
      real(4) :: cov8_j !!
      real(4) :: cov9_k !!
      real(4) :: nou1_ip1 !!
      real(4) :: nou2_jp1 !!
      real(4) :: nou3_kp1 !!
      real(4) :: nou4_ip1 !!
      real(4) :: nou5_jp1 !!
      real(4) :: nou6_kp1 !!
      real(4) :: nou7_ip1 !!
      real(4) :: nou8_jp1 !!
      real(4) :: nou9_kp1 !!
      real(4) :: diu1_ip1 !!
      real(4) :: diu2_jp1 !!
      real(4) :: diu3_kp1 !!
      real(4) :: diu4_ip1 !!
      real(4) :: diu5_jp1 !!
      real(4) :: diu6_kp1 !!
      real(4) :: diu7_ip1 !!
      real(4) :: diu8_jp1 !!
      real(4) :: diu9_kp1 !!
      real(4) :: cov1_ip1 !!
      real(4) :: cov2_jp1 !!
      real(4) :: cov3_kp1 !!
      real(4) :: cov4_ip1 !!
      real(4) :: cov5_jp1 !!
      real(4) :: cov6_kp1 !!
      real(4) :: cov7_ip1 !!
      real(4) :: cov8_jp1 !!
      real(4) :: cov9_kp1 !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(in) :: u !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(in) :: v !!
      real(4), dimension(0:ip+1,-1:jp+1,-1:kp+1), intent(in) :: w !!
      real(4), dimension(0:ip+1,0:jp+1), intent(out) :: uspd !!
      real(4), dimension(0:ip+1,0:jp+1), intent(out) :: vspd !!
      integer :: i !!
      integer :: j !!
      integer :: k !!
      real(4) :: covc !!
      real(4) :: covx1 !!
      real(4) :: covy1 !!
      real(4) :: covz1 !!

            ! Extra declarations
            real (kind=4) :: exectime

            ! Buffer declarations
            integer(8) :: state_ptr_buf
            integer(8) :: u_buf
            integer(8) :: v_buf
            integer(8) :: dx1_buf
            integer(8) :: dy1_buf
            integer(8) :: cov1_buf
            integer(8) :: cov2_buf
            integer(8) :: cov3_buf
            integer(8) :: cov4_buf
            integer(8) :: cov5_buf
            integer(8) :: cov6_buf
            integer(8) :: cov7_buf
            integer(8) :: cov8_buf
            integer(8) :: dzn_buf
            integer(8) :: cov9_buf

            integer, dimension(1) :: state_ptr

            ! Size declarations
            integer, dimension(1) :: state_ptr_sz
            integer, dimension(3) :: u_sz
            integer, dimension(3) :: v_sz
            integer, dimension(1) :: dx1_sz
            integer, dimension(1) :: dy1_sz
            integer, dimension(3) :: cov1_sz
            integer, dimension(3) :: cov2_sz
            integer, dimension(3) :: cov3_sz
            integer, dimension(3) :: cov4_sz
            integer, dimension(3) :: cov5_sz
            integer, dimension(3) :: cov6_sz
            integer, dimension(3) :: cov7_sz
            integer, dimension(3) :: cov8_sz
            integer, dimension(1) :: dzn_sz
            integer, dimension(3) :: cov9_sz

            
! Size assignments
            state_ptr_sz = shape(state_ptr)
            u_sz = shape(u)
            v_sz = shape(v)
            dx1_sz = shape(dx1)
            dy1_sz = shape(dy1)
            cov1_sz = shape(cov1)
            cov2_sz = shape(cov2)
            cov3_sz = shape(cov3)
            cov4_sz = shape(cov4)
            cov5_sz = shape(cov5)
            cov6_sz = shape(cov6)
            cov7_sz = shape(cov7)
            cov8_sz = shape(cov8)
            dzn_sz = shape(dzn)
            cov9_sz = shape(cov9)

            ! Buffer loads
            call oclLoadBuffer(STATE_PTR_BUF_IDX, state_ptr_buf)
            call oclLoadBuffer(U_BUF_IDX, u_buf)
            call oclLoadBuffer(V_BUF_IDX, v_buf)
            call oclLoadBuffer(DX1_BUF_IDX, dx1_buf)
            call oclLoadBuffer(DY1_BUF_IDX, dy1_buf)
            call oclLoadBuffer(COV1_BUF_IDX, cov1_buf)
            call oclLoadBuffer(COV2_BUF_IDX, cov2_buf)
            call oclLoadBuffer(COV3_BUF_IDX, cov3_buf)
            call oclLoadBuffer(COV4_BUF_IDX, cov4_buf)
            call oclLoadBuffer(COV5_BUF_IDX, cov5_buf)
            call oclLoadBuffer(COV6_BUF_IDX, cov6_buf)
            call oclLoadBuffer(COV7_BUF_IDX, cov7_buf)
            call oclLoadBuffer(COV8_BUF_IDX, cov8_buf)
            call oclLoadBuffer(DZN_BUF_IDX, dzn_buf)
            call oclLoadBuffer(COV9_BUF_IDX, cov9_buf)

            ! Original code with buffer writes and reads
      call vel2( &
            nou1,nou5,nou9,nou2,nou3,nou4,nou6,nou7,nou8,&
            diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,&
            cov1,cov2,cov3,cov4,cov5,cov6,cov7,cov8,cov9,&
            u,v,w,dx1,dy1,dzn,dzs,uspd,vspd)
! ---- BEGIN velFG_map_135 ----------------------------------------------------------------------------------------------------
            oclGlobalRange = (300 * 300)
            oclLocalRange = 0
            state_ptr(1) = ST_VELFG_MAP_135
            
            call oclWrite3DFloatArrayBuffer(u_buf,u_sz,u)
            call oclWrite3DFloatArrayBuffer(v_buf,v_sz,v)
            call oclWrite1DFloatArrayBuffer(dx1_buf,dx1_sz,dx1)
            call oclWrite1DFloatArrayBuffer(dy1_buf,dy1_sz,dy1)
            call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
            call runOcl(oclGlobalRange,oclLocalRange,exectime)
            ! call velFG_map_135

! ---- END --------------------------------------------------------------------------------------------------------------------
        write(6,*) 'CHK_uspd_vspd=',uspd(ip/2,jp/2),vspd(ip/2,jp/2)
! ---- BEGIN velFG_map_148 ----------------------------------------------------------------------------------------------------
            oclGlobalRange = (80 * (300 * 300))
            oclLocalRange = 0
            state_ptr(1) = ST_VELFG_MAP_148
            
            call oclWrite3DFloatArrayBuffer(cov1_buf,cov1_sz,cov1)
            call oclWrite3DFloatArrayBuffer(cov2_buf,cov2_sz,cov2)
            call oclWrite3DFloatArrayBuffer(cov3_buf,cov3_sz,cov3)
            call oclWrite3DFloatArrayBuffer(cov4_buf,cov4_sz,cov4)
            call oclWrite3DFloatArrayBuffer(cov5_buf,cov5_sz,cov5)
            call oclWrite3DFloatArrayBuffer(cov6_buf,cov6_sz,cov6)
            call oclWrite3DFloatArrayBuffer(cov7_buf,cov7_sz,cov7)
            call oclWrite3DFloatArrayBuffer(cov8_buf,cov8_sz,cov8)
            call oclWrite1DFloatArrayBuffer(dzn_buf,dzn_sz,dzn)
            call oclWrite3DFloatArrayBuffer(cov9_buf,cov9_sz,cov9)
            call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
            call runOcl(oclGlobalRange,oclLocalRange,exectime)
            ! call velFG_map_148

! ---- END --------------------------------------------------------------------------------------------------------------------
      return
      end subroutine velFG
! Footer (produceCode_progUnit c)
end module module_velFG
