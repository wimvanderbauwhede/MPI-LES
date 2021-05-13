module module_velFG
 contains 
      subroutine velfg(dx1,dy1,dzn,f, &
      g,dzs,h,u,v,w)

      use module_adam_feedbf_les_press_velfg_ve_etc_superkernel_init
      use oclWrapper
!      implicit none
      use params_common_sn

      implicit none

! otherStatements

! remainingDecls
      real(4), dimension(-1:ip+1), intent(in) :: dx1 !!
      real(4), dimension(0:jp+1), intent(in) :: dy1 !!
      real(4), dimension(-1:kp+2), intent(in) :: dzn !!
      real(4), dimension(-1:kp+2), intent(in) :: dzs !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(out) :: f !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(out) :: g !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(out) :: h !!
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
      real(4) :: uspd !!
      real(4) :: vspd !!
      integer :: i !!
      integer :: j !!
      integer :: k !!
      real(4) :: covc !!
      real(4) :: covx1 !!
      real(4) :: covy1 !!
      real(4) :: covz1 !!
      integer, parameter :: u0 = 0  !!

      ! Extra declarations
      real (kind=4) :: exectime

      ! Buffer declarations
      integer(8) :: state_ptr_buf

      integer, dimension(1) :: state_ptr

      ! Size declarations
      integer, dimension(1) :: state_ptr_sz

      
! Size assignments
      state_ptr_sz = shape(state_ptr)

      ! Buffer loads
      call oclLoadBuffer(STATE_PTR_BUF_IDX, state_ptr_buf)

      ! Original code with buffer writes and reads
! ---- BEGIN velFG_map_95 -----------------------------------------------------------------------------------------------------
      oclGlobalRange = (80 * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_VELFG_MAP_95
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call velFG_map_95

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN velFG_map_135 ----------------------------------------------------------------------------------------------------
      oclGlobalRange = (80 * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_VELFG_MAP_135
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call velFG_map_135

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN velFG_map_175 ----------------------------------------------------------------------------------------------------
      oclGlobalRange = ((80 - 1) * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_VELFG_MAP_175
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call velFG_map_175

! ---- END --------------------------------------------------------------------------------------------------------------------
      return
      end subroutine velFG
! Footer (produceCode_progUnit c)
end module module_velFG
