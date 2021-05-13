module module_les
      use module_boundsm 
 contains 
        subroutine les(u,v,w,f,g,h,sm,dx1,dy1,dzn,dzs,dxs,dys)

      use module_adam_feedbf_les_press_velfg_ve_etc_superkernel_init
      use oclWrapper
!      implicit none
      use params_common_sn

      implicit none

! otherStatements

! remainingDecls
      real(4), dimension(kp) :: delx1 !!
      real(4) :: diu1_i_j_k !!
      real(4) :: diu2_i_j_k !!
      real(4) :: diu3_i_j_k !!
      real(4) :: diu4_i_j_k !!
      real(4) :: diu5_i_j_k !!
      real(4) :: diu6_i_j_k !!
      real(4) :: diu7_i_j_k !!
      real(4) :: diu8_i_j_k !!
      real(4) :: diu9_i_j_k !!
      real(4) :: diu1_ip1_j_k !!
      real(4) :: diu2_im1_j_k !!
      real(4) :: diu2_im1_jp1_k !!
      real(4) :: diu2_i_jp1_k !!
      real(4) :: diu3_im1_j_k !!
      real(4) :: diu3_im1_j_kp1 !!
      real(4) :: diu3_i_j_kp1 !!
      real(4) :: diu4_i_jm1_k !!
      real(4) :: diu4_ip1_j_k !!
      real(4) :: diu4_ip1_jm1_k !!
      real(4) :: diu5_i_jp1_k !!
      real(4) :: diu6_i_jm1_k !!
      real(4) :: diu6_i_jm1_kp1 !!
      real(4) :: diu6_i_j_kp1 !!
      real(4) :: diu7_i_j_km1 !!
      real(4) :: diu7_ip1_j_k !!
      real(4) :: diu7_ip1_j_km1 !!
      real(4) :: diu8_i_j_km1 !!
      real(4) :: diu8_i_jp1_k !!
      real(4) :: diu8_i_jp1_km1 !!
      real(4) :: diu9_i_j_kp1 !!
      real(4) :: diu1_i_j_1 !!
      real(4) :: diu1_ip1_j_1 !!
      real(4) :: diu2_i_j_1 !!
      real(4) :: diu2_i_jp1_1 !!
      real(4) :: diu4_ip1_j_1 !!
      real(4) :: diu2_im1_jp1_1 !!
      real(4) :: diu3_i_j_2 !!
      real(4) :: diu4_i_j_1 !!
      real(4) :: diu4_ip1_jm1_1 !!
      real(4) :: diu5_i_jp1_1 !!
      real(4) :: diu5_i_j_1 !!
      real(4) :: diu6_i_j_2 !!
      real(4) :: diu7_ip1_j_1 !!
      real(4) :: diu8_i_jp1_1 !!
      real(4), dimension(-1:kp+2) :: dzs !!
      real(4), dimension(-1:ip+1), intent(in) :: dx1 !!
      real(4), dimension(0:jp+1), intent(in) :: dy1 !!
      real(4), dimension(-1:kp+2), intent(in) :: dzn !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: f !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: g !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: h !!
      real(4), dimension(-1:ip+1,-1:jp+1,0:kp+1), intent(out) :: sm !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(in) :: u !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(in) :: v !!
      real(4), dimension(0:ip+1,-1:jp+1,-1:kp+1), intent(in) :: w !!
      real(4), dimension(0:ip), intent(in) :: dxs !!
      real(4), dimension(0:jp), intent(in) :: dys !!
      integer :: i !!
      integer :: j !!
      integer :: k !!
      real(4) :: csx1 !!
      real(4) :: dudxx1 !!
      real(4) :: dudyx1 !!
      real(4) :: dudzx1 !!
      real(4) :: dvdxx1 !!
      real(4) :: dvdyx1 !!
      real(4) :: dvdzx1 !!
      real(4) :: dwdxx1 !!
      real(4) :: dwdyx1 !!
      real(4) :: dwdzx1 !!
      real(4) :: visux2 !!
      real(4) :: visux1 !!
      real(4) :: visuy2 !!
      real(4) :: visuy1 !!
      real(4) :: visuz2 !!
      real(4) :: visuz1 !!
      real(4) :: visvx2 !!
      real(4) :: visvx1 !!
      real(4) :: visvy2 !!
      real(4) :: visvy1 !!
      real(4) :: visvz2 !!
      real(4) :: visvz1 !!
      real(4) :: viswx2 !!
      real(4) :: viswx1 !!
      real(4) :: viswy2 !!
      real(4) :: viswy1 !!
      real(4) :: viswz2 !!
      real(4) :: viswz1 !!
      real(4) :: evsx2 !!
      real(4) :: evsx1 !!
      real(4) :: evsy2 !!
      real(4) :: evsy1 !!
      real(4) :: evsz2 !!
      real(4) :: evsz1 !!
      real(4) :: vfu !!
      real(4) :: vfv !!
      real(4) :: vfw !!

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
! ---- BEGIN les_map_119 ------------------------------------------------------------------------------------------------------
      oclGlobalRange = (80 * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_LES_MAP_119
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call les_map_119

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN les_map_164 ------------------------------------------------------------------------------------------------------
      oclGlobalRange = (((80 - 2) + 1) * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_LES_MAP_164
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call les_map_164

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN les_map_202 ------------------------------------------------------------------------------------------------------
      oclGlobalRange = (((80 - 2) + 1) * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_LES_MAP_202
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call les_map_202

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN les_map_240 ------------------------------------------------------------------------------------------------------
      oclGlobalRange = (80 * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_LES_MAP_240
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call les_map_240

! ---- END --------------------------------------------------------------------------------------------------------------------
      return
      end subroutine les
! Footer (produceCode_progUnit c)
end module module_les
