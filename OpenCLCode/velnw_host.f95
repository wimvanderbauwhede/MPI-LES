module module_velnw
 contains 
      subroutine velnw(p,ro,dxs,u,dt,f,dys,v,g,dzs,w,h)

      use module_adam_bondv1_feedbf_les_press_v_etc_superkernel_init
      use oclWrapper
!      implicit none
      use params_common_sn

      implicit none

! otherStatements

! remainingDecls
      real(4), intent(in) :: dt !!
      real(4), dimension(0:ip), intent(in) :: dxs !!
      real(4), dimension(0:jp), intent(in) :: dys !!
      real(4), dimension(-1:kp+2), intent(in) :: dzs !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(in) :: f !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(in) :: g !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(in) :: h !!
      real(4), dimension(0:1,0:ip+2,0:jp+2,0:kp+1), intent(in) :: p !!
      real(4), intent(in) :: ro !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(inout) :: u !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(inout) :: v !!
      real(4), dimension(0:ip+1,-1:jp+1,-1:kp+1), intent(inout) :: w !!
      integer :: i !!
      integer :: j !!
      integer :: k !!
      real(4) :: pz !!

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
! ---- BEGIN velnw_map_36 -----------------------------------------------------------------------------------------------------
      oclGlobalRange = (80 * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_VELNW_MAP_36
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call velnw_map_36

! ---- END --------------------------------------------------------------------------------------------------------------------
      return
      end subroutine velnw
! Footer (produceCode_progUnit c)
end module module_velnw
