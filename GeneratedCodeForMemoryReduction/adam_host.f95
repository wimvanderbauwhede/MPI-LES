module module_adam
 contains 
subroutine adam(n,nmax,data21,fold,gold,hold,&
f,g,h)

        use module_adam_feedbf_les_press_velfg_ve_etc_superkernel_init
        use oclWrapper
!      implicit none
      use params_common_sn

      implicit none

! otherStatements

! remainingDecls
      character*70, intent(in) :: data21 !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: f !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: g !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: h !!
      real(4), dimension(ip,jp,kp), intent(inout) :: fold !!
      real(4), dimension(ip,jp,kp), intent(inout) :: gold !!
      real(4), dimension(ip,jp,kp), intent(inout) :: hold !!
      integer, intent(in) :: n !!
      integer, intent(in) :: nmax !!
      integer :: i !!
      integer :: j !!
      integer :: k !!
      real(4) :: fd !!
      real(4) :: gd !!
      real(4) :: hd !!

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
! ---- BEGIN adam_map_36 ------------------------------------------------------------------------------------------------------
        oclGlobalRange = (80 * (300 * 300))
        oclLocalRange = 0
        state_ptr(1) = ST_ADAM_MAP_36
        
        call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
        call runOcl(oclGlobalRange,oclLocalRange,exectime)
        ! call adam_map_36

! ---- END --------------------------------------------------------------------------------------------------------------------
end subroutine adam
! Footer (produceCode_progUnit c)
end module module_adam
