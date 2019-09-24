module module_bondv1
 contains 
subroutine bondv1(u,z2,dzn,v,w,n,n0,dt,dxs)

            use module_adam_bondv1_feedbf_les_press_v_etc_superkernel_init
            use oclWrapper
!      implicit none
      use params_common_sn

      implicit none

! otherStatements

! remainingDecls
      real(4), intent(in) :: dt !!
      real(4), dimension(0:ip), intent(in) :: dxs !!
      real(4), dimension(-1:kp+2), intent(in) :: dzn !!
      integer, intent(in) :: n !!
      integer, intent(in) :: n0 !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(inout) :: u !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(inout) :: v !!
      real(4), dimension(0:ip+1,-1:jp+1,-1:kp+1), intent(inout) :: w !!
      real(4), dimension(0:kp+2), intent(in) :: z2 !!
      real(4) :: u_val !!
      integer :: i !!
      integer :: j !!
      integer :: k !!
      real(4) :: aaa !!
      real(4) :: bbb !!
      real(4) :: uout !!
      real(4) :: gaaa !!
      real(4) :: gbbb !!

            ! Extra declarations
            real (kind=4) :: exectime
            real(kind=4), dimension(1:NUNITS) :: global_bbb_array
            real(kind=4), dimension(1:NUNITS) :: global_aaa_array

            ! Buffer declarations
            integer(8) :: state_ptr_buf
            integer(8) :: global_aaa_array_buf
            integer(8) :: global_bbb_array_buf
            integer(8) :: uout_buf

            integer, dimension(1) :: state_ptr

            ! Size declarations
            integer, dimension(1) :: state_ptr_sz
            integer, dimension(1) :: uout_ptr_sz
            integer, dimension(1) :: global_aaa_array_sz
            integer, dimension(1) :: global_bbb_array_sz
            real(kind=4), dimension(1) :: uout_ptr

            Integer :: r_iter
            
! Size assignments
            state_ptr_sz = shape(state_ptr)
            uout_ptr_sz = shape(uout_ptr)
            global_aaa_array_sz = shape(global_aaa_array)
            global_bbb_array_sz = shape(global_bbb_array)

            ! Buffer loads
            call oclLoadBuffer(STATE_PTR_BUF_IDX, state_ptr_buf)
            call oclLoadBuffer(GLOBAL_AAA_ARRAY_BUF_IDX, global_aaa_array_buf)
            call oclLoadBuffer(GLOBAL_BBB_ARRAY_BUF_IDX, global_bbb_array_buf)
            call oclLoadBuffer(UOUT_BUF_IDX, uout_buf)

            ! Original code with buffer writes and reads
! ---- BEGIN bondv1_map_38 ----------------------------------------------------------------------------------------------------
            oclGlobalRange = (((1 - 0) + 1) * (78 * 300))
            oclLocalRange = 0
            state_ptr(1) = ST_BONDV1_MAP_38
            
            call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
            call runOcl(oclGlobalRange,oclLocalRange,exectime)
            ! call bondv1_map_38

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN bondv1_map_48 ----------------------------------------------------------------------------------------------------
            oclGlobalRange = (((1 - 0) + 1) * (((80 - 79) + 1) * 300))
            oclLocalRange = 0
            state_ptr(1) = ST_BONDV1_MAP_48
            
            call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
            call runOcl(oclGlobalRange,oclLocalRange,exectime)
            ! call bondv1_map_48

! ---- END --------------------------------------------------------------------------------------------------------------------
            if ((n == n0)) then
! ---- BEGIN bondv1_map_58 ----------------------------------------------------------------------------------------------------
                oclGlobalRange = (80 * (300 * ((300 - 2) + 1)))
                oclLocalRange = 0
                state_ptr(1) = ST_BONDV1_MAP_58
                
                call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
                call runOcl(oclGlobalRange,oclLocalRange,exectime)
                ! call bondv1_map_58

! ---- END --------------------------------------------------------------------------------------------------------------------
            end if
    aaa = 0.0
! ---- BEGIN bondv1_reduce_69 -------------------------------------------------------------------------------------------------
            oclGlobalRange = 1
            oclLocalRange = 1
            state_ptr(1) = ST_BONDV1_REDUCE_69
            
            call oclWrite1DFloatArrayBuffer(global_aaa_array_buf,global_aaa_array_sz,global_aaa_array)
            call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)

            call runOcl(oclGlobalRange,oclLocalRange,exectime)
            ! call bondv1_reduce_69
            call oclRead1DFloatArrayBuffer(global_aaa_array_buf,global_aaa_array_sz,global_aaa_array)

            aaa = amax1(aaa,global_aaa_array(1))

! ---- END --------------------------------------------------------------------------------------------------------------------
    gaaa = aaa
    bbb = 1e38 
! ---- BEGIN bondv1_reduce_76 -------------------------------------------------------------------------------------------------
            oclGlobalRange = 1
            oclLocalRange = 1
            state_ptr(1) = ST_BONDV1_REDUCE_76
            
            call oclWrite1DFloatArrayBuffer(global_bbb_array_buf,global_bbb_array_sz,global_bbb_array)
            call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)

            call runOcl(oclGlobalRange,oclLocalRange,exectime)
            ! call bondv1_reduce_76
            call oclRead1DFloatArrayBuffer(global_bbb_array_buf,global_bbb_array_sz,global_bbb_array)

            bbb = amin1(bbb,global_bbb_array(1))

! ---- END --------------------------------------------------------------------------------------------------------------------
    gbbb=bbb
    uout = (gaaa+gbbb)/2.
! ---- BEGIN bondv1_map_83 ----------------------------------------------------------------------------------------------------
            oclGlobalRange = (80 * 300)
            oclLocalRange = 0
            state_ptr(1) = ST_BONDV1_MAP_83
            
            uout_ptr(1) = uout
            call oclWrite1DFloatArrayBuffer(uout_buf,uout_ptr_sz,uout_ptr)! Automatic conversion to array
            call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
            call runOcl(oclGlobalRange,oclLocalRange,exectime)
            ! call bondv1_map_83

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN bondv1_map_98 ----------------------------------------------------------------------------------------------------
            oclGlobalRange = ((((80 + 1) - 0) + 1) * (((300 + 1) - 0) + 1))
            oclLocalRange = 0
            state_ptr(1) = ST_BONDV1_MAP_98
            
            call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
            call runOcl(oclGlobalRange,oclLocalRange,exectime)
            ! call bondv1_map_98

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN bondv1_map_116 ---------------------------------------------------------------------------------------------------
            oclGlobalRange = ((((300 + 1) - 0) + 1) * (((300 + 1) - 0) + 1))
            oclLocalRange = 0
            state_ptr(1) = ST_BONDV1_MAP_116
            
            call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
            call runOcl(oclGlobalRange,oclLocalRange,exectime)
            ! call bondv1_map_116

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN bondv1_map_128 ---------------------------------------------------------------------------------------------------
            oclGlobalRange = ((((300 + 1) - (-1)) + 1) * (((300 + 1) - 0) + 1))
            oclLocalRange = 0
            state_ptr(1) = ST_BONDV1_MAP_128
            
            call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
            call runOcl(oclGlobalRange,oclLocalRange,exectime)
            ! call bondv1_map_128

! ---- END --------------------------------------------------------------------------------------------------------------------
end subroutine bondv1
! Footer (produceCode_progUnit c)
end module module_bondv1
