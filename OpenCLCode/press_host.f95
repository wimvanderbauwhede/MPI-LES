module module_press
    use module_bondFG 
    use module_boundp 
 contains 
subroutine press(u,v,w,p,rhs,f,g,h,dx1,dy1,dzn,dxs,dys,dzs,dt,n,nmax &
)

      use module_adam_bondv1_feedbf_les_press_v_etc_superkernel_init
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
      real(4), dimension(0:1,0:ip+2,0:jp+2,0:kp+1) :: p !!
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
      integer(8) :: w_buf
      integer(8) :: global_rhsav_array_buf
      integer(8) :: global_area_array_buf
      integer(8) :: rhsav_buf
      integer(8) :: dzs_buf
      integer(8) :: nrd_buf
      integer(8) :: global_pav_array_buf
      integer(8) :: global_pco_array_buf
      integer(8) :: pav_buf

      integer, dimension(1) :: state_ptr

      ! Size declarations
      integer, dimension(1) :: state_ptr_sz
      integer, dimension(3) :: w_sz
      integer, dimension(1) :: rhsav_ptr_sz
      integer, dimension(1) :: dzs_sz
      integer, dimension(1) :: nrd_ptr_sz
      integer, dimension(1) :: pav_ptr_sz
      integer, dimension(1) :: global_rhsav_array_sz
      integer, dimension(1) :: global_area_array_sz
      integer, dimension(1) :: global_pav_array_sz
      integer, dimension(1) :: global_pco_array_sz
      real(kind=4), dimension(1) :: rhsav_ptr
      integer, dimension(1) :: nrd_ptr
      real(kind=4), dimension(1) :: pav_ptr

      Integer :: r_iter
      
! Size assignments
      state_ptr_sz = shape(state_ptr)
      w_sz = shape(w)
      rhsav_ptr_sz = shape(rhsav_ptr)
      dzs_sz = shape(dzs)
      nrd_ptr_sz = shape(nrd_ptr)
      pav_ptr_sz = shape(pav_ptr)
      global_rhsav_array_sz = shape(global_rhsav_array)
      global_area_array_sz = shape(global_area_array)
      global_pav_array_sz = shape(global_pav_array)
      global_pco_array_sz = shape(global_pco_array)

      ! Buffer loads
      call oclLoadBuffer(STATE_PTR_BUF_IDX, state_ptr_buf)
      call oclLoadBuffer(W_BUF_IDX, w_buf)
      call oclLoadBuffer(GLOBAL_RHSAV_ARRAY_BUF_IDX, global_rhsav_array_buf)
      call oclLoadBuffer(GLOBAL_AREA_ARRAY_BUF_IDX, global_area_array_buf)
      call oclLoadBuffer(RHSAV_BUF_IDX, rhsav_buf)
      call oclLoadBuffer(DZS_BUF_IDX, dzs_buf)
      call oclLoadBuffer(NRD_BUF_IDX, nrd_buf)
      call oclLoadBuffer(GLOBAL_PAV_ARRAY_BUF_IDX, global_pav_array_buf)
      call oclLoadBuffer(GLOBAL_PCO_ARRAY_BUF_IDX, global_pco_array_buf)
      call oclLoadBuffer(PAV_BUF_IDX, pav_buf)

      ! Original code with buffer writes and reads
! ---- BEGIN press_map_64 -----------------------------------------------------------------------------------------------------
      oclGlobalRange = (80 * 300)
      oclLocalRange = 0
      state_ptr(1) = ST_PRESS_MAP_64
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call press_map_64

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN press_map_69 -----------------------------------------------------------------------------------------------------
      oclGlobalRange = (80 * 300)
      oclLocalRange = 0
      state_ptr(1) = ST_PRESS_MAP_69
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call press_map_69

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN press_map_74 -----------------------------------------------------------------------------------------------------
      oclGlobalRange = (300 * 300)
      oclLocalRange = 0
      state_ptr(1) = ST_PRESS_MAP_74
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call press_map_74

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN press_map_80 -----------------------------------------------------------------------------------------------------
      oclGlobalRange = (80 * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_PRESS_MAP_80
      
      call oclWrite3DFloatArrayBuffer(w_buf,w_sz,w)
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call press_map_80

! ---- END --------------------------------------------------------------------------------------------------------------------
    rhsav = 0.0
    area = 0.0
! ---- BEGIN press_reduce_93 --------------------------------------------------------------------------------------------------
      oclGlobalRange = 1
      oclLocalRange = 1
      state_ptr(1) = ST_PRESS_REDUCE_93
      
      call oclWrite1DFloatArrayBuffer(global_rhsav_array_buf,global_rhsav_array_sz,global_rhsav_array)
      call oclWrite1DFloatArrayBuffer(global_area_array_buf,global_area_array_sz,global_area_array)
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)

      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call press_reduce_93
      call oclRead1DFloatArrayBuffer(global_rhsav_array_buf,global_rhsav_array_sz,global_rhsav_array)
      call oclRead1DFloatArrayBuffer(global_area_array_buf,global_area_array_sz,global_area_array)

      rhsav = (rhsav + global_rhsav_array(1))
      area = (area + global_area_array(1))

! ---- END --------------------------------------------------------------------------------------------------------------------
    rhsav = rhsav/area
! ---- BEGIN press_map_102 ----------------------------------------------------------------------------------------------------
      oclGlobalRange = (80 * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_PRESS_MAP_102
      
      rhsav_ptr(1) = rhsav
      call oclWrite1DFloatArrayBuffer(rhsav_buf,rhsav_ptr_sz,rhsav_ptr)! Automatic conversion to array
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call press_map_102

! ---- END --------------------------------------------------------------------------------------------------------------------
      do l=1, 50
        sor = 0.0
          do nrd=0, 1
! ---- BEGIN press_map_112 ----------------------------------------------------------------------------------------------------
              oclGlobalRange = (80 * (300 * 300))
              oclLocalRange = 0
              state_ptr(1) = ST_PRESS_MAP_112
              
              call oclWrite1DFloatArrayBuffer(dzs_buf,dzs_sz,dzs)
              nrd_ptr(1) = nrd
              call oclWrite1DIntArrayBuffer(nrd_buf,nrd_ptr_sz,nrd_ptr)! Automatic conversion to array
              call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
              call runOcl(oclGlobalRange,oclLocalRange,exectime)
              ! call press_map_112

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN press_map_140 ----------------------------------------------------------------------------------------------------
              oclGlobalRange = ((((80 + 1) - 0) + 1) * (((300 + 1) - 0) + 1))
              oclLocalRange = 0
              state_ptr(1) = ST_PRESS_MAP_140
              
              call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
              call runOcl(oclGlobalRange,oclLocalRange,exectime)
              ! call press_map_140

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN press_map_146 ----------------------------------------------------------------------------------------------------
              oclGlobalRange = ((((80 + 1) - 0) + 1) * (((300 + 1) - 0) + 1))
              oclLocalRange = 0
              state_ptr(1) = ST_PRESS_MAP_146
              
              call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
              call runOcl(oclGlobalRange,oclLocalRange,exectime)
              ! call press_map_146

! ---- END --------------------------------------------------------------------------------------------------------------------
          end do
! ---- BEGIN press_map_153 ----------------------------------------------------------------------------------------------------
          oclGlobalRange = ((((300 + 1) - 0) + 1) * (((300 + 1) - 0) + 1))
          oclLocalRange = 0
          state_ptr(1) = ST_PRESS_MAP_153
          
          call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
          call runOcl(oclGlobalRange,oclLocalRange,exectime)
          ! call press_map_153

! ---- END --------------------------------------------------------------------------------------------------------------------
      end do
    pav = 0.0
    pco = 0.0
! ---- BEGIN press_reduce_162 -------------------------------------------------------------------------------------------------
      oclGlobalRange = 1
      oclLocalRange = 1
      state_ptr(1) = ST_PRESS_REDUCE_162
      
      call oclWrite1DFloatArrayBuffer(global_pav_array_buf,global_pav_array_sz,global_pav_array)
      call oclWrite1DFloatArrayBuffer(global_pco_array_buf,global_pco_array_sz,global_pco_array)
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)

      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call press_reduce_162
      call oclRead1DFloatArrayBuffer(global_pav_array_buf,global_pav_array_sz,global_pav_array)
      call oclRead1DFloatArrayBuffer(global_pco_array_buf,global_pco_array_sz,global_pco_array)

      pav = (pav + global_pav_array(1))
      pco = (pco + global_pco_array(1))

! ---- END --------------------------------------------------------------------------------------------------------------------
    pav = pav/pco
! ---- BEGIN press_map_171 ----------------------------------------------------------------------------------------------------
      oclGlobalRange = (80 * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_PRESS_MAP_171
      
      pav_ptr(1) = pav
      call oclWrite1DFloatArrayBuffer(pav_buf,pav_ptr_sz,pav_ptr)! Automatic conversion to array
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call press_map_171

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN press_map_178 ----------------------------------------------------------------------------------------------------
      oclGlobalRange = ((((80 + 1) - 0) + 1) * (((300 + 1) - 0) + 1))
      oclLocalRange = 0
      state_ptr(1) = ST_PRESS_MAP_178
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call press_map_178

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN press_map_184 ----------------------------------------------------------------------------------------------------
      oclGlobalRange = ((((80 + 1) - 0) + 1) * (((300 + 1) - 0) + 1))
      oclLocalRange = 0
      state_ptr(1) = ST_PRESS_MAP_184
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call press_map_184

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN press_map_190 ----------------------------------------------------------------------------------------------------
      oclGlobalRange = ((((300 + 1) - 0) + 1) * (((300 + 1) - 0) + 1))
      oclLocalRange = 0
      state_ptr(1) = ST_PRESS_MAP_190
      
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call press_map_190

! ---- END --------------------------------------------------------------------------------------------------------------------
end subroutine press
! Footer (produceCode_progUnit c)
end module module_press
