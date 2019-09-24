module module_les
      use module_boundsm 
 contains 
      subroutine les(delx1,dx1,dy1,dzn,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,sm,f,g, &
      h,u,v,uspd,vspd,dxs,dys,n)

          use module_adam_bondv1_feedbf_les_press_v_etc_superkernel_init
          use oclWrapper
!      implicit none
      use params_common_sn

      implicit none

! otherStatements

! remainingDecls
      real(4), dimension(kp), intent(out) :: delx1 !!
      real(4), dimension(-1:ip+2,0:jp+2,0:kp+2), intent(in) :: diu1 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(in) :: diu2 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(in) :: diu3 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(in) :: diu4 !!
      real(4), dimension(-1:ip+2,0:jp+2,0:kp+2), intent(in) :: diu5 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(in) :: diu6 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(in) :: diu7 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(in) :: diu8 !!
      real(4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(in) :: diu9 !!
      integer, intent(in) :: n !!
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
      real(4), dimension(0:ip+1,0:jp+1), intent(in) :: uspd !!
      real(4), dimension(0:ip+1,0:jp+1), intent(in) :: vspd !!
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
          integer(8) :: diu1_buf
          integer(8) :: diu2_buf
          integer(8) :: diu3_buf
          integer(8) :: diu4_buf
          integer(8) :: diu5_buf
          integer(8) :: diu6_buf
          integer(8) :: diu7_buf
          integer(8) :: diu8_buf
          integer(8) :: diu9_buf

          integer, dimension(1) :: state_ptr

          ! Size declarations
          integer, dimension(1) :: state_ptr_sz
          integer, dimension(3) :: diu1_sz
          integer, dimension(3) :: diu2_sz
          integer, dimension(3) :: diu3_sz
          integer, dimension(3) :: diu4_sz
          integer, dimension(3) :: diu5_sz
          integer, dimension(3) :: diu6_sz
          integer, dimension(3) :: diu7_sz
          integer, dimension(3) :: diu8_sz
          integer, dimension(3) :: diu9_sz

          
! Size assignments
          state_ptr_sz = shape(state_ptr)
          diu1_sz = shape(diu1)
          diu2_sz = shape(diu2)
          diu3_sz = shape(diu3)
          diu4_sz = shape(diu4)
          diu5_sz = shape(diu5)
          diu6_sz = shape(diu6)
          diu7_sz = shape(diu7)
          diu8_sz = shape(diu8)
          diu9_sz = shape(diu9)

          ! Buffer loads
          call oclLoadBuffer(STATE_PTR_BUF_IDX, state_ptr_buf)
          call oclLoadBuffer(DIU1_BUF_IDX, diu1_buf)
          call oclLoadBuffer(DIU2_BUF_IDX, diu2_buf)
          call oclLoadBuffer(DIU3_BUF_IDX, diu3_buf)
          call oclLoadBuffer(DIU4_BUF_IDX, diu4_buf)
          call oclLoadBuffer(DIU5_BUF_IDX, diu5_buf)
          call oclLoadBuffer(DIU6_BUF_IDX, diu6_buf)
          call oclLoadBuffer(DIU7_BUF_IDX, diu7_buf)
          call oclLoadBuffer(DIU8_BUF_IDX, diu8_buf)
          call oclLoadBuffer(DIU9_BUF_IDX, diu9_buf)

          ! Original code with buffer writes and reads
! ---- BEGIN les_map_88 -------------------------------------------------------------------------------------------------------
          oclGlobalRange = 80
          oclLocalRange = 0
          state_ptr(1) = ST_LES_MAP_88
          
          call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
          call runOcl(oclGlobalRange,oclLocalRange,exectime)
          ! call les_map_88

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN les_map_91 -------------------------------------------------------------------------------------------------------
          oclGlobalRange = (80 * (300 * 300))
          oclLocalRange = 0
          state_ptr(1) = ST_LES_MAP_91
          
          call oclWrite3DFloatArrayBuffer(diu1_buf,diu1_sz,diu1)
          call oclWrite3DFloatArrayBuffer(diu2_buf,diu2_sz,diu2)
          call oclWrite3DFloatArrayBuffer(diu3_buf,diu3_sz,diu3)
          call oclWrite3DFloatArrayBuffer(diu4_buf,diu4_sz,diu4)
          call oclWrite3DFloatArrayBuffer(diu5_buf,diu5_sz,diu5)
          call oclWrite3DFloatArrayBuffer(diu6_buf,diu6_sz,diu6)
          call oclWrite3DFloatArrayBuffer(diu7_buf,diu7_sz,diu7)
          call oclWrite3DFloatArrayBuffer(diu8_buf,diu8_sz,diu8)
          call oclWrite3DFloatArrayBuffer(diu9_buf,diu9_sz,diu9)
          call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
          call runOcl(oclGlobalRange,oclLocalRange,exectime)
          ! call les_map_91

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN les_map_109 ------------------------------------------------------------------------------------------------------
          oclGlobalRange = ((((80 + 1) - 0) + 1) * (((300 + 1) - (-1)) + 1))
          oclLocalRange = 0
          state_ptr(1) = ST_LES_MAP_109
          
          call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
          call runOcl(oclGlobalRange,oclLocalRange,exectime)
          ! call les_map_109

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN les_map_115 ------------------------------------------------------------------------------------------------------
          oclGlobalRange = ((((80 + 1) - 0) + 1) * (((300 + 1) - 0) + 1))
          oclLocalRange = 0
          state_ptr(1) = ST_LES_MAP_115
          
          call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
          call runOcl(oclGlobalRange,oclLocalRange,exectime)
          ! call les_map_115

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN les_map_121 ------------------------------------------------------------------------------------------------------
          oclGlobalRange = ((((300 + 1) - (-1)) + 1) * (((300 + 1) - 0) + 1))
          oclLocalRange = 0
          state_ptr(1) = ST_LES_MAP_121
          
          call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
          call runOcl(oclGlobalRange,oclLocalRange,exectime)
          ! call les_map_121

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN les_map_127 ------------------------------------------------------------------------------------------------------
          oclGlobalRange = (((80 - 2) + 1) * (300 * 300))
          oclLocalRange = 0
          state_ptr(1) = ST_LES_MAP_127
          
          call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
          call runOcl(oclGlobalRange,oclLocalRange,exectime)
          ! call les_map_127

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN les_map_155 ------------------------------------------------------------------------------------------------------
          oclGlobalRange = (300 * 300)
          oclLocalRange = 0
          state_ptr(1) = ST_LES_MAP_155
          
          call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
          call runOcl(oclGlobalRange,oclLocalRange,exectime)
          ! call les_map_155

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN les_map_175 ------------------------------------------------------------------------------------------------------
          oclGlobalRange = (((80 - 2) + 1) * (300 * 300))
          oclLocalRange = 0
          state_ptr(1) = ST_LES_MAP_175
          
          call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
          call runOcl(oclGlobalRange,oclLocalRange,exectime)
          ! call les_map_175

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN les_map_203 ------------------------------------------------------------------------------------------------------
          oclGlobalRange = (300 * 300)
          oclLocalRange = 0
          state_ptr(1) = ST_LES_MAP_203
          
          call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
          call runOcl(oclGlobalRange,oclLocalRange,exectime)
          ! call les_map_203

! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN les_map_223 ------------------------------------------------------------------------------------------------------
          oclGlobalRange = (80 * (300 * 300))
          oclLocalRange = 0
          state_ptr(1) = ST_LES_MAP_223
          
          call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
          call runOcl(oclGlobalRange,oclLocalRange,exectime)
          ! call les_map_223

! ---- END --------------------------------------------------------------------------------------------------------------------
      return
      end subroutine les
! Footer (produceCode_progUnit c)
end module module_les
