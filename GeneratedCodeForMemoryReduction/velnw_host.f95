module module_velnw
 contains 
      subroutine velnw(p0,ro,dxs,u,dt,f,dys,v,g,dzs,w,h)

      use module_adam_feedbf_les_press_velfg_ve_etc_superkernel_init
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
      real(4), dimension(0:ip+2,0:jp+2,0:kp+1), intent(in) :: p0 !!
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
      integer(8) :: p0_buf
      integer(8) :: ro_buf
      integer(8) :: dxs_buf
      integer(8) :: u_buf
      integer(8) :: dt_buf
      integer(8) :: f_buf
      integer(8) :: dys_buf
      integer(8) :: v_buf
      integer(8) :: g_buf
      integer(8) :: dzs_buf
      integer(8) :: w_buf
      integer(8) :: h_buf

      integer, dimension(1) :: state_ptr

      ! Size declarations
      integer, dimension(1) :: state_ptr_sz
      integer, dimension(3) :: p0_sz
      integer, dimension(1) :: ro_ptr_sz
      integer, dimension(1) :: dxs_sz
      integer, dimension(3) :: u_sz
      integer, dimension(1) :: dt_ptr_sz
      integer, dimension(3) :: f_sz
      integer, dimension(1) :: dys_sz
      integer, dimension(3) :: v_sz
      integer, dimension(3) :: g_sz
      integer, dimension(1) :: dzs_sz
      integer, dimension(3) :: w_sz
      integer, dimension(3) :: h_sz
      real(kind=4), dimension(1) :: ro_ptr
      real(kind=4), dimension(1) :: dt_ptr

      
! Size assignments
      state_ptr_sz = shape(state_ptr)
      p0_sz = shape(p0)
      ro_ptr_sz = shape(ro_ptr)
      dxs_sz = shape(dxs)
      u_sz = shape(u)
      dt_ptr_sz = shape(dt_ptr)
      f_sz = shape(f)
      dys_sz = shape(dys)
      v_sz = shape(v)
      g_sz = shape(g)
      dzs_sz = shape(dzs)
      w_sz = shape(w)
      h_sz = shape(h)

      ! Buffer loads
      call oclLoadBuffer(STATE_PTR_BUF_IDX, state_ptr_buf)
      call oclLoadBuffer(P0_BUF_IDX, p0_buf)
      call oclLoadBuffer(RO_BUF_IDX, ro_buf)
      call oclLoadBuffer(DXS_BUF_IDX, dxs_buf)
      call oclLoadBuffer(U_BUF_IDX, u_buf)
      call oclLoadBuffer(DT_BUF_IDX, dt_buf)
      call oclLoadBuffer(F_BUF_IDX, f_buf)
      call oclLoadBuffer(DYS_BUF_IDX, dys_buf)
      call oclLoadBuffer(V_BUF_IDX, v_buf)
      call oclLoadBuffer(G_BUF_IDX, g_buf)
      call oclLoadBuffer(DZS_BUF_IDX, dzs_buf)
      call oclLoadBuffer(W_BUF_IDX, w_buf)
      call oclLoadBuffer(H_BUF_IDX, h_buf)

      ! Original code with buffer writes and reads
! ---- BEGIN velnw_map_36 -----------------------------------------------------------------------------------------------------
      oclGlobalRange = (80 * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_VELNW_MAP_36
      
      call oclWrite3DFloatArrayBuffer(p0_buf,p0_sz,p0)
      ro_ptr(1) = ro
      call oclWrite1DFloatArrayBuffer(ro_buf,ro_ptr_sz,ro_ptr)! Automatic conversion to array
      call oclWrite1DFloatArrayBuffer(dxs_buf,dxs_sz,dxs)
      dt_ptr(1) = dt
      call oclWrite1DFloatArrayBuffer(dt_buf,dt_ptr_sz,dt_ptr)! Automatic conversion to array
      call oclWrite3DFloatArrayBuffer(f_buf,f_sz,f)
      call oclWrite3DFloatArrayBuffer(u_buf,u_sz,u)
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call velnw_map_36

      call oclRead3DFloatArrayBuffer(u_buf,u_sz,u)
! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN velnw_map_44 -----------------------------------------------------------------------------------------------------
      oclGlobalRange = (80 * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_VELNW_MAP_44
      
      call oclWrite3DFloatArrayBuffer(p0_buf,p0_sz,p0)
      ro_ptr(1) = ro
      call oclWrite1DFloatArrayBuffer(ro_buf,ro_ptr_sz,ro_ptr)! Automatic conversion to array
      call oclWrite1DFloatArrayBuffer(dys_buf,dys_sz,dys)
      dt_ptr(1) = dt
      call oclWrite1DFloatArrayBuffer(dt_buf,dt_ptr_sz,dt_ptr)! Automatic conversion to array
      call oclWrite3DFloatArrayBuffer(g_buf,g_sz,g)
      call oclWrite3DFloatArrayBuffer(v_buf,v_sz,v)
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call velnw_map_44

      call oclRead3DFloatArrayBuffer(v_buf,v_sz,v)
! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN velnw_map_52 -----------------------------------------------------------------------------------------------------
      oclGlobalRange = ((80 - 1) * (300 * 300))
      oclLocalRange = 0
      state_ptr(1) = ST_VELNW_MAP_52
      
      call oclWrite3DFloatArrayBuffer(p0_buf,p0_sz,p0)
      ro_ptr(1) = ro
      call oclWrite1DFloatArrayBuffer(ro_buf,ro_ptr_sz,ro_ptr)! Automatic conversion to array
      call oclWrite1DFloatArrayBuffer(dzs_buf,dzs_sz,dzs)
      dt_ptr(1) = dt
      call oclWrite1DFloatArrayBuffer(dt_buf,dt_ptr_sz,dt_ptr)! Automatic conversion to array
      call oclWrite3DFloatArrayBuffer(h_buf,h_sz,h)
      call oclWrite3DFloatArrayBuffer(w_buf,w_sz,w)
      call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
      call runOcl(oclGlobalRange,oclLocalRange,exectime)
      ! call velnw_map_52

      call oclRead3DFloatArrayBuffer(w_buf,w_sz,w)
! ---- END --------------------------------------------------------------------------------------------------------------------
      return
      end subroutine velnw
! Footer (produceCode_progUnit c)
end module module_velnw
