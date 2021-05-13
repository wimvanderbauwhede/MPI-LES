module module_feedbf
 contains 
subroutine feedbf(usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,&
                  dt,beta,fx,fy,fz,f,g,h,n)

        use module_adam_feedbf_les_press_velfg_ve_etc_superkernel_init
        use oclWrapper
!      implicit none
      use params_common_sn

      implicit none

! otherStatements

! remainingDecls
      real(4), intent(in) :: alpha !!
      real(4), intent(in) :: beta !!
      real(4), dimension(-1:ip+1,0:jp+1,0:kp+1), intent(in) :: bmask1 !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(in) :: cmask1 !!
      real(4), dimension(0:ip+1,0:jp+1,0:kp+1), intent(in) :: dmask1 !!
      real(4), intent(in) :: dt !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: f !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(out) :: fx !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(out) :: fy !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(out) :: fz !!
      integer, intent(in) :: n !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: g !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: h !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(in) :: u !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: usum !!
      real(4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(in) :: v !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: vsum !!
      real(4), dimension(0:ip+1,-1:jp+1,-1:kp+1), intent(in) :: w !!
      real(4), dimension(0:ip,0:jp,0:kp), intent(inout) :: wsum !!
      integer :: i !!
      integer :: j !!
      integer :: k !!
      real(4) :: f1x !!
      real(4) :: f1y !!
      real(4) :: f1z !!
      real(4) :: f2x !!
      real(4) :: f2y !!
      real(4) :: f2z !!

        ! Extra declarations
        real (kind=4) :: exectime

        ! Buffer declarations
        integer(8) :: state_ptr_buf
        integer(8) :: usum_buf
        integer(8) :: u_buf
        integer(8) :: bmask1_buf
        integer(8) :: vsum_buf
        integer(8) :: v_buf
        integer(8) :: cmask1_buf
        integer(8) :: wsum_buf
        integer(8) :: w_buf
        integer(8) :: dmask1_buf
        integer(8) :: alpha_buf
        integer(8) :: dt_buf
        integer(8) :: beta_buf
        integer(8) :: fx_buf
        integer(8) :: fy_buf
        integer(8) :: fz_buf
        integer(8) :: f_buf
        integer(8) :: g_buf
        integer(8) :: h_buf

        integer, dimension(1) :: state_ptr

        ! Size declarations
        integer, dimension(1) :: state_ptr_sz
        integer, dimension(3) :: usum_sz
        integer, dimension(3) :: u_sz
        integer, dimension(3) :: bmask1_sz
        integer, dimension(3) :: vsum_sz
        integer, dimension(3) :: v_sz
        integer, dimension(3) :: cmask1_sz
        integer, dimension(3) :: wsum_sz
        integer, dimension(3) :: w_sz
        integer, dimension(3) :: dmask1_sz
        integer, dimension(1) :: alpha_ptr_sz
        integer, dimension(1) :: dt_ptr_sz
        integer, dimension(1) :: beta_ptr_sz
        integer, dimension(3) :: fx_sz
        integer, dimension(3) :: fy_sz
        integer, dimension(3) :: fz_sz
        integer, dimension(3) :: f_sz
        integer, dimension(3) :: g_sz
        integer, dimension(3) :: h_sz
        real(kind=4), dimension(1) :: alpha_ptr
        real(kind=4), dimension(1) :: dt_ptr
        real(kind=4), dimension(1) :: beta_ptr

        
! Size assignments
        state_ptr_sz = shape(state_ptr)
        usum_sz = shape(usum)
        u_sz = shape(u)
        bmask1_sz = shape(bmask1)
        vsum_sz = shape(vsum)
        v_sz = shape(v)
        cmask1_sz = shape(cmask1)
        wsum_sz = shape(wsum)
        w_sz = shape(w)
        dmask1_sz = shape(dmask1)
        alpha_ptr_sz = shape(alpha_ptr)
        dt_ptr_sz = shape(dt_ptr)
        beta_ptr_sz = shape(beta_ptr)
        fx_sz = shape(fx)
        fy_sz = shape(fy)
        fz_sz = shape(fz)
        f_sz = shape(f)
        g_sz = shape(g)
        h_sz = shape(h)

        ! Buffer loads
        call oclLoadBuffer(STATE_PTR_BUF_IDX, state_ptr_buf)
        call oclLoadBuffer(USUM_BUF_IDX, usum_buf)
        call oclLoadBuffer(U_BUF_IDX, u_buf)
        call oclLoadBuffer(BMASK1_BUF_IDX, bmask1_buf)
        call oclLoadBuffer(VSUM_BUF_IDX, vsum_buf)
        call oclLoadBuffer(V_BUF_IDX, v_buf)
        call oclLoadBuffer(CMASK1_BUF_IDX, cmask1_buf)
        call oclLoadBuffer(WSUM_BUF_IDX, wsum_buf)
        call oclLoadBuffer(W_BUF_IDX, w_buf)
        call oclLoadBuffer(DMASK1_BUF_IDX, dmask1_buf)
        call oclLoadBuffer(ALPHA_BUF_IDX, alpha_buf)
        call oclLoadBuffer(DT_BUF_IDX, dt_buf)
        call oclLoadBuffer(BETA_BUF_IDX, beta_buf)
        call oclLoadBuffer(FX_BUF_IDX, fx_buf)
        call oclLoadBuffer(FY_BUF_IDX, fy_buf)
        call oclLoadBuffer(FZ_BUF_IDX, fz_buf)
        call oclLoadBuffer(F_BUF_IDX, f_buf)
        call oclLoadBuffer(G_BUF_IDX, g_buf)
        call oclLoadBuffer(H_BUF_IDX, h_buf)

        ! Original code with buffer writes and reads
! ---- BEGIN feedbf_map_49 ----------------------------------------------------------------------------------------------------
        oclGlobalRange = (80 * (300 * 300))
        oclLocalRange = 0
        state_ptr(1) = ST_FEEDBF_MAP_49
        
        call oclWrite3DFloatArrayBuffer(u_buf,u_sz,u)
        call oclWrite3DFloatArrayBuffer(bmask1_buf,bmask1_sz,bmask1)
        call oclWrite3DFloatArrayBuffer(v_buf,v_sz,v)
        call oclWrite3DFloatArrayBuffer(cmask1_buf,cmask1_sz,cmask1)
        call oclWrite3DFloatArrayBuffer(w_buf,w_sz,w)
        call oclWrite3DFloatArrayBuffer(dmask1_buf,dmask1_sz,dmask1)
        alpha_ptr(1) = alpha
        call oclWrite1DFloatArrayBuffer(alpha_buf,alpha_ptr_sz,alpha_ptr)! Automatic conversion to array
        dt_ptr(1) = dt
        call oclWrite1DFloatArrayBuffer(dt_buf,dt_ptr_sz,dt_ptr)! Automatic conversion to array
        beta_ptr(1) = beta
        call oclWrite1DFloatArrayBuffer(beta_buf,beta_ptr_sz,beta_ptr)! Automatic conversion to array
        call oclWrite3DFloatArrayBuffer(usum_buf,usum_sz,usum)
        call oclWrite3DFloatArrayBuffer(vsum_buf,vsum_sz,vsum)
        call oclWrite3DFloatArrayBuffer(wsum_buf,wsum_sz,wsum)
        call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
        call runOcl(oclGlobalRange,oclLocalRange,exectime)
        ! call feedbf_map_49

        call oclRead3DFloatArrayBuffer(fx_buf,fx_sz,fx)
        call oclRead3DFloatArrayBuffer(fy_buf,fy_sz,fy)
        call oclRead3DFloatArrayBuffer(fz_buf,fz_sz,fz)
        call oclRead3DFloatArrayBuffer(usum_buf,usum_sz,usum)
        call oclRead3DFloatArrayBuffer(vsum_buf,vsum_sz,vsum)
        call oclRead3DFloatArrayBuffer(wsum_buf,wsum_sz,wsum)
! ---- END --------------------------------------------------------------------------------------------------------------------
! ---- BEGIN feedbf_map_67 ----------------------------------------------------------------------------------------------------
        oclGlobalRange = (80 * (300 * 300))
        oclLocalRange = 0
        state_ptr(1) = ST_FEEDBF_MAP_67
        
        call oclWrite3DFloatArrayBuffer(fx_buf,fx_sz,fx)
        call oclWrite3DFloatArrayBuffer(fy_buf,fy_sz,fy)
        call oclWrite3DFloatArrayBuffer(fz_buf,fz_sz,fz)
        call oclWrite3DFloatArrayBuffer(f_buf,f_sz,f)
        call oclWrite3DFloatArrayBuffer(g_buf,g_sz,g)
        call oclWrite3DFloatArrayBuffer(h_buf,h_sz,h)
        call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
        call runOcl(oclGlobalRange,oclLocalRange,exectime)
        ! call feedbf_map_67

        call oclRead3DFloatArrayBuffer(f_buf,f_sz,f)
        call oclRead3DFloatArrayBuffer(g_buf,g_sz,g)
        call oclRead3DFloatArrayBuffer(h_buf,h_sz,h)
! ---- END --------------------------------------------------------------------------------------------------------------------
end subroutine feedbf
! Footer (produceCode_progUnit c)
end module module_feedbf
