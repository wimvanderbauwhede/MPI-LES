module singleton_module_feedbf_map_49

contains

subroutine feedbf_map_49_scal(usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,fx,fy,fz)
    ! local vars: i,j,k,f1x,f1y,f1z,f2x,f2y,f2z
    real(kind=4) :: f1x
    real(kind=4) :: f1y
    real(kind=4) :: f1z
    real(kind=4) :: f2x
    real(kind=4) :: f2y
    real(kind=4) :: f2z
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: u
    real, dimension((-1):(ip+1),0:(jp+1),0:(kp+1)), intent(in) :: bmask1
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: v
    real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: cmask1
    real, dimension(0:(ip+1),(-1):(jp+1),(-1):(kp+1)), intent(in) :: w
    real, dimension(0:(ip+1),0:(jp+1),0:(kp+1)), intent(in) :: dmask1
    real(kind=4), intent(in) :: alpha
    real(kind=4), intent(in) :: dt
    real(kind=4), intent(in) :: beta
! WRITTEN
    real, dimension(0:ip,0:jp,0:kp), intent(out) :: fx
    real, dimension(0:ip,0:jp,0:kp), intent(out) :: fy
    real, dimension(0:ip,0:jp,0:kp), intent(out) :: fz
! READ & WRITTEN
    real, dimension(0:ip,0:jp,0:kp), intent(inout) :: usum
    real, dimension(0:ip,0:jp,0:kp), intent(inout) :: vsum
    real, dimension(0:ip,0:jp,0:kp), intent(inout) :: wsum
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
                usum(i,j,k) = (usum(i,j,k) + u(i,j,k)) * bmask1(i,j,k)
                vsum(i,j,k) = (vsum(i,j,k) + v(i,j,k)) * cmask1(i,j,k)
                wsum(i,j,k) = (wsum(i,j,k) + w(i,j,k)) * dmask1(i,j,k)
                f1x = alpha * usum(i,j,k) * dt
                f1y = alpha * vsum(i,j,k) * dt
                f1z = alpha * wsum(i,j,k) * dt
                f2x = beta * u(i,j,k) * bmask1(i,j,k)
                f2y = beta * v(i,j,k) * cmask1(i,j,k)
                f2z = beta * w(i,j,k) * dmask1(i,j,k)
                fx(i,j,k) = f1x + f2x
                fy(i,j,k) = f1y + f2y
                fz(i,j,k) = f1z + f2z
end subroutine feedbf_map_49

end module singleton_module_feedbf_map_49

