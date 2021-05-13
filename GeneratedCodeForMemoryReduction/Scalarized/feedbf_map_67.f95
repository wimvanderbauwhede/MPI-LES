module singleton_module_feedbf_map_67

contains

subroutine feedbf_map_67_scal(f,fx,g,fy,h,fz)
    ! local vars: i,j,k
    ! parallelfortran: synthesised loop variable decls
! READ
    real, dimension(0:ip,0:jp,0:kp), intent(in) :: fx
    real, dimension(0:ip,0:jp,0:kp), intent(in) :: fy
    real, dimension(0:ip,0:jp,0:kp), intent(in) :: fz
! WRITTEN
! READ & WRITTEN
    real, dimension(0:ip,0:jp,0:kp), intent(inout) :: f
    real, dimension(0:ip,0:jp,0:kp), intent(inout) :: g
    real, dimension(0:ip,0:jp,0:kp), intent(inout) :: h
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
                f(i,j,k) = f(i,j,k) + fx(i,j,k)
                g(i,j,k) = g(i,j,k) + fy(i,j,k)
                h(i,j,k) = h(i,j,k) + fz(i,j,k)
end subroutine feedbf_map_67

end module singleton_module_feedbf_map_67

