module singleton_module_adam_map_36

contains

subroutine adam_map_36_scal(f,g,h,fold,gold,hold)
    ! local vars: fd,i,j,k,gd,hd
    real(kind=4) :: fd
    real(kind=4) :: gd
    real(kind=4) :: hd
    ! parallelfortran: synthesised loop variable decls
! READ
! WRITTEN
! READ & WRITTEN
    real, dimension(0:ip,0:jp,0:kp), intent(inout) :: f
    real, dimension(0:ip,0:jp,0:kp), intent(inout) :: g
    real, dimension(0:ip,0:jp,0:kp), intent(inout) :: h
    real, dimension(1:ip,1:jp,1:kp), intent(inout) :: fold
    real, dimension(1:ip,1:jp,1:kp), intent(inout) :: gold
    real, dimension(1:ip,1:jp,1:kp), intent(inout) :: hold
! globalIdDeclaration
! globalIdInitialisation
! ptrAssignments_fseq
    ! parallelfortran: synthesised loop variables
    ! parallelfortran: original code
                fd = f(i,j,k)
                gd = g(i,j,k)
                hd = h(i,j,k)
                f(i,j,k) = 1.5 * f(i,j,k) - 0.5 * fold(i,j,k)
                g(i,j,k) = 1.5 * g(i,j,k) - 0.5 * gold(i,j,k)
                h(i,j,k) = 1.5 * h(i,j,k) - 0.5 * hold(i,j,k)
                fold(i,j,k) = fd
                gold(i,j,k) = gd
                hold(i,j,k) = hd
end subroutine adam_map_36

end module singleton_module_adam_map_36

