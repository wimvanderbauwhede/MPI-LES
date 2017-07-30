module module_boundp
#ifdef MPI
    use communication_helper_real
#else
    use params_common_sn
#endif
implicit none

contains
#ifdef NESTED_LES
subroutine boundp2(n,p)
    use common_sn ! create_new_include_statements() line 102
    integer, intent(In) :: n
#else
subroutine boundp2(p)
    use common_sn ! create_new_include_statements() line 102
#endif

    real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1) , intent(InOut) :: p
    integer :: i, j
!
! --computational boundary(neumann condition)
    do j = 0,jp+1
        do i = 0,ip+1
            p(i,j,   0) = p(i,j,1)
            p(i,j,kp+1) = p(i,j,kp)
        end do
    end do
#ifdef MPI
#ifdef NESTED_LES
!   if (syncTicks == 0  .and. n > 2) then
   if (syncTicks == 0) then
#endif
! --halo exchanges
    call exchangeRealHalos(p, procPerRow, neighbours, 1, 2, 1, 2)
#ifdef NESTED_LES
   end if
#endif
#endif
end subroutine boundp2

#ifdef NESTED_LES
subroutine boundp1(n,p)
    integer, intent(In) :: n
#else
subroutine boundp1(p)
    use common_sn ! create_new_include_statements() line 102
#endif

    real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1) , intent(InOut) :: p
#if !defined(MPI) || (PROC_PER_ROW==1)
    integer :: i, j, k
#else
    integer :: j, k
#endif
!
! --computational boundary(neumann condition)
#ifdef MPI
    if (isTopRow(procPerRow) .or. isBottomRow(procPerRow)) then
#endif
        do k = 0,kp+1
            do j = 0,jp+1
#ifdef MPI
                if (isTopRow(procPerRow)) then
#endif
                    p(   0,j,k) = p(1 ,j,k)
#ifdef MPI
                else
#endif
                    p(ip+1,j,k) = p(ip,j,k)
#ifdef MPI
                end if
#endif
            end do
        end do
#ifdef MPI
    end if
#endif
! --side flow exchanges
#if !defined(MPI) || (PROC_PER_ROW==1)
    do k = 0,kp+1
        do i = 0,ip+1
            p(i,   0,k) = p(i,jp,k) ! right to left
            p(i,jp+1,k) = p(i, 1,k) ! left to right
        end do
    end do
#else
    call sideflowRightLeft(p, procPerRow, jp+1, 1, 0, 1, 0, 0)
    call sideflowLeftRight(p, procPerRow, 2, jp+2, 0, 1, 0, 0)
#endif
#ifdef MPI
! --halo exchanges
#ifdef NESTED_LES
!   if (syncTicks == 0  .and. n > 2) then
   if (syncTicks == 0) then
#endif
    call exchangeRealHalos(p, procPerRow, neighbours, 1, 2, 1, 2)
#ifdef NESTED_LES
   end if
#endif
#endif
end subroutine boundp1

end module module_boundp
