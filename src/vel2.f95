module module_vel2
#ifdef MPI
    use communication_helper_real
#endif
    implicit none
contains

#ifdef NESTED_LES
      subroutine vel2(n,nou1,u,diu1,dx1,nou5,v,diu5,dy1,nou9,w,diu9,dzn,cov1,cov5,cov9,nou2, &
      diu2,cov2,nou3,diu3,dzs,cov3,nou4,diu4,cov4,nou6,diu6,cov6,nou7,diu7,cov7,nou8,diu8,cov8,uspd,vspd)
      use common_sn ! create_new_include_statements() line 102

        integer, intent(In) :: n
#else
      subroutine vel2(nou1,u,diu1,dx1,nou5,v,diu5,dy1,nou9,w,diu9,dzn,cov1,cov5,cov9,nou2, &
      diu2,cov2,nou3,diu3,dzs,cov3,nou4,diu4,cov4,nou6,diu6,cov6,nou7,diu7,cov7,nou8,diu8,cov8,uspd,vspd)
      use common_sn ! create_new_include_statements() line 102
#endif
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov9
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu9
        real(kind=4), dimension(-1:ip+1) , intent(In) :: dx1
        real(kind=4), dimension(0:jp+1) , intent(In) :: dy1
        real(kind=4), dimension(-1:kp+2) , intent(In) :: dzn
        real(kind=4), dimension(-1:kp+2) , intent(In) :: dzs
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou9
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: u
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(In) :: w
!
!wall function
        real(kind=4), dimension(0:ip+1,0:jp+1) , intent(out) :: uspd
        real(kind=4), dimension(0:ip+1,0:jp+1) , intent(out) :: vspd
        integer :: i,j,k
        integer, parameter  :: u0 = 0

        
      do j=1,jp
        do i=1,ip
         uspd(i,j)=(u(i,j,1)**2+((0.5*(v(i,j-1,1)+v(i,j,1))*dx1(i+1)&
     +0.5*(v(i+1,j-1,1)+v(i+1,j,1))*dx1(i))/(dx1(i)+dx1(i+1)))**2)**0.5
        end do
        end do
        do j=1,jp
        do i=1,ip
         vspd(i,j)=(v(i,j,1)**2+((0.5*(u(i-1,j,1)+u(i,j,1))*dy1(j+1)&
     +0.5*(u(i-1,j+1,1)+u(i,j+1,1))*dy1(j))/(dy1(j)+dy1(j+1)))**2)**0.5
        end do
        end do
! WV: to get a point somewhere near the middle of the domain
#ifdef MPI
       if (rank == mpi_size / 2 + procPerRow / 2 - 1 ) then
#endif
        write(6,*) 'CHK_uspd=',uspd(ip/2,jp/2),vspd(ip/2,jp/2)
#ifdef MPI
       end if
#endif
!
      do k = 1,kp
      do j = 1,jp
      do i = 1,ip
        nou1(i,j,k) = ( u(i-1,j,k)+u(i,j,k))/2. ! i.e 2-point average of u in i-direction
        diu1(i,j,k) = (-u(i-1,j,k)+u(i,j,k))/dx1(i) ! i.e. du/dx
        nou5(i,j,k) = ( v(i,j-1,k)+v(i,j,k))/2.
        diu5(i,j,k) = (-v(i,j-1,k)+v(i,j,k))/dy1(j)
        nou9(i,j,k) = ( w(i,j,k-1)+w(i,j,k))/2.
        diu9(i,j,k) = (-w(i,j,k-1)+w(i,j,k))/dzn(k)
!
        cov1(i,j,k) = nou1(i,j,k)*diu1(i,j,k) !
        cov5(i,j,k) = nou5(i,j,k)*diu5(i,j,k)
        cov9(i,j,k) = nou9(i,j,k)*diu9(i,j,k)
      end do
      end do
      end do

#ifdef GR_DEBUG
    print*, 'GR: SUM(nou1) = ', sum(nou1)
    print*, 'GR: SUM(diu1) = ', sum(diu1)
    print*, 'GR: SUM(nou5) = ', sum(nou5)
    print*, 'GR: SUM(diu5) = ', sum(diu5)
    print*, 'GR: SUM(nou9) = ', sum(nou9)
    print*, 'GR: SUM(diu9) = ', sum(diu9)
    print*, 'GR: SUM(cov1) = ', sum(cov1)
    print*, 'GR: SUM(cov5) = ', sum(cov5)
    print*, 'GR: SUM(cov9) = ', sum(cov9)
#endif

      do k = 1,kp
      do j = 1,jp
      do i = 1,ip
        nou2(i,j,k) = (dx1(i+1)*v(i,j-1,k)+dx1(i)*v(i+1,j-1,k)) /(dx1(i)+dx1(i+1))
        diu2(i,j,k) = 2.*(-u(i,j-1,k)+u(i,j,k))/(dy1(j-1)+dy1(j))
        cov2(i,j,k) = nou2(i,j,k)*diu2(i,j,k)
      end do
      end do
      end do

#ifdef GR_DEBUG
    print*, 'GR: SUM(nou2) = ', sum(nou2)
    print*, 'GR: SUM(diu2) = ', sum(diu2)
    print*, 'GR: SUM(cov2) = ', sum(cov2)
#endif
!
      do k = 2,kp+1
      do j = 1,jp
      do i = 1,ip
        nou3(i,j,k) = (dx1(i+1)*w(i,j,k-1)+dx1(i)*w(i+1,j,k-1)) /(dx1(i)+dx1(i+1))
        diu3(i,j,k) = (-u(i,j,k-1)+u(i,j,k))/dzs(k-1)
        cov3(i,j,k) = nou3(i,j,k)*diu3(i,j,k)
      end do
      end do
      end do

      do j=1,jp
      do i=1,ip
       nou3(i,j,1) = 0.5*(dx1(i+1)*w(i,j,1)+dx1(i)*w(i+1,j,1))/(dx1(i)+dx1(i+1))
       diu3(i,j,1)=uspd(i,j)*0.4/alog(0.5*dzn(1)/0.1)/(0.5*dzn(1))/0.4*u(i,j,1)/uspd(i,j)
       cov3(i,j,1) = nou3(i,j,1)*diu3(i,j,1)
      end do
      end do


      do k = 1,kp
      do j = 1,jp
      do i = 1,ip
        nou4(i,j,k) = (dy1(j+1)*u(i-1,j,k)+dy1(j)*u(i-1,j+1,k)) /(dy1(j)+dy1(j+1))
        diu4(i,j,k) = 2.*(-v(i-1,j,k)+v(i,j,k))/(dx1(i-1)+dx1(i))
        cov4(i,j,k) = (nou4(i,j,k)-u0)*diu4(i,j,k)
      end do
      end do
      end do
!
      do k = 2,kp+1
      do j = 1,jp
      do i = 1,ip
        nou6(i,j,k) = (dy1(j+1)*w(i,j,k-1)+dy1(j)*w(i,j+1,k-1)) /(dy1(j)+dy1(j+1))
        diu6(i,j,k) = (-v(i,j,k-1)+v(i,j,k))/dzs(k-1)
        cov6(i,j,k) = nou6(i,j,k)*diu6(i,j,k)
      end do
      end do
      end do
!
      do j=1,jp
      do i=1,ip
       nou6(i,j,1) = 0.5*(dy1(j+1)*w(i,j,1)+dy1(j)*w(i,j+1,1))/(dy1(j)+dy1(j+1))
       diu6(i,j,1)=vspd(i,j)*0.4/alog(0.5*dzn(1)/0.1)/(0.5*dzn(1))/0.4*v(i,j,1)/vspd(i,j)
       cov6(i,j,1) = nou6(i,j,1)*diu6(i,j,1)
      end do
      end do


      do k = 1,kp-1
      do j = 1,jp
      do i = 1,ip
        nou7(i,j,k) = (dzn(k+1)*u(i-1,j,k)+dzn(k)*u(i-1,j,k+1)) /(dzn(k)+dzn(k+1))
        diu7(i,j,k) = 2.*(-w(i-1,j,k)+w(i,j,k))/(dx1(i-1)+dx1(i))
        cov7(i,j,k) = (nou7(i,j,k)-u0)*diu7(i,j,k)
      end do
      end do
      end do
!
      do k = 1,kp-1
      do j = 1,jp
      do i = 1,ip
        nou8(i,j,k) = (dzn(k+1)*v(i,j-1,k)+dzn(k)*v(i,j-1,k+1)) /(dzn(k)+dzn(k+1))
        diu8(i,j,k) = 2.*(-w(i,j-1,k)+w(i,j,k))/(dy1(j-1)+dy1(j))
        cov8(i,j,k) = nou8(i,j,k)*diu8(i,j,k)
      end do
      end do
      end do
! ====================================
#ifdef MPI
    if (isBottomRow(procPerRow)) then
#endif
      do k = 1,kp
      do j = 1,jp
        nou1(ip+1,j,k) = nou1(ip,j,k)
        diu1(ip+1,j,k) = diu1(ip,j,k)
        cov1(ip+1,j,k) = cov1(ip,j,k)
      end do
      end do
#ifdef MPI
    end if
#endif
#if !defined(MPI) || (PROC_PER_ROW==1)
      do k = 1,kp
      do i = 1,ip
        nou2(i,0,k) = nou2(i,jp,k)
        diu2(i,0,k) = diu2(i,jp,k)
        cov2(i,0,k) = cov2(i,jp,k)
        nou2(i,jp+1,k) = nou2(i,1,k)
        diu2(i,jp+1,k) = diu2(i,1,k)
        cov2(i,jp+1,k) = cov2(i,1,k)
      end do
      end do
#else
    call sideflowRightLeft(nou2, procPerRow, jp+1, 1, 1, 2, 1, 2)
    call sideflowRightLeft(diu2, procPerRow, jp+1, 1, 1, 2, 1, 2)
    call sideflowRightLeft(cov2, procPerRow, jp+1, 1, 1, 2, 1, 2)
    call sideflowLeftRight(nou2, procPerRow, 2, jp+2, 1, 2, 1, 2)
    call sideflowLeftRight(diu2, procPerRow, 2, jp+2, 1, 2, 1, 2)
    call sideflowLeftRight(cov2, procPerRow, 2, jp+2, 1, 2, 1, 2)
#endif
#ifdef MPI
    if (isBottomRow(procPerRow)) then
#endif
      do k = 1,kp
      do j = 1,jp
        nou4(ip+1,j,k) = nou4(ip,j,k)
        diu4(ip+1,j,k) = diu4(ip,j,k)
        cov4(ip+1,j,k) = cov4(ip,j,k)
      end do
      end do
#ifdef MPI
    end if
#endif
#if !defined(MPI) || (PROC_PER_ROW==1)
      do k = 1,kp
      do i = 1,ip
        nou5(i,0,k) = nou5(i,jp,k)
        diu5(i,0,k) = diu5(i,jp,k)
        cov5(i,0,k) = cov5(i,jp,k)
        nou5(i,jp+1,k) = nou5(i,1,k)
        diu5(i,jp+1,k) = diu5(i,1,k)
        cov5(i,jp+1,k) = cov5(i,1,k)
      end do
      end do
#else
    call sideflowRightLeft(nou5, procPerRow, jp+2, 2, 2, 2, 1, 2)
    call sideflowRightLeft(diu5, procPerRow, jp+2, 2, 2, 2, 1, 2)
    call sideflowRightLeft(cov5, procPerRow, jp+2, 2, 2, 2, 1, 2)
    call sideflowLeftRight(nou5, procPerRow, 3, jp+3, 2, 2, 1, 2)
    call sideflowLeftRight(diu5, procPerRow, 3, jp+3, 2, 2, 1, 2)
    call sideflowLeftRight(cov5, procPerRow, 3, jp+3, 2, 2, 1, 2)
#endif
#ifdef MPI
    if (isBottomRow(procPerRow)) then
#endif
      do k = 1,kp-1
      do j = 1,jp
        nou7(ip+1,j,k) = nou7(ip,j,k)
        diu7(ip+1,j,k) = diu7(ip,j,k)
        cov7(ip+1,j,k) = cov7(ip,j,k)
      end do
      end do
#ifdef MPI
    end if
#endif
#if !defined(MPI) || (PROC_PER_ROW==1)
      do k = 1,kp-1
      do i = 1,ip
        nou8(i,0,k) = nou8(i,jp,k)
        diu8(i,0,k) = diu8(i,jp,k)
        cov8(i,0,k) = cov8(i,jp,k)
        nou8(i,jp+1,k) = nou8(i,1,k)
        diu8(i,jp+1,k) = diu8(i,1,k)
        cov8(i,jp+1,k) = cov8(i,1,k)
      end do
      end do
#else
    call sideflowRightLeft(nou8, procPerRow, jp+1, 1, 1, 2, 1, 3)
    call sideflowRightLeft(diu8, procPerRow, jp+1, 1, 1, 2, 1, 3)
    call sideflowRightLeft(cov8, procPerRow, jp+1, 1, 1, 2, 1, 3)
    call sideflowLeftRight(nou8, procPerRow, 2, jp+2, 1, 2, 1, 3)
    call sideflowLeftRight(diu8, procPerRow, 2, jp+2, 1, 2, 1, 3)
    call sideflowLeftRight(cov8, procPerRow, 2, jp+2, 1, 2, 1, 3)
#endif
! --les
#ifdef MPI
    if (isBottomRow(procPerRow)) then
#endif
      do k = 1,kp+1
      do j = 1,jp+1
        diu2(ip+1,j,k) = diu2(ip,j,k)
        diu3(ip+1,j,k) = diu3(ip,j,k)
      end do
      end do
#ifdef MPI
    end if
    if (isTopRow(procPerRow)) then
#endif
      do k = 1,kp+1
      do j = 1,jp+1
        diu2(0,j,k) = diu2(1,j,k)
        diu3(0,j,k) = diu3(1,j,k)
      end do
      end do
#ifdef MPI
    end if
#endif
#if !defined(MPI) || (PROC_PER_ROW==1)
      do k = 1,kp+1
      do i = 1,ip+1
        diu4(i,0,k) = diu4(i,jp,k)
        diu6(i,0,k) = diu6(i,jp,k)
      end do
      end do
#else
    call sideflowRightLeft(diu4, procPerRow, jp+1, 1, 1, 1, 1, 1)
    call sideflowRightLeft(diu6, procPerRow, jp+1, 1, 1, 1, 1, 1)
#endif

#ifdef MPI
#ifdef NESTED_LES
!   if (syncTicks == 0  .and. n > 2) then
   if (syncTicks == 0) then
#endif
    call exchangeRealHalos(nou1, procPerRow, neighbours, 1, 2, 2, 2)
    call exchangeRealHalos(diu1, procPerRow, neighbours, 1, 2, 2, 2)
    call exchangeRealHalos(cov1, procPerRow, neighbours, 1, 2, 2, 2)
    call exchangeRealHalos(nou2, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(diu2, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov2, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(nou3, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(diu3, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov3, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(nou4, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(diu4, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov4, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(nou5, procPerRow, neighbours, 1, 2, 2, 2)
    call exchangeRealHalos(diu5, procPerRow, neighbours, 1, 2, 2, 2)
    call exchangeRealHalos(cov5, procPerRow, neighbours, 1, 2, 2, 2)
    call exchangeRealHalos(nou6, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(diu6, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov6, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(nou7, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(diu7, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov7, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(nou8, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(diu8, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov8, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(nou9, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(diu9, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov9, procPerRow, neighbours, 1, 2, 1, 2)
#ifdef NESTED_LES
   end if
#endif
#endif

!uspd,vspd
!        do j=1,jp
!        do i=1,ip
!         uspd(i,j)=u(i,j,1)/uspd(i,j)
!         vspd(i,j)=v(i,j,1)/vspd(i,j)
!        end do
!        end do



#ifdef WV_DEBUG
    print *, 'F95 DIU SUMS:',sum(diu1),sum(diu2),sum(diu3),sum(diu4),sum(diu5),sum(diu6),sum(diu7),sum(diu8),sum(diu9)
#endif

end subroutine vel2

#if 0
#ifdef NESTED_LES
      subroutine vel2_scalarized(n,u,diu1,dx1,v,diu5,dy1,w,diu9,dzn,cov1,cov5,cov9, &
      diu2,cov2,diu3,dzs,cov3,diu4,cov4,diu6,cov6,diu7,cov7,diu8,cov8,uspd,vspd)
      use common_sn ! create_new_include_statements() line 102

        integer, intent(In) :: n
#else
      subroutine vel2_scalarized(u,diu1,dx1,v,diu5,dy1,w,diu9,dzn,cov1,cov5,cov9, &
      diu2,cov2,diu3,dzs,cov3,diu4,cov4,diu6,cov6,diu7,cov7,diu8,cov8,uspd,vspd)
      use common_sn ! create_new_include_statements() line 102
#endif
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov9
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu9
        real(kind=4), dimension(-1:ip+1) , intent(In) :: dx1
        real(kind=4), dimension(0:jp+1) , intent(In) :: dy1
        real(kind=4), dimension(-1:kp+2) , intent(In) :: dzn
        real(kind=4), dimension(-1:kp+2) , intent(In) :: dzs
        real(kind=4),  :: nou1
        real(kind=4),  :: nou2
        real(kind=4),  :: nou3
        real(kind=4),  :: nou4
        real(kind=4),  :: nou5
        real(kind=4),  :: nou6
        real(kind=4),  :: nou7
        real(kind=4),  :: nou8
        real(kind=4),  :: nou9
        integer :: kp1
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: u
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(In) :: w
!
!wall function
        real(kind=4), dimension(0:ip+1,0:jp+1) , intent(out) :: uspd
        real(kind=4), dimension(0:ip+1,0:jp+1) , intent(out) :: vspd
        integer :: i,j,k
        integer, parameter  :: u0 = 0


      do j=1,jp
        do i=1,ip
         uspd(i,j)=(u(i,j,1)**2+((0.5*(v(i,j-1,1)+v(i,j,1))*dx1(i+1)&
     +0.5*(v(i+1,j-1,1)+v(i+1,j,1))*dx1(i))/(dx1(i)+dx1(i+1)))**2)**0.5
        end do
        end do
        do j=1,jp
        do i=1,ip
         vspd(i,j)=(v(i,j,1)**2+((0.5*(u(i-1,j,1)+u(i,j,1))*dy1(j+1)&
     +0.5*(u(i-1,j+1,1)+u(i,j+1,1))*dy1(j))/(dy1(j)+dy1(j+1)))**2)**0.5
        end do
        end do
! WV: to get a point somewhere near the middle of the domain
#ifdef MPI
       if (rank == mpi_size / 2 + procPerRow / 2 - 1 ) then
#endif
        write(6,*) 'CHK_uspd=',uspd(ip/2,jp/2),vspd(ip/2,jp/2)
#ifdef MPI
       end if
#endif
!
      do k = 1,kp
      do j = 1,jp
      do i = 1,ip
        nou1 = (u(i-1,j,k)+u(i,j,k))/2.
        diu1(i,j,k) = (-u(i-1,j,k)+u(i,j,k))/dx1(i)
        nou5 = (v(i,j-1,k)+v(i,j,k))/2.
        diu5(i,j,k) = (-v(i,j-1,k)+v(i,j,k))/dy1(j)
        nou9 = (w(i,j,k-1)+w(i,j,k))/2.
        diu9(i,j,k) = (-w(i,j,k-1)+w(i,j,k))/dzn(k)
!
        cov1(i,j,k) = nou1*diu1(i,j,k)
        cov5(i,j,k) = nou5*diu5(i,j,k)
        cov9(i,j,k) = nou9*diu9(i,j,k)

        nou2 = (dx1(i+1)*v(i,j-1,k)+dx1(i)*v(i+1,j-1,k)) /(dx1(i)+dx1(i+1))
        diu2(i,j,k) = 2.*(-u(i,j-1,k)+u(i,j,k))/(dy1(j-1)+dy1(j))
        cov2(i,j,k) = nou2*diu2(i,j,k)

        kp1 = k+1
        nou3 = (dx1(i+1)*w(i,j,k)+dx1(i)*w(i+1,j,k)) /(dx1(i)+dx1(i+1))
        diu3(i,j,kp1) = (-u(i,j,k)+u(i,j,kp1))/dzs(k)
        cov3(i,j,kp1) = nou3*diu3(i,j,kp1)

       if (k==1) then
       nou3 = 0.5*(dx1(i+1)*w(i,j,1)+dx1(i)*w(i+1,j,1))/(dx1(i)+dx1(i+1))
       diu3(i,j,1)=uspd(i,j)*0.4/alog(0.5*dzn(1)/0.1)/(0.5*dzn(1))/0.4*u(i,j,1)/uspd(i,j)
       cov3(i,j,1) = nou3*diu3(i,j,1)
       end if

        nou4 = (dy1(j+1)*u(i-1,j,k)+dy1(j)*u(i-1,j+1,k)) /(dy1(j)+dy1(j+1))
        diu4(i,j,k) = 2.*(-v(i-1,j,k)+v(i,j,k))/(dx1(i-1)+dx1(i))
        cov4(i,j,k) = (nou4-u0)*diu4(i,j,k)
        kp1 = k+1
        nou6 = (dy1(j+1)*w(i,j,k)+dy1(j)*w(i,j+1,k)) /(dy1(j)+dy1(j+1))
        diu6(i,j,kp1) = (-v(i,j,k)+v(i,j,kp1))/dzs(k)
        cov6(i,j,kp1) = nou6*diu6(i,j,kp1)
      if (k==1) then
       nou6 = 0.5*(dy1(j+1)*w(i,j,1)+dy1(j)*w(i,j+1,1))/(dy1(j)+dy1(j+1))
       diu6(i,j,1)=vspd(i,j)*0.4/alog(0.5*dzn(1)/0.1)/(0.5*dzn(1))/0.4*v(i,j,1)/vspd(i,j)
       cov6(i,j,1) = nou6*diu6(i,j,1)
      end if
      if (k<kp) then
        nou7 = (dzn(k+1)*u(i-1,j,k)+dzn(k)*u(i-1,j,k+1)) /(dzn(k)+dzn(k+1))
        diu7(i,j,k) = 2.*(-w(i-1,j,k)+w(i,j,k))/(dx1(i-1)+dx1(i))
        cov7(i,j,k) = (nou7-u0)*diu7(i,j,k)
        end if
      if (k<kp) then
        nou8 = (dzn(k+1)*v(i,j-1,k)+dzn(k)*v(i,j-1,k+1)) /(dzn(k)+dzn(k+1))
        diu8(i,j,k) = 2.*(-w(i,j-1,k)+w(i,j,k))/(dy1(j-1)+dy1(j))
        cov8(i,j,k) = nou8*diu8(i,j,k)
        end if
      end do
      end do
      end do
! ====================================
#ifdef MPI
    if (isBottomRow(procPerRow)) then
#endif
! WV: nou1 is only ever accessed between 1 and {i/j/k}p
      do k = 1,kp
      do j = 1,jp
      if (i==ip) then
!        nou1(ip+1,j,k) = nou1(ip,j,k)
        diu1(ip+1,j,k) = diu1(ip,j,k)
        cov1(ip+1,j,k) = cov1(ip,j,k)
      end if
      end do
      end do
#ifdef MPI
    end if
#endif
#if !defined(MPI) || (PROC_PER_ROW==1)
! WV: nou2 is only ever accessed between 1 and {i/j/k}p
      do k = 1,kp
      do i = 1,ip
        if (j==jp) then
!        nou2(i,0,k) = nou2(i,jp,k)
        diu2(i,0,k) = diu2(i,jp,k)
        cov2(i,0,k) = cov2(i,jp,k)
!        nou2(i,jp+1,k) = nou2(i,1,k)
        diu2(i,jp+1,k) = diu2(i,1,k)
        cov2(i,jp+1,k) = cov2(i,1,k)
        end if
      end do
      end do
#else
!    call sideflowRightLeft(nou2, procPerRow, jp+1, 1, 1, 2, 1, 2)
    call sideflowRightLeft(diu2, procPerRow, jp+1, 1, 1, 2, 1, 2)
    call sideflowRightLeft(cov2, procPerRow, jp+1, 1, 1, 2, 1, 2)
!    call sideflowLeftRight(nou2, procPerRow, 2, jp+2, 1, 2, 1, 2)
    call sideflowLeftRight(diu2, procPerRow, 2, jp+2, 1, 2, 1, 2)
    call sideflowLeftRight(cov2, procPerRow, 2, jp+2, 1, 2, 1, 2)
#endif
#ifdef MPI
    if (isBottomRow(procPerRow)) then
#endif
! WV: nou4 is only ever accessed between 1 and {i/j/k}p
      do k = 1,kp
      do j = 1,jp
!        nou4(ip+1,j,k) = nou4(ip,j,k)
        diu4(ip+1,j,k) = diu4(ip,j,k)
        cov4(ip+1,j,k) = cov4(ip,j,k)
      end do
      end do
#ifdef MPI
    end if
#endif
#if !defined(MPI) || (PROC_PER_ROW==1)
! WV: nou4 is only ever accessed between 1 and {i/j/k}p
      do k = 1,kp
      do i = 1,ip
!        nou5(i,0,k) = nou5(i,jp,k)
        diu5(i,0,k) = diu5(i,jp,k)
        cov5(i,0,k) = cov5(i,jp,k)
!        nou5(i,jp+1,k) = nou5(i,1,k)
        diu5(i,jp+1,k) = diu5(i,1,k)
        cov5(i,jp+1,k) = cov5(i,1,k)
      end do
      end do
#else
!    call sideflowRightLeft(nou5, procPerRow, jp+2, 2, 2, 2, 1, 2)
    call sideflowRightLeft(diu5, procPerRow, jp+2, 2, 2, 2, 1, 2)
    call sideflowRightLeft(cov5, procPerRow, jp+2, 2, 2, 2, 1, 2)
!    call sideflowLeftRight(nou5, procPerRow, 3, jp+3, 2, 2, 1, 2)
    call sideflowLeftRight(diu5, procPerRow, 3, jp+3, 2, 2, 1, 2)
    call sideflowLeftRight(cov5, procPerRow, 3, jp+3, 2, 2, 1, 2)
#endif
#ifdef MPI
    if (isBottomRow(procPerRow)) then
#endif
! WV: nou7 is only ever accessed between 1 and {i/j/k}p, even kp-1
      do k = 1,kp-1
      do j = 1,jp
!        nou7(ip+1,j,k) = nou7(ip,j,k)
        diu7(ip+1,j,k) = diu7(ip,j,k)
        cov7(ip+1,j,k) = cov7(ip,j,k)
      end do
      end do
#ifdef MPI
    end if
#endif
#if !defined(MPI) || (PROC_PER_ROW==1)
! WV: nou8 is only ever accessed between 1 and {i/j/k}p, even kp-1
      do k = 1,kp-1
      do i = 1,ip
!        nou8(i,0,k) = nou8(i,jp,k)
        diu8(i,0,k) = diu8(i,jp,k)
        cov8(i,0,k) = cov8(i,jp,k)
!        nou8(i,jp+1,k) = nou8(i,1,k)
        diu8(i,jp+1,k) = diu8(i,1,k)
        cov8(i,jp+1,k) = cov8(i,1,k)
      end do
      end do
#else
!    call sideflowRightLeft(nou8, procPerRow, jp+1, 1, 1, 2, 1, 3)
    call sideflowRightLeft(diu8, procPerRow, jp+1, 1, 1, 2, 1, 3)
    call sideflowRightLeft(cov8, procPerRow, jp+1, 1, 1, 2, 1, 3)
!    call sideflowLeftRight(nou8, procPerRow, 2, jp+2, 1, 2, 1, 3)
    call sideflowLeftRight(diu8, procPerRow, 2, jp+2, 1, 2, 1, 3)
    call sideflowLeftRight(cov8, procPerRow, 2, jp+2, 1, 2, 1, 3)
#endif
! --les
#ifdef MPI
    if (isBottomRow(procPerRow)) then
#endif
      do k = 1,kp+1
      do j = 1,jp+1
        diu2(ip+1,j,k) = diu2(ip,j,k)
        diu3(ip+1,j,k) = diu3(ip,j,k)
      end do
      end do
#ifdef MPI
    end if
    if (isTopRow(procPerRow)) then
#endif
      do k = 1,kp+1
      do j = 1,jp+1
        diu2(0,j,k) = diu2(1,j,k)
        diu3(0,j,k) = diu3(1,j,k)
      end do
      end do
#ifdef MPI
    end if
#endif
#if !defined(MPI) || (PROC_PER_ROW==1)
      do k = 1,kp+1
      do i = 1,ip+1
        diu4(i,0,k) = diu4(i,jp,k)
        diu6(i,0,k) = diu6(i,jp,k)
      end do
      end do
#else
    call sideflowRightLeft(diu4, procPerRow, jp+1, 1, 1, 1, 1, 1)
    call sideflowRightLeft(diu6, procPerRow, jp+1, 1, 1, 1, 1, 1)
#endif

#ifdef MPI
#ifdef NESTED_LES
!   if (syncTicks == 0  .and. n > 2) then
   if (syncTicks == 0) then
#endif

    call exchangeRealHalos(diu1, procPerRow, neighbours, 1, 2, 2, 2)
    call exchangeRealHalos(cov1, procPerRow, neighbours, 1, 2, 2, 2)

    call exchangeRealHalos(diu2, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov2, procPerRow, neighbours, 1, 2, 1, 2)

    call exchangeRealHalos(diu3, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov3, procPerRow, neighbours, 1, 2, 1, 2)
!    call exchangeRealHalos(nou4, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(diu4, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov4, procPerRow, neighbours, 1, 2, 1, 2)

    call exchangeRealHalos(diu5, procPerRow, neighbours, 1, 2, 2, 2)
    call exchangeRealHalos(cov5, procPerRow, neighbours, 1, 2, 2, 2)
!    call exchangeRealHalos(nou6, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(diu6, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov6, procPerRow, neighbours, 1, 2, 1, 2)
!    call exchangeRealHalos(nou7, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(diu7, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov7, procPerRow, neighbours, 1, 2, 1, 2)
!    call exchangeRealHalos(nou8, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(diu8, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov8, procPerRow, neighbours, 1, 2, 1, 2)

    call exchangeRealHalos(diu9, procPerRow, neighbours, 1, 2, 1, 2)
    call exchangeRealHalos(cov9, procPerRow, neighbours, 1, 2, 1, 2)
#ifdef NESTED_LES
   end if
#endif
#endif

!uspd,vspd
!        do j=1,jp
!        do i=1,ip
!         uspd(i,j)=u(i,j,1)/uspd(i,j)
!         vspd(i,j)=v(i,j,1)/vspd(i,j)
!        end do
!        end do



#ifdef WV_DEBUG
    print *, 'F95 DIU SUMS:',sum(diu1),sum(diu2),sum(diu3),sum(diu4),sum(diu5),sum(diu6),sum(diu7),sum(diu8),sum(diu9)
#endif

end subroutine vel2_scalarized
#endif

end module module_vel2
