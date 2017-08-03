module module_ifdata

      use module_bondv1
      use module_velFG
      use module_les
      use module_boundp
#if IFBF == 1
      use module_feedbfm
      use module_feedbf
#endif
#if IADAM == 1
      use module_adam
#endif
#ifdef WV_NEW
    implicit none
#endif

contains

#ifdef WV_NEW
      subroutine zero_arrays(dfu1, dfv1, dfw1)
#else
      subroutine zero_arrays(cov1,cov2,cov3,cov4,cov5,cov6,cov7,cov8,cov9, &
              dfu1, dfv1, dfw1, &
              diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,  &
!              f,g,h, &
              nou1,nou2,nou3,nou4,nou5,nou6,nou7,nou8,nou9 &
              )
#endif
#ifdef WV_NEW
    use params_common_sn
#else
    use common_sn ! create_new_include_statements() line 102
#endif
#ifndef WV_NEW
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov9
#endif
        real(kind=4), dimension(0:ip,jp,kp) , intent(Out) :: dfu1
        real(kind=4), dimension(ip,0:jp,kp) , intent(Out) :: dfv1
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: dfw1
#ifndef WV_NEW
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu9
!        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: f
!        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: g
!        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: h
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou9


              cov1 = 0.0
              cov5 = 0.0
              diu1 = 0.0
              diu5 = 0.0
              nou1 = 0.0
              nou5 = 0.0
              cov2 = 0.0
              cov3 = 0.0
              cov4 = 0.0
              cov6 = 0.0
              cov7 = 0.0
              cov8 = 0.0
              cov9 = 0.0
              diu2 = 0.0
              diu3 = 0.0
              diu4 = 0.0
              diu6 = 0.0
              diu7 = 0.0
              diu8 = 0.0
              diu9 = 0.0
              nou2 = 0.0
              nou3 = 0.0
              nou4 = 0.0
              nou6 = 0.0
              nou7 = 0.0
              nou8 = 0.0
              nou9 = 0.0
#endif
              dfu1 = 0.0
              dfv1 = 0.0
              dfw1 = 0.0
!              f = 0.0
!              g = 0.0
!              h = 0.0

      end subroutine

#ifdef WV_NEW
      subroutine ifdata( &
!#if ICAL == 1
      fold,gold,hold, time, &
!#endif
      n,u,v,w,p,usum,vsum,wsum, &
      f,g,h, &
      ical,nif)
#else
      subroutine ifdata( &
!#if ICAL == 1
      fold,gold,hold,fghold, time, &
!#endif
      n,u,v,w,p,usum,vsum,wsum, &
      delx1,dx1,dy1,dzn,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,sm,f,g,h,z2,dt, &
      dxs,cov1,cov2,cov3,dfu1,vn,cov4,cov5,cov6,dfv1,cov7,cov8,cov9,dfw1,dzs,nou1,nou5,nou9,nou2, &
      nou3,nou4,nou6,nou7,nou8,bmask1,cmask1,dmask1,alpha,beta,fx,fy,fz,amask1,zbm,ical,nif)
#endif

#ifdef WV_NEW
    use params_common_sn
#else
    use common_sn ! create_new_include_statements() line 102
#endif
!#if ICAL == 1

        real(kind=4), dimension(ip,jp,kp) , intent(InOut) :: fold
        real(kind=4), dimension(ip,jp,kp) , intent(InOut) :: gold
        real(kind=4), dimension(ip,jp,kp) , intent(InOut) :: hold
        real(kind=4), intent(InOut) :: time
        integer, intent(In) :: ical

!#endif
#ifndef WV_NEW
        real(kind=4), dimension(ip,jp,kp) , intent(InOut) :: fghold
        real(kind=4), intent(In) :: alpha
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1) , intent(Out) :: amask1
        real(kind=4), intent(In) :: beta
        real(kind=4), dimension(-1:ip+1,0:jp+1,0:kp+1) , intent(InOut) :: bmask1
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(InOut) :: cmask1
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov9
        real(kind=4), dimension(kp) , intent(Out) :: delx1

        real(kind=4), dimension(0:ip,jp,kp) , intent(Out) :: dfu1
        real(kind=4), dimension(ip,0:jp,kp) , intent(Out) :: dfv1
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: dfw1

        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu9
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1) , intent(InOut) :: dmask1

        real(kind=4), intent(In) :: dt
        real(kind=4), dimension(-1:ip+1) , intent(In) :: dx1
        real(kind=4), dimension(0:ip) , intent(In) :: dxs
        real(kind=4), dimension(0:jp+1) , intent(In) :: dy1
        real(kind=4), dimension(-1:kp+2) , intent(In) :: dzn
        real(kind=4), dimension(-1:kp+2) , intent(In) :: dzs
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: fx
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: fy
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: fz
#endif
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: f
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: g
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: h
        integer, intent(InOut) :: n
        integer, intent(in) :: nif   
#ifndef WV_NEW
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou9
#endif
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1) , intent(InOut) :: p
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(InOut) :: u
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: usum
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(InOut) :: v
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: vsum
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(InOut) :: w
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: wsum
#ifndef WV_NEW
        real(kind=4), intent(In) :: vn
        real(kind=4), dimension(-1:ip+1,-1:jp+1,0:kp+1) , intent(Out) :: sm
        real(kind=4), dimension(0:kp+2) , intent(In) :: z2
        real(kind=4), dimension(-1:ipmax+1,-1:jpmax+1) , intent(InOut) :: zbm
#endif
        real(kind=4),allocatable :: ua(:,:,:), ur(:,:,:)
        real(kind=4),allocatable :: va(:,:,:), vr(:,:,:)
        real(kind=4),allocatable :: wa(:,:,:), wr(:,:,:)
        real(kind=4),allocatable :: usuma(:,:,:)
        real(kind=4),allocatable :: vsuma(:,:,:)
        real(kind=4),allocatable :: wsuma(:,:,:)
        real(kind=4),allocatable :: pa(:,:,:),pr(:,:,:)
        real(kind=4),allocatable :: fa(:,:,:),fr(:,:,:)
        real(kind=4),allocatable :: ga(:,:,:),gr(:,:,:)
        real(kind=4),allocatable :: ha(:,:,:),hr(:,:,:)
        real(kind=4),allocatable :: folda(:,:,:)
        real(kind=4),allocatable :: golda(:,:,:)
        real(kind=4),allocatable :: holda(:,:,:)

#if defined( MPI_NEW_WV ) || defined(WV_NEW)
    character(len=70) :: filename
    integer :: i,j,k
#endif

#if IADAM == 1
        character(len=70) :: data21dummy
        integer :: n1,n2
#endif

!    if ((im/=ip) .or. (jp/=jp) .or. (kp/=kp)) then
!            print *, "ip,kp,kp is different from ip,jp,kp, aborting!"
!            call exit(-1)
!    end if


!
!
!
!      ical =  0;initial start(grid write)
!           =  1;continuous data read,start
!           = 10;continuous data write
!
!#if ICAL == 1
       if(ical == 1) then
#ifdef MPI
        if (isMaster()) then
            write(filename, '("../data/data30",i6.6, ".dat")') nif
            open(unit=30,file=filename,form='unformatted',status='unknown')
            read(30) n,time
        end if

#ifdef MPI_NEW_WV2
        if (isMaster()) then
            allocate(ua(ipmax,jpmax,kp))
            read(30) (((ua(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if
        allocate(ur(ip,jp,kp))
        call MPI_Scatter(ur, ip*jp*kp, MPI_REAL, ua, ipmax*jpmax*kp, MPI_REAL, 0, communicator, ierror)
        u(1:ip,1:jp,1:kp)=ur
        if (isMaster()) deallocate(ua)
        deallocate(ur)
#else
        allocate(ua(0:ipmax+1,-1:jpmax+1,0:kp+1))
         
        if (isMaster()) then
        read(30) (((ua(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if

        call distributeifu(ua, ip, jp, kp, ipmax, jpmax, procPerRow)
        do k=1,kp
        do j=1,jp
        do i=1,ip
         u(i,j,k)=ua(i,j,k)
        end do
        end do
        end do
        deallocate(ua)
#endif

#ifdef MPI_NEW_WV2
        if (isMaster()) then
            allocate(va(ipmax,jpmax,kp))
            read(30) (((va(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if
        allocate(vr(ip,jp,kp))
        call MPI_Scatter(vr, ip*jp*(kp+2), MPI_REAL, va, ipmax*jpmax*kp, MPI_REAL, 0, communicator, ierror)
        v(1:ip,1:jp,1:kp)=vr
        if (isMaster()) deallocate(va)
        deallocate(vr)
#else
        allocate(va(0:ipmax+1,-1:jpmax+1,0:kp+1))               
        if (isMaster()) then
        read(30) (((va(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if

        call distributeifu(va, ip, jp, kp, ipmax, jpmax, procPerRow)
        do k=1,kp
        do j=1,jp
        do i=1,ip
         v(i,j,k)=va(i,j,k)
        end do
        end do
        end do
        deallocate(va)
#endif

#ifdef MPI_NEW_WV2
        if (isMaster()) then
            allocate(wa(ipmax,jpmax,kp))
            read(30) (((wa(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if
        allocate(wr(ip,jp,kp))
        call MPI_Scatter(wr, ip*jp*(kp+2), MPI_REAL, wa, ipmax*jpmax*kp, MPI_REAL, 0, communicator, ierror)
        w(1:ip,1:jp,1:kp)=wr
        if (isMaster()) deallocate(wa)
        deallocate(wr)
#else

        allocate(wa(0:ipmax+1,-1:jpmax+1,-1:kp+1))
        if (isMaster()) then
        read(30) (((wa(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if

        call distributeifw(wa, ip, jp, kp, ipmax, jpmax, procPerRow)
        do k=1,kp
        do j=1,jp
        do i=1,ip
         w(i,j,k)=wa(i,j,k)
        end do
        end do
        end do
        deallocate(wa)
#endif

#ifdef MPI_NEW_WV2
        if (isMaster()) then
            allocate(pa(ipmax,jpmax,kp))
            read(30) (((pa(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if
        allocate(pr(ip,jp,kp))
        call MPI_Scatter(pr, ip*jp*(kp+2), MPI_REAL, pa, ipmax*jpmax*kp, MPI_REAL, 0, communicator, ierror)
        p(1:ip,1:jp,1:kp)=pr
        if (isMaster()) deallocate(pa)
        deallocate(pr)
#else
        allocate(pa(0:ipmax+2,0:jpmax+2,0:kp+1))
        if (isMaster()) then
        read(30) (((pa(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if

        call distributeifp(pa, ip, jp, kp, ipmax, jpmax, procPerRow)
        do k=1,kp
        do j=1,jp
        do i=1,ip
         p(i,j,k)=pa(i,j,k)
        end do
        end do
        end do
        deallocate(pa)
#endif

        allocate(usuma(0:ipmax,0:jpmax,0:kp))
        if (isMaster()) then
        read(30) (((usuma(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if

        call distributeifusum(usuma, ip, jp, kp, ipmax, jpmax, procPerRow)
        do k=1,kp
        do j=1,jp
        do i=1,ip
         usum(i,j,k)=usuma(i,j,k)
        end do
        end do
        end do
        deallocate(usuma)


        allocate(vsuma(0:ipmax,0:jpmax,0:kp))
        if (isMaster()) then
        read(30) (((vsuma(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if

        call distributeifusum(vsuma, ip, jp, kp, ipmax, jpmax, procPerRow)
        do k=1,kp
        do j=1,jp
        do i=1,ip
         vsum(i,j,k)=vsuma(i,j,k)
        end do
        end do
        end do
        deallocate(vsuma)


        allocate(wsuma(0:ipmax,0:jpmax,0:kp))
        if (isMaster()) then
        read(30) (((wsuma(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if

        call distributeifusum(wsuma, ip, jp, kp, ipmax, jpmax, procPerRow)
        do k=1,kp
        do j=1,jp
        do i=1,ip
         wsum(i,j,k)=wsuma(i,j,k)
        end do
        end do
        end do
        deallocate(wsuma)

        if (isMaster()) then
        close(30)
        end if


!31
        if (isMaster()) then
        write(filename, '("../data/data31",i6.6, ".dat")') nif
        open(unit=31,file=filename,form='unformatted',status='unknown')
        end if

        allocate(fa(0:ipmax,0:jpmax,0:kp))

        if (isMaster()) then
        read(31) (((fa(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if

        call distributeiff(fa, ip, jp, kp, ipmax, jpmax, procPerRow)
        do k=1,kp
        do j=1,jp
        do i=1,ip
         f(i,j,k)=fa(i,j,k)
        end do
        end do
        end do
        deallocate(fa)

        allocate(ga(0:ipmax,0:jpmax,0:kp))
        if (isMaster()) then
        read(31) (((ga(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if

        call distributeiff(ga, ip, jp, kp, ipmax, jpmax, procPerRow)
        do k=1,kp
        do j=1,jp
        do i=1,ip
         g(i,j,k)=ga(i,j,k)
        end do
        end do
        end do
        deallocate(ga)


        allocate(ha(0:ipmax,0:jpmax,0:kp))
        if (isMaster()) then
        read(31) (((ha(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if

        call distributeiff(ha, ip, jp, kp, ipmax, jpmax, procPerRow)
        do k=1,kp
        do j=1,jp
        do i=1,ip
         h(i,j,k)=ha(i,j,k)
        end do
        end do
        end do
        deallocate(ha)


        allocate(folda(ipmax,jpmax,kp))

        if (isMaster()) then
        read(31) (((folda(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if

        call distributeiffold(folda, ip, jp, kp, ipmax, jpmax, procPerRow)
        do k=1,kp
        do j=1,jp
        do i=1,ip
         fold(i,j,k)=folda(i,j,k)
        end do
        end do
        end do
        deallocate(folda)


        allocate(golda(ipmax,jpmax,kp))

        if (isMaster()) then
        read(31) (((golda(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if

        call distributeiffold(golda, ip, jp, kp, ipmax, jpmax, procPerRow)
        do k=1,kp
        do j=1,jp
        do i=1,ip
         gold(i,j,k)=golda(i,j,k)
        end do
        end do
        end do
        deallocate(golda)


        allocate(holda(ipmax,jpmax,kp))

        if (isMaster()) then
        read(31) (((holda(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        end if

        call distributeiffold(holda, ip, jp, kp, ipmax, jpmax, procPerRow)
        do k=1,kp
        do j=1,jp
        do i=1,ip
         hold(i,j,k)=holda(i,j,k)
        end do
        end do
        end do
        deallocate(holda)

        if (isMaster()) then
        close(31)
        end if

!#endif


            call boundp2(p)

#else
! NO MPI

!30
        write(filename, '("../data/data30",i6.6, ".dat")') nif
        open(unit=30,file=filename,form='unformatted',status='unknown')
        read(30) n,time

        read(30) (((u(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        read(30) (((v(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        read(30) (((w(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        read(30) (((p(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        read(30) (((usum(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        read(30) (((vsum(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        read(30) (((wsum(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        close(30)

!31
        write(filename, '("../data/data31",i6.6, ".dat")') nif
        open(unit=31,file=filename,form='unformatted',status='unknown')
        read(31) (((f(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        read(31) (((g(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        read(31) (((h(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        read(31) (((fold(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        read(31) (((gold(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        read(31) (((hold(i,j,k),i=1,ipmax),j=1,jpmax),k=1,kp)
        close(31)

            call boundp2(p)

#endif
! WV: I added this routine to explicitly set all arrays to zero
#ifdef WV_NEW
!        call zero_arrays(dfu1, dfv1, dfw1)
#else
        call zero_arrays( &
              cov1,cov2,cov3,cov4,cov5,cov6,cov7,cov8,cov9, &
              dfu1, dfv1, dfw1, &
              diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,  &
!              f,g,h, &
              nou1,nou2,nou3,nou4,nou5,nou6,nou7,nou8,nou9 &
              )
#endif

     end if

            call boundp1(p)
            call boundp2(p)

!        call velfg(kp,jp,ip,dx1,cov1,cov2,cov3,dfu1,diu1,diu2,dy1,diu3,dzn,vn,f,cov4,cov5,cov6,dfv1, &
!      diu4,diu5,diu6,g,cov7,cov8,cov9,dfw1,diu7,diu8,diu9,dzs,h,nou1,u,nou5,v,nou9,w,nou2,nou3, &
!      nou4,nou6,nou7,nou8)
#if IFBF == 1
!        if(ifbf == 1) then
!        call feedbf(kp,jp,ip,usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,fx,fy,fz,f,g, &
!      h)
!        call feedbfm(kp,jp,ip,amask1,bmask1,cmask1,dmask1,zbm,z2,dzn)
!        endif
#endif
!        call les(kp,delx1,dx1,dy1,dzn,jp,ip,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,sm,f,g,h)

! --adam
! WV iadam is not defined!
#if IADAM == 1
! WV        if(iadam.eq.1) then
!            n1=1
!            n2=2
!            data21dummy=""
!          call adam(n1,n2,data21dummy,fold,ip,jp,kp,gold,hold,fghold,f,g,h)
! WV        end if
#endif
!
      end subroutine ifdata

end module module_ifdata
