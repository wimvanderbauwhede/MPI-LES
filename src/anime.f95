module module_anime
#ifdef MPI
    use communication_helper_real
#ifdef NESTED_LES
    use nesting_support
#endif
#ifdef MPI_NEW_WV
    use params_common_sn
    implicit none
#endif
#endif

contains

! Added by WV for use with MPI_NEW_WV2
#ifdef MPI_NEW_WV2
real function calc_avg_ua(ua,i,j,k) result(avg)
    real, dimension(:,:,:), intent(In) :: ua
    integer, intent(In) :: i,j,k
    real :: uam1

    if (i==1) then
        uam1=ua(i,j,k)
    else
        uam1=ua(i-1,j,k)
    end if
    avg =real(0.5*(uam1+ua(i,j,k)))
end function  calc_avg_ua

real function calc_avg_va(va,i,j,k) result(avg)
    real, dimension(:,:,:), intent(In) :: va
    integer, intent(In) :: i,j,k
    real :: va_jm1

    if (j==1) then
        va_jm1=va(i,jpmax,k)
    else
        va_jm1=va(i,j-1,k)
    end if
    avg =real(0.5*(va_jm1+va(i,j,k)))
end function  calc_avg_va

real function calc_avg_wa(wa,i,j,k) result(avg)
    real, dimension(:,:,:), intent(In) :: wa
    integer, intent(In) :: i,j,k
    real :: wa_km1

    if (k==1) then
        wa_km1=0.0
    else
        wa_km1=wa(i,j,k-1)
    end if
    avg =real(0.5*(wa_km1+wa(i,j,k)))
end function  calc_avg_wa
#endif

subroutine anime(n,n0,n1,nmax,km,jm,im,dxl,dx1,dyl,dy1,z2,data22,data23,u,w,v,p,amask1,zbm)

    use common_sn ! create_new_include_statements() line 102
    real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1) , intent(InOut) :: amask1
    character(len=70), intent(In) :: data22
    character(len=70), intent(In) :: data23
    real(kind=4), dimension(-1:ip+1) , intent(In) :: dx1
    real(kind=4), dimension(0:ip) , intent(In) :: dxl
    real(kind=4), dimension(0:jp+1) , intent(In) :: dy1
    real(kind=4), dimension(0:jp) , intent(In) :: dyl
    integer, intent(In) :: im
    integer, intent(In) :: jm
    integer, intent(In) :: km
    integer, intent(In) :: n
    integer, intent(In) :: n0
    integer, intent(In) :: n1
    integer, intent(In) :: nmax
    real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: u
    real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
    real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(In) :: w
    real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1) , intent(In) :: p
    real(kind=4), dimension(kp+2) , intent(In) :: z2
    real(kind=4), dimension(-1:ipmax+1,-1:jpmax+1) , intent(In)  :: zbm
!average_out
#ifdef MPI_NEW_WV
    integer :: irec, i,j,k
    character(len=70) :: filename
#endif
#ifdef MPI_NEW_WV2
    real(kind=4), dimension(ip,jp,kp)  :: uani
    real(kind=4), dimension(ip,jp,kp)  :: vani
    real(kind=4), dimension(ip,jp,kp) :: wani
    real(kind=4), dimension(ip,jp,kp) :: pani
#else
    real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: uani
    real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: vani
    real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) :: wani
    real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1) :: pani

#endif
!mpi_out
!    real(kind=4), dimension(0:ipmax+1,-1:jpmax+1,0:kp+1)  :: ua
!    real(kind=4), dimension(0:ipmax+1,-1:jpmax+1,0:kp+1)  :: va
!    real(kind=4), dimension(0:ipmax+1,-1:jpmax+1,-1:kp+1)  :: wa
!    real(kind=4), dimension(0:ipmax+1,0:jpmax+1,0:kp+1)   :: amask1a
!    real(kind=4), dimension(-1:ipmax+1,-1:jpmax+1,1) , intent(in) :: zbm1
!    character(len=70) :: filename

    real(kind=4),allocatable :: ua(:,:,:)
    real(kind=4),allocatable :: va(:,:,:)
    real(kind=4),allocatable :: wa(:,:,:)
    real(kind=4),allocatable :: pa(:,:,:)
    real(kind=4),allocatable :: amask1a(:,:,:)
#ifdef NESTED_LES
    integer :: nn
#endif
!    if(n == n0.or.n == nmax.or.mod(n,1000) == 0.) then
!        do k = 1,km
!            do j = 1,jm
!                do i = 1,im
!                    a1(i,j,k) = real(dxl(i-1)+dx1(i))
!                    a2(i,j,k) = real(dyl(j-1)+dy1(j))
!                    a3(i,j,k) = real(z2(k))
!                end do
!            end do
!        end do
!        open(unit=22,file=data22,form='unformatted',status='unknown')
!        write(22) im,jm,km
!        write(22) (((real(a1(i,j,k)),i=1,im),j=1,jm),k=1,km), &
!                  (((real(a3(i,j,k)),i=1,im),j=1,jm),k=1,km), &
!                  (((real(a2(i,j,k)),i=1,im),j=1,jm),k=1,km)
!        close(unit=22)
!   end if
#ifdef MPI

!reset
!      if(n.eq.n1) then  
!      do k=0,km
!      do j=0,jm
!      do i=0,im
!      uani(i,j,k)=0.
!      vani(i,j,k)=0.
!      wani(i,j,k)=0.
!      pani(i,j,k)=0.
!      end do
!      end do
!      end do
!      end if

#ifdef NESTED_LES
    if (syncTicks == 0) then
#endif
!average_out
#ifndef MPI_NEW_WV2
      do k=0,km
      do j=0,jm
      do i=0,im
      uani(i,j,k)=uani(i,j,k)+u(i,j,k)
      vani(i,j,k)=vani(i,j,k)+v(i,j,k)
      wani(i,j,k)=wani(i,j,k)+w(i,j,k)
      pani(i,j,k)=pani(i,j,k)+p(i,j,k)
      end do
      end do
      end do
#else
      do k=1,km
          do j=1,jm
              do i=1,im
                  uani(i,j,k)=uani(i,j,k)+u(i,j,k)
                  vani(i,j,k)=vani(i,j,k)+v(i,j,k)
                  wani(i,j,k)=wani(i,j,k)+w(i,j,k)
                  pani(i,j,k)=pani(i,j,k)+p(i,j,k)
              end do
          end do
      end do
#endif
!#ifdef NESTED_LES
!    end if
!#endif
!
!#ifdef NESTED_LES
!    if (syncTicks == 0) then
!#endif

    ! Normalizing the time step
#ifdef NESTED_LES
    nn = n
    if (inNestedGrid()) then
        nn = n0+(n-n0)/int(dt_orig/dt_nest)
    end if
    if(n>=n1 .and. mod(nn,avetime) == 0) then !default
#else
    if(n.ge.n1.and.mod(n,avetime).eq.0) then !default
#endif
!        print *, n,nn,rank,uani(ip/2,jp/2,kp/2)
!       if(mod(n,avetime).eq.0) then !default

           if (isMaster()) then
               write(filename, '("../out/data23",i6.6, ".dat")') n
               open(unit=23,file=filename,form='unformatted',access='direct',recl=4*ipmax*jpmax)
           end if

! ---- ua ----

#ifdef MPI_NEW_WV2
        if (isMaster()) then
            allocate(ua(ipmax,jpmax,kp))
        end if
        call MPI_Gather(uani, ip*jp*kp, MPI_REAL, ua, ip*jp*kp, MPI_REAL, 0, communicator, ierror)
        call checkMPIError()

       if (isMaster()) then
           ua=ua/real(avetime)
!           print *,ua(ipmax/2,jpmax/2,kp/2)
           irec = 1
           do  k=1,km_sl
                write(23,rec=irec) ((calc_avg_ua(ua,i,j,k),i=1,ipmax),j=1,jpmax)
                irec = irec + 1
           end do
           deallocate(ua)
       end if
#else
       allocate(ua(0:ipmax+1,-1:jpmax+1,0:kp+1))
       call distributeu(ua, uani, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
            do k = 1,km
                do j = 1,jm
                    do i = 1,im
                     ua(i,j,k) = uani(i,j,k)
                    end do
                end do
            end do

            do k=1,km
                do j=1,jpmax
                 do i=1,ipmax
                    ua(i,j,k)=ua(i,j,k)/real(avetime)
                 end do
                end do
            end do
            !       print *,ua(ipmax/2,jpmax/2,kp/2)
            !boundary
            do k = 1,km
                 do j = 1,jpmax
                    ua(0,j,k) = ua(1,j,k)
                 end do
            end do

            irec = 1
            do  k=1,km_sl
              write(23,rec=irec) ((real(0.5*(ua(i-1,j,k)+ua(i,j,k))),i=1,ipmax),j=1,jpmax)
              irec = irec + 1
            end do
       end if
       deallocate(ua)
#endif
! ---- wa -----
#ifdef MPI_NEW_WV2
        if (isMaster()) then
            allocate(wa(ipmax,jpmax,kp))
        end if
        call MPI_Gather(wani, ip*jp*kp, MPI_REAL, wa, ip*jp*kp, MPI_REAL, 0, communicator, ierror)
        call checkMPIError()

       if (isMaster()) then
           wa=wa/real(avetime)
!           irec = 1
           do  k=1,km_sl
                write(23,rec=irec) ((calc_avg_wa(wa,i,j,k),i=1,ipmax),j=1,jpmax)
                irec = irec + 1
           end do
           deallocate(wa)
       end if
#else
        allocate(wa(0:ipmax+1,-1:jpmax+1,-1:kp+1))
        call distributew(wa, wani, ip, jp, kp, ipmax, jpmax, procPerRow)
        if (isMaster()) then
            do k = 1,km
                do j = 1,jm
                    do i = 1,im
                     wa(i,j,k) = wani(i,j,k)
                    end do
                end do
            end do

            do k=1,km
                do j=1,jpmax
                 do i=1,ipmax
                    wa(i,j,k)=wa(i,j,k)/real(avetime)
                 end do
                end do
            end do

        !boundary
            do j = 1,jpmax
                do i = 1,ipmax
                    wa(i,j,0) = 0.0
                end do
            end do

        !       do  k=1,km
            do  k=1,km_sl
                write(23,rec=irec) ((real(0.5*(wa(i,j,k-1)+wa(i,j,k))),i=1,ipmax),j=1,jpmax)
                irec = irec + 1
            end do
        end if
        deallocate(wa)
#endif


! ---- va ----
#ifdef MPI_NEW_WV2
        if (isMaster()) allocate(va(ipmax,jpmax,kp))

        call MPI_Gather(vani, ip*jp*kp, MPI_REAL, va, ip*jp*kp, MPI_REAL, 0, communicator, ierror)
        call checkMPIError()

       if (isMaster()) then
           va=va/real(avetime)
!           irec = 1
           do  k=1,km_sl
                write(23,rec=irec) ((calc_avg_va(va,i,j,k),i=1,ipmax),j=1,jpmax)
                irec = irec + 1
           end do
           deallocate(va)
       end if

#else
       allocate(va(0:ipmax+1,-1:jpmax+1,0:kp+1))
        call distributev(va, vani, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
          do k = 1,km
            do j = 1,jm
                do i = 1,im
                 va(i,j,k) = vani(i,j,k)
                end do
            end do
          end do

       do k=1,km
        do j=1,jpmax
         do i=1,ipmax
            va(i,j,k)=va(i,j,k)/real(avetime)
         end do
        end do
       end do

!boundary
            do k = 1,km
                do i = 1,ipmax
                  va(i,0,k) = va(i,jpmax,k)
                end do
            end do


       do  k=1,km_sl
       write(23,rec=irec) ((real(0.5*(va(i,j-1,k)+va(i,j,k))),i=1,ipmax),j=1,jpmax)
       irec = irec + 1
       end do
!------ if you output p, have to comment out this close
!       close(23)
!----------------------
       end if
       deallocate(va)
#endif


!--------pressure-----------
#ifdef MPI_NEW_WV2
       if (isMaster()) allocate(pa(ipmax,jpmax,kp))
       call MPI_Gather(pani, ip*jp*kp, MPI_REAL, pa, ip*jp*kp, MPI_REAL, 0, communicator, ierror)
       call checkMPIError()

       if (isMaster()) then
           pa=pa/real(avetime)
!           irec = 1
           do  k=1,km_sl
                write(23,rec=irec) ((pa(i,j,k),i=1,ipmax),j=1,jpmax)
                irec = irec + 1
           end do
           close(23)
           deallocate(pa)
       end if
#else
       allocate(pa(0:ipmax+2,0:jpmax+2,0:kp+1))
       call distributep(pa, pani, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
         do k = 1,km
            do j = 1,jm
                do i = 1,im
                 pa(i,j,k) = pani(i,j,k)
                end do
            end do
          end do

       do k=1,km
        do j=1,jpmax
         do i=1,ipmax
            pa(i,j,k)=pa(i,j,k)/real(avetime)
         end do
        end do
       end do


!       do  k=1,km
       do  k=1,km_sl
       write(23,rec=irec) ((pa(i,j,k),i=1,ipmax),j=1,jpmax)
       irec = irec + 1
       end do
       close(23)
       end if
       deallocate(pa)
#endif


#ifdef MPI_NEW_WV2
      uani=0.
      vani=0.
      wani=0.
      pani=0.
#else
      do k=0,km
      do j=0,jm
      do i=0,im
      uani(i,j,k)=0.
      vani(i,j,k)=0.
      wani(i,j,k)=0.
      pani(i,j,k)=0.
      end do
      end do
      end do
#endif
      end if ! timestep condition
#ifdef NESTED_LES
    end if
#endif
#endif

end subroutine anime


! WV: TODO
!data30,31
subroutine ifdata_out(n,n0,n1,nmax,time,km,jm,im,u,w,v,p,usum,vsum,wsum,f,g,h,fold,gold,hold)


    use common_sn ! create_new_include_statements() line 102
    integer, intent(In) :: im
    integer, intent(In) :: jm
    integer, intent(In) :: km
    integer, intent(In) :: n
    integer, intent(In) :: n0
    integer, intent(In) :: nmax
    real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: u
    real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
    real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(In) :: w
    integer, intent(In) :: n1
    real, intent(In) :: time
!ifdata
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In)  :: usum
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In)  :: vsum
    real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In)  :: wsum

    real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1),intent(In)  :: p

    real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In)  :: f
    real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In)  :: g
    real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In)  :: h
!    real(kind=4), dimension(0:ipmax,0:jpmax,0:kp) :: fa
!    real(kind=4), dimension(0:ipmax,0:jpmax,0:kp) :: ga
!    real(kind=4), dimension(0:ipmax,0:jpmax,0:kp) :: ha

    real(kind=4), dimension(ip,jp,kp) , intent(In)  :: fold
    real(kind=4), dimension(ip,jp,kp) , intent(In)  :: gold
    real(kind=4), dimension(ip,jp,kp) , intent(In)  :: hold
!    real(kind=4), dimension(ipmax,jpmax,kp) :: folda
!    real(kind=4), dimension(ipmax,jpmax,kp) :: golda
!    real(kind=4), dimension(ipmax,jpmax,kp) :: holda

!mpi_out
!    real(kind=4), dimension(0:ipmax+1,-1:jpmax+1,0:kp+1)  :: ua
!    real(kind=4), dimension(0:ipmax+1,-1:jpmax+1,0:kp+1)  :: va
!    real(kind=4), dimension(0:ipmax+1,-1:jpmax+1,-1:kp+1)  :: wa
!    real(kind=4), dimension(0:ipmax+1,0:jpmax+1,0:kp+1)   :: amask1a
!    real(kind=4), dimension(-1:ipmax+1,-1:jpmax+1,1) , intent(in) :: zbm1
!    character(len=70) :: filename

    real(kind=4),allocatable :: ua(:,:,:)
    real(kind=4),allocatable :: va(:,:,:)
    real(kind=4),allocatable :: wa(:,:,:)
    real(kind=4),allocatable :: usuma(:,:,:)
    real(kind=4),allocatable :: vsuma(:,:,:)
    real(kind=4),allocatable :: wsuma(:,:,:)
    real(kind=4),allocatable :: pa(:,:,:)
    real(kind=4),allocatable :: fa(:,:,:)
    real(kind=4),allocatable :: ga(:,:,:)
    real(kind=4),allocatable :: ha(:,:,:)
    real(kind=4),allocatable :: folda(:,:,:)
    real(kind=4),allocatable :: golda(:,:,:)
    real(kind=4),allocatable :: holda(:,:,:)
#ifdef MPI_NEW_WV
    integer :: irec, i,j,k
    character(len=70) :: filename
#endif

! WV: considering  do n = n0,nmax, and n1 = 1, n != n1-1
#ifdef NESTED_LES
!    if (syncTicks == 0 .and. n > 2) then
    if (syncTicks == 0) then
       if(  n == nmax ) then !default
#else
       if((n.eq.n1-1).or.(n.eq.nmax))  then      
#endif
        if (isMaster()) then
        write(filename, '("../data/data30",i6.6, ".dat")') n

        open(unit=30,file=filename,form='unformatted',status='replace')

        write(30) n,time
        end if

       allocate(ua(0:ipmax+1,-1:jpmax+1,0:kp+1))
        call distributeu(ua, u, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
          do k = 1,km
            do j = 1,jm
                do i = 1,im
                    ua(i,j,k) = u(i,j,k)
                end do
            end do
          end do
        write(30) (((ua(i,j,k),i=1,ipmax),j=1,jpmax),k=1,km)
        end if
        deallocate(ua)


       allocate(va(0:ipmax+1,-1:jpmax+1,0:kp+1))
        call distributev(va, v, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
          do k = 1,km
            do j = 1,jm
                do i = 1,im 
                    va(i,j,k) = v(i,j,k)
                end do
            end do
          end do
        write(30) (((va(i,j,k),i=1,ipmax),j=1,jpmax),k=1,km)
        end if
        deallocate(va)


       allocate(wa(0:ipmax+1,-1:jpmax+1,-1:kp+1))
        call distributew(wa, w, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
          do k = 1,km
            do j = 1,jm
                do i = 1,im 
                    wa(i,j,k) = w(i,j,k)
                end do
            end do
          end do
        write(30) (((wa(i,j,k),i=1,ipmax),j=1,jpmax),k=1,km)
        end if
        deallocate(wa)


       allocate(pa(0:ipmax+2,0:jpmax+2,0:kp+1))
        call distributep(pa, p, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
          do k = 1,km
            do j = 1,jm
                do i = 1,im 
                    pa(i,j,k) = p(i,j,k)
                end do
            end do
          end do
        write(30) (((pa(i,j,k),i=1,ipmax),j=1,jpmax),k=1,km)
        end if
        deallocate(pa)


       allocate(usuma(0:ipmax,0:jpmax,0:kp))
        call distributeusum(usuma, usum, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
              do k = 1,km
                do j = 1,jm
                    do i = 1,im
                        usuma(i,j,k) = usum(i,j,k)
                    end do
                end do
              end do
              write(30) (((usuma(i,j,k),i=1,ipmax),j=1,jpmax),k=1,km)
        end if
        deallocate(usuma)


      allocate(vsuma(0:ipmax,0:jpmax,0:kp))
        call distributeusum(vsuma, vsum, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
              do k = 1,km
                do j = 1,jm
                    do i = 1,im
                        vsuma(i,j,k) = vsum(i,j,k)
                    end do
                end do
              end do
              write(30) (((vsuma(i,j,k),i=1,ipmax),j=1,jpmax),k=1,km)
        end if
        deallocate(vsuma)

       allocate(wsuma(0:ipmax,0:jpmax,0:kp))
       call distributeusum(wsuma, wsum, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
              do k = 1,km
                do j = 1,jm
                    do i = 1,im
                        wsuma(i,j,k) = wsum(i,j,k)
                    end do
                end do
              end do
            write(30) (((wsuma(i,j,k),i=1,ipmax),j=1,jpmax),k=1,km)
            close(30)
        end if
        deallocate(wsuma)

        if (isMaster()) then
             write(filename, '("../data/data31",i6.6, ".dat")') n
             open(unit=31,file=filename,form='unformatted',status='replace')
        end if


     allocate(fa(0:ipmax,0:jpmax,0:kp))
        call distributef(fa, f, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
          do k = 1,km
            do j = 1,jm
                do i = 1,im
                    fa(i,j,k) = f(i,j,k)
                end do
            end do
          end do
        write(31) (((fa(i,j,k),i=1,ipmax),j=1,jpmax),k=1,km)
        end if
        deallocate(fa)


     allocate(ga(0:ipmax,0:jpmax,0:kp))
        call distributef(ga, g, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
          do k = 1,km
            do j = 1,jm
                do i = 1,im
                    ga(i,j,k) = g(i,j,k)
                end do
            end do
          end do
        write(31) (((ga(i,j,k),i=1,ipmax),j=1,jpmax),k=1,km)
        end if
        deallocate(ga)


     allocate(ha(0:ipmax,0:jpmax,0:kp))
        call distributef(ha, h, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
          do k = 1,km
            do j = 1,jm
                do i = 1,im
                    ha(i,j,k) = h(i,j,k)
                end do
            end do
          end do
        write(31) (((ha(i,j,k),i=1,ipmax),j=1,jpmax),k=1,km)
        end if
        deallocate(ha)


     allocate(folda(ipmax,jpmax,kp))
        call distributefold(folda, fold, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
          do k = 1,km
            do j = 1,jm
                do i = 1,im
                    folda(i,j,k) = fold(i,j,k)
                end do
            end do
          end do
        write(31) (((folda(i,j,k),i=1,ipmax),j=1,jpmax),k=1,km)
        end if
        deallocate(folda)


     allocate(golda(ipmax,jpmax,kp))
        call distributefold(golda, gold, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
          do k = 1,km
            do j = 1,jm
                do i = 1,im
                    golda(i,j,k) = gold(i,j,k)
                end do
            end do
          end do
        write(31) (((golda(i,j,k),i=1,ipmax),j=1,jpmax),k=1,km)
        end if
        deallocate(golda)


     allocate(holda(ipmax,jpmax,kp))
        call distributefold(holda, fold, ip, jp, kp, ipmax, jpmax, procPerRow)
       if (isMaster()) then
          do k = 1,km
            do j = 1,jm
                do i = 1,im
                    holda(i,j,k) = hold(i,j,k)
                end do
            end do
          end do
        write(31) (((holda(i,j,k),i=1,ipmax),j=1,jpmax),k=1,km)
        close(31)

        end if
        deallocate(holda)

        end if
#ifdef NESTED_LES
    end if
#endif
end subroutine ifdata_out




end module module_anime
