
module params_common_sn
! WV I think it would be better to have
! use params_mpi
! which would contain

#ifdef MPI
    use communication_helper
    integer, parameter :: procPerRow = PROC_PER_ROW, procPerCol = PROC_PER_COL, dimensions = 2
    integer :: dimensionSizes(dimensions)
    logical :: periodicDimensions(dimensions)
    integer :: coordinates(dimensions), neighbours(2*dimensions)
    logical :: reorder
    data dimensionSizes /procPerCol,procPerRow/, periodicDimensions /.false.,.false./, &
    reorder /.false./
#endif
    integer, parameter :: ipmax = 300, jpmax = 300 !
#ifndef TEST_SMALL_DOMAIN
#ifdef MPI
    integer, parameter :: ip = ipmax/PROC_PER_COL ! rows per process
    integer, parameter :: jp = jpmax/PROC_PER_ROW ! columns per process
    integer, parameter :: kp=80
#else
    integer, parameter :: ip = 300, jp = 300, kp = 80
#endif
#else
    integer, parameter :: ip = 25, jp = 25, kp = 80
#endif

    integer, parameter :: bipmax = 300, bjpmax = 300, bx = 0, by = 0
! WV: not nice
    character(300) :: datafile = '../GIS/Kyoto_1km2_4m_with_buffer.txt'

!-- grid
    real, parameter :: dxgrid = 4. ! meters
    real, parameter :: dygrid = 4.

!-- les
    real, parameter :: cs0 = 0.14 !smagorinsky constant
 
!-- parameter for anime
    integer, parameter :: i_anime=0
    integer, parameter :: avetime=20, km_sl=80

!-- parameter for aveflow
    integer, parameter :: i_aveflow=0

!-- for output-if
    integer, parameter :: i_ifdata_out=0


end module params_common_sn

