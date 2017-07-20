
module params_common_sn
! WV I think it would be better to have
! use params_mpi
! which would contain

#ifdef MPI
!    use communication_helper
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
    integer, parameter :: ip = 300, jp = 300, kp = 80 ! so @4m, max is 1200m x 1200m
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

#ifdef MPI
#ifdef NESTED_LES
! Original grid size
    integer, parameter :: orig_grid_y = 150 !750 ! 3km
    integer, parameter :: orig_grid_x = 150 !3000 ! 12km
! Nested grid location and size
    integer, parameter :: nested_grid_x = 150 !2000 ! 4km
    integer, parameter :: nested_grid_start_x = 75 !250
    integer, parameter :: nested_grid_end_x  = nested_grid_x + nested_grid_start_x

    integer, parameter :: nested_grid_y = 150 !500 ! 1km
    integer, parameter :: nested_grid_start_y = 75 !275
    integer, parameter :: nested_grid_end_y  = nested_grid_y + nested_grid_start_y

! Nest grid resolution
    real, parameter :: dxgrid_nest = 4.0 !2.0
    real, parameter :: dygrid_nest = 4.0 !2.0
    real, parameter :: dxgrid_orig = 4.0
    real, parameter :: dygrid_orig = 4.0

! Subgrid coordinates for nest
    integer, parameter :: i_s_nest_start =  nested_grid_start_x / ip ! 76 / (300 / 4) = 1
    integer, parameter :: i_s_nest_end =  i_s_nest_start  nested_grid_x / ip - 1 ! 1 + (150 / 75)  - 1 = 2
    integer, parameter :: j_s_nest_start =  nested_grid_start_y / ip
    integer, parameter :: j_s_nest_end =  j_s_nest_start + nested_grid_y / ip

! Time steps
    real, parameter :: dt_nest = 0.05 !0.025 ! seconds
    real, parameter :: dt_orig = 0.05 ! seconds
#endif
#endif
end module params_common_sn

