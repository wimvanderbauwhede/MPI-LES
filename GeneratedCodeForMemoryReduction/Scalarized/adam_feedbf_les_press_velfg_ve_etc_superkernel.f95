module singleton_module_adam_feedbf_les_press_velfg_ve_etc_superkernel

contains

subroutine adam_feedbf_les_press_velfg_ve_etc_superkernel_scal(f,g,h,fold,gold,hold,usum,u,bmask1, &
      vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,fx,fy,fz,dx1,dy1,dzs,dzn,delx1,sm,dxs,dys,rhs, &
      global_rhsav_array,global_area_array,rhsav,nrd,p0,p1,global_pav_array,global_pco_array,pav, &
      ro,state_ptr)
use module_adam_feedbf_les_press_velfg_ve_etc_superkernel_init
    integer :: global_id
    integer :: i
    integer :: i_range
    integer :: i_rel
    integer :: j
    integer :: j_range
    integer :: j_rel
    integer :: k
    integer :: k_rel
  real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: u
  real, dimension((-1):(ip+1),0:(jp+1),0:(kp+1)), intent(in) :: bmask1
  real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: v
  real, dimension(0:(ip+1),(-1):(jp+1),0:(kp+1)), intent(in) :: cmask1
  real, dimension(0:(ip+1),(-1):(jp+1),(-1):(kp+1)), intent(in) :: w
  real, dimension(0:(ip+1),0:(jp+1),0:(kp+1)), intent(in) :: dmask1
  real, dimension(1:1), intent(in) :: alpha
  real, dimension(1:1), intent(in) :: dt
  real, dimension(1:1), intent(in) :: beta
  real, dimension(0:ip,0:jp,0:kp), intent(out) :: fx
  real, dimension(0:ip,0:jp,0:kp), intent(out) :: fy
  real, dimension(0:ip,0:jp,0:kp), intent(out) :: fz
  real, dimension((-1):(ip+1)), intent(in) :: dx1
  real, dimension(0:(jp+1)), intent(in) :: dy1
  real, dimension((-1):(kp+2)), intent(in) :: dzs
  real, dimension((-1):(kp+2)), intent(in) :: dzn
  real, dimension(1:kp), intent(in) :: delx1
  real, dimension((-1):(ip+1),(-1):(jp+1),0:(kp+1)), intent(out) :: sm
  real, dimension(0:ip), intent(in) :: dxs
  real, dimension(0:jp), intent(in) :: dys
  real, dimension(0:ip,0:jp,0:kp), intent(inout) :: f
  real, dimension(0:ip,0:jp,0:kp), intent(inout) :: g
  real, dimension(0:ip,0:jp,0:kp), intent(inout) :: h
  real, dimension(0:(ip+1),0:(jp+1),0:(kp+1)), intent(out) :: rhs
  real, dimension(1:1), intent(in) :: rhsav
  integer, dimension(1:1), intent(in) :: nrd
  real, dimension(0:(ip+2),0:(jp+2),0:(kp+1)), intent(in) :: p0
  real, dimension(1:1), intent(in) :: pav
  real, dimension(1:1), intent(in) :: ro
  real, dimension(1:ip,1:jp,1:kp), intent(in) :: fold
  real, dimension(1:ip,1:jp,1:kp), intent(in) :: gold
  real, dimension(1:ip,1:jp,1:kp), intent(in) :: hold
  real, dimension(0:ip,0:jp,0:kp), intent(in) :: usum
  real, dimension(0:ip,0:jp,0:kp), intent(in) :: vsum
  real, dimension(0:ip,0:jp,0:kp), intent(in) :: wsum
  real, dimension(1:nunits), intent(out) :: global_rhsav_array
  real, dimension(1:nunits), intent(out) :: global_area_array
  real, dimension(0:(ip+2),0:(jp+2),0:(kp+1)), intent(out) :: p1
  real, dimension(1:nunits), intent(out) :: global_pav_array
  real, dimension(1:nunits), intent(out) :: global_pco_array
  integer :: state
  integer, dimension(1:1), intent(In) :: state_ptr
  parameter(st_adam_map_36=0)
  parameter(st_feedbf_map_49=1)
  parameter(st_feedbf_map_67=2)
  parameter(st_les_map_119=3)
  parameter(st_les_map_164=4)
  parameter(st_les_map_202=5)
  parameter(st_les_map_240=6)
  parameter(st_press_map_65=7)
  parameter(st_press_reduce_78=8)
  parameter(st_press_map_87=9)
  parameter(st_press_map_97=10)
  parameter(st_press_reduce_129=11)
  parameter(st_press_map_138=12)
  parameter(st_velfg_map_95=13)
  parameter(st_velfg_map_135=14)
  parameter(st_velfg_map_175=15)
  parameter(st_velnw_map_36=16)
  parameter(st_velnw_map_44=17)
  parameter(st_velnw_map_52=18)
  state = state_ptr(1) 
! SUPERKERNEL BODY
  select case(state)
    case (st_adam_map_36)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call adam_map_36(f,g,h,fold,gold,hold)

    case (st_feedbf_map_49)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call feedbf_map_49(usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,fx,fy,fz)

    case (st_feedbf_map_67)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call feedbf_map_67(f,fx,g,fy,h,fz)

    case (st_les_map_119)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call les_map_119(u,dx1,dy1,dzs,v,w,dzn,delx1,sm)

    case (st_les_map_164)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 2)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call les_map_164(sm,dy1,dx1,dzn,u,v,dzs,w,dxs,f)

    case (st_les_map_202)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 2)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call les_map_202(sm,dy1,dx1,dzn,u,v,dzs,w,dys,g)

    case (st_les_map_240)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call les_map_240(sm,dzn,dx1,dy1,u,dzs,w,v,h)

    case (st_press_map_65)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call press_map_65(u,dx1,v,dy1,w,dzn,f,g,h,rhs,dt)

    case (st_press_reduce_78)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call press_reduce_78(dx1,dy1,dzn,rhs,global_rhsav_array,global_area_array)

    case (st_press_map_87)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call press_map_87(rhs,rhsav)

    case (st_press_map_97)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call press_map_97(dzs,dys,dxs,nrd,p0,rhs,p1)

    case (st_press_reduce_129)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call press_reduce_129(p0,dx1,dy1,dzn,global_pav_array,global_pco_array)

    case (st_press_map_138)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call press_map_138(p0,pav)

    case (st_velfg_map_95)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call velfg_map_95(u,dx1,v,dy1,w,dzs,dzn,f)

    case (st_velfg_map_135)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call velfg_map_135(dy1,u,v,dx1,w,dzs,dzn,g)

    case (st_velfg_map_175)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call velfg_map_175(dzn,u,w,dx1,v,dy1,h)

    case (st_velnw_map_36)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call velnw_map_36(p0,ro,dxs,u,dt,f)

    case (st_velnw_map_44)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call velnw_map_44(p0,ro,dys,v,dt,g)

    case (st_velnw_map_52)
    call get_global_id(global_id,0)
    j_range = ((300 - 1) + 1)
    i_range = ((300 - 1) + 1)
    k_rel = (global_id / (j_range * i_range))
    k = (k_rel + 1)
    j_rel = ((global_id - (k_rel * (j_range * i_range))) / i_range)
    j = (j_rel + 1)
    i_rel = ((global_id - (k_rel * (j_range * i_range))) - (j_rel * i_range))
    i = (i_rel + 1)
call velnw_map_52(p0,ro,dzs,w,dt,h)

  end select
end subroutine adam_feedbf_les_press_velfg_ve_etc_superkernel

end module singleton_module_adam_feedbf_les_press_velfg_ve_etc_superkernel

