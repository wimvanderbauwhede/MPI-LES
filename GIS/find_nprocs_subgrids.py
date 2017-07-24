#!/usr/bin/env python
from math import sqrt,floor

def find_nprocs_subgrids(grid_setup,np_max):
    nested_grid_x = grid_setup['nested_grid_x']
    nested_grid_y = grid_setup['nested_grid_y']
    orig_grid_x = grid_setup['orig_grid_x']
    orig_grid_y = grid_setup['orig_grid_y']
    dxgrid_orig = grid_setup['dxgrid_orig']
    dygrid_orig = grid_setup['dygrid_orig']
    dxgrid_nest = grid_setup['dxgrid_nest']
    dygrid_nest = grid_setup['dygrid_nest']
    nested_grid_start_x = grid_setup['nested_grid_start_x']
    nested_grid_start_y = grid_setup['nested_grid_start_y']    
    print('Possible solutions with number of processes (nprocs) and grid sizes (ip,jp)')
    print('nprocs, ip x jp, proc_col x proc_row') #, npx, npy')        
    for np in range(np_max,0,-1):
        b = nested_grid_y / nested_grid_x
        a_x = (nested_grid_x*(1-dxgrid_nest/dxgrid_orig) + orig_grid_x)        
        a_y = (nested_grid_y*(1-dygrid_nest/dygrid_orig) + orig_grid_y)
        a_xy = a_x*a_y
        a2 = a_xy/np
#         print(np,a_x,a_y,a2)
        if (a2==floor(a2)):   
#             print('Candidate:',a2, ' for np=',np)
            b_n = nested_grid_y / nested_grid_x     
            b_o = orig_grid_y / orig_grid_x
            if b_n < 1:
                b_n = 1/b_n
            if b_o < 1: 
                b_o = 1/b_o
            b = b_n
            if b_o > b:
                b = b_o
            if b<2:
                b=2                
            dx_guess = sqrt(a2)
            dx_b = floor(dx_guess/b)
            dx_e = int(dx_guess*b)
#             print(dx_b, dx_e,b)
            for dx in range(dx_b,dx_e):
                dy = a2/dx
                nx = nested_grid_start_x/dx
                ny = nested_grid_start_y/dy
                rx = (orig_grid_x-(nested_grid_x*dxgrid_nest)/dxgrid_orig - nested_grid_start_x)/dx
                ry = (orig_grid_y-(nested_grid_y*dygrid_nest)/dygrid_orig - nested_grid_start_y)/dy
                if (dy==floor(dy)):
#                     print(np,' Valid dy:',dx,int(dy),';',nx,ny,rx,ry)
                        p_col = (orig_grid_x + nested_grid_x*(1 - dxgrid_nest/dxgrid_orig) ) / dx
                        p_row = (orig_grid_y + nested_grid_y*(1 - dygrid_nest/dygrid_orig) ) / dy
                        if (nx==floor(nx)) and (ny==floor(ny)) and \
                            (rx==floor(rx)) and (ry==floor(ry)) and \
                            (p_row==floor(p_row)) and (p_col==floor(p_col)) :
                            print(np,dx,'x',int(dy),int(p_col),'x',int(p_row))
#                               int(nx),int(ny),rx,ry)


grid_setup = {
'nested_grid_x' : 200,
'nested_grid_y' : 200,
'nested_grid_start_x' : 100,
'nested_grid_start_y' : 100,
'orig_grid_x' : 300,
'orig_grid_y' : 300,
'dxgrid_orig' : 4,
'dygrid_orig' : 4,
'dxgrid_nest' : 2,
'dygrid_nest' : 2
}

np_max  = 24 
                   
find_nprocs_subgrids(grid_setup, np_max)                            


grid_setup_Kyoto = {
'nested_grid_x' : 4000/2,
'nested_grid_y' : 1000/2,
'nested_grid_start_x' : 1000/4,
'nested_grid_start_y' : 800/4,
'orig_grid_x' : 12000/4,
'orig_grid_y' : 2600/4,
'dxgrid_orig' : 4,
'dygrid_orig' : 4,
'dxgrid_nest' : 2,
'dygrid_nest' : 2
}
np_max = 256

find_nprocs_subgrids(grid_setup_Kyoto, np_max)  
