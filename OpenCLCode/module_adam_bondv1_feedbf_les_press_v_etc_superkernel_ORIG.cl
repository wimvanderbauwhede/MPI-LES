


#include "array_index_f2c1d.h"

void adam_map_36(__global float *f,__global float *g,__global float *h,__global float *fold,__global float *gold,__global float *hold) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: fd,i,j,k,gd,hd
    float fd;
    int i;
    int j;
    int k;
    float gd;
    float hd;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        k_rel = (global_id*1.0/(j_range*i_range));
        k = k_rel+1;
        j_rel = ((global_id-(k_rel*(j_range*i_range)))*1.0/i_range);
        j = j_rel+1;
        i_rel = ((global_id-(k_rel*(j_range*i_range)))-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                                fd = f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)];
                                gd = g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)];
                                hd = h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)];
                                f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = 1.5*f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]-0.5*fold[F3D2C(((ip - 1 )+1),((jp - 1 )+1) , 1,1,1 , i,j,k)];
                                g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = 1.5*g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]-0.5*gold[F3D2C(((ip - 1 )+1),((jp - 1 )+1) , 1,1,1 , i,j,k)];
                                h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = 1.5*h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]-0.5*hold[F3D2C(((ip - 1 )+1),((jp - 1 )+1) , 1,1,1 , i,j,k)];
                                fold[F3D2C(((ip - 1 )+1),((jp - 1 )+1) , 1,1,1 , i,j,k)] = fd;
                                gold[F3D2C(((ip - 1 )+1),((jp - 1 )+1) , 1,1,1 , i,j,k)] = gd;
                                hold[F3D2C(((ip - 1 )+1),((jp - 1 )+1) , 1,1,1 , i,j,k)] = hd;
                }
void bondv1_map_38(__global float *z2,__global float *dzn,__global float *u,__global float *v,__global float *w) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: u_val,k,i,j
    float u_val;
    int k;
    int i;
    int j;
        // parallelfortran: synthesised loop variable decls
    int i_range;
    int k_range;
    int j_range;
    int i_rel;
    int k_rel;
    int j_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        i_range = ((1-0)+1);
        k_range = ((78-1)+1);
        j_range = ((300-1)+1);
        i_rel = (global_id*1.0/(k_range*j_range));
        i = i_rel+0;
        k_rel = ((global_id-(i_rel*(k_range*j_range)))*1.0/j_range);
        k = k_rel+1;
        j_rel = ((global_id-(i_rel*(k_range*j_range)))-(k_rel*j_range));
        j = j_rel+1;
        // parallelfortran: original code
                                        u_val = 5.*pow((float)(((z2[F1D2C(0 , k)]+0.5*dzn[F1D2C((-1) , k)])*1.0/600.)));
                                        u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)] = u_val;
                                        v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)] = 0.0;
                                        w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)] = 0.0;
                    }
void bondv1_map_48(__global float *u,__global float *v,__global float *w) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: k,i,j
    int k;
    int i;
    int j;
        // parallelfortran: synthesised loop variable decls
    int i_range;
    int k_range;
    int j_range;
    int i_rel;
    int k_rel;
    int j_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        i_range = ((1-0)+1);
        k_range = ((80-79)+1);
        j_range = ((300-1)+1);
        i_rel = (global_id*1.0/(k_range*j_range));
        i = i_rel+0;
        k_rel = ((global_id-(i_rel*(k_range*j_range)))*1.0/j_range);
        k = k_rel+79;
        j_rel = ((global_id-(i_rel*(k_range*j_range)))-(k_rel*j_range));
        j = j_rel+1;
        // parallelfortran: original code
                                        u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)] = u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)];
                                        v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)] = 0.0;
                                        w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)] = 0.0;
                    }
void bondv1_map_58(__global float *u,__global float *v,__global float *w) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: k,i,j
    int k;
    int i;
    int j;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        i_range = ((300-2)+1);
        k_rel = (global_id*1.0/(j_range*i_range));
        k = k_rel+1;
        j_rel = ((global_id-(k_rel*(j_range*i_range)))*1.0/i_range);
        j = j_rel+1;
        i_rel = ((global_id-(k_rel*(j_range*i_range)))-(j_rel*i_range));
        i = i_rel+2;
        // parallelfortran: original code
                                        u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)] = u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , j,k)];
                                        v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)] = v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , j,k)];
                                        w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)] = w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , j,k)];
                    }
void bondv1_reduce_69(__global float *u,__global float *global_aaa_array) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // missing args: 
        // local vars: k,j
    int k;
    int j;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int k_rel;
    int j_rel;
    int chunk_size;
    int local_id;
    int local_id_fortran;
    int group_id;
    int group_id_fortran;
    int global_id;
    int r_iter;
    int local_chunk_size;
    int start_position;
    #if NTH > 1
        // arrays prefixed with __PH0__ should be declared using the __PH1__ modifier in c kernel version
    __local float local_aaa_array[(NTH - 1 +1)];
    #endif
    float clk_local_mem_fence;
    float local_aaa;
        group_id = get_group_id(0);
        global_id = get_global_id(0);
    #if NTH > 1
        local_id = get_local_id(0);
        // local_id_fortran and group_id_fortran are used to reconcile the fact that fortran arrays are referenced from 1
        // not 0 like other opencl supporting languages
        local_id_fortran = local_id+1;
    #endif
        group_id_fortran = group_id+1;
    #if NTH > 1
        local_chunk_size = (((80*300)*1.0/NTH)*1.0/NUNITS);
    #else
        local_chunk_size = ((80*300)*1.0/NUNITS);
    #endif
        start_position = local_chunk_size*global_id;
        local_aaa = 0;
    for (r_iter = start_position;r_iter <= ((start_position + local_chunk_size) - 1);r_iter += 1) {
                k_range = ((80-1)+1);
                j_range = ((300-1)+1);
                k_rel = r_iter*1.0/j_range;
                k = k_rel+1;
                j_rel = (r_iter-(k_rel*j_range));
                j = j_rel+1;
                local_aaa = (float)fmax(local_aaa,u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , j,k)]);
    }
    #if NTH > 1
        local_aaa_array[F1D2C(1 , local_id_fortran)] = local_aaa;
    #ifdef BARRIER_OK
    barrier(CLK_LOCAL_MEM_FENCE);
    #endif
        local_aaa = 0;
    for (r_iter = 1;r_iter <= NTH;r_iter += 1) {
                local_aaa = (float)fmax(local_aaa,local_aaa_array[F1D2C(1 , r_iter)]);
    }
    #endif
        global_aaa_array[F1D2C(1 , group_id_fortran)] = local_aaa;
    }
void bondv1_reduce_76(__global float *u,__global float *global_bbb_array) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // missing args: 
        // local vars: k,j
    int k;
    int j;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int k_rel;
    int j_rel;
    int chunk_size;
    int local_id;
    int local_id_fortran;
    int group_id;
    int group_id_fortran;
    int global_id;
    int r_iter;
    int local_chunk_size;
    int start_position;
    #if NTH > 1
        // arrays prefixed with __PH0__ should be declared using the __PH1__ modifier in c kernel version
    __local float local_bbb_array[(NTH - 1 +1)];
    #endif
    float clk_local_mem_fence;
    float local_bbb;
        group_id = get_group_id(0);
        global_id = get_global_id(0);
    #if NTH > 1
        local_id = get_local_id(0);
        // local_id_fortran and group_id_fortran are used to reconcile the fact that fortran arrays are referenced from 1
        // not 0 like other opencl supporting languages
        local_id_fortran = local_id+1;
    #endif
        group_id_fortran = group_id+1;
    #if NTH > 1
        local_chunk_size = (((80*300)*1.0/NTH)*1.0/NUNITS);
    #else
        local_chunk_size = ((80*300)*1.0/NUNITS);
    #endif
        start_position = local_chunk_size*global_id;
        local_bbb = 1e38;
    for (r_iter = start_position;r_iter <= ((start_position + local_chunk_size) - 1);r_iter += 1) {
                k_range = ((80-1)+1);
                j_range = ((300-1)+1);
                k_rel = r_iter*1.0/j_range;
                k = k_rel+1;
                j_rel = (r_iter-(k_rel*j_range));
                j = j_rel+1;
                local_bbb = (float)fmin(local_bbb,u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , j,k)]);
    }
    #if NTH > 1
        local_bbb_array[F1D2C(1 , local_id_fortran)] = local_bbb;
    #ifdef BARRIER_OK
    barrier(CLK_LOCAL_MEM_FENCE);
    #endif
        local_bbb = 1e38;
    for (r_iter = 1;r_iter <= NTH;r_iter += 1) {
                local_bbb = (float)fmin(local_bbb,local_bbb_array[F1D2C(1 , r_iter)]);
    }
    #endif
        global_bbb_array[F1D2C(1 , group_id_fortran)] = local_bbb;
    }
void bondv1_map_83(__global float *u,__global float *dt,__global float *uout,__global float *dxs,__global float *v,__global float *w) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: k,j
    int k;
    int j;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int k_rel;
    int j_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        k_rel = global_id*1.0/j_range;
        k = k_rel+1;
        j_rel = (global_id-(k_rel*j_range));
        j = j_rel+1;
        // parallelfortran: original code
                            u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , ip,j,k)] = u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , ip,j,k)]-dt*uout*(u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , ip,j,k)]-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , ip-1,j,k)])*1.0/dxs[F1D2C(0 , ip)];
                            v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , ip+1,j,k)] = v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , ip+1,j,k)]-dt*uout*(v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , ip+1,j,k)]-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , ip,j,k)])*1.0/dxs[F1D2C(0 , ip)];
                            w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , ip+1,j,k)] = w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , ip+1,j,k)]-dt*uout*(w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , ip+1,j,k)]-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , ip,j,k)])*1.0/dxs[F1D2C(0 , ip)];
              }
void bondv1_map_98(__global float *u,__global float *v,__global float *w) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: k,i
    int k;
    int i;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int i_range;
    int k_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = (((80+1)-0)+1);
        i_range = (((300+1)-0)+1);
        k_rel = global_id*1.0/i_range;
        k = k_rel+0;
        i_rel = (global_id-(k_rel*i_range));
        i = i_rel+0;
        // parallelfortran: original code
                        u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,k)] = u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,jp,k)];
                        u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,jp+1,k)] = u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,k)];
                        v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,k)] = v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,jp,k)];
                        v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,jp+1,k)] = v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,k)];
    if ((k<80)) {
                        w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,k)] = w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,jp,k)];
                        w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,jp+1,k)] = w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,k)];
    }
    }
void bondv1_map_116(__global float *u,__global float *v) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: i,j
    int i;
    int j;
        // parallelfortran: synthesised loop variable decls
    int j_range;
    int i_range;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        j_range = (((300+1)-0)+1);
        i_range = (((300+1)-0)+1);
        j_rel = global_id*1.0/i_range;
        j = j_rel+0;
        i_rel = (global_id-(j_rel*i_range));
        i = i_rel+0;
        // parallelfortran: original code
                        u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)] = -u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)];
                        u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,kp+1)] = u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,kp)];
                        v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)] = -v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)];
                        v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,kp+1)] = v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,kp)];
            }
void bondv1_map_128(__global float *w) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: i,j
    int i;
    int j;
        // parallelfortran: synthesised loop variable decls
    int j_range;
    int i_range;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        j_range = (((300+1)-(-1))+1);
        i_range = (((300+1)-0)+1);
        j_rel = global_id*1.0/i_range;
        j = (j_rel+(-1));
        i_rel = (global_id-(j_rel*i_range));
        i = i_rel+0;
        // parallelfortran: original code
                        w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j)] = 0.0;
                        w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,kp)] = 0.0;
            }
void feedbf_map_48(__global float *zbm,__global float *z2,__global float *dzn,__global float *usum,__global float *u,__global float *vsum,__global float *v,__global float *wsum,__global float *w,__global float *alpha,__global float *dt,__global float *beta,__global float *f,__global float *g,__global float *h) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: abcd_mask,i,j,k,f1x,f1y,f1z,f2x,f2y,f2z,fx,fy,fz
    float *abcd_mask;
    int i;
    int j;
    int k;
    float f1x;
    float f1y;
    float f1z;
    float f2x;
    float f2y;
    float f2z;
    float fx;
    float fy;
    float fz;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        k_rel = (global_id*1.0/(j_range*i_range));
        k = k_rel+1;
        j_rel = ((global_id-(k_rel*(j_range*i_range)))*1.0/i_range);
        j = j_rel+1;
        i_rel = ((global_id-(k_rel*(j_range*i_range)))-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                                abcd_mask[F1D2C(0 , )] = 1.;
                                *abcd_mask = 0.;
                                abcd_mask[F1D2C(0 , )] = 0.;
                                abcd_mask[F1D2C(0 , )] = 0.;
                if (zbm[F2D2C((((ipmax+1) - (-1) )+1) , (-1),(-1) , i,j)]>z2[F1D2C(0 , k)]+0.5*dzn[F1D2C((-1) , k)]) {
                                        abcd_mask[F1D2C(0 , )] = 0.0;
                                        *abcd_mask = 1.0;
                                        abcd_mask[F1D2C(0 , )] = 1.0;
                                        abcd_mask[F1D2C(0 , )] = 1.0;
                }
                                usum[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = (usum[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*;
                                vsum[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = (vsum[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*abcd_mask[F1D2C(0 , )];
                                wsum[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = (wsum[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)])*abcd_mask[F1D2C(0 , )];
                                f1x = alpha*usum[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]*dt;
                                f1y = alpha*vsum[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]*dt;
                                f1z = alpha*wsum[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]*dt;
                                f2x = beta*u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]*;
                                f2y = beta*v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]*abcd_mask[F1D2C(0 , )];
                                f2z = beta*w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]*abcd_mask[F1D2C(0 , )];
                                fx = f1x+f2x;
                                fy = f1y+f2y;
                                fz = f1z+f2z;
                                f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]+fx;
                                g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]+fy;
                                h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]+fz;
                }
void les_map_121(__global float *dx1,__global float *dy1,__global float *dzn,__global float *delx1) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: k
    int k;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int k_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_rel = global_id;
        k = k_rel+1;
        // parallelfortran: original code
                    delx1[F1D2C(1 , k)] = pow((float)((dx1[F1D2C((-1) , )]*dy1[F1D2C(0 , )]*dzn[F1D2C((-1) , k)])),(float)(1.*1.0/3.));
          }
void les_map_124(__global float *u,__global float *dx1,__global float *dy1,__global float *dzs,__global float *v,__global float *w,__global float *dzn,__global float *delx1,__global float *sm) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: k,diu1_i_j_k,i,j,dudxx1,diu2_i_j_k,diu2_im1_j_k,diu2_im1_jp1_k,diu2_i_jp1_k,dudyx1,diu3_i_j_k,diu3_im1_j_k,diu3_im1_j_kp1,diu3_i_j_kp1,dudzx1,diu4_i_j_k,diu4_i_jm1_k,diu4_ip1_j_k,diu4_ip1_jm1_k,dvdxx1,diu5_i_j_k,dvdyx1,diu6_i_j_k,diu6_i_jm1_k,diu6_i_jm1_kp1,diu6_i_j_kp1,dvdzx1,diu7_i_j_k,diu7_i_j_km1,diu7_ip1_j_k,diu7_ip1_j_km1,dwdxx1,diu8_i_j_k,diu8_i_j_km1,diu8_i_jp1_k,diu8_i_jp1_km1,dwdyx1,diu9_i_j_k,dwdzx1,csx1
    int k;
    float diu1_i_j_k;
    int i;
    int j;
    float dudxx1;
    float diu2_i_j_k;
    float diu2_im1_j_k;
    float diu2_im1_jp1_k;
    float diu2_i_jp1_k;
    float dudyx1;
    float diu3_i_j_k;
    float diu3_im1_j_k;
    float diu3_im1_j_kp1;
    float diu3_i_j_kp1;
    float dudzx1;
    float diu4_i_j_k;
    float diu4_i_jm1_k;
    float diu4_ip1_j_k;
    float diu4_ip1_jm1_k;
    float dvdxx1;
    float diu5_i_j_k;
    float dvdyx1;
    float diu6_i_j_k;
    float diu6_i_jm1_k;
    float diu6_i_jm1_kp1;
    float diu6_i_j_kp1;
    float dvdzx1;
    float diu7_i_j_k;
    float diu7_i_j_km1;
    float diu7_ip1_j_k;
    float diu7_ip1_j_km1;
    float dwdxx1;
    float diu8_i_j_k;
    float diu8_i_j_km1;
    float diu8_i_jp1_k;
    float diu8_i_jp1_km1;
    float dwdyx1;
    float diu9_i_j_k;
    float dwdzx1;
    float csx1;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        k_rel = (global_id*1.0/(j_range*i_range));
        k = k_rel+1;
        j_rel = ((global_id-(k_rel*(j_range*i_range)))*1.0/i_range);
        j = j_rel+1;
        i_rel = ((global_id-(k_rel*(j_range*i_range)))-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                diu1_i_j_k = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dx1[F1D2C((-1) , i)];
                dudxx1 = diu1_i_j_k;
                diu2_i_j_k = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/(dy1[F1D2C(0 , j-1)]+dy1[F1D2C(0 , j)]);
                diu2_im1_j_k = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j-1,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)])*1.0/(dy1[F1D2C(0 , j-1)]+dy1[F1D2C(0 , j)]);
                diu2_im1_jp1_k = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
                diu2_i_jp1_k = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
                dudyx1 = (diu2_im1_j_k+diu2_im1_jp1_k+diu2_i_j_k+diu2_i_jp1_k)*.25;
                diu3_i_j_k = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k-1)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dzs[F1D2C((-1) , k-1)];
                diu3_im1_j_k = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k-1)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)])*1.0/dzs[F1D2C((-1) , k-1)];
                diu3_im1_j_kp1 = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k+1)])*1.0/dzs[F1D2C((-1) , k)];
                diu3_i_j_kp1 = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k+1)])*1.0/dzs[F1D2C((-1) , k)];
                dudzx1 = (diu3_im1_j_k+diu3_im1_j_kp1+diu3_i_j_k+diu3_i_j_kp1)*.25;
                diu4_i_j_k = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)]);
                diu4_i_jm1_k = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j-1,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)]);
                diu4_ip1_j_k = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
                diu4_ip1_jm1_k = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j-1,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
                dvdxx1 = (diu4_i_j_k+diu4_i_jm1_k+diu4_ip1_j_k+diu4_ip1_jm1_k)*.25;
                diu5_i_j_k = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dy1[F1D2C(0 , j)];
                dvdyx1 = diu5_i_j_k;
                diu6_i_j_k = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k-1)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dzs[F1D2C((-1) , k-1)];
                diu6_i_jm1_k = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k-1)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)])*1.0/dzs[F1D2C((-1) , k-1)];
                diu6_i_jm1_kp1 = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k+1)])*1.0/dzs[F1D2C((-1) , k)];
                diu6_i_j_kp1 = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k+1)])*1.0/dzs[F1D2C((-1) , k)];
                dvdzx1 = (diu6_i_jm1_k+diu6_i_jm1_kp1+diu6_i_j_k+diu6_i_j_kp1)*.25;
                diu7_i_j_k = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i-1,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)]);
                diu7_i_j_km1 = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i-1,j,k-1)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k-1)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)]);
                diu7_ip1_j_k = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
                diu7_ip1_j_km1 = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , j,k-1)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i+1,j,k-1)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
                dwdxx1 = (diu7_i_j_k+diu7_i_j_km1+diu7_ip1_j_k+diu7_ip1_j_km1)*.25;
                diu8_i_j_k = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j-1,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)])*1.0/(dy1[F1D2C(0 , j-1)]+dy1[F1D2C(0 , j)]);
                diu8_i_j_km1 = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j-1,k-1)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k-1)])*1.0/(dy1[F1D2C(0 , j-1)]+dy1[F1D2C(0 , j)]);
                diu8_i_jp1_k = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
                diu8_i_jp1_km1 = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k-1)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j+1,k-1)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
                dwdyx1 = (diu8_i_j_k+diu8_i_j_km1+diu8_i_jp1_k+diu8_i_jp1_km1)*.25;
                diu9_i_j_k = (-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k-1)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)])*1.0/dzn[F1D2C((-1) , k)];
                dwdzx1 = diu9_i_j_k;
            csx1 = cs0;
            sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)] = pow((float)((csx1*delx1[F1D2C(1 , k)])))*sqrt(2.*(pow((float)(dudxx1))+pow((float)(dvdyx1))+pow((float)(dwdzx1)))+pow((float)(dudyx1+dvdxx1))+pow((float)(dwdyx1+dvdzx1))+pow((float)(dudzx1+dwdxx1)));
      }
void les_map_169(__global float *sm) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: k,j
    int k;
    int j;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int k_rel;
    int j_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = (((80+1)-0)+1);
        j_range = (((300+1)-(-1))+1);
        k_rel = global_id*1.0/j_range;
        k = k_rel+0;
        j_rel = (global_id-(k_rel*j_range));
        j = (j_rel+(-1));
        // parallelfortran: original code
                                        sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , j,k)] = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , j,k)];
                                        sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , ip+1,j,k)] = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , ip,j,k)];
                    }
void les_map_175(__global float *sm) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: k,i
    int k;
    int i;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int i_range;
    int k_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = (((80+1)-0)+1);
        i_range = (((300+1)-0)+1);
        k_rel = global_id*1.0/i_range;
        k = k_rel+0;
        i_rel = (global_id-(k_rel*i_range));
        i = i_rel+0;
        // parallelfortran: original code
                                        sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,jp+1,k)] = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,jp,k)];
                                        sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,k)] = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,k)];
                    }
void les_map_181(__global float *sm) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: i,j
    int i;
    int j;
        // parallelfortran: synthesised loop variable decls
    int j_range;
    int i_range;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        j_range = (((300+1)-(-1))+1);
        i_range = (((300+1)-0)+1);
        j_rel = global_id*1.0/i_range;
        j = (j_rel+(-1));
        i_rel = (global_id-(j_rel*i_range));
        i = i_rel+0;
        // parallelfortran: original code
                        sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j)] = -sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j)];
                        sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,kp+1)] = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,kp)];
            }
void les_map_187(__global float *sm,__global float *dy1,__global float *dx1,__global float *dzn,__global float *u,__global float *v,__global float *dzs,__global float *w,__global float *dxs,__global float *f) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: k,diu1_i_j_k,i,j,diu2_i_j_k,diu2_i_jp1_k,diu3_i_j_k,diu3_i_j_kp1,diu4_ip1_j_k,diu4_ip1_jm1_k,diu7_ip1_j_k,diu7_ip1_j_km1,evsx2,evsx1,evsy2,evsy1,evsz2,evsz1,diu1_ip1_j_k,visux2,visux1,visuy2,visuy1,visuz2,visuz1,vfu
    int k;
    float diu1_i_j_k;
    int i;
    int j;
    float diu2_i_j_k;
    float diu2_i_jp1_k;
    float diu3_i_j_k;
    float diu3_i_j_kp1;
    float diu4_ip1_j_k;
    float diu4_ip1_jm1_k;
    float diu7_ip1_j_k;
    float diu7_ip1_j_km1;
    float evsx2;
    float evsx1;
    float evsy2;
    float evsy1;
    float evsz2;
    float evsz1;
    float diu1_ip1_j_k;
    float visux2;
    float visux1;
    float visuy2;
    float visuy1;
    float visuz2;
    float visuz1;
    float vfu;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-2)+1);
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        k_rel = (global_id*1.0/(j_range*i_range));
        k = k_rel+2;
        j_rel = ((global_id-(k_rel*(j_range*i_range)))*1.0/i_range);
        j = j_rel+1;
        i_rel = ((global_id-(k_rel*(j_range*i_range)))-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
            evsx2 = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k)];
            evsx1 = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)];
            evsy2 = (dy1[F1D2C(0 , j+1)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]))+dy1[F1D2C(0 , j)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j+1,k)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j+1,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            evsy1 = (dy1[F1D2C(0 , j+1)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j-1,k)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j-1,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]))+dy1[F1D2C(0 , j)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            evsz2 = (dzn[F1D2C((-1) , k+1)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]))+dzn[F1D2C((-1) , k)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k+1)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k+1)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))*1.0/(dzn[F1D2C((-1) , k)]+dzn[F1D2C((-1) , k+1)]);
            evsz1 = (dzn[F1D2C((-1) , k)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k-1)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k-1)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]))+dzn[F1D2C((-1) , k-1)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))*1.0/(dzn[F1D2C((-1) , k-1)]+dzn[F1D2C((-1) , k)]);
            diu1_i_j_k = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dx1[F1D2C((-1) , i)];
            diu1_ip1_j_k = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j,k)])*1.0/dx1[F1D2C((-1) , i+1)];
            visux2 = (evsx2)*2.*diu1_ip1_j_k;
            visux1 = (evsx1)*2.*diu1_i_j_k;
            diu2_i_jp1_k = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            diu4_ip1_j_k = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
            visuy2 = (evsy2)*(diu2_i_jp1_k+diu4_ip1_j_k);
            diu2_i_j_k = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/(dy1[F1D2C(0 , j-1)]+dy1[F1D2C(0 , j)]);
            diu4_ip1_jm1_k = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j-1,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
            visuy1 = (evsy1)*(diu2_i_j_k+diu4_ip1_jm1_k);
            diu3_i_j_kp1 = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k+1)])*1.0/dzs[F1D2C((-1) , k)];
            diu7_ip1_j_k = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
            visuz2 = (evsz2)*(diu3_i_j_kp1+diu7_ip1_j_k);
            diu3_i_j_k = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k-1)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dzs[F1D2C((-1) , k-1)];
            diu7_ip1_j_km1 = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , j,k-1)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i+1,j,k-1)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
            visuz1 = (evsz1)*(diu3_i_j_k+diu7_ip1_j_km1);
            vfu = (visux2-visux1)*1.0/dxs[F1D2C(0 , i)]+(visuy2-visuy1)*1.0/dy1[F1D2C(0 , j)]+(visuz2-visuz1)*1.0/dzn[F1D2C((-1) , k)];
            f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = (f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]+vfu);
      }
void les_map_225(__global float *sm,__global float *dy1,__global float *dx1,__global float *dzn,__global float *u,__global float *v,__global float *dzs,__global float *w,__global float *uspd,__global float *dxs,__global float *f) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: i,j,evsx2,evsx1,evsy2,evsy1,evsz2,visux2,visux1,visuy2,visuy1,visuz2,visuz1,vfu,diu1_i_j_1,diu1_ip1_j_1,diu2_i_jp1_1,diu4_ip1_j_1,diu2_i_j_1,diu4_ip1_jm1_1,diu3_i_j_2,diu7_ip1_j_1
    int i;
    int j;
    float evsx2;
    float evsx1;
    float evsy2;
    float evsy1;
    float evsz2;
    float visux2;
    float visux1;
    float visuy2;
    float visuy1;
    float visuz2;
    float visuz1;
    float vfu;
    float diu1_i_j_1;
    float diu1_ip1_j_1;
    float diu2_i_jp1_1;
    float diu4_ip1_j_1;
    float diu2_i_j_1;
    float diu4_ip1_jm1_1;
    float diu3_i_j_2;
    float diu7_ip1_j_1;
        // parallelfortran: synthesised loop variable decls
    int j_range;
    int i_range;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        j_rel = global_id*1.0/i_range;
        j = j_rel+1;
        i_rel = (global_id-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
            evsx2 = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j)];
            evsx1 = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j)];
            evsy2 = (dy1[F1D2C(0 , j+1)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]))+dy1[F1D2C(0 , j)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j+1)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j+1)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            evsy1 = (dy1[F1D2C(0 , j+1)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j-1)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j-1)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]))+dy1[F1D2C(0 , j)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            evsz2 = (dzn[F1D2C((-1) , )]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]))+*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))*1.0/(+dzn[F1D2C((-1) , )]);
            diu1_i_j_1 = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)])*1.0/dx1[F1D2C((-1) , i)];
            diu1_ip1_j_1 = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j)])*1.0/dx1[F1D2C((-1) , i+1)];
            visux2 = (evsx2)*2.*diu1_ip1_j_1;
            visux1 = (evsx1)*2.*diu1_i_j_1;
            diu2_i_jp1_1 = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j+1)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            diu4_ip1_j_1 = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
            visuy2 = (evsy2)*(diu2_i_jp1_1+diu4_ip1_j_1);
            diu2_i_j_1 = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)])*1.0/(dy1[F1D2C(0 , j-1)]+dy1[F1D2C(0 , j)]);
            diu4_ip1_jm1_1 = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j-1)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
            visuy1 = (evsy1)*(diu2_i_j_1+diu4_ip1_jm1_1);
            diu3_i_j_2 = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)])*1.0/;
            diu7_ip1_j_1 = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i+1,j)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
            visuz2 = (evsz2)*(diu3_i_j_2+diu7_ip1_j_1);
            visuz1 = pow((float)((0.4*uspd[F2D2C((((ip+1) - 0 )+1) , 0,0 , i,j)]*1.0/(float)log(0.5**1.0/0.1))))*(u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)]*1.0/uspd[F2D2C((((ip+1) - 0 )+1) , 0,0 , i,j)]);
            vfu = (visux2-visux1)*1.0/dxs[F1D2C(0 , i)]+(visuy2-visuy1)*1.0/dy1[F1D2C(0 , j)]+(visuz2-visuz1)*1.0/;
            f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j)] = (f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j)]+vfu);
      }
void les_map_253(__global float *sm,__global float *dy1,__global float *dx1,__global float *dzn,__global float *u,__global float *v,__global float *dzs,__global float *w,__global float *dys,__global float *g) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: k,i,j,diu2_im1_jp1_k,diu2_i_jp1_k,diu4_i_j_k,diu4_ip1_j_k,diu5_i_j_k,diu6_i_j_k,diu6_i_j_kp1,diu8_i_jp1_k,diu8_i_jp1_km1,evsx2,evsx1,evsy2,evsy1,evsz2,evsz1,visvx2,visvx1,diu5_i_jp1_k,visvy2,visvy1,visvz2,visvz1,vfv
    int k;
    int i;
    int j;
    float diu2_im1_jp1_k;
    float diu2_i_jp1_k;
    float diu4_i_j_k;
    float diu4_ip1_j_k;
    float diu5_i_j_k;
    float diu6_i_j_k;
    float diu6_i_j_kp1;
    float diu8_i_jp1_k;
    float diu8_i_jp1_km1;
    float evsx2;
    float evsx1;
    float evsy2;
    float evsy1;
    float evsz2;
    float evsz1;
    float visvx2;
    float visvx1;
    float diu5_i_jp1_k;
    float visvy2;
    float visvy1;
    float visvz2;
    float visvz1;
    float vfv;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-2)+1);
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        k_rel = (global_id*1.0/(j_range*i_range));
        k = k_rel+2;
        j_rel = ((global_id-(k_rel*(j_range*i_range)))*1.0/i_range);
        j = j_rel+1;
        i_rel = ((global_id-(k_rel*(j_range*i_range)))-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
            evsy2 = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j+1,k)];
            evsy1 = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)];
            evsx2 = (dy1[F1D2C(0 , j+1)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]))+dy1[F1D2C(0 , j)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j+1,k)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j+1,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            evsx1 = (dy1[F1D2C(0 , j+1)]*((dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i-1,j,k)]+dx1[F1D2C((-1) , i-1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)]))+dy1[F1D2C(0 , j)]*((dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i-1,j+1,k)]+dx1[F1D2C((-1) , i-1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j+1,k)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)])))*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            evsz2 = (dzn[F1D2C((-1) , k+1)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]))+dzn[F1D2C((-1) , k)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k+1)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k+1)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))*1.0/(dzn[F1D2C((-1) , k)]+dzn[F1D2C((-1) , k+1)]);
            evsz1 = (dzn[F1D2C((-1) , k)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k-1)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k-1)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]))+dzn[F1D2C((-1) , k-1)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))*1.0/(dzn[F1D2C((-1) , k-1)]+dzn[F1D2C((-1) , k)]);
            diu2_i_jp1_k = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            diu4_ip1_j_k = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
            visvx2 = (evsx2)*(diu2_i_jp1_k+diu4_ip1_j_k);
            diu2_im1_jp1_k = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            diu4_i_j_k = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)]);
            visvx1 = (evsx1)*(diu2_im1_jp1_k+diu4_i_j_k);
            diu5_i_jp1_k = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j+1,k)])*1.0/dy1[F1D2C(0 , j+1)];
            visvy2 = (evsy2)*2.*diu5_i_jp1_k;
            diu5_i_j_k = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dy1[F1D2C(0 , j)];
            visvy1 = (evsy1)*2.*diu5_i_j_k;
            diu6_i_j_kp1 = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k+1)])*1.0/dzs[F1D2C((-1) , k)];
            diu8_i_jp1_k = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            visvz2 = (evsz2)*(diu6_i_j_kp1+diu8_i_jp1_k);
            diu6_i_j_k = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k-1)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dzs[F1D2C((-1) , k-1)];
            diu8_i_jp1_km1 = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k-1)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j+1,k-1)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            visvz1 = (evsz1)*(diu6_i_j_k+diu8_i_jp1_km1);
            vfv = (visvx2-visvx1)*1.0/dx1[F1D2C((-1) , i)]+(visvy2-visvy1)*1.0/dys[F1D2C(0 , j)]+(visvz2-visvz1)*1.0/dzn[F1D2C((-1) , k)];
            g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = (g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]+vfv);
      }
void les_map_291(__global float *sm,__global float *dy1,__global float *dx1,__global float *dzn,__global float *u,__global float *v,__global float *dzs,__global float *w,__global float *vspd,__global float *dys,__global float *g) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: i,j,evsx2,evsx1,evsy2,evsy1,evsz2,diu2_i_jp1_1,diu4_ip1_j_1,visvx2,visvx1,visvy2,visvy1,visvz2,visvz1,vfv,diu2_im1_jp1_1,diu4_i_j_1,diu5_i_jp1_1,diu5_i_j_1,diu6_i_j_2,diu8_i_jp1_1
    int i;
    int j;
    float evsx2;
    float evsx1;
    float evsy2;
    float evsy1;
    float evsz2;
    float diu2_i_jp1_1;
    float diu4_ip1_j_1;
    float visvx2;
    float visvx1;
    float visvy2;
    float visvy1;
    float visvz2;
    float visvz1;
    float vfv;
    float diu2_im1_jp1_1;
    float diu4_i_j_1;
    float diu5_i_jp1_1;
    float diu5_i_j_1;
    float diu6_i_j_2;
    float diu8_i_jp1_1;
        // parallelfortran: synthesised loop variable decls
    int j_range;
    int i_range;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        j_rel = global_id*1.0/i_range;
        j = j_rel+1;
        i_rel = (global_id-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
            evsy2 = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j+1)];
            evsy1 = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j)];
            evsx2 = (dy1[F1D2C(0 , j+1)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]))+dy1[F1D2C(0 , j)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j+1)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j+1)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            evsx1 = (dy1[F1D2C(0 , j+1)]*((dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i-1,j)]+dx1[F1D2C((-1) , i-1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)]))+dy1[F1D2C(0 , j)]*((dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i-1,j+1)]+dx1[F1D2C((-1) , i-1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j+1)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)])))*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            evsz2 = (dzn[F1D2C((-1) , )]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]))+*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))*1.0/(+dzn[F1D2C((-1) , )]);
            diu2_i_jp1_1 = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j+1)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            diu4_ip1_j_1 = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
            visvx2 = (evsx2)*(diu2_i_jp1_1+diu4_ip1_j_1);
            diu2_im1_jp1_1 = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j+1)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            diu4_i_j_1 = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)]);
            visvx1 = (evsx1)*(diu2_im1_jp1_1+diu4_i_j_1);
            diu5_i_jp1_1 = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j+1)])*1.0/dy1[F1D2C(0 , j+1)];
            visvy2 = (evsy2)*2.*diu5_i_jp1_1;
            diu5_i_j_1 = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)])*1.0/dy1[F1D2C(0 , j)];
            visvy1 = (evsy1)*2.*diu5_i_j_1;
            diu6_i_j_2 = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)])*1.0/;
            diu8_i_jp1_1 = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j+1)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            visvz2 = (evsz2)*(diu6_i_j_2+diu8_i_jp1_1);
            visvz1 = pow((float)((0.4*vspd[F2D2C((((ip+1) - 0 )+1) , 0,0 , i,j)]*1.0/(float)log(0.5**1.0/0.1))))*(v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)]*1.0/vspd[F2D2C((((ip+1) - 0 )+1) , 0,0 , i,j)]);
            vfv = (visvx2-visvx1)*1.0/dx1[F1D2C((-1) , i)]+(visvy2-visvy1)*1.0/dys[F1D2C(0 , j)]+(visvz2-visvz1)*1.0/;
            g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j)] = (g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j)]+vfv);
      }
void les_map_319(__global float *sm,__global float *dzn,__global float *dx1,__global float *dy1,__global float *u,__global float *dzs,__global float *w,__global float *v,__global float *h) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: k,i,j,diu3_im1_j_kp1,diu3_i_j_kp1,diu6_i_jm1_kp1,diu6_i_j_kp1,diu7_i_j_k,diu7_ip1_j_k,diu8_i_j_k,diu8_i_jp1_k,diu9_i_j_k,evsx2,evsx1,evsy2,evsy1,evsz2,evsz1,viswx2,viswx1,viswy2,viswy1,diu9_i_j_kp1,viswz2,viswz1,vfw
    int k;
    int i;
    int j;
    float diu3_im1_j_kp1;
    float diu3_i_j_kp1;
    float diu6_i_jm1_kp1;
    float diu6_i_j_kp1;
    float diu7_i_j_k;
    float diu7_ip1_j_k;
    float diu8_i_j_k;
    float diu8_i_jp1_k;
    float diu9_i_j_k;
    float evsx2;
    float evsx1;
    float evsy2;
    float evsy1;
    float evsz2;
    float evsz1;
    float viswx2;
    float viswx1;
    float viswy2;
    float viswy1;
    float diu9_i_j_kp1;
    float viswz2;
    float viswz1;
    float vfw;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        k_rel = (global_id*1.0/(j_range*i_range));
        k = k_rel+1;
        j_rel = ((global_id-(k_rel*(j_range*i_range)))*1.0/i_range);
        j = j_rel+1;
        i_rel = ((global_id-(k_rel*(j_range*i_range)))-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
            evsz2 = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k+1)];
            evsz1 = sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)];
            evsx2 = (dzn[F1D2C((-1) , k+1)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]))+dzn[F1D2C((-1) , k)]*((dx1[F1D2C((-1) , i+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k+1)]+dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i+1,j,k+1)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))*1.0/(dzn[F1D2C((-1) , k)]+dzn[F1D2C((-1) , k+1)]);
            evsx1 = (dzn[F1D2C((-1) , k+1)]*((dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i-1,j,k)]+dx1[F1D2C((-1) , i-1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)]))+dzn[F1D2C((-1) , k)]*((dx1[F1D2C((-1) , i)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i-1,j,k+1)]+dx1[F1D2C((-1) , i-1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k+1)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)])))*1.0/(dzn[F1D2C((-1) , k)]+dzn[F1D2C((-1) , k+1)]);
            evsy2 = (dzn[F1D2C((-1) , k+1)]*((dy1[F1D2C(0 , j+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)]+dy1[F1D2C(0 , j)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]))+dzn[F1D2C((-1) , k)]*((dy1[F1D2C(0 , j+1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k+1)]+dy1[F1D2C(0 , j)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j+1,k+1)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)])))*1.0/(dzn[F1D2C((-1) , k)]+dzn[F1D2C((-1) , k+1)]);
            evsy1 = (dzn[F1D2C((-1) , k+1)]*((dy1[F1D2C(0 , j)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j-1,k)]+dy1[F1D2C(0 , j-1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k)])*1.0/(dy1[F1D2C(0 , j-1)]+dy1[F1D2C(0 , j)]))+dzn[F1D2C((-1) , k)]*((dy1[F1D2C(0 , j)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j-1,k+1)]+dy1[F1D2C(0 , j-1)]*sm[F3D2C((((ip+1) - (-1) )+1),(((jp+1) - (-1) )+1) , (-1),(-1),0 , i,j,k+1)])*1.0/(dy1[F1D2C(0 , j-1)]+dy1[F1D2C(0 , j)])))*1.0/(dzn[F1D2C((-1) , k)]+dzn[F1D2C((-1) , k+1)]);
            diu3_i_j_kp1 = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k+1)])*1.0/dzs[F1D2C((-1) , k)];
            diu7_ip1_j_k = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
            viswx2 = (evsx2)*(diu3_i_j_kp1+diu7_ip1_j_k);
            diu3_im1_j_kp1 = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k+1)])*1.0/dzs[F1D2C((-1) , k)];
            diu7_i_j_k = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i-1,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)]);
            viswx1 = (evsx1)*(diu3_im1_j_kp1+diu7_i_j_k);
            diu6_i_j_kp1 = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k+1)])*1.0/dzs[F1D2C((-1) , k)];
            diu8_i_jp1_k = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
            viswy2 = (evsy2)*(diu6_i_j_kp1+diu8_i_jp1_k);
            diu6_i_jm1_kp1 = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k+1)])*1.0/dzs[F1D2C((-1) , k)];
            diu8_i_j_k = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j-1,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)])*1.0/(dy1[F1D2C(0 , j-1)]+dy1[F1D2C(0 , j)]);
            viswy1 = (evsy1)*(diu6_i_jm1_kp1+diu8_i_j_k);
            diu9_i_j_kp1 = (-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k+1)])*1.0/dzn[F1D2C((-1) , k+1)];
            viswz2 = (evsz2)*2.*diu9_i_j_kp1;
            diu9_i_j_k = (-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k-1)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)])*1.0/dzn[F1D2C((-1) , k)];
            viswz1 = (evsz1)*2.*diu9_i_j_k;
            vfw = (viswx2-viswx1)*1.0/dx1[F1D2C((-1) , i)]+(viswy2-viswy1)*1.0/dy1[F1D2C(0 , j)]+(viswz2-viswz1)*1.0/dzn[F1D2C((-1) , k)];
            h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = (h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]+vfw);
      }
void press_map_64(__global float *f) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // local vars: j,k
    int j;
    int k;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int k_rel;
    int j_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        k_rel = global_id*1.0/j_range;
        k = k_rel+1;
        j_rel = (global_id-(k_rel*j_range));
        j = j_rel+1;
        // parallelfortran: original code
                f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , j,k)] = f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , j,k)];
        }
void press_map_69(__global float *g) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // local vars: k,i
    int k;
    int i;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int i_range;
    int k_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        i_range = ((300-1)+1);
        k_rel = global_id*1.0/i_range;
        k = k_rel+1;
        i_rel = (global_id-(k_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,k)] = g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,jp,k)];
        }
void press_map_74(__global float *h) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // local vars: j,i
    int j;
    int i;
        // parallelfortran: synthesised loop variable decls
    int j_range;
    int i_range;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        j_rel = global_id*1.0/i_range;
        j = j_rel+1;
        i_rel = (global_id-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j)] = 0.0;
                h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,kp)] = 0.0;
        }
void press_map_80(__global float *u,__global float *dx1,__global float *v,__global float *dy1,__global float *w,__global float *dzn,__global float *f,__global float *g,__global float *h,__global float *rhs,__global float *dt) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // local vars: j,k,i
    int j;
    int k;
    int i;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        k_rel = (global_id*1.0/(j_range*i_range));
        k = k_rel+1;
        j_rel = ((global_id-(k_rel*(j_range*i_range)))*1.0/i_range);
        j = j_rel+1;
        i_rel = ((global_id-(k_rel*(j_range*i_range)))-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                                rhs[F3D2C((((ip+1) - 0 )+1),(((jp+1) - 0 )+1) , 0,0,0 , i,j,k)] = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dx1[F1D2C((-1) , i)]+(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dy1[F1D2C(0 , j)]+(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k-1)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)])*1.0/dzn[F1D2C((-1) , k)];
                                rhs[F3D2C((((ip+1) - 0 )+1),(((jp+1) - 0 )+1) , 0,0,0 , i,j,k)] = (f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]-f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i-1,j,k)])*1.0/dx1[F1D2C((-1) , i)]+(g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]-g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j-1,k)])*1.0/dy1[F1D2C(0 , j)]+(h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]-h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k-1)])*1.0/dzn[F1D2C((-1) , k)]+rhs[F3D2C((((ip+1) - 0 )+1),(((jp+1) - 0 )+1) , 0,0,0 , i,j,k)]*1.0/dt;
                }
void press_reduce_93(__global float *dx1,__global float *dy1,__global float *dzn,__global float *rhs,__global float *global_rhsav_array,__global float *global_area_array) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // missing args: 
        // local vars: j,k,i
    int j;
    int k;
    int i;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    int chunk_size;
    int local_id;
    int local_id_fortran;
    int group_id;
    int group_id_fortran;
    int global_id;
    int r_iter;
    int local_chunk_size;
    int start_position;
    #if NTH > 1
        // arrays prefixed with __PH0__ should be declared using the __PH1__ modifier in c kernel version
    __local float local_rhsav_array[(NTH - 1 +1)];
    __local float local_area_array[(NTH - 1 +1)];
    #endif
    float local_rhsav;
    float clk_local_mem_fence;
    float local_area;
        group_id = get_group_id(0);
        global_id = get_global_id(0);
    #if NTH > 1
        local_id = get_local_id(0);
        // local_id_fortran and group_id_fortran are used to reconcile the fact that fortran arrays are referenced from 1
        // not 0 like other opencl supporting languages
        local_id_fortran = local_id+1;
    #endif
        group_id_fortran = group_id+1;
    #if NTH > 1
        local_chunk_size = (((80*(300*300))*1.0/NTH)*1.0/NUNITS);
    #else
        local_chunk_size = ((80*(300*300))*1.0/NUNITS);
    #endif
        start_position = local_chunk_size*global_id;
        local_rhsav = 0;
        local_area = 0;
    for (r_iter = start_position;r_iter <= ((start_position + local_chunk_size) - 1);r_iter += 1) {
                k_range = ((80-1)+1);
                j_range = ((300-1)+1);
                i_range = ((300-1)+1);
                k_rel = (r_iter*1.0/(j_range*i_range));
                k = k_rel+1;
                j_rel = ((r_iter-(k_rel*(j_range*i_range)))*1.0/i_range);
                j = j_rel+1;
                i_rel = ((r_iter-(k_rel*(j_range*i_range)))-(j_rel*i_range));
                i = i_rel+1;
                local_rhsav = (local_rhsav+(((dx1[F1D2C((-1) , i)]*dy1[F1D2C(0 , j)])*dzn[F1D2C((-1) , k)])*rhs[F3D2C((((ip+1) - 0 )+1),(((jp+1) - 0 )+1) , 0,0,0 , i,j,k)]));
                local_area = (local_area+((dx1[F1D2C((-1) , i)]*dy1[F1D2C(0 , j)])*dzn[F1D2C((-1) , k)]));
    }
    #if NTH > 1
        local_rhsav_array[F1D2C(1 , local_id_fortran)] = local_rhsav;
        local_area_array[F1D2C(1 , local_id_fortran)] = local_area;
    #ifdef BARRIER_OK
    barrier(CLK_LOCAL_MEM_FENCE);
    #endif
        local_rhsav = 0;
        local_area = 0;
    for (r_iter = 1;r_iter <= NTH;r_iter += 1) {
                local_rhsav = (local_rhsav+local_rhsav_array[F1D2C(1 , r_iter)]);
                local_area = (local_area+local_area_array[F1D2C(1 , r_iter)]);
    }
    #endif
        global_rhsav_array[F1D2C(1 , group_id_fortran)] = local_rhsav;
        global_area_array[F1D2C(1 , group_id_fortran)] = local_area;
    }
void press_map_102(__global float *rhs,__global float *rhsav) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // local vars: j,k,i
    int j;
    int k;
    int i;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        k_rel = (global_id*1.0/(j_range*i_range));
        k = k_rel+1;
        j_rel = ((global_id-(k_rel*(j_range*i_range)))*1.0/i_range);
        j = j_rel+1;
        i_rel = ((global_id-(k_rel*(j_range*i_range)))-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                                rhs[F3D2C((((ip+1) - 0 )+1),(((jp+1) - 0 )+1) , 0,0,0 , i,j,k)] = rhs[F3D2C((((ip+1) - 0 )+1),(((jp+1) - 0 )+1) , 0,0,0 , i,j,k)]-rhsav;
                }
void press_map_112(__global float *dzs,__global float *dys,__global float *dxs,__global int *nrd,__global float *p,__global float *rhs) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // local vars: j,k,i,dz1,dz2,cn4s,cn4l,cn3s,cn3l,cn2s,cn2l,cn1,reltmp
    int j;
    int k;
    int i;
    float dz1;
    float dz2;
    float cn4s;
    float cn4l;
    float cn3s;
    float cn3l;
    float cn2s;
    float cn2l;
    float cn1;
    float reltmp;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        k_rel = (global_id*1.0/(j_range*i_range));
        k = k_rel+1;
        j_rel = ((global_id-(k_rel*(j_range*i_range)))*1.0/i_range);
        j = j_rel+1;
        i_rel = ((global_id-(k_rel*(j_range*i_range)))-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                                dz1 = dzs[F1D2C((-1) , k-1)];
                                dz2 = dzs[F1D2C((-1) , k)];
                                cn4s = 2.*1.0/(dz1*(dz1+dz2));
                                cn4l = 2.*1.0/(dz2*(dz1+dz2));
                                        cn3s = 2.*1.0/(dys[F1D2C(0 , j-1)]*(dys[F1D2C(0 , j-1)]+dys[F1D2C(0 , j)]));
                                        cn3l = 2.*1.0/(dys[F1D2C(0 , j)]*(dys[F1D2C(0 , j-1)]+dys[F1D2C(0 , j)]));
                                                cn2s = 2.*1.0/(dxs[F1D2C(0 , i-1)]*(dxs[F1D2C(0 , i-1)]+dxs[F1D2C(0 , i)]));
                                                cn2l = 2.*1.0/(dxs[F1D2C(0 , i)]*(dxs[F1D2C(0 , i-1)]+dxs[F1D2C(0 , i)]));
                                                cn1 = 1.*1.0/(2.*1.0/(dxs[F1D2C(0 , i-1)]*dxs[F1D2C(0 , i)])+2.*1.0/(dys[F1D2C(0 , j-1)]*dys[F1D2C(0 , j)])+2.*1.0/(dz1*dz2));
                      if (nrd==0) {
                                                reltmp = omega*(cn1*(cn2l*p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i+1,j,k)]+cn2s*p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i-1,j,k)]+cn3l*p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j+1,k)]+cn3s*p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j-1,k)]+cn4l*p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k+1)]+cn4s*p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k-1)]-rhs[F3D2C((((ip+1) - 0 )+1),(((jp+1) - 0 )+1) , 0,0,0 , i,j,k)])-p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k)]);
                                                p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k)]+reltmp;
                       } else {
                                                reltmp = omega*(cn1*(cn2l*p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i+1,j,k)]+cn2s*p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i-1,j,k)]+cn3l*p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j+1,k)]+cn3s*p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j-1,k)]+cn4l*p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k+1)]+cn4s*p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k-1)]-rhs[F3D2C((((ip+1) - 0 )+1),(((jp+1) - 0 )+1) , 0,0,0 , i,j,k)])-p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k)]);
                                                p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k)]+reltmp;
                      }
                      }
void press_map_140(__global float *p) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // local vars: j,k
    int j;
    int k;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int k_rel;
    int j_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = (((80+1)-0)+1);
        j_range = (((300+1)-0)+1);
        k_rel = global_id*1.0/j_range;
        k = k_rel+0;
        j_rel = (global_id-(k_rel*j_range));
        j = j_rel+0;
        // parallelfortran: original code
                        p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , j,k)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , j,k)];
                        p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , ip+1,j,k)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , ip,j,k)];
            }
void press_map_146(__global float *p) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // local vars: k,i
    int k;
    int i;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int i_range;
    int k_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = (((80+1)-0)+1);
        i_range = (((300+1)-0)+1);
        k_rel = global_id*1.0/i_range;
        k = k_rel+0;
        i_rel = (global_id-(k_rel*i_range));
        i = i_rel+0;
        // parallelfortran: original code
                        p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,k)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,jp,k)];
                        p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,jp+1,k)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,k)];
            }
void press_map_153(__global float *p) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // local vars: j,i
    int j;
    int i;
        // parallelfortran: synthesised loop variable decls
    int j_range;
    int i_range;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        j_range = (((300+1)-0)+1);
        i_range = (((300+1)-0)+1);
        j_rel = global_id*1.0/i_range;
        j = j_rel+0;
        i_rel = (global_id-(j_rel*i_range));
        i = i_rel+0;
        // parallelfortran: original code
                p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j)];
                p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,kp+1)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,kp)];
        }
void press_reduce_162(__global float *p,__global float *dx1,__global float *dy1,__global float *dzn,__global float *global_pav_array,__global float *global_pco_array) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // missing args: 
        // local vars: j,k,i
    int j;
    int k;
    int i;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    int chunk_size;
    int local_id;
    int local_id_fortran;
    int group_id;
    int group_id_fortran;
    int global_id;
    int r_iter;
    int local_chunk_size;
    int start_position;
    #if NTH > 1
        // arrays prefixed with __PH0__ should be declared using the __PH1__ modifier in c kernel version
    __local float local_pav_array[(NTH - 1 +1)];
    __local float local_pco_array[(NTH - 1 +1)];
    #endif
    float local_pav;
    float clk_local_mem_fence;
    float local_pco;
        group_id = get_group_id(0);
        global_id = get_global_id(0);
    #if NTH > 1
        local_id = get_local_id(0);
        // local_id_fortran and group_id_fortran are used to reconcile the fact that fortran arrays are referenced from 1
        // not 0 like other opencl supporting languages
        local_id_fortran = local_id+1;
    #endif
        group_id_fortran = group_id+1;
    #if NTH > 1
        local_chunk_size = (((80*(300*300))*1.0/NTH)*1.0/NUNITS);
    #else
        local_chunk_size = ((80*(300*300))*1.0/NUNITS);
    #endif
        start_position = local_chunk_size*global_id;
        local_pav = 0;
        local_pco = 0;
    for (r_iter = start_position;r_iter <= ((start_position + local_chunk_size) - 1);r_iter += 1) {
                k_range = ((80-1)+1);
                j_range = ((300-1)+1);
                i_range = ((300-1)+1);
                k_rel = (r_iter*1.0/(j_range*i_range));
                k = k_rel+1;
                j_rel = ((r_iter-(k_rel*(j_range*i_range)))*1.0/i_range);
                j = j_rel+1;
                i_rel = ((r_iter-(k_rel*(j_range*i_range)))-(j_rel*i_range));
                i = i_rel+1;
                local_pav = (local_pav+(((p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k)]*dx1[F1D2C((-1) , i)])*dy1[F1D2C(0 , j)])*dzn[F1D2C((-1) , k)]));
                local_pco = (local_pco+((dx1[F1D2C((-1) , i)]*dy1[F1D2C(0 , j)])*dzn[F1D2C((-1) , k)]));
    }
    #if NTH > 1
        local_pav_array[F1D2C(1 , local_id_fortran)] = local_pav;
        local_pco_array[F1D2C(1 , local_id_fortran)] = local_pco;
    #ifdef BARRIER_OK
    barrier(CLK_LOCAL_MEM_FENCE);
    #endif
        local_pav = 0;
        local_pco = 0;
    for (r_iter = 1;r_iter <= NTH;r_iter += 1) {
                local_pav = (local_pav+local_pav_array[F1D2C(1 , r_iter)]);
                local_pco = (local_pco+local_pco_array[F1D2C(1 , r_iter)]);
    }
    #endif
        global_pav_array[F1D2C(1 , group_id_fortran)] = local_pav;
        global_pco_array[F1D2C(1 , group_id_fortran)] = local_pco;
    }
void press_map_171(__global float *p,__global float *pav) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // local vars: j,k,i
    int j;
    int k;
    int i;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        k_rel = (global_id*1.0/(j_range*i_range));
        k = k_rel+1;
        j_rel = ((global_id-(k_rel*(j_range*i_range)))*1.0/i_range);
        j = j_rel+1;
        i_rel = ((global_id-(k_rel*(j_range*i_range)))-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                                p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k)]-pav;
                }
void press_map_178(__global float *p) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // local vars: j,k
    int j;
    int k;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int k_rel;
    int j_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = (((80+1)-0)+1);
        j_range = (((300+1)-0)+1);
        k_rel = global_id*1.0/j_range;
        k = k_rel+0;
        j_rel = (global_id-(k_rel*j_range));
        j = j_rel+0;
        // parallelfortran: original code
                        p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , j,k)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , j,k)];
                        p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , ip+1,j,k)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , ip,j,k)];
            }
void press_map_184(__global float *p) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // local vars: k,i
    int k;
    int i;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int i_range;
    int k_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = (((80+1)-0)+1);
        i_range = (((300+1)-0)+1);
        k_rel = global_id*1.0/i_range;
        k = k_rel+0;
        i_rel = (global_id-(k_rel*i_range));
        i = i_rel+0;
        // parallelfortran: original code
                        p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,k)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,jp,k)];
                        p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,jp+1,k)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,k)];
            }
void press_map_190(__global float *p) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        const float pjuge = 0.0001;
        const int nmaxp = 50;
        const float omega = 1.;
        // local vars: j,i
    int j;
    int i;
        // parallelfortran: synthesised loop variable decls
    int j_range;
    int i_range;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        j_range = (((300+1)-0)+1);
        i_range = (((300+1)-0)+1);
        j_rel = global_id*1.0/i_range;
        j = j_rel+0;
        i_rel = (global_id-(j_rel*i_range));
        i = i_rel+0;
        // parallelfortran: original code
                p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j)];
                p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,kp+1)] = p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,kp)];
        }
void velfg_map_97(__global float *u,__global float *v,__global float *dx1,__global float *dy1,__global float *uspd,__global float *vspd) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
            const int u0 = 0;
        // local vars: i,j
    int i;
    int j;
        // parallelfortran: synthesised loop variable decls
    int j_range;
    int i_range;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        j_rel = global_id*1.0/i_range;
        j = j_rel+1;
        i_rel = (global_id-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                  uspd[F2D2C((((ip+1) - 0 )+1) , 0,0 , i,j)] = pow((float)((pow((float)(u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)]))+pow((float)(((0.5*(v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)])*dx1[F1D2C((-1) , i+1)]+0.5*(v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j-1)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j)])*dx1[F1D2C((-1) , i)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)])))))));
                  vspd[F2D2C((((ip+1) - 0 )+1) , 0,0 , i,j)] = pow((float)((pow((float)(v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)]))+pow((float)(((0.5*(u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j)])*dy1[F1D2C(0 , j+1)]+0.5*(u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j+1)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j+1)])*dy1[F1D2C(0 , j)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)])))))));
         }
void velfg_map_110(__global float *u,__global float *dx1,__global float *v,__global float *dy1,__global float *w,__global float *dzs,__global float *dzn,__global float *f,__global float *g,__global float *h) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
            const int u0 = 0;
        // local vars: i,j,nou1_i,k,diu1_i,cov1_i,nou1_ip1,diu1_ip1,cov1_ip1,nou2_j,diu2_j,cov2_j,nou2_jp1,diu2_jp1,cov2_jp1,nou3_k,diu3_k,cov3_k,nou3_kp1,diu3_kp1,cov3_kp1,covx1,covy1,covz1,covc,nou4_i,diu4_i,cov4_i,nou4_ip1,diu4_ip1,cov4_ip1,nou5_j,diu5_j,cov5_j,nou5_jp1,diu5_jp1,cov5_jp1,nou6_k,diu6_k,cov6_k,nou6_kp1,diu6_kp1,cov6_kp1,nou7_i,diu7_i,cov7_i,nou7_ip1,diu7_ip1,cov7_ip1,nou8_j,diu8_j,cov8_j,nou8_jp1,diu8_jp1,cov8_jp1,nou9_k,diu9_k,cov9_k,nou9_kp1,diu9_kp1,cov9_kp1
    int i;
    int j;
    float nou1_i;
    int k;
    float diu1_i;
    float cov1_i;
    float nou1_ip1;
    float diu1_ip1;
    float cov1_ip1;
    float nou2_j;
    float diu2_j;
    float cov2_j;
    float nou2_jp1;
    float diu2_jp1;
    float cov2_jp1;
    float nou3_k;
    float diu3_k;
    float cov3_k;
    float nou3_kp1;
    float diu3_kp1;
    float cov3_kp1;
    float covx1;
    float covy1;
    float covz1;
    float covc;
    float nou4_i;
    float diu4_i;
    float cov4_i;
    float nou4_ip1;
    float diu4_ip1;
    float cov4_ip1;
    float nou5_j;
    float diu5_j;
    float cov5_j;
    float nou5_jp1;
    float diu5_jp1;
    float cov5_jp1;
    float nou6_k;
    float diu6_k;
    float cov6_k;
    float nou6_kp1;
    float diu6_kp1;
    float cov6_kp1;
    float nou7_i;
    float diu7_i;
    float cov7_i;
    float nou7_ip1;
    float diu7_ip1;
    float cov7_ip1;
    float nou8_j;
    float diu8_j;
    float cov8_j;
    float nou8_jp1;
    float diu8_jp1;
    float cov8_jp1;
    float nou9_k;
    float diu9_k;
    float cov9_k;
    float nou9_kp1;
    float diu9_kp1;
    float cov9_kp1;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        k_rel = (global_id*1.0/(j_range*i_range));
        k = k_rel+1;
        j_rel = ((global_id-(k_rel*(j_range*i_range)))*1.0/i_range);
        j = j_rel+1;
        i_rel = ((global_id-(k_rel*(j_range*i_range)))-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                nou1_i = (u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/2.;
                diu1_i = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dx1[F1D2C((-1) , i)];
                cov1_i = nou1_i*diu1_i;
                nou1_ip1 = (u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j,k)])*1.0/2.;
                diu1_ip1 = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j,k)])*1.0/dx1[F1D2C((-1) , i+1)];
                cov1_ip1 = nou1_ip1*diu1_ip1;
                nou2_j = (dx1[F1D2C((-1) , i+1)]*v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+dx1[F1D2C((-1) , i)]*v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j-1,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
                diu2_j = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/(dy1[F1D2C(0 , j-1)]+dy1[F1D2C(0 , j)]);
                cov2_j = nou2_j*diu2_j;
                nou2_jp1 = (dx1[F1D2C((-1) , i+1)]*v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+dx1[F1D2C((-1) , i)]*v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
                diu2_jp1 = 2.*(-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
                cov2_jp1 = nou2_jp1*diu2_jp1;
                nou3_k = (dx1[F1D2C((-1) , i+1)]*w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k-1)]+dx1[F1D2C((-1) , i)]*w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i+1,j,k-1)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
                diu3_k = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k-1)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dzs[F1D2C((-1) , k-1)];
                cov3_k = nou3_k*diu3_k;
                nou3_kp1 = (dx1[F1D2C((-1) , i+1)]*w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+dx1[F1D2C((-1) , i)]*w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
                diu3_kp1 = (-u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k+1)])*1.0/dzs[F1D2C((-1) , k)];
                cov3_kp1 = nou3_kp1*diu3_kp1;
                covx1 = (dx1[F1D2C((-1) , i+1)]*cov1_i+dx1[F1D2C((-1) , i)]*cov1_ip1)*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
                covy1 = (cov2_j+cov2_jp1)*1.0/2.;
                covz1 = (cov3_k+cov3_kp1)*1.0/2.;
                covc = covx1+covy1+covz1;
                f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = -covc;
                nou4_i = (dy1[F1D2C(0 , j+1)]*u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+dy1[F1D2C(0 , j)]*u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
                diu4_i = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)]);
                cov4_i = (nou4_i-u0)*diu4_i;
                nou4_ip1 = (dy1[F1D2C(0 , j+1)]*u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+dy1[F1D2C(0 , j)]*u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
                diu4_ip1 = 2.*(-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
                cov4_ip1 = (nou4_ip1-u0)*diu4_ip1;
                nou5_j = (v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/2.;
                diu5_j = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dy1[F1D2C(0 , j)];
                cov5_j = nou5_j*diu5_j;
                nou5_jp1 = (v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j+1,k)])*1.0/2.;
                diu5_jp1 = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j+1,k)])*1.0/dy1[F1D2C(0 , j+1)];
                cov5_jp1 = nou5_jp1*diu5_jp1;
                nou6_k = (dy1[F1D2C(0 , j+1)]*w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k-1)]+dy1[F1D2C(0 , j)]*w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j+1,k-1)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
                diu6_k = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k-1)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)])*1.0/dzs[F1D2C((-1) , k-1)];
                cov6_k = nou6_k*diu6_k;
                nou6_kp1 = (dy1[F1D2C(0 , j+1)]*w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+dy1[F1D2C(0 , j)]*w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
                diu6_kp1 = (-v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k+1)])*1.0/dzs[F1D2C((-1) , k)];
                cov6_kp1 = nou6_kp1*diu6_kp1;
                covx1 = (cov4_i+cov4_ip1)*1.0/2.;
                covy1 = (dy1[F1D2C(0 , j+1)]*cov5_j+dy1[F1D2C(0 , j)]*cov5_jp1)*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
                covz1 = (cov6_k+cov6_kp1)*1.0/2.;
                covc = covx1+covy1+covz1;
                g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = -covc;
    if ((k<80-1)) {
                nou7_i = (dzn[F1D2C((-1) , k+1)]*u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k)]+dzn[F1D2C((-1) , k)]*u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i-1,j,k+1)])*1.0/(dzn[F1D2C((-1) , k)]+dzn[F1D2C((-1) , k+1)]);
                diu7_i = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i-1,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)])*1.0/(dx1[F1D2C((-1) , i-1)]+dx1[F1D2C((-1) , i)]);
                cov7_i = (nou7_i-u0)*diu7_i;
                nou7_ip1 = (dzn[F1D2C((-1) , k+1)]*u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+dzn[F1D2C((-1) , k)]*u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k+1)])*1.0/(dzn[F1D2C((-1) , k)]+dzn[F1D2C((-1) , k+1)]);
                diu7_ip1 = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i+1,j,k)])*1.0/(dx1[F1D2C((-1) , i)]+dx1[F1D2C((-1) , i+1)]);
                cov7_ip1 = (nou7_ip1-u0)*diu7_ip1;
                nou8_j = (dzn[F1D2C((-1) , k+1)]*v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k)]+dzn[F1D2C((-1) , k)]*v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j-1,k+1)])*1.0/(dzn[F1D2C((-1) , k)]+dzn[F1D2C((-1) , k+1)]);
                diu8_j = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j-1,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)])*1.0/(dy1[F1D2C(0 , j-1)]+dy1[F1D2C(0 , j)]);
                cov8_j = nou8_j*diu8_j;
                nou8_jp1 = (dzn[F1D2C((-1) , k+1)]*v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+dzn[F1D2C((-1) , k)]*v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k+1)])*1.0/(dzn[F1D2C((-1) , k)]+dzn[F1D2C((-1) , k+1)]);
                diu8_jp1 = 2.*(-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j+1,k)])*1.0/(dy1[F1D2C(0 , j)]+dy1[F1D2C(0 , j+1)]);
                cov8_jp1 = nou8_jp1*diu8_jp1;
                nou9_k = (w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k-1)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)])*1.0/2.;
                diu9_k = (-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k-1)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)])*1.0/dzn[F1D2C((-1) , k)];
                cov9_k = nou9_k*diu9_k;
                nou9_kp1 = (w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k+1)])*1.0/2.;
                diu9_kp1 = (-w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k+1)])*1.0/dzn[F1D2C((-1) , k+1)];
                cov9_kp1 = nou9_kp1*diu9_kp1;
              covx1 = (cov7_i+cov7_ip1)*1.0/2.;
              covy1 = (cov8_j+cov8_jp1)*1.0/2.;
              covz1 = (dzn[F1D2C((-1) , k+1)]*cov9_k+dzn[F1D2C((-1) , k)]*cov9_kp1)*1.0/(dzn[F1D2C((-1) , k)]+dzn[F1D2C((-1) , k+1)]);
              covc = covx1+covy1+covz1;
                h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)] = -covc;
    }
    }
void velfg_map_197(__global float *f) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
            const int u0 = 0;
        // local vars: j,k
    int j;
    int k;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int k_rel;
    int j_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        k_rel = global_id*1.0/j_range;
        k = k_rel+1;
        j_rel = (global_id-(k_rel*j_range));
        j = j_rel+1;
        // parallelfortran: original code
                                f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , j,k)] = f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , j,k)];
                }
void velfg_map_202(__global float *g) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
            const int u0 = 0;
        // local vars: i,k
    int i;
    int k;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int i_range;
    int k_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        i_range = ((300-1)+1);
        k_rel = global_id*1.0/i_range;
        k = k_rel+1;
        i_rel = (global_id-(k_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                        g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,k)] = g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,jp,k)];
            }
void velfg_map_207(__global float *h) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
            const int u0 = 0;
        // local vars: i,j
    int i;
    int j;
        // parallelfortran: synthesised loop variable decls
    int j_range;
    int i_range;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        j_rel = global_id*1.0/i_range;
        j = j_rel+1;
        i_rel = (global_id-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                        h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j)] = 0.0;
                        h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,kp)] = 0.0;
            }
void velnw_map_36(__global float *p,__global float *ro,__global float *dxs,__global float *u,__global float *dt,__global float *f,__global float *dys,__global float *v,__global float *g,__global float *dzs,__global float *w,__global float *h) {

    const int kp = 80;
        const int ip = 300;
        const int jp = 300;
        const int ipmax = ip;
        const int jpmax = jp;
        const float dxgrid = 4.;
        const float dygrid = 4.;
        const float cs0 = 0.14;
        const int i_anime = 1;
        const int avetime = 2;
        const int km_sl = 80;
        const int i_aveflow = 0;
        const int i_ifdata_out = 0;
        const float dt_orig = 0.05;
        // local vars: pz,i,j,k
    float pz;
    int i;
    int j;
    int k;
        // parallelfortran: synthesised loop variable decls
    int k_range;
    int j_range;
    int i_range;
    int k_rel;
    int j_rel;
    int i_rel;
    // READ
    // WRITTEN
    // READ & WRITTEN
    // globalIdDeclaration
    int global_id;
    // globalIdInitialisation
        global_id = get_global_id(0);
    // ptrAssignments_fseq
        // parallelfortran: synthesised loop variables
        k_range = ((80-1)+1);
        j_range = ((300-1)+1);
        i_range = ((300-1)+1);
        k_rel = (global_id*1.0/(j_range*i_range));
        k = k_rel+1;
        j_rel = ((global_id-(k_rel*(j_range*i_range)))*1.0/i_range);
        j = j_rel+1;
        i_rel = ((global_id-(k_rel*(j_range*i_range)))-(j_rel*i_range));
        i = i_rel+1;
        // parallelfortran: original code
                pz = (-p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k)]+p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i+1,j,k)])*1.0/ro*1.0/dxs[F1D2C(0 , i)];
                u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)] = u[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+dt*(f[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]-pz);
                pz = (-p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k)]+p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j+1,k)])*1.0/ro*1.0/dys[F1D2C(0 , j)];
                v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)] = v[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),0 , i,j,k)]+dt*(g[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]-pz);
    if ((k<80-1)) {
                pz = (-p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k)]+p[F4D2C(((1 - 0 )+1),(((ip+2) - 0 )+1),(((jp+2) - 0 )+1) , 0,0,0,0 , i,j,k+1)])*1.0/ro*1.0/dzs[F1D2C((-1) , k)];
                w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)] = w[F3D2C((((ip+1) - 0 )+1),(((jp+1) - (-1) )+1) , 0,(-1),(-1) , i,j,k)]+dt*(h[F3D2C(((ip - 0 )+1),((jp - 0 )+1) , 0,0,0 , i,j,k)]-pz);
    }
    }
__kernel void adam_bondv1_feedbf_les_press_v_etc_superkernel(__global float *f,__global float *g,__global float *h,__global float *fold,__global float *gold,__global float *hold,__global float *z2,__global float *dzn,__global float *u,__global float *v,__global float *w,__global float *global_aaa_array,__global float *global_bbb_array,__global float *dt,__global float *uout,__global float *dxs,__global float *zbm,__global float *usum,__global float *vsum,__global float *wsum,__global float *alpha,__global float *beta,__global float *dx1,__global float *dy1,__global float *delx1,__global float *dzs,__global float *sm,__global float *uspd,__global float *dys,__global float *vspd,__global float *rhs,__global float *global_rhsav_array,__global float *global_area_array,__global float *rhsav,__global int *nrd,__global float *p,__global float *global_pav_array,__global float *global_pco_array,__global float *pav,__global float *ro,__global int *state_ptr) {

  int state;
  const int st_adam_map_36 = 0;
  const int st_bondv1_map_38 = 1;
  const int st_bondv1_map_48 = 2;
  const int st_bondv1_map_58 = 3;
  const int st_bondv1_reduce_69 = 4;
  const int st_bondv1_reduce_76 = 5;
  const int st_bondv1_map_83 = 6;
  const int st_bondv1_map_98 = 7;
  const int st_bondv1_map_116 = 8;
  const int st_bondv1_map_128 = 9;
  const int st_feedbf_map_48 = 10;
  const int st_les_map_121 = 11;
  const int st_les_map_124 = 12;
  const int st_les_map_169 = 13;
  const int st_les_map_175 = 14;
  const int st_les_map_181 = 15;
  const int st_les_map_187 = 16;
  const int st_les_map_225 = 17;
  const int st_les_map_253 = 18;
  const int st_les_map_291 = 19;
  const int st_les_map_319 = 20;
  const int st_press_map_64 = 21;
  const int st_press_map_69 = 22;
  const int st_press_map_74 = 23;
  const int st_press_map_80 = 24;
  const int st_press_reduce_93 = 25;
  const int st_press_map_102 = 26;
  const int st_press_map_112 = 27;
  const int st_press_map_140 = 28;
  const int st_press_map_146 = 29;
  const int st_press_map_153 = 30;
  const int st_press_reduce_162 = 31;
  const int st_press_map_171 = 32;
  const int st_press_map_178 = 33;
  const int st_press_map_184 = 34;
  const int st_press_map_190 = 35;
  const int st_velfg_map_97 = 36;
  const int st_velfg_map_110 = 37;
  const int st_velfg_map_197 = 38;
  const int st_velfg_map_202 = 39;
  const int st_velfg_map_207 = 40;
  const int st_velnw_map_36 = 41;
    state = *state_ptr;
  // SUPERKERNEL BODY
  switch ( state ) {
        case (st_adam_map_36): {
            adam_map_36(f,g,h,fold,gold,hold);
        } break;
        case (st_bondv1_map_38): {
            bondv1_map_38(z2,dzn,u,v,w);
        } break;
        case (st_bondv1_map_48): {
            bondv1_map_48(u,v,w);
        } break;
        case (st_bondv1_map_58): {
            bondv1_map_58(u,v,w);
        } break;
        case (st_bondv1_reduce_69): {
            bondv1_reduce_69(u,global_aaa_array);
        } break;
        case (st_bondv1_reduce_76): {
            bondv1_reduce_76(u,global_bbb_array);
        } break;
        case (st_bondv1_map_83): {
            bondv1_map_83(u,dt,uout,dxs,v,w);
        } break;
        case (st_bondv1_map_98): {
            bondv1_map_98(u,v,w);
        } break;
        case (st_bondv1_map_116): {
            bondv1_map_116(u,v);
        } break;
        case (st_bondv1_map_128): {
            bondv1_map_128(w);
        } break;
        case (st_feedbf_map_48): {
            feedbf_map_48(zbm,z2,dzn,usum,u,vsum,v,wsum,w,alpha,dt,beta,f,g,h);
        } break;
        case (st_les_map_121): {
            les_map_121(dx1,dy1,dzn,delx1);
        } break;
        case (st_les_map_124): {
            les_map_124(u,dx1,dy1,dzs,v,w,dzn,delx1,sm);
        } break;
        case (st_les_map_169): {
            les_map_169(sm);
        } break;
        case (st_les_map_175): {
            les_map_175(sm);
        } break;
        case (st_les_map_181): {
            les_map_181(sm);
        } break;
        case (st_les_map_187): {
            les_map_187(sm,dy1,dx1,dzn,u,v,dzs,w,dxs,f);
        } break;
        case (st_les_map_225): {
            les_map_225(sm,dy1,dx1,dzn,u,v,dzs,w,uspd,dxs,f);
        } break;
        case (st_les_map_253): {
            les_map_253(sm,dy1,dx1,dzn,u,v,dzs,w,dys,g);
        } break;
        case (st_les_map_291): {
            les_map_291(sm,dy1,dx1,dzn,u,v,dzs,w,vspd,dys,g);
        } break;
        case (st_les_map_319): {
            les_map_319(sm,dzn,dx1,dy1,u,dzs,w,v,h);
        } break;
        case (st_press_map_64): {
            press_map_64(f);
        } break;
        case (st_press_map_69): {
            press_map_69(g);
        } break;
        case (st_press_map_74): {
            press_map_74(h);
        } break;
        case (st_press_map_80): {
            press_map_80(u,dx1,v,dy1,w,dzn,f,g,h,rhs,dt);
        } break;
        case (st_press_reduce_93): {
            press_reduce_93(dx1,dy1,dzn,rhs,global_rhsav_array,global_area_array);
        } break;
        case (st_press_map_102): {
            press_map_102(rhs,rhsav);
        } break;
        case (st_press_map_112): {
            press_map_112(dzs,dys,dxs,nrd,p,rhs);
        } break;
        case (st_press_map_140): {
            press_map_140(p);
        } break;
        case (st_press_map_146): {
            press_map_146(p);
        } break;
        case (st_press_map_153): {
            press_map_153(p);
        } break;
        case (st_press_reduce_162): {
            press_reduce_162(p,dx1,dy1,dzn,global_pav_array,global_pco_array);
        } break;
        case (st_press_map_171): {
            press_map_171(p,pav);
        } break;
        case (st_press_map_178): {
            press_map_178(p);
        } break;
        case (st_press_map_184): {
            press_map_184(p);
        } break;
        case (st_press_map_190): {
            press_map_190(p);
        } break;
        case (st_velfg_map_97): {
            velfg_map_97(u,v,dx1,dy1,uspd,vspd);
        } break;
        case (st_velfg_map_110): {
            velfg_map_110(u,dx1,v,dy1,w,dzs,dzn,f,g,h);
        } break;
        case (st_velfg_map_197): {
            velfg_map_197(f);
        } break;
        case (st_velfg_map_202): {
            velfg_map_202(g);
        } break;
        case (st_velfg_map_207): {
            velfg_map_207(h);
        } break;
        case (st_velnw_map_36): {
            velnw_map_36(p,ro,dxs,u,dt,f,dys,v,g,dzs,w,h);
      }
  }
  }
  
