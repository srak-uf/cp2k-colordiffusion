/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2018  CP2K developers group                         *
 *****************************************************************************/

//  Authors: Andreas Gloess <andreas.gloess@chem.uzh.ch>

#if defined (__ACC)

// kernel (input) parameters

/*
// Kernel optimized for CUDA on K20x running on K20c, with __ldg() and cudaFuncSetSharedMemConfig()
// Kernel_dnt_largeDB2(m=23, n=23, k=23, tile_m=3, tile_n=2, w=8, v=12, threads=96, grouping=16, minblocks=12) , # 362.853 GFlop/s (K20x)
// Results LIBTEST(6440):  OpenCL = 176.8 GFlops, CUDA = 209.7 GFlops
#define m 23
#define n 23
#define k 23
#define M 3
#define N 2
#define w 8
#define v 12
#define blockdim 96
#define grouping 16
#define minblocks 12
*/

/*
// Kernel optimized for CUDA on K20c running on K20c, with __ldg() and cudaFuncSetSharedMemConfig()
// Kernel_dnt_largeDB2(m=23, n=23, k=23, tile_m=2, tile_n=3, w=8, v=20, threads=96, grouping=16, minblocks=12) , # 315.641 GFlop/s (K20c)
// Results LIBTEST(6440):  OpenCL = 162.8 GFlops, CUDA = 241.5 GFlops
#define m 23
#define n 23
#define k 23
#define M 2
#define N 3
#define w 8
#define v 20
#define blockdim 96
#define grouping 16
#define minblocks 12
*/

// Kernel optimized for CUDA on K20c running on K20c, w/o __ldg() and cudaFuncSetSharedMemConfig()
// Kernel_dnt_largeDB2(m=23, n=23, k=23, tile_m=3, tile_n=3, w=4, v=22, threads=96, grouping=16, minblocks=1), # 240.196 GFlop/s
// Results LIBTEST(6440): OpenCL = 187.8 GFlops, CUDA = 193.0 GFlops
#define m 23
#define n 23
#define k 23
#define M 3
#define N 3
#define w 4
#define v 22
#define blockdim 96
#define grouping 16
#define minblocks 1



// kernel (input) dependent parameters
#define mya_size ((w * m + blockdim - 1) / blockdim)
#define myb_size ((w * n + blockdim - 1) / blockdim)
#define buff_size  MAX(m * w + w * n, v * m)
#define wa (k - (k / w) * w)

//**************************************************************************//
__kernel __attribute__ ((reqd_work_group_size(blockdim, 1, 1)))
  void clsmm_dnt_largeDB2_16_23_23_12_23_96_3_2_12_8 (
                __global int    *param_stack,
                         int    careful,
                         int    nruns,
                __global double *a_data,
                __global double *b_data,
                __global double *c_data)
{
    // registers to store thread's result tile
    __private double myc[N * M];

    // registers to store input slabs during double buffering
    // If there are too few thread, each thread has to store
    // multiple elements of the input slabs in it's registers.
    __private double mya[mya_size];
    __private double myb[myb_size];

     // initialize the thread's result tile to zero
    for (int i = 0; i < N * M; i++)
        myc[i] = 0.0;

    // buffer needs to hold input and output slabs (not both simultaneously).
    __local double buff[buff_size];

    // conveniece pointers
    __local double *buff_l = buff;
    __local double *buff_r = &(buff[m * w]);

    // first stack entry to be processed by this thread-block
    int psp = 7 * (get_group_id(0) * grouping);

    // grouping is the number of stack entries process by each thread-block
    // careful is the number of launched thread-blocks.
    // nruns is the number of stack entries process by the last thread-block
    int nrun = (get_group_id(0) == careful) ? nruns : grouping;

    // all stack entries relavant for this thread-block are loaded at once
    // allows to look ahead and and flush result tile only when really needed
    __local int param_stack_s[4 * grouping];

    // load parameter stack, might read beyond
    for (int i = get_local_id(0); i < 7 * nrun; i += blockdim) {
        //int p_tmp = __ldg(param_stack + psp + i);
        int p_tmp = *(param_stack + psp + i);
        if (i % 7 > 2)
            param_stack_s[(i / 7) * 4 + i % 7 - 3] = p_tmp - 1;
    }

    psp = 0;

    barrier(CLK_LOCAL_MEM_FENCE);

   // get the offsets for the a-block and the b-block from the stack
    int srcA = param_stack_s[psp];
    int srcB = param_stack_s[psp + 1];
    // start off double buffering by loading the first data
    load_gmem_into_regs(a_data + srcA, mya, m * w, blockdim);
    load_gmem_into_regs(b_data + srcB, myb, n * w, blockdim);

    barrier(CLK_LOCAL_MEM_FENCE);
    // in each run we process one stack entry
    for (int run = 0; run < nrun; run++) {

        // load the first slab for multiplication into the smem
        load_regs_into_smem(mya, buff_l, m * w, blockdim);
        load_regs_into_smem(myb, buff_r, n * w, blockdim);
        barrier(CLK_LOCAL_MEM_FENCE);

        // this is the actual double buffering loop
        for (int t = 0; t < (k / w -1) * w ; t += w) {
            // load next input slab from global memory into registers
            srcA += m * w;
            srcB += n * w;
            load_gmem_into_regs(a_data + srcA, mya, m * w, blockdim);
            load_gmem_into_regs(b_data + srcB, myb, n * w, blockdim);
            // multiply previous slab, which is stored in shared memory,
            // and accumulate the results in the registers myc
            multiply(buff_l, buff_r, myc, w, m, n, M, N, blockdim);
            barrier(CLK_LOCAL_MEM_FENCE);
            // copy next slab from registers to shared memory
            load_regs_into_smem(mya, buff_l, m * w, blockdim);
            load_regs_into_smem(myb, buff_r, n * w, blockdim);
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        if (wa != 0) { // is there a tail-slab?
            // If the input slab witdh w is not a divisor of k,
            // a smaller tail-slab of width wa has to be process
            // load tail-slab into registers
            srcA += m * w;
            srcB += n * w;
            load_gmem_into_regs(a_data + srcA, mya, m * wa, blockdim);
            load_gmem_into_regs(b_data + srcB, myb, n * wa, blockdim);
            // multiply last regular slab, which the loop left in shared memory
            multiply(buff_l, buff_r, myc, w, m, n, M, N, blockdim);
            barrier(CLK_LOCAL_MEM_FENCE);
            // copy tail-slab from register into shared mem
            load_regs_into_smem(mya, buff_l, m * wa, blockdim);
            load_regs_into_smem(myb, buff_r, n * wa, blockdim);
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        psp = 4 * run + 4;

        if(run < nrun-1){
            // get the offsets for the a-block and the b-block from the stack
            srcA = param_stack_s[psp];
            srcB = param_stack_s[psp + 1];
            // load the data for the next iteration of the loop
            load_gmem_into_regs(a_data + srcA, mya, m * w, blockdim);
            load_gmem_into_regs(b_data + srcB, myb, n * w, blockdim);
        }

        if (wa != 0) { // is there a tail-slab?
            // multiply the tail-slab
            multiply(buff_l, buff_r, myc, wa, m, n, M, N, blockdim);
        }else{
            // multiply last regular slab, which the loop left in shared memory
            multiply(buff_l, buff_r, myc, w, m, n, M, N, blockdim);
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        // multiplication for this run done
        // do we have to flush the result tile?
        if(run == nrun-1 || (param_stack_s[psp - 1] != param_stack_s[psp + 3])) {
            int c_loc = param_stack_s[psp - 2];
            writeback_results(myc, c_data + c_loc, buff, m, n, M, N, v, blockdim);
        }
    }
}

#endif
//EOF
