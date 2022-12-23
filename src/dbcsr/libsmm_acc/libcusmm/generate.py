#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from os import path
from glob import glob
from itertools import chain
from optparse import OptionParser

from kernels.cusmm_dnt_largeDB1 import Kernel_dnt_largeDB1
from kernels.cusmm_dnt_largeDB2 import Kernel_dnt_largeDB2
from kernels.cusmm_dnt_medium   import Kernel_dnt_medium
from kernels.cusmm_dnt_small    import Kernel_dnt_small
from kernels.cusmm_dnt_tiny     import Kernel_dnt_tiny


#===============================================================================
def main(argv):
    triples  = combinations(23)                 # blocked H2O (benchmark)
    triples += combinations(6)                  # idem min basis
    triples += combinations(14,16,29)           # RPA water
    triples += combinations(14,32,29)
    triples += combinations(5, 32, 13, 24, 26)
    triples += combinations(9, 32, 22)
    triples += combinations(64)
    triples += combinations(78)
    triples += combinations(16,29,55)
    triples += combinations(32,29,55)
    triples += combinations(12)
    triples += combinations(4,5,7,9,13,25,26,28,32,45)
    triples += combinations(4,10,15)
    #triples += combinations(5,12,13,26,32)
    #triples += combinations(12,25)
    #triples += combinations(9,12)
    #triples += combinations(5,13,16,32)
    triples += combinations(6,7,8)
    triples += combinations(13,14,25,26,32)

    usage = "Generator of LibCuSMM. The Library for Cuda Small Matrix Multiplications."
    parser = OptionParser(usage)
    parser.add_option("-p", "--params", metavar="filename.txt",
        default="parameters_P100.txt",
        help="Default: %default")

    (options, args) = parser.parse_args(argv)
    assert(len(args)==0)
    param_fn = options.params

    plan = make_plan(triples, param_fn)
    gen_library(plan)

    print("Libcusmm: Generated %d kernels."%len(plan))


#===============================================================================
def make_plan(triples, param_fn):
    all_kernels = eval(open(param_fn).read())
    print("Libcusmm: Found %d parameter sets."%len(all_kernels))

    triples = list(set(triples))

    plan = {}
    for (m, n, k) in triples:
        possible_kernels = [kern for kern in all_kernels if kern.can_handle(m,n,k)]
        if(len(possible_kernels) == 1):
            plan[(m,n,k)] = possible_kernels[0]
        elif(len(possible_kernels) > 1):
            raise(Exception("found more than one kernel for %dx%dx%d"%(m,n,k)))
        else:
            raise(Exception("missing kernel parameters for %dx%dx%d"%(m,n,k)))

    return(plan)



#===============================================================================
def gen_library(plan):
    header  = "/*****************************************************************************\n"
    header += " *  CP2K: A general program to perform molecular dynamics simulations        *\n"
    header += " *  Copyright (C) 2000 - 2018  CP2K developers group                         *\n"
    header += " *****************************************************************************/\n"

    for i in get_includes(plan):
        header += '#include "%s"\n'%i
    header += "\n\n"

    all_kernels = [plan[k] for k in sorted(plan.keys())]

    # write part files
    LAUNCHERS_PER_OBJ = 100
    for i in range(len(all_kernels)//LAUNCHERS_PER_OBJ + 1):
        a = i*LAUNCHERS_PER_OBJ
        b = min((i+1)*LAUNCHERS_PER_OBJ, len(all_kernels))
        output = header
        for j in range(a, b):
            output += all_kernels[j].launcher_code() + "\n\n"
        output += "//EOF\n"
        writefile("libcusmm_part%.2d.cu"%i, output)

    # write main file
    output = header
    for kern in all_kernels:
        l = "launch_"+kern.name
        output += "int "+ l + "(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data);\n"
    output += "\n\n"

    output += gen_process(plan)
    output += "\n\n"

    left_sizes = [(k,m) for (m,n,k) in plan.keys()]
    output += gen_transpose(left_sizes)
    output += "\n\n"

    output += gen_list(plan)
    output += "//EOF\n"
    writefile("libcusmm.cu", output)


#===============================================================================
def gen_process(plan):
    output  = "int libcusmm_process_d(int *param_stack, int stack_size,"
    output += "cudaStream_t stream, int m, int n, int k, "
    output += "double *a_data, double *b_data, double *c_data){\n"


    #generate jump table -------------------------------------------------------
    m_vals = sorted(list(set([m for (m,n,k) in plan.keys()])))
    n_vals = sorted(list(set([n for (m,n,k) in plan.keys()])))
    k_vals = sorted(list(set([k for (m,n,k) in plan.keys()])))
    assert(len(m_vals) * len(n_vals) * len(k_vals) < pow(2,16))

    output += "int idx = 0;\n"
    output += "bool missing = false;\n\n"
    output += "switch(m){\n"
    for i, m in enumerate(m_vals):
        output += "case %d: idx = %d; break;\n"%(m, i)
    output += "default: missing = true;\n"
    output += "}\n\n"

    output += "idx *= %d;\n"%len(n_vals)
    output += "switch(n){\n"
    for i, n in enumerate(n_vals):
        output += "case %d: idx += %d; break;\n"%(n, i)
    output += "default: missing = true;\n"
    output += "}\n\n"

    output += "idx *= %d;\n"%len(k_vals)
    output += "switch(k){\n"
    for i, k in enumerate(k_vals):
        output += "case %d: idx += %d; break;\n"%(k, i)
    output += "default: missing = true;\n"
    output += "}\n\n"

    output += "if(missing) return -1;\n"

    idx_map = dict()
    for (m,n,k) in plan.keys():
        idx = (m_vals.index(m)*len(n_vals) + n_vals.index(n))*len(k_vals) + k_vals.index(k)
        idx_map[idx] = (m,n,k)

    output += "switch(idx){\n"
    for idx in sorted(idx_map.keys()):
        mnk = idx_map[idx]
        output += "case %d:\n"%idx
        output += "// m=%d, n=%d, k=%d\n"%mnk
        output += "return launch_"+plan[mnk].name
        output += "(param_stack, stack_size, stream, %d, %d, %d, "%mnk
        output += "a_data, b_data, c_data);\n\n"
    output += "}\n\n"

    output += "return -1; // should never happen\n"
    output += "}\n\n\n"

    output += '#define dbcsr_type_real_4     1\n'
    output += '#define dbcsr_type_real_8     3\n'
    output += '#define dbcsr_type_complex_4  5\n'
    output += '#define dbcsr_type_complex_8  7\n\n'

    output += 'extern "C" int libsmm_acc_process '
    output += '(void *param_stack, int stack_size, int nparams, int datatype,'
    output += ' void *a_data, void *b_data, void *c_data,'
    output += ' int m_max, int n_max, int k_max, int def_mnk, void *stream){\n'
    output += 'cudaStream_t* custream = (cudaStream_t*) stream;\n'
    output += 'if(def_mnk!=1)\n'
    output += '  return(-1); // inhomogenous stacks not supported\n'
    output += 'if(datatype==dbcsr_type_real_8)\n'
    output += '  return(libcusmm_process_d ((int *) param_stack, stack_size, *custream,'
    output += ' m_max, n_max, k_max,'
    output += '(double *) a_data, (double *) b_data, (double *) c_data));\n\n'
    output += 'return(-1); // datatype not supported\n'
    output += '};\n'

    return(output)


#===============================================================================
def gen_list(plan):
    output  = "void libcusmm_list_blocksizes_d(const int **list, int *length) {\n"
    output += "static const int blocksizes_d[] = { \n"
    for mnk in sorted(plan.keys()):
        output += " %d, %d, %d,\n"%mnk
    output += "};\n\n"

    output += "*list = blocksizes_d;\n"
    output += "*length = %d;\n"%len(plan)
    output += "}\n"

    return(output)

#===============================================================================
def gen_transpose(sizes):
    output  = "int libcusmm_transpose_d(int *trs_stack, int offset, int nblks,\n"
    output += "double *buffer, int m, int n, cudaStream_t * stream) {\n"

    m_vals = sorted(list(set([m for (m,n) in sizes])))
    n_vals = sorted(list(set([n for (m,n) in sizes])))
    assert(len(m_vals) * len(n_vals) < pow(2,16))

    output += "int idx = 0;\n"
    output += "bool missing = false;\n\n"
    output += "switch(m){\n"
    for i, m in enumerate(m_vals):
        output += "case %d: idx = %d; break;\n"%(m, i)
    output += "default: missing = true;\n"
    output += "}\n\n"

    output += "idx *= %d;\n"%len(n_vals)
    output += "switch(n){\n"
    for i, n in enumerate(n_vals):
        output += "case %d: idx += %d; break;\n"%(n, i)
    output += "default: missing = true;\n"
    output += "}\n\n"

    output += "// If there is no kernel for these blocks, we don't need to transpose them.\n"
    output += "if(missing) return 0;\n\n"

    idx_map = dict()
    for (m,n) in sizes:
        idx = m_vals.index(m)*len(n_vals) + n_vals.index(n)
        idx_map[idx] = (m,n)

    output += "typedef void (*kernel)(int *, int , double*);\n"
    for idx in sorted(idx_map.keys()):
        mn = idx_map[idx]
        output += "// m=%d, n=%d\n"%mn
        output += "static kernel kern_func_%d = transpose_d<%d,%d>;\n"%(idx, mn[0], mn[1])
        output += "static bool configured_%d = false;\n"%idx
    output += "\n\n"

    output += "switch(idx){\n"
    for idx in sorted(idx_map.keys()):
        output += "case %d:\n"%idx
        output += "if(configured_%d == false){\n"%idx
        output += "  cudaError_t err = cudaFuncSetSharedMemConfig(kern_func_%d, cudaSharedMemBankSizeEightByte);\n"%idx
        output += "  if(err != cudaSuccess) return(-1);\n"
        output += "  configured_%d = true;\n"%idx
        output += "}\n"
        output += "kern_func_%d<<<nblks, 128, 0, *stream>>>(trs_stack+offset, nblks, buffer);\n"%idx
        output += "break;\n"

    output += "// If there is no kernel for these blocks, we don't need to transpose them.\n"
    output += "default: return 0;\n"
    output += "}\n\n"

    output += "return(cudaGetLastError());\n"

    output += "}\n\n\n"

    output += 'extern "C" int libsmm_acc_transpose '
    output += '(void *trs_stack, int offset, int nblks, void *buffer,'
    output += 'int datatype, int m, int n, void* stream) {\n'
    output += 'cudaStream_t* custream = (cudaStream_t*) stream;\n'
    output += 'if(datatype != dbcsr_type_real_8)\n'
    output += '  return 0; //transpose not needed\n'
    output += 'return libcusmm_transpose_d((int *) trs_stack, offset, nblks, (double *) buffer, m, n, custream);\n'
    output += '};\n'

    return(output)


#===============================================================================
def get_includes(plan):
    includes = sorted(set(["./kernels/"+kern.include() for kern in plan.values() ]))
    includes += ["./kernels/cusmm_common.h", "./kernels/cusmm_transpose.h", "../include/libsmm_acc.h"]
    return(includes)


#===============================================================================
def writefile(fn, content):
    if(path.exists(fn)):
        f = open(fn, "r")
        old_content = f.read()
        f.close()
        if(old_content == content):
            return

    f = open(fn, "w")
    f.write(content)
    f.close()

    #if(fn.endswith(".h") or fn.endswith(".cu")):
    #    check_call(["indent",fn])
    #print("Wrote: "+fn)

#===============================================================================
def combinations(*sizes):
     return(list(product(sizes, sizes, sizes)))

#===============================================================================
def product(*args):
    """ Python 2.4 workaround """
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    pools = map(tuple, args)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

#===============================================================================
if(len(sys.argv)==2 and sys.argv[-1]=="--selftest"):
    main(argv=[])
else:
    main(argv=sys.argv[1:])

#EOF

