#
# Author: Alfio Lazzaro, alfio.lazzaro@mat.ethz.ch (2013-2015)
# Library for the generate script used in LIBSMM library    
#

#
# Header of the makefile for small and tiny compilation parts
#
write_makefile_header() {
    #
    # generate the list of indices files
    #
    prefix_file=$1 ; shift

    printf "LIBSMM_DIMS         = $@\n"
    printf "LIBSMM_DIMS_INDICES = \$(foreach m,\$(LIBSMM_DIMS),\$(foreach n,\$(LIBSMM_DIMS),\$(foreach k,\$(LIBSMM_DIMS),\$m_\$n_\$k)))"
    printf "\n\n"

    #
    # consider only a sub-range of all indices
    #
    printf "LIBSMM_SI = 1\n"
    printf "LIBSMM_EI = \$(words \$(LIBSMM_DIMS_INDICES))\n"
    printf "LIBSMM_INDICES = \$(wordlist \$(LIBSMM_SI),\$(LIBSMM_EI),\$(LIBSMM_DIMS_INDICES))\n\n"

    #
    # output directory for compiled and results files 
    #
    printf "LIBSMM_WORKDIR=${work_dir}\n\n"

    #
    # list of source files
    #
    printf "LIBSMM_SRCFILES=\$(patsubst %%,${prefix_file}_find_%%.f90,\$(LIBSMM_INDICES)) \n"
    
    #
    # list of executables
    #
    printf "LIBSMM_OBJFILES=\$(patsubst %%,\$(LIBSMM_WORKDIR)/${prefix_file}_find_%%.o,\$(LIBSMM_INDICES)) \n"
    
    #
    # list of output files
    #
    printf "LIBSMM_OUTFILES=\$(LIBSMM_OBJFILES:.o=.out) \n\n"

    #
    # name of the executable
    #
    printf "LIBSMM_EXE=\$(LIBSMM_WORKDIR)/${prefix_file}_find_\$(firstword \$(LIBSMM_INDICES))__\$(lastword \$(LIBSMM_INDICES)).x \n\n"

    #
    # main target
    #
    printf ".PHONY: all_libsmm compile_libsmm source_libsmm \$(LIBSMM_EXE:.x=.f90) \n"
    printf "all_libsmm: \n\n"

    #
    # include makefile for source master code generation
    #
    printf "LIBSMM_DATATYPE=${strdat}\n"
    if [ -n "${target_compile_offload}" ]; then
	printf "LIBSMM_MIC_OFFLOAD=1\n"
    fi
    printf "include ../make.gen\n\n"

    #
    # include Makefile for libxsmm kernels
    #
    if [[ ("$prefix_file" = "small") && ( -n "${libxsmm_dir}") ]]; then
	printf "LIBXSMM_DIR=${libxsmm_dir}\n"
	printf "LIBXSMM_WORKDIR=${libxsmm_work_dir}\n"
	printf "ROW_MAJOR=0\n"
	printf "${SIMD_libxsmm_target}\n"
	printf "include \$(LIBXSMM_DIR)/Makefile\n"
	printf "LIBXSMM_OBJFILES=\$(patsubst %%,\$(LIBXSMM_WORKDIR)/libxsmm${type_label}_${SIMD_libxsmm}_${prefix_file}_find_%%.o,\$(LIBSMM_INDICES))\n\n"
    fi

    #
    # write general targets
    #
    printf "all_libsmm: \$(LIBSMM_EXE) \n"
    printf "\t rm -f \$(LIBSMM_OUTFILES) \n"
    printf "\t export OMP_NUM_THREADS="
    if [ -n "${MIC_OMP_NUM_THREADS}" ]; then
	printf "${MIC_OMP_NUM_THREADS}"
    else
	printf "${ntasks}"
    fi
    printf " ; ./\$< \n\n"

    printf "\$(LIBSMM_EXE): \$(LIBSMM_OBJFILES) \$(LIBSMM_EXE:.x=.f90)"
    if [[ ("$prefix_file" = "small") && ( -n "${libxsmm_dir}") ]]; then
	printf " \$(LIBXSMM_OBJFILES)"
    fi
    printf "\n"
    if [ -n "${target_compile_offload}" ]; then
	printf "\t ${target_compile_offload} "
    else
	printf "\t ${target_compile} "
    fi
    printf "\$^ -o \$@ ${blas_linking} \n\n"
}


#
# Function to run the make command
# Execution in several jobs if requested
#
run_make() {
    cd ${run_dir}
    mkdir -p ${work_dir}
    if [ -n "${target_compile_offload}" ]; then
	chmod o+w ${work_dir}
    fi

    ndims=$#
    nelements=$((ndims*ndims*ndims))
    nelements_in=$((nelements / jobs))
    nelements_out=$((nelements % jobs))
    element_start=1
    element_end=0

    echo "# Elements = ${ndims}^3 = ${nelements} split in ${jobs} jobs."

    for (( ijob=1 ; ijob<=jobs; ijob++ )); do

        element_end=$(( element_end + nelements_in ))
        if [ ${ijob} -le ${nelements_out} ]; then
            element_end=$(( element_end + 1))
        fi

        test_name=${run_dir}_job${ijob}
        echo "Launching elements ${element_start} --> ${element_end} (${test_name})"
        ${run_cmd} make -j ${ntasks} -f ../${make_file} ${target} LIBSMM_SI=${element_start} LIBSMM_EI=${element_end}

        element_start=$(( element_end + 1))
    done

    cd ..
}


#
# Macro to collected results (small and tiny parts)
#
collect_results() {

    #
    # Stop the execution if the output files were not produced or it is running in batch mode
    #
    if [ "$target" = "source" -o "$target" = "compile" -o "$run_cmd" = "batch_cmd" ]; then
	return
    fi

    echo "Collecting results (it can take several minutes)..."

    out_file=$1 ; shift
    suffix=$1 ; shift

    #
    # analyse results finding optimal tiny or small mults
    #
    (
	for m in $@  ; do
	    for n in $@  ; do
		for k in $@  ; do
		    file=${suffix}_find_${m}_${n}_${k}.out
		    res=`tail -n 1 ${run_dir}/${work_dir}/$file`
		    if [ -z "${res}" ]; then
			res=`tail -n 2 ${run_dir}/${work_dir}/$file`
		    fi
		    echo "$m $n $k $res"
		done ; done ; done
    ) > ${out_file}
}


do_generate_tiny() {
    # 
    # skip the compilation part if it needs only to collect the results
    #
    if [ "$run_cmd" != "true" ]; then
        #
        # compile the generator of tiny mults
        #
	${host_compile} -c mults.f90 
	${host_compile} mults.o tiny_gen.f90 -o tiny_gen.x
	
        #
        # for easy parallelism go via a Makefile
        #
	rm -f ${make_file}
	(
	    write_makefile_header tiny "${dims_tiny}"

	    #
	    # Write specific targets
	    #
	    printf "compile_libsmm: \$(LIBSMM_OBJFILES) \n"
	    printf "\$(LIBSMM_WORKDIR)/%%.o: %%.f90 \n"
	    printf "\t ${target_compile} -c \$< -o \$@ \n\n"

	    printf "source_libsmm: \$(LIBSMM_SRCFILES) \n"
	    printf "%%.f90: \n"
	    printf '\t .././tiny_gen.x `echo $* | awk -F_ '\''{ print $$3" "$$4" "$$5 }'\''`'
	    printf " ${transpose_flavor} ${data_type} > \$@ \n\n"
	) > ${make_file}

	run_make ${dims_tiny}

    fi # close if [ "$run_cmd" != "true" ]

    collect_results ${tiny_file} "tiny" "${dims_tiny}"
}


do_generate_small() {
    #
    # Check if tiny file exists
    #
    if [ ! -s ${tiny_file} ]; then
	echo "Tiny file ${tiny_file} doesn't exist. Have you ran tiny1/tiny2 phase?"
	echo "Abort execution."
	echo
	exit
    fi

    #
    # skip the compilation part if it needs only to collect the results
    #
    if [ "$run_cmd" != "true" ]; then
	#
	# compile the generator of small mults
	#
	${host_compile} -c mults.f90
	${host_compile} -c multrec_gen.f90
	${host_compile} mults.o multrec_gen.o small_gen.f90 -o small_gen.x
	#
	# compile libxsmm generator
	#
	if [[ ( -n "${libxsmm_dir}") ]]; then
	    libxsmm_work_dir=libxsmm
	    mkdir -p ${run_dir}/${libxsmm_work_dir}
	    ( cd ${libxsmm_dir} && make -j ${ntasks} generator )
	fi
	#
	# for easy parallelism go via a Makefile
	#
	rm -f ${make_file}
	(
	    write_makefile_header small "${dims_small}"
    	    #
	    # Write specific targets
	    #
	    printf "compile_libsmm: \$(LIBSMM_OBJFILES)"
	    if [[ ( -n "${libxsmm_dir}") ]]; then
		printf " \$(LIBXSMM_OBJFILES)"
	    fi
	    printf "\n"
	    printf "\$(LIBSMM_WORKDIR)/%%.o: \$(LIBSMM_WORKDIR)/%%.f90 \n"
	    printf "\t ${target_compile} -c \$< -o \$@ \n"
    	    #
	    # compile libxsmm kernel
	    #
	    if [[ ( -n "${libxsmm_dir}") ]]; then
		printf "\$(LIBXSMM_WORKDIR)/%%.o: \$(LIBXSMM_WORKDIR)/%%.c \n"
		printf "\t \$(CC) \$(CFLAGS) \$(DFLAGS) \$(TARGET) -c \$< -o \$@ \n"
	    fi
	    printf "\n"

	    printf "source_libsmm: \$(addprefix \$(LIBSMM_WORKDIR)/,\$(LIBSMM_SRCFILES))"
	    if [[ ( -n "${libxsmm_dir}") ]]; then
		printf " \$(LIBXSMM_OBJFILES:.o=.c)"
	    fi
	    printf "\n"
	    printf "\$(LIBSMM_WORKDIR)/%%.f90: ../${tiny_file} \n"
	    printf '\t .././small_gen.x `echo $* | awk -F_ '\''{ print $$3" "$$4" "$$5 }'\''`'
	    printf " ${transpose_flavor} ${data_type} ${SIMD_size} ../${tiny_file} "
	    if [ -n "${libxsmm_dir}" ]; then
		printf "1 "
	    fi
	    printf "> \$@\n"
	    #
	    # generate libxsmm kernel
	    #
	    if [[ ( -n "${libxsmm_dir}") ]]; then
		printf "\$(LIBXSMM_WORKDIR)/%%.c:\n"
		printf "\t rm -f \$@\n"
		printf "\t \$(LIBXSMM_DIR)/bin/./generator dense"
		printf " \$@"
		printf " libxsmm_"
		printf '`echo $* | awk -F_ '\''{ print $$6"_"$$7"_"$$8" "$$6" "$$7" "$$8" "$$6" "$$8" "$$6 }'\''` '
		printf "1 1 0 0 ${SIMD_libxsmm} nopf ${data_libxsmm}" 
	    fi
	    printf "\n\n"
	) > ${make_file}

	run_make ${dims_small}

    fi # close if [ "$run_cmd" != "true" ]

    collect_results ${small_file} "small" "${dims_small}"

}

do_generate_lib() {
    #
    # Check if tiny file exists
    #
    if [ ! -s ${tiny_file} ]; then
	echo "Tiny file ${tiny_file} doesn't exist. Have you ran tiny1/tiny2 phase?"
	echo "Abort execution."
	echo
	exit
    fi    
    
    #
    # Check if small file exists
    #
    if [ ! -s ${small_file} ]; then
	echo "Small file ${small_file} doesn't exist. Have you ran small1/small2 phase?"
	echo "Abort execution."
	echo
	exit
    fi    

    #
    # compile the generator of small mults
    #
    ${host_compile} -c mults.f90 
    ${host_compile} -c multrec_gen.f90 
    ${host_compile} mults.o multrec_gen.o lib_gen.f90 -o lib_gen.x

    #
    # directory for the library
    #
    mkdir -p lib

    #
    # generate the generic caller
    #
    echo "Generate the generic caller... it can take a while..."

    maxsize=-1
    numsize=0

    for myn in ${dims_small}
    do
	numsize=$((numsize+1))
	maxsize=`echo "$myn $maxsize" | awk '{if ($1>$2) { print $1 } else { print $2 } }'`
    done

    #
    # generate a translation array for the jump table
    #
    count=0
    eles="(/0"
    for i in `seq 1 $maxsize`
    do
	
	found=0
	for myn in ${dims_small}
	do
	    if [ "$myn" == "$i" ]; then
		found=1
	    fi
	done
	if [ "$found" == 1 ]; then
	    count=$((count+1))
	    ele=$count
	else
	    ele=0 
	fi
	eles="$eles,$ele"
    done
    eles="$eles/)"
   
    cd ${run_dir}

    file="smm${type_label}.f90"
    rm -f ${file}

    write_routine() {

	if [ $# -eq 0 ]; then
	    printf "SUBROUTINE smm_vec${type_label}(M,N,K,A,B,C,stack_size,dbcsr_ps_width,params,p_a_first,p_b_first,p_c_first)\n" >> ${file}
	else
	    printf "SUBROUTINE smm${type_label}(M,N,K,A,B,C)\n" >> ${file}
	fi
	printf " INTEGER :: M,N,K,LDA,LDB\n ${strdat} :: C(*), B(*), A(*)\n" >> ${file}
	if [ $# -eq 0 ]; then
	    printf " INTEGER :: stack_size, dbcsr_ps_width \n" >> ${file}
	    printf " INTEGER :: params(dbcsr_ps_width,stack_size) \n" >> ${file}
	    printf " INTEGER :: p_a_first,p_b_first,p_c_first \n" >> ${file}
	    printf " INTEGER :: sp\n" >> ${file}
	fi
	printf " INTEGER, PARAMETER :: indx(0:$maxsize)=&\n $eles\n" >> ${file}
	printf " INTEGER :: im,in,ik,itot\n" >> ${file}
	printf " $strdat, PARAMETER :: one=1\n" >> ${file}
	printf " IF (M<=$maxsize) THEN\n   im=indx(M)\n ELSE\n   im=0\n ENDIF\n" >> ${file}
	printf " IF (N<=$maxsize) THEN\n   in=indx(N)\n ELSE\n   in=0\n ENDIF\n" >> ${file}
	printf " IF (K<=$maxsize) THEN\n   ik=indx(K)\n ELSE\n   ik=0\n ENDIF\n" >> ${file}
	printf " itot=(ik*($numsize+1)+in)*($numsize+1)+im\n" >> ${file}
	
	count=0
	printf " SELECT CASE(itot)\n" >> ${file}
	for myk in 0 ${dims_small}
	do
	    for myn in 0 ${dims_small}
	    do
		for mym in 0 ${dims_small}
		do
		    printf " CASE($count)\n " >> ${file}
		    prod=$((myk*myn*mym))
		    if [[ "$prod" == "0" ]]; then
			printf '   GOTO 999\n' >> ${file}
		    else
			if [ -n "${target_compile_offload}" ]; then
			    printf '!dir$ attributes offload:mic ::' >> ${file}
			    printf "smm${type_label}_${mym}_${myn}_${myk}" >> ${file}
			    if [ $# -eq 0 ]; then
				printf "_stack" >> ${file}
			    fi
			    printf "\n" >> ${file}
			fi
			if [ $# -eq 0 ]; then
			    printf "   CALL smm${type_label}_${mym}_${myn}_${myk}_stack(A,B,C,stack_size,dbcsr_ps_width,params,p_a_first,p_b_first,p_c_first)\n" >> ${file}
			else
			    printf "   CALL smm${type_label}_${mym}_${myn}_${myk}(A,B,C)\n" >> ${file}
			fi
		    fi
		    count=$((count+1))
		done
	    done
	done
	printf " END SELECT\n" >> ${file}
	printf " RETURN\n" >> ${file}
	printf "999 CONTINUE \n" >> ${file}
	printf " ${lds}\n" >> ${file}
	if [ -n "${target_compile_offload}" ]; then
	    printf '!dir$ attributes offload:mic ::' >> ${file}
	    printf "${gemm}\n" >> ${file}
	fi
	if [ $# -eq 0 ]; then
	    printf " DO sp = 1, stack_size\n" >> ${file}
	    printf "   CALL ${gemm}('%s','%s',M,N,K,one,A(params(p_a_first,sp)),LDA,B(params(p_b_first,sp)),LDB,one,C(params(p_c_first,sp)),M)\n" $ta $tb >> ${file}
	    printf " ENDDO\n" >> ${file}
	    printf "END SUBROUTINE smm_vec${type_label}\n\n" >> ${file}
	else
	    printf " CALL ${gemm}('%s','%s',M,N,K,one,A,LDA,B,LDB,one,C,M)\n" $ta $tb >> ${file}
	    printf "END SUBROUTINE smm${type_label}\n\n" >> ${file}
	fi
    }
    
    write_routine
    write_routine 0

    cd ..

    #
    # compile libxsmm generator
    #
    if [[ ( -n "${libxsmm_dir}") ]]; then
	libxsmm_work_dir=libxsmm
	mkdir -p ${run_dir}/${libxsmm_work_dir}
	( cd ${libxsmm_dir} && make -j ${ntasks} generator )
    fi

    #
    # for easy parallelism go via a Makefile
    #
    rm -f ${make_file}
    rm -f ${archive}

    (
	printf "LIBSMM_DIMS = ${dims_small}\n"
	printf "LIBSMM_DIMS_INDICES = \$(foreach m,\$(LIBSMM_DIMS),\$(foreach n,\$(LIBSMM_DIMS),\$(foreach k,\$(LIBSMM_DIMS),\$m_\$n_\$k)))"
	printf "\n\n"

	#
	# output directory for compiled and results files
	#
	printf "LIBSMM_WORKDIR=${work_dir}\n\n"

	#
	# driver file
	#
	printf "LIBSMM_DRIVER=${file%.*}.o \n\n"

	#
	# list of source files
	#
	printf "LIBSMM_SRCFILES=\$(patsubst %%,\$(LIBSMM_WORKDIR)/smm${type_label}_%%.f90,\$(LIBSMM_DIMS_INDICES)) \n"

	#
	# list of obj files
	#
	printf "LIBSMM_OBJFILES=\$(patsubst %%,\$(LIBSMM_WORKDIR)/smm${type_label}_%%.o,\$(LIBSMM_DIMS_INDICES)) \n\n"

	printf ".PHONY: \$(LIBSMM_WORKDIR)/\$(LIBSMM_DRIVER) all_libsmm \n\n"
	
	printf "all_libsmm: \n\n"

        #
        # include Makefile for libxsmm kernels
        #
	if [[ ( -n "${libxsmm_dir}") ]]; then
            printf "LIBXSMM_DIR=${libxsmm_dir}\n"
            printf "LIBXSMM_WORKDIR=${libxsmm_work_dir}\n"
            printf "ROW_MAJOR=0\n"
            printf "${SIMD_libxsmm_target}\n"
            printf "include \$(LIBXSMM_DIR)/Makefile\n"
            printf "LIBXSMM_OBJFILES=\$(patsubst %%,\$(LIBXSMM_WORKDIR)/libxsmm${type_label}_${SIMD_libxsmm}_%%.o,\$(LIBSMM_DIMS_INDICES))\n\n"
	fi

	#
	# generation source rule
	#
	printf "source_libsmm: \$(LIBSMM_SRCFILES) \n"
	printf "\$(LIBSMM_WORKDIR)/%%.f90: ../${small_file} ../${tiny_file} \n"
	printf '\t .././lib_gen.x `echo $* | awk -F_ '\''{ print $$3" "$$4" "$$5 }'\''`'
	printf " ${transpose_flavor} ${data_type} ${SIMD_size} ../${small_file} ../${tiny_file} > \$@\n"
        #
        # generate libxsmm kernel
        #
        if [[ ( -n "${libxsmm_dir}") ]]; then
            printf "\$(LIBXSMM_WORKDIR)/%%.c:\n"
            printf "\t rm -f \$@\n"
            printf "\t \$(LIBXSMM_DIR)/bin/./generator dense"
            printf " \$@"
            printf " libxsmm_"
            printf '`echo $* | awk -F_ '\''{ print $$4"_"$$5"_"$$6" "$$4" "$$5" "$$6" "$$4" "$$6" "$$4 }'\''` '
            printf "1 1 0 0 ${SIMD_libxsmm} nopf ${data_libxsmm}\n" 
        fi
	printf "\n"

	#
	# compile rule
	#
	printf "compile_libsmm: \$(LIBSMM_OBJFILES) \n"
	printf "\$(LIBSMM_WORKDIR)/%%.o: \$(LIBSMM_WORKDIR)/%%.f90 \n"
	printf "\t ${target_compile} -c \$< -o \$@ \n\n"
        #
        # compile libxsmm kernel
        #
        if [[ ( -n "${libxsmm_dir}") ]]; then
            printf "\$(LIBXSMM_WORKDIR)/%%.o: \$(LIBXSMM_WORKDIR)/%%.c \n"
            printf "\t \$(CC) \$(CFLAGS) \$(DFLAGS) \$(TARGET) -c \$< -o \$@ \n"
        fi
        printf "\n"

	printf "\$(LIBSMM_WORKDIR)/\$(LIBSMM_DRIVER): \n"
	printf "\t ${target_compile} -c \$(notdir \$*).f90 -o \$@ \n\n"

	printf "all_libsmm: ${archive} \n"
	printf "${archive}: \$(LIBSMM_OBJFILES) \$(LIBSMM_WORKDIR)/\$(LIBSMM_DRIVER)"
	if [ -n "${libxsmm_dir}" ]; then
	    printf " \$(LIBXSMM_OBJFILES)"
	fi
	printf "\n"
	if [ -n "${target_compile_offload}" ]; then
	    printf "\t xiar -rs -qoffload-build \$@ \$^ \n"
	else
	    printf "\t ar -rs \$@ \$^ \n"
	fi
	printf "\t @echo 'Library produced at `pwd`/lib/lib${library}.a'\n\n"

    ) > ${make_file}

    cd ${run_dir}

    #
    # make the output dir
    #
    mkdir -p ${work_dir}

    #
    # execute makefile compiling all variants and executing them
    #
    test_name="${run_dir}_job${jobs}"
    ${run_cmd} make -j ${ntasks} -f ../${make_file} ${target}

}

do_check() {
    #
    # check if the library exists
    #
    archive="lib/lib${library}.a"
    if [ ! -s ${archive} ]; then
	echo "ERROR: the library $archive doesn't exist."
	exit
    fi

    # 
    # skip the compilation part if it needs only to collect the results
    #
    if [ "$run_cmd" != "true" ]; then

	cd ${run_dir}

	echo "Split in ${jobs} jobs."

	ndims=`echo ${dims_small} | wc -w`
	nelements=$((ndims*ndims*ndims))
	nelements_in=$((nelements / jobs))
	nelements_out=$((nelements % jobs))
	element=1
	element_end=0
	ijob=0

        #
        # make the output dir
        #
        mkdir -p ${work_dir}

	for m in ${dims_small}  ; do
	    for n in ${dims_small}  ; do
		for k in ${dims_small}  ; do

		    if [ $element_end -eq 0 -o $element -gt $element_end ]; then
			element_end=$(( element_end + nelements_in ))
			ijob=$(( ijob + 1))
			
			if [ ${ijob} -le ${nelements_out} ]; then
			    element_end=$(( element_end + 1))
			fi
			
			echo "Preparing test program for job #$ijob..."
			filename=${test_file}_job$ijob

cat << EOF > ${work_dir}/${filename}.f90
MODULE WTF_job$ijob
  INTERFACE MYRAND
    MODULE PROCEDURE SMYRAND, DMYRAND, CMYRAND, ZMYRAND
  END INTERFACE
CONTAINS
EOF
if [ -n "${target_compile_offload}" ]; then
    printf '!dir$ attributes offload:mic :: DMYRAND \n' >> ${work_dir}/${filename}.f90
fi
cat << EOF >> ${work_dir}/${filename}.f90
  SUBROUTINE DMYRAND(A)
    REAL(KIND=KIND(0.0D0)), DIMENSION(:,:) :: A
    REAL(KIND=KIND(0.0)), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Aeq
    CALL RANDOM_NUMBER(Aeq)
    A=Aeq
  END SUBROUTINE
EOF
if [ -n "${target_compile_offload}" ]; then
    printf '!dir$ attributes offload:mic :: SMYRAND \n' >> ${work_dir}/${filename}.f90
fi
cat << EOF >> ${work_dir}/${filename}.f90
  SUBROUTINE SMYRAND(A)
    REAL(KIND=KIND(0.0)), DIMENSION(:,:) :: A
    REAL(KIND=KIND(0.0)), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Aeq
    CALL RANDOM_NUMBER(Aeq)
    A=Aeq
  END SUBROUTINE
EOF
if [ -n "${target_compile_offload}" ]; then
    printf '!dir$ attributes offload:mic :: CMYRAND \n' >> ${work_dir}/${filename}.f90
fi
cat << EOF >> ${work_dir}/${filename}.f90
  SUBROUTINE CMYRAND(A)
    COMPLEX(KIND=KIND(0.0)), DIMENSION(:,:) :: A
    REAL(KIND=KIND(0.0)), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Aeq,Beq
    CALL RANDOM_NUMBER(Aeq)
    CALL RANDOM_NUMBER(Beq)
    A=CMPLX(Aeq,Beq,KIND=KIND(0.0))
  END SUBROUTINE
EOF
if [ -n "${target_compile_offload}" ]; then
    printf '!dir$ attributes offload:mic :: ZMYRAND \n' >> ${work_dir}/${filename}.f90
fi
cat << EOF >> ${work_dir}/${filename}.f90
  SUBROUTINE ZMYRAND(A)
    COMPLEX(KIND=KIND(0.0D0)), DIMENSION(:,:) :: A
    REAL(KIND=KIND(0.0)), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Aeq,Beq
    CALL RANDOM_NUMBER(Aeq)
    CALL RANDOM_NUMBER(Beq)
    A=CMPLX(Aeq,Beq,KIND=KIND(0.0D0))
  END SUBROUTINE
END MODULE

EOF
if [ -n "${target_compile_offload}" ]; then
    printf '!dir$ attributes offload:mic :: testit, ' >> ${work_dir}/${filename}.f90
    printf "${gemm}, smm${type_label} \n" >> ${work_dir}/${filename}.f90
fi
cat << EOF >> ${work_dir}/${filename}.f90
SUBROUTINE testit(M,N,K)
  USE WTF_job$ijob
  IMPLICIT NONE
  INTEGER :: M,N,K

  $strdat :: C1(M,N), C2(M,N)
  $strdat :: ${decl}
  $strdat, PARAMETER :: one=1
  INTEGER :: i,LDA,LDB

  REAL(KIND=KIND(0.0D0)) :: flops,gflop
  REAL :: t1,t2,t3,t4
  INTEGER :: Niter

  flops=2*REAL(M,KIND=KIND(0.0D0))*N*K
  gflop=1000.0D0*1000.0D0*1000.0D0
  ! assume we would like to do 1 Gflop for testing a subroutine
  Niter=MAX(1,CEILING(MIN(10000000.0D0,1*gflop/flops)))
  ${lds}

  DO i=1,3
     CALL MYRAND(A)
     CALL MYRAND(B)
     CALL MYRAND(C1)
     C2=C1

     CALL ${gemm}("$ta","$tb",M,N,K,one,A,LDA,B,LDB,one,C1,M) 
     CALL smm${type_label}(M,N,K,A,B,C2)

     IF (MAXVAL(ABS(C2-C1))>100*EPSILON(REAL(1.0,KIND=KIND(A(1,1))))) THEN
        write(6,*) "Matrix size",M,N,K
        write(6,*) "Max diff=",MAXVAL(ABS(C2-C1))
        write(6,*) "A=",A
        write(6,*) "B=",B
        write(6,*) "C1=",C1
        write(6,*) "C2=",C2
        write(6,*) "BLAS and smm yield different results : possible compiler bug... do not use the library ;-)"
        ERROR STOP
     ENDIF
  ENDDO

  A=0; B=0; C1=0 ; C2=0

  CALL CPU_TIME(t1) 
  DO i=1,Niter
     CALL ${gemm}("$ta","$tb",M,N,K,one,A,LDA,B,LDB,one,C1,M) 
  ENDDO
  CALL CPU_TIME(t2) 

  CALL CPU_TIME(t3)
  DO i=1,Niter
     CALL smm${type_label}(M,N,K,A,B,C2)
  ENDDO
  CALL CPU_TIME(t4)

  WRITE(6,'(A,I5,I5,I5,A,F6.3,A,F6.3,A,F12.3,A)') "Matrix size ",M,N,K, &
        " smm: ",Niter*flops/(t4-t3)/gflop," Gflops. Linked blas: ",Niter*flops/(t2-t1)/gflop,&
        " Gflops. Performance ratio: ",((t2-t1)/(t4-t3))*100,"%"

END SUBROUTINE 
EOF
if [ -n "${target_compile_offload}" ]; then
    printf '!dir$ attributes offload:mic :: testit \n ' >> ${work_dir}/${filename}.f90
fi
cat << EOF >> ${work_dir}/${filename}.f90
PROGRAM tester
  IMPLICIT NONE
EOF
if [ -n "${target_compile_offload}" ]; then
    printf '!dir$ offload begin target(mic) \n' >> ${work_dir}/${filename}.f90
fi
cat << EOF >> ${work_dir}/${filename}.f90
  !\$omp parallel
EOF

		    fi

		    echo "   CALL testit(${m},${n},${k})" >> ${work_dir}/${filename}.f90

		    # last entry for a job
		    if [ $element -eq $element_end ]; then

			# last job
			if [ $ijob -eq $jobs -a ! -n "${target_compile_offload}" ]; then
cat << EOF >> ${work_dir}/${filename}.f90
  ! checking 'peak' performance (and a size likely outside of the library)
  CALL testit(1000,1000,1000)
EOF
			fi

cat << EOF >> ${work_dir}/${filename}.f90
  !\$omp end parallel
EOF
if [ -n "${target_compile_offload}" ]; then
    printf '!dir$ end offload \n' >> ${work_dir}/${filename}.f90
fi
cat << EOF >> ${work_dir}/${filename}.f90
END PROGRAM
EOF

			echo "Launching job #$ijob"
			test_name="${run_dir}_job${ijob}"

			exe=${filename}.x

			rm -f ${work_dir}/${filename}.sh
			#
			# Prepare the script for compile the benchmarking 
			# and testing program for the smm library
			#
			(
			    printf "#!/bin/bash -e \n\n"
			    printf "set -o pipefail\n\n"
			    printf "cd ${work_dir}\n"
			    if [ -n "${target_compile_offload}" ]; then
				printf "${target_compile_offload}"
			    else
				printf "${target_compile}"
			    fi
			    printf " ${filename}.f90 -o ${exe} ../../${archive} ${blas_linking}\n"
			    printf "export OMP_NUM_THREADS=1 ; ./${exe} | tee ${filename}.out\n"
			) > ${work_dir}/${filename}.sh
			chmod +x ${work_dir}/${filename}.sh

			${run_cmd} ./${work_dir}/${filename}.sh
		    fi

		    element=$(( element + 1 ))  

		done ; done ; done

	#
	# Stop the execution if it is running in batch mode
	#
	if [ "$run_cmd" = "batch_cmd" ]; then
	    return
	fi

	cd ..

    fi # close if [ "$run_cmd" != "true" ]

    #
    # Merge the output files in the case of $run_cmd" = "true"
    #
    rm -f ${test_file}_${config_file_name}.out

    for (( ijob=1 ; ijob<=jobs; ijob++ )); do
	if [ ! -s ${run_dir}/${work_dir}/${test_file}"_job$ijob".out ]; then
	    echo "ERROR: Empty check file \"${run_dir}/${work_dir}/${test_file}_job$ijob.out\""
	    exit
	fi
	cat ${run_dir}/${work_dir}/${test_file}"_job$ijob".out | grep "Matrix size" >> ${test_file}_${config_file_name}.out
    done

    #
    # We're done... protect the user from bad compilers
    #
    set +e
    grep "BLAS and smm yield different results" ${test_file}_${config_file_name}.out >& /dev/null
    if [ "$?" == "0" ]; then
	echo "Library is miscompiled... removing ${archive}"
	echo
	rm -f ${archive}
    else
	pathhere=`pwd -P`
	echo "Done... check performance looking at ${test_file}_${config_file_name}.out"
	echo "Final library can be found at ${archive}"
	echo
    fi
    
}
