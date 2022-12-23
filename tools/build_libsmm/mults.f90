MODULE mults

  IMPLICIT NONE

CONTAINS
  FUNCTION trdat(data_type,write_target,in_intent_label)
    INTEGER :: data_type
    LOGICAL :: write_target
    CHARACTER(LEN=*), OPTIONAL :: in_intent_label
    CHARACTER(LEN=50) :: options
    CHARACTER(LEN=50) :: trdat

    options=""    
    IF (PRESENT(in_intent_label)) THEN
       IF (in_intent_label/="") THEN
          options=", INTENT("//TRIM(in_intent_label)//")"
       ENDIF
    ENDIF
    IF (write_target) THEN
       options=TRIM(options)//", TARGET"
    ENDIF
    SELECT CASE(data_type)
    CASE(1)
      trdat="REAL(KIND=KIND(0.0D0))"//TRIM(options)
    CASE(2)
      trdat="REAL(KIND=KIND(0.0))"//TRIM(options)
    CASE(3)
      trdat="COMPLEX(KIND=KIND(0.0D0))"//TRIM(options)
    CASE(4)
      trdat="COMPLEX(KIND=KIND(0.0))"//TRIM(options)
    END SELECT
  END FUNCTION

  FUNCTION trgemm(data_type)
    INTEGER :: data_type
    CHARACTER(LEN=5) :: trgemm
    SELECT CASE(data_type)
    CASE(1)
      trgemm="DGEMM"
    CASE(2)
      trgemm="SGEMM"
    CASE(3)
      trgemm="ZGEMM"
    CASE(4)
      trgemm="CGEMM"
    END SELECT
  END FUNCTION

  FUNCTION trstr(transpose_flavor,data_type)
    INTEGER :: transpose_flavor, data_type
    CHARACTER(LEN=3) :: trstr
    CHARACTER(LEN=1) :: dstr
    SELECT CASE(data_type)
    CASE(1)
     dstr="d"
    CASE(2)
     dstr="s"
    CASE(3)
     dstr="z"
    CASE(4)
     dstr="c"
    END SELECT
    SELECT CASE(transpose_flavor)
    CASE(1)
     trstr=dstr//"nn"
    CASE(2)
     trstr=dstr//"tn"
    CASE(3)
     trstr=dstr//"nt"
    CASE(4)
     trstr=dstr//"tt"
    END SELECT
  END FUNCTION trstr

  FUNCTION trparam(stack_size_label)
    CHARACTER(LEN=*), OPTIONAL :: stack_size_label
    CHARACTER(LEN=128) :: trparam
    if (PRESENT(stack_size_label)) THEN
       trparam = "A,B,C,"//TRIM(stack_size_label)//",dbcsr_ps_width,params,p_a_first,p_b_first,p_c_first"
    ELSE
       trparam = "A,B,C"
    ENDIF
  END FUNCTION trparam

  SUBROUTINE write_stack_params(data_type,stack_size_label)
    INTEGER          :: data_type
    CHARACTER(LEN=*), OPTIONAL :: stack_size_label
    CALL write_matrix_defs(data_type=data_type,write_intent=.TRUE.,write_target=.FALSE.)
    IF (PRESENT(stack_size_label)) THEN
       write(6,'(A)')                    "        INTEGER, INTENT(IN) :: "//TRIM(stack_size_label)//", dbcsr_ps_width"
       write(6,'(A)')                    "        INTEGER, INTENT(IN) :: params(dbcsr_ps_width, "//TRIM(stack_size_label)//")"
       write(6,'(A)')                    "        INTEGER, INTENT(IN) :: p_a_first, p_b_first, p_c_first"
    ENDIF
  END SUBROUTINE write_stack_params

  SUBROUTINE write_matrix_defs(M,N,K,transpose_flavor,data_type,write_intent,&
       write_target,stack_size_label,padding)
   INTEGER, OPTIONAL          :: M,N,K,transpose_flavor
   INTEGER                    :: data_type
   LOGICAL                    :: write_intent, write_target
   CHARACTER(LEN=*), OPTIONAL :: stack_size_label
   LOGICAL, OPTIONAL          :: padding
   CHARACTER(LEN=50)          :: intent_label   
   LOGICAL                    :: do_padding

   IF (PRESENT(M).AND.PRESENT(N).AND.PRESENT(K).AND.PRESENT(transpose_flavor)) THEN
      IF (PRESENT(stack_size_label)) THEN
         ! +8 ... the buffered routines need to be able to read past the last 'used' elements of the C array.
         !        the array therefore needs to be padded appropriately.
         write(6,'(A)') "      "//trdat(data_type,write_target)// &
              " :: C(M*N*"//TRIM(stack_size_label)// &
              "+8), B(K*N*"//TRIM(stack_size_label)// &
              "), A(M*K*"//TRIM(stack_size_label)//")"
      ELSE
         IF (write_intent) THEN
            write(6,'(A,I0,A,I0,A)') &
                 "      "//trdat(data_type,write_target,"INOUT")//" :: C(",M,",",N,")"
            intent_label="IN"
         ELSE
            do_padding=.FALSE.
            IF (PRESENT(padding)) THEN
               IF (padding) do_padding=.TRUE.
            ENDIF
            IF (do_padding) THEN
               write(6,'(A)') &
                    "      "//trdat(data_type,write_target)//" :: C(M*N+8)"
            ELSE
               write(6,'(A,I0,A,I0,A)') &
                    "      "//trdat(data_type,write_target)//" :: C(",M,",",N,")"
            ENDIF
            intent_label=""
         ENDIF
         SELECT CASE(transpose_flavor)
         CASE(1)
            write(6,'(A,I0,A,I0,A,I0,A,I0,A)') &
                 "      "//trdat(data_type,write_target,intent_label)//" :: B(",K,",",N,"), A(",M,",",K,")"                 
         CASE(2)
            write(6,'(A,I0,A,I0,A,I0,A,I0,A)') &
                 "      "//trdat(data_type,write_target,intent_label)//" :: B(",K,",",N,"), A(",K,",",M,")"
         CASE(3)
            write(6,'(A,I0,A,I0,A,I0,A,I0,A)') &
                 "      "//trdat(data_type,write_target,intent_label)//" :: B(",N,",",K,"), A(",M,",",K,")"
         CASE(4)
            write(6,'(A,I0,A,I0,A,I0,A,I0,A)') & 
                 "      "//trdat(data_type,write_target,intent_label)//" :: B(",N,",",K,"), A(",K,",",M,")"
         END SELECT
      ENDIF
   ELSE
      IF (write_intent) THEN
         write(6,'(A)') "      "//trdat(data_type,write_target,"INOUT")//" :: C(*)"
         write(6,'(A)') "      "//trdat(data_type,write_target,"IN")//" :: B(*), A(*)"
      ELSE
         write(6,'(A)') "      "//trdat(data_type,write_target)//" :: C(*)"
         write(6,'(A)') "      "//trdat(data_type,write_target)//" :: B(*), A(*)"
      ENDIF
   ENDIF
  END SUBROUTINE write_matrix_defs

  SUBROUTINE smm_inner(mi,mf,ni,nf,ki,kf,iloop,mu,nu,ku,transpose_flavor,data_type)
     INTEGER :: mi,mf,ni,nf,ki,kf,iloop,mu,nu,ku,transpose_flavor,data_type
     INTEGER :: im,in,ik,ido
     INTEGER :: loop_order(3,6),have_loops

     loop_order(:,1)=(/1,2,3/)
     loop_order(:,2)=(/2,1,3/)
     loop_order(:,3)=(/2,3,1/)
     loop_order(:,4)=(/1,3,2/)
     loop_order(:,5)=(/3,1,2/)
     loop_order(:,6)=(/3,2,1/)
     have_loops=0
     CALL out_loop(mi,mf,ni,nf,ki,kf,mu,nu,ku,loop_order(1,iloop),have_loops)
     CALL out_loop(mi,mf,ni,nf,ki,kf,mu,nu,ku,loop_order(2,iloop),have_loops)
     CALL out_loop(mi,mf,ni,nf,ki,kf,mu,nu,ku,loop_order(3,iloop),have_loops) 
     ! what is the fastest order for these loops ? Does it matter ?
     DO im=0,mu-1
     DO in=0,nu-1
     DO ik=0,ku-1
        SELECT CASE(transpose_flavor) 
        CASE(1)
          write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
           "        C(i+",im,",j+",in,")=C(i+",im,",j+",in,")+A(i+",im,",l+",ik,")*B(l+",ik,",j+",in,")"
        CASE(2)
          write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
           "        C(i+",im,",j+",in,")=C(i+",im,",j+",in,")+A(l+",ik,",i+",im,")*B(l+",ik,",j+",in,")"
        CASE(3)
          write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
           "        C(i+",im,",j+",in,")=C(i+",im,",j+",in,")+A(i+",im,",l+",ik,")*B(j+",in,",l+",ik,")"
        CASE(4)
          write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
           "        C(i+",im,",j+",in,")=C(i+",im,",j+",in,")+A(l+",ik,",i+",im,")*B(j+",in,",l+",ik,")"
        END SELECT
     ENDDO
     ENDDO
     ENDDO
     DO ido=1,have_loops
     write(6,'(A)') "     ENDDO "
     ENDDO
  END SUBROUTINE smm_inner

  SUBROUTINE out_loop(mi,mf,ni,nf,ki,kf,mu,nu,ku,ichoice,have_loops)
     INTEGER :: mi,mf,ni,nf,ki,kf,ichoice,mu,nu,ku,have_loops
     IF (ichoice==1) THEN
        IF (nf-ni+1>nu) THEN
           write(6,'(A,I0,A,I0,A,I0)') "     DO j=",ni,",",nf,",",nu
           have_loops=have_loops+1
        ELSE
           write(6,'(A,I0)') "     j=",ni 
        ENDIF
     ENDIF
     IF (ichoice==2) THEN
        IF (mf-mi+1>mu) THEN
           write(6,'(A,I0,A,I0,A,I0)') "     DO i=",mi,",",mf,",",mu
           have_loops=have_loops+1
        ELSE
           write(6,'(A,I0)') "     i=",mi 
        ENDIF
     ENDIF
     IF (ichoice==3) THEN
        IF (kf-ki+1>ku) THEN
           write(6,'(A,I0,A,I0,A,I0)') "     DO l=",ki,",",kf,",",ku
           have_loops=have_loops+1
        ELSE
           write(6,'(A,I0)') "     l=",ki 
        ENDIF
     ENDIF
  END SUBROUTINE

END MODULE mults
