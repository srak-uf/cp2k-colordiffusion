!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE colordiffusion_types
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set

  USE cp_linked_list_input,            ONLY: cp_sll_val_next,&
                                             cp_sll_val_type
                                             
  USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&  
                                             cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_constants,                 ONLY: do_ne_ccd, &
                                             do_ne_const_curr
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_list_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE simpar_types,                    ONLY: simpar_type
  USE particle_list_types,             ONLY: particle_list_type
  !!! Read &INDEX_CC
  USE input_val_types,                 ONLY: val_get,&
                                             val_type
  USE cp_parser_methods,               ONLY: read_float_object,&
                                             read_integer_object

#include "../base/base_uses.f90"
    
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: colordiffusion_type, create_colordiffusion

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'colordiffusion_methods'
  INTEGER, PRIVATE, SAVE               :: last_colordiffusion_id_nr = 0

  TYPE colordiffusion_type
    INTEGER                                 :: id_nr, ref_count, id_cd
    CHARACTER(LEN=default_string_length)    :: label
    REAL(KIND=dp)                           :: FX, FY, FZ, JX, JY, JZ
    REAL(KIND=dp), DIMENSION(:), POINTER    :: REGX_PLUS, REGX_MINUS
  END TYPE colordiffusion_type


CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param nemd ...
!> \param nemd_section ...
! **************************************************************************************************
  SUBROUTINE create_colordiffusion(colordiffusion, md_section, simpar, force_env)
    TYPE(colordiffusion_type), POINTER              :: colordiffusion
    TYPE(simpar_type), POINTER                      :: simpar
    TYPE(section_vals_type), POINTER                :: md_section
    TYPE(section_vals_type), POINTER                :: color_section
    TYPE(particle_list_type), POINTER               :: particles  ! This derived type name has not been declared.
    INTEGER, DIMENSION(:), POINTER                  :: tmplist, tmplist_plus,&
                                                       tmplist_minus
    REAL(KIND=dp)                                   :: TOTCHG
    TYPE(force_env_type), POINTER                   :: force_env

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colordiffusion'
    CHARACTER(len=*), PARAMETER :: routineP = moduleN//':'//routineN

    LOGICAL                                         :: check,explicit
    TYPE(cp_subsys_type), POINTER                   :: subsys
    INTEGER                                         :: i, il, ipart, j
    TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: atomic_kind_set
    LOGICAL                                         :: shell_present
    TYPE(atomic_kind_list_type), POINTER            :: atomic_kinds

    !! Read
    !! &INDEX_CC
    !!  [Index]   [Color charge]
    !! &END &NDEX_CC
    TYPE(cp_sll_val_type), POINTER                  :: list
    TYPE(val_type), POINTER                         :: val
    LOGICAL                                         :: is_ok
    CHARACTER(len=2*default_string_length)          :: line_att
    INTEGER                                         :: start_c, end_c, iatom, wrd,&
                                                       n_rep_val
    CHARACTER(len=default_string_length)            :: error_message
    TYPE(section_vals_type), POINTER                :: index_cc_section


    check = .NOT. ASSOCIATED(colordiffusion)
    simpar%do_cdnemd =.FALSE.
    CPASSERT(check)

    color_section => section_vals_get_subs_vals(md_section, "COLOR")
    CALL section_vals_val_get(color_section, "CDTYPE", explicit=explicit)
    IF (explicit) THEN
        ALLOCATE(colordiffusion)
        last_colordiffusion_id_nr = last_colordiffusion_id_nr+1
        colordiffusion%id_nr = last_colordiffusion_id_nr
        colordiffusion%ref_count = 1
        simpar%do_cdnemd = .TRUE.
        CALL section_vals_val_get(color_section, "CDTYPE", i_val=colordiffusion%id_cd)
        CALL section_vals_val_get(color_section, "FX", r_val=colordiffusion%FX)
        CALL section_vals_val_get(color_section, "FY", r_val=colordiffusion%FY)
        CALL section_vals_val_get(color_section, "FZ", r_val=colordiffusion%FZ)
        ! colordiffusion%FX = 0.019446905434120372 * colordiffusion%FX
        ! colordiffusion%FY = 0.019446905434120372 * colordiffusion%FY
        ! colordiffusion%FZ = 0.019446905434120372 * colordiffusion%FZ
        CALL section_vals_val_get(color_section, "JX", r_val=colordiffusion%JX)
        CALL section_vals_val_get(color_section, "JY", r_val=colordiffusion%JY)
        CALL section_vals_val_get(color_section, "JZ", r_val=colordiffusion%JZ)

        ! force_env_get => subsys => atomic_kinds koredeiinoka???? nidodema ppoi
        CALL force_env_get(force_env, subsys=subsys)
        CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds,particles=particles)

        atomic_kind_set => atomic_kinds%els
        CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, shell_present=shell_present)
        IF (shell_present) THEN
          CALL cp_abort(__LOCATION__, "CDTYPE not compatible with shell_present: particles only!")
        END IF
        particles%els(:)%ccharge = 0.0_dp
        particles%els(:)%ctag = 0


        SELECT CASE (colordiffusion%id_cd)
          CASE(do_ne_ccd, do_ne_const_curr)
            index_cc_section => section_vals_get_subs_vals(color_section,"INDEX_CC")
            CALL section_vals_val_get(index_cc_section, "_DEFAULT_KEYWORD_", n_rep_val=n_rep_val)
            CALL section_vals_list_get(index_cc_section, "_DEFAULT_KEYWORD_", list=list)
            DO iatom = 1, n_rep_val
              is_ok = cp_sll_val_next(list, val)
              CALL val_get(val, c_val=line_att)
              start_c = 1
              DO wrd = 1, 2
                DO j = start_c, LEN(line_att)
                  IF (line_att(j:j) /= ' ') THEN
                    start_c = j
                    EXIT
                  END IF
                END DO

                ! IF (LEN_TRIM(line_att(start_c:end_c-1)) == 0) &
                !     CALL cp_abort(__LOCATION__, &
                !                   "Incorrectly formatted input line for atom "// &
                !                   TRIM(ADJUSTL(cp_to_string(iatom)))// &
                !                   " found in INDEX_CC section. Input line: <"// &
                !                   TRIM(line_att)//"> ")

                end_c = LEN(line_att)+1
                DO j = start_c, LEN(line_att)
                  IF (line_att(j:j) == ' ') THEN
                    end_c = j
                    EXIT
                  END IF
                END DO
                SELECT CASE (wrd)
                  CASE (1)
                    CALL read_integer_object(line_att(start_c:end_c-1), &
                                            ipart, error_message)
                  CASE (2)
                    CALL read_float_object(line_att(start_c:end_c-1), &
                                          particles%els(ipart)%ccharge, error_message)
                    IF (LEN_TRIM(error_message) /= 0) &
                      CALL cp_abort(__LOCATION__, &
                                    "Incorrectly formatted input line for atom "// &
                                    TRIM(ADJUSTL(cp_to_string(iatom)))// &
                                    " found in INDEX_CC section. "//TRIM(error_message)// &
                                    " Input line: <"//TRIM(line_att)//"> ")
                END SELECT
                start_c = end_c
                IF (start_c > LEN_TRIM(line_att)) EXIT
              END DO
            END DO
        END SELECT

        SELECT CASE(colordiffusion%id_cd)
          CASE(do_ne_ccd, do_ne_const_curr)
            IF (abs(sum(particles%els(:)%ccharge)) > 1e-10) THEN
              CALL cp_abort(__LOCATION__, "Color charge not balanced!!")
            ELSE
              TOTCHG = sum(particles%els(:)%ccharge) !!!!!
              DO j = 1, size(particles%els(:))
                IF (particles%els(j)%ccharge > 0.1 ) THEN
                  particles%els(j)%ccharge = particles%els(j)%ccharge - TOTCHG * 1.0_dp
                  EXIT
                END IF
              END DO
            END IF 
        END SELECT

        CALL CD_info(colordiffusion, particles, color_section)

    END IF
  END SUBROUTINE create_colordiffusion


! *****************************************************************************
!> \brief ...
! *****************************************************************************
  SUBROUTINE CD_info(colordiffusion, particles, section)
    TYPE(colordiffusion_type), POINTER       :: colordiffusion
    TYPE(section_vals_type), POINTER         :: section
    TYPE(particle_list_type), POINTER        :: particles

    CHARACTER(len=*), PARAMETER :: routineN = 'CD_info', routineP = moduleN//':'//routineN

    INTEGER                                  :: iw
    LOGICAL                                  :: new_file
    TYPE(cp_logger_type), POINTER            :: logger
    integer                                  :: i


    NULLIFY(logger)
    logger => cp_get_default_logger()
    iw = cp_print_key_unit_nr(logger,section,"PRINT%CD_INFO",extension=".log")
    IF (iw > 0) THEN
       WRITE ( iw , '( /,A )' ) ' COLOR-DIFFUSION | NE Info for MD'
       IF (colordiffusion%id_cd == do_ne_ccd) THEN
          WRITE ( iw, '( A,T69,A )' ) ' COLOR-DIFFUSION | Type of algorithm','const.F-CCD'
       ELSE IF (colordiffusion%id_cd == do_ne_const_curr) THEN
          WRITE ( iw, '( A,A,A,T69,F12.9 )' ) &
            ' COLOR-DIFFUSION | Field along X [', ' hartree*bohr^-1 ' ,']',colordiffusion%FX
          WRITE ( iw, '( A,A,A,T69,F12.9 )' ) &
            ' COLOR-DIFFUSION | Force along Y [', ' hartree*bohr^-1 ' ,']',colordiffusion%FY
          WRITE ( iw, '( A,A,A,T69,F12.9 )' ) &
            ' COLOR-DIFFUSION | Force along Z [', ' hartree*bohr^-1 ' ,']',colordiffusion%FZ
          WRITE ( iw, '( A,A,A,T69,F12.9 )' ) &
            ' COLOR-DIFFUSION | Preset flux along X [', ' a.u. ' ,']',colordiffusion%JX
          WRITE ( iw, '( A,A,A,T69,F12.9 )' ) &
            ' COLOR-DIFFUSION | Preset flux along Y [', ' a.u. ' ,']',colordiffusion%JY
          WRITE ( iw, '( A,A,A,T69,F12.9 )' ) &
            ' COLOR-DIFFUSION | Preset flux along Z [', ' a.u. ' ,']',colordiffusion%JZ
       END IF
       
       WRITE ( iw, '( A,T69,A,/ )' ) ' Atom', 'Color charge'
       DO i=1,size(particles%els(:))
          WRITE ( iw, '( I6,T69,F12.9 )' ) i, particles%els(i)%ccharge
       END DO
       write(iw, '(/,/)')
    END IF
    CALL cp_print_key_finished_output(iw,logger,section,"PRINT%CD_INFO")

  END SUBROUTINE CD_info

END MODULE colordiffusion_types
