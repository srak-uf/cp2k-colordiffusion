!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> 
!> 
!> 
! **************************************************************************************************
MODULE colordiffusion_utils
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set

  USE cp_linked_list_input,            ONLY: cp_sll_val_next,&
                                             cp_sll_val_type
                                             
  USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_create,&
                                              cp_logger_set,&
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
  USE particle_types,                  ONLY: particle_type
  USE simpar_types,                    ONLY: simpar_type
  USE particle_list_types,             ONLY: particle_list_type
  USE colordiffusion_types,            ONLY: colordiffusion_type
  USE md_environment_types,            ONLY: get_md_env,&
                                              md_environment_type

#include "../base/base_uses.f90"
    
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: setup_flux

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'colordiffusion_utils'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param nemd ...
!> \param nemd_section ...
! **************************************************************************************************
  SUBROUTINE setup_flux(force_env, cdvar, md_section, rstvelflag)
    TYPE(force_env_type), POINTER                      :: force_env
    TYPE(cp_subsys_type), POINTER                      :: subsys
    TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
    TYPE(particle_list_type), POINTER                  :: particles
    TYPE(atomic_kind_type), POINTER                    :: atomic_kind
    INTEGER                                            :: natoms, i
    REAL(KIND=dp)                                      :: mass
    TYPE(colordiffusion_type), POINTER                 :: cdvar
    REAL(KIND=dp)                                      :: scale, j_cx, j_cy, j_cz, &
                                                          diff_x, diff_y, diff_z
    TYPE(section_vals_type), POINTER                   :: md_section
    TYPE(section_vals_type), POINTER                   :: print_section
    TYPE(cp_logger_type), POINTER                      :: logger
    integer                                            :: num_x, iw
    CHARACTER(LEN=default_string_length)               :: my_format
    logical                                            :: rstvelflag

    NULLIFY (atomic_kind, particles)
    NULLIFY (particle_set)
    CALL force_env_get(force_env, subsys=subsys)
    CALL cp_subsys_get(subsys,&
                       particles=particles)

    
    particle_set => particles%els
    natoms = SIZE(particle_set)

    ! get scale factor
    scale = 0
    DO i = 1, natoms
      atomic_kind => particle_set(i)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
      scale = scale + particle_set(i)%ccharge * particle_set(i)%ccharge / mass
    END DO
    scale = 1/scale

    IF (rstvelflag == .TRUE.) THEN
      j_cx = 0.0_dp
      j_cy = 0.0_dp
      j_cz = 0.0_dp
      DO i = 1, natoms
        IF (cdvar%JX > 0.0_dp) THEN
          j_cx = j_cx + particle_set(i)%v(1) * particle_set(i)%ccharge
        ELSE IF (cdvar%JY > 0.0_dp) THEN
          j_cy = j_cy + particle_set(i)%v(2) * particle_set(i)%ccharge
        ELSE IF (cdvar%JZ > 0.0_dp) THEN
          j_cz = j_cz + particle_set(i)%v(3) * particle_set(i)%ccharge
        END IF
      END DO
      IF (cdvar%JX > 0.0_dp) THEN
        diff_x = cdvar%JX - j_cx 
      ELSE IF (cdvar%JY > 0.0_dp) THEN
        diff_y = cdvar%JY - j_cy
      ELSE IF (cdvar%JZ > 0.0_dp) THEN
        diff_z = cdvar%JZ - j_cz
      END IF
      DO i = 1, natoms
        atomic_kind => particle_set(i)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
        IF (cdvar%JX > 0.0_dp) THEN
          particle_set(i)%v(1) = particle_set(i)%v(1) + scale * diff_x * particle_set(i)%ccharge/mass
        ELSE IF (cdvar%JY > 0.0_dp) THEN
          particle_set(i)%v(2) = particle_set(i)%v(2) + scale * diff_y * particle_set(i)%ccharge/mass
        ELSE IF (cdvar%JZ > 0.0_dp) THEN
          particle_set(i)%v(3) = particle_set(i)%v(3) + scale * diff_z * particle_set(i)%ccharge/mass
        END IF
      END DO
    ELSE
      DO i = 1, natoms
        atomic_kind => particle_set(i)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
        IF (cdvar%JX > 0.0_dp) THEN
          particle_set(i)%v(1) = scale*particle_set(i)%ccharge/mass*cdvar%JX
        ELSE IF (cdvar%JY > 0.0_dp) THEN
          particle_set(i)%v(2) = scale*particle_set(i)%ccharge/mass*cdvar%JY
        ELSE IF (cdvar%JZ > 0.0_dp) THEN
          particle_set(i)%v(3) = scale*particle_set(i)%ccharge/mass*cdvar%JZ
        END IF
      END DO
    END IF
    
   

  END SUBROUTINE setup_flux

END MODULE colordiffusion_utils
