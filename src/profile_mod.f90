      module profile_mod
      
!--------------------------------------------------------------------------

    use swim_global_data_mod, only : &
            & rspec, ispec, &               ! int: kind specification for real and integer
            & swim_error                    ! error routine
!            & swim_string_length, &         ! length of strings for names, files, etc.
    
      implicit none


!
!   PLASMA STATE DATA that will be given to aorsa via a namelist
!
!--------------------------------------------------------------------------
   
    !-----------------------------------
    ! Time at beginning and end of time step
    !-----------------------------------
    REAL (KIND = rspec) ::  &
        S_t0,                  &   ! time at beginning of step [msec]
      & S_t1                       ! time at end of step [msec]
   
    !-----------------------------------
    ! Basic Geometry
    !-----------------------------------
    REAL (KIND = rspec) ::  &
        S_r_axis,              & ! major radius of magnetic axis [m]
        S_z_axis,              & ! Z of magnetic axis [m]
        S_r0_mach,             & ! Z of machine center [m]
        S_z0_mach,             & ! major radius of machine center [m]
        S_r_min,               & ! major radius of inside of bounding box [m]
        S_r_max,               & ! major radius of outside of bounding box [m]
        S_z_min,               & ! Z of bottom of bounding box [m]
        S_z_max                  ! Z of top of bounding box [m]
            
    !-----------------------------------
    ! Particle Species
    !-----------------------------------
    
    INTEGER :: S_nspec       ! number of ion species = nspec_th + nspec_nonMax

    !-----------------------------------
    ! Main (thermal) Plasma Species
    !-----------------------------------m
    
    integer, parameter :: nrho_max = 181
    integer, parameter :: n_spec_th_max = 6
    integer, parameter :: n_spec_max = 7
    integer, parameter :: n_spec_nm_max = 2
    integer, parameter :: swim_string_length = 128

    
    
    INTEGER :: S_nspec_th                    ! number of thermal ion species
    character(len = 32) ::  &
        S_s_name(0:n_spec_max)              ! names of main species, (0:nspec_th)
    REAL (KIND = rspec) :: &
        S_q_s(0:n_spec_th_max),              & ! charge of species s [C], (0:nspec_th)
     &  S_m_s(0:n_spec_th_max)                 ! mass of species s [kg], (0:nspec_th)
    
    INTEGER :: S_nrho_n          ! number of rho values in thermal species density grid
    REAL (KIND = rspec) :: &
        S_rho_n_grid(nrho_max),       & ! rho values in density grid, (1:nrho_n)
     &  S_n_s(nrho_max, 0:n_spec_th_max),           & ! density profile of species s, (1:nrho_n, 0:nspec_th)
     &  S_zeff(nrho_max),           & ! effective impurity charge profile, (1:nrho_n)
     &  S_m_impurity(nrho_max)              ! effective impurity mass profile, (1:nrho_n)
 
    INTEGER :: S_nrho_T          ! number of rho values in temperature grid
    REAL (KIND = rspec) :: &
        S_rho_T_grid(nrho_max),       & ! rho values in temperature grid, (1:nrho_T)
      & S_T_s(nrho_max, 0:n_spec_th_max)              ! Temperature profile of species s, (1:nrho_T, 0:nspec_th)
 
    INTEGER :: S_nrho_v_par      ! number of main rho values in parallel velocity grid
!    REAL (KIND = rspec), ALLOCATABLE :: &
!        PS_rho_v_par_grid(:),   & ! rho values in parallel velocity grid, (1:nrho_v_par)
!      & PS_v_par_s(:, :)        & ! v parallel profile of species s, 
                                  ! (1:nrho_v_par, 0:nspec_th)
 
    !-----------------------------------
    ! Non-Maxwellian Species
    !-----------------------------------
    
    INTEGER :: S_nspec_nonMax    ! number of non-Maxwellian species
    character(len=32), dimension(n_spec_nm_max ) :: &
        S_nonMax_name         ! names of non-Maxwellian species, (1:nspec_nonMax)
    
    REAL (KIND = rspec), dimension(n_spec_nm_max ) :: &
        S_q_nonMax_s,       & ! charge of species s [C], (1:nspec_nonMax)
        S_m_nonMaX_s          ! mass of species s [kg], (1:nspec_nonMax)
    
    INTEGER :: S_ntheta_n        ! number of theta values in 2D density grid

    REAL (KIND = rspec), ALLOCATABLE :: &
        S_n_nonMax2D_s(:, :,:)  ! 2D density profile of non-Maxwellian species s,
                                 ! (1:nrho_n, 1:ntheta_n, 1:nspec_nonMax)

    REAL (KIND = rspec), ALLOCATABLE :: &
        S_n_nonMax_s(:, :)   ! Flux surface average density profile of 
                              ! non-Maxwellian species s, (1:nrho_n, 1:nspec_nonMax)
 
    character(len = swim_string_length) :: &
        S_dist_fun_s         ! distribution function of non-Maxwellian  
                              ! species s, (1:nspec_nonMax) N.B. For now a distribution
                              ! function type is a file name
 
 
    !-----------------------------------
    ! Magnetics
    !
    ! magnetics: B(x), magnetic field.  Like distribution_fn there is a
    ! user-defined type that contains a file with the data.
    ! AORSA and TORIC get all their magntics data by reading an eqdisk file.
    ! Eventually all the magnetics data will appear separately in the Plasma
    ! state.
    !-----------------------------------

    character(len = swim_string_length) :: S_eqdsk_file   ! eqdisk file
        
    REAL (KIND = rspec) ::  &
        S_B_axis               ! Field at magnetic axis [T]

    !--------------------------------------------------------------------------
    !
    ! RF Data
    !
    ! Allow multiple RF sources, ICF, ICRH. So RF_frequency, power, etc may
    ! not be a scalar in that case.
    !   Assumption 1: the RF component invocation will
    !       involve a loop over each RF source, each of which can have its own
    !       RF_frequency, etc.  (vs. adding another dimension to the arrays).
    !   Assumption 2 : each source involves invoking another executable.
    !--------------------------------------------------------------------------
   
    
    INTEGER :: S_nrf_src         ! number of RF sources
       ! names of rf sources, (1:nrf_src)
        

    
 
    character(len = swim_string_length) :: S_ant_model_src !file name for antenna model 
    !---------------------------------------------------------------------------
    ! Note:
    ! Antenna model is currently defined in a file. The PREPARE_CODE_INPUT program
    ! should extract from it the data needed to define the geometry and operation
    ! i.e. phasing or mode number spectrum
    !
    ! For these, see the example namelist attached to the end
    ! of the aorsa.doc file sent by Fred to Swim list.  For now the data in this
    ! file includes:
    !   nphi =toroidal mode number (find another name not conflicting toroidal angle)
    !   antlen = vertical height of antenna [m]
    !   dpsiant0 = radial thickness of antenna in rho
    !   rant = radial location of antenna in major radius [m]
    !   yant = vertical location of antenna center [m]
    !
    !  N.B. In this scheme the toroidal mode number comes in through the antenna model
    !  The antenna geometry should eventually come from one of the standard machine
    !  definition files. 
    !
    !---------------------------------------------------------------------------
  
  
    ! RF Outputs that go back into the Plasma State.  Profiles are flux surface averages.
    
    ! N.B. We will want to put in 2D power deposition profiles, but I don't think they
    ! are needed for our initial coupling
        
    INTEGER ::  &
        S_nrho_prf,    &   ! number rho values for RF power deposition grid
        S_ntheta_prf   ! number of theta values in 2D RF power dep grid
        
    REAL (KIND = rspec), ALLOCATABLE :: &
        S_rho_prf_grid(:), & ! rho values in RF power deposition grid, (1:nrho__prf)
        S_prf2D_src_s(:,:,:,:),   & ! 2D Power deposition from each source into each 
                                                        ! species, (1:nrho__prf, 1:nrf_src, 0:nspec)
        S_prf_src_s(:,:,:)  ! Power deposition profile from each source into each 
                              ! species, (1:nrho__prf, 1:nrf_src, 0:nspec)
    real(kind = rspec) ::   S_prf_total_s(nrho_max,0:n_spec_max)   ! Total rf power deposition profile into each species
                              ! summed over sources, (1:nrho__prf, 0:nspec)
    
    integer :: S_nrho_cdrf(n_spec_max) ! # of rho values for RF current drive grid for each species
        
    REAL (KIND = rspec), ALLOCATABLE :: &
        S_rho_cdrf_grid(:),    & ! rho values in RF current drive grid, (1:nrho__cdrf)
        S_cdrf_src_s(:,:,:),   & ! Driven current profile from each source, in each species
                                  ! (1:nrho__cdrf, 1:nrf_src, 1:nspec_nonMax)
        S_cdrf_total_s(:,:)      ! Total current driven by all sources in each species
    
 
    character(len = swim_string_length), dimension(n_spec_nm_max)  :: &
        S_ql_operator          ! quasilinear operator for each non Maxwellian--file name 
                                ! species, (1:nspec_nonMax)
    character(len = swim_string_length), dimension(n_spec_nm_max) :: &
        S_distribution_fun      ! distribution function for each non Maxwellian species



    namelist /state/ S_t0, S_t1, S_r_axis, S_z_axis, S_r0_mach, &
       S_z0_mach, S_r_min, S_r_max, S_z_min, S_z_max, &
       S_nspec,S_nspec_th,S_s_name, S_q_s, S_m_s,  &
       S_nrho_n, S_rho_n_grid, S_n_s, S_zeff, S_m_impurity, &
       S_nrho_T, S_rho_T_grid   , S_T_S , S_ant_model_src, S_ql_operator, &
       S_distribution_fun     


    end module profile_mod
      

