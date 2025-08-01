!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: isotope_mod.F90
!
! !DESCRIPTION: Module ISOTOPE_MOD contains variables and routines
! for simulating CH4, CO, CO2, and OCS with an online calculation of the
! chemistry between them using KPP. It was adapted directly from
! the module CH4_CO_CO2_MOD.F provided by Beata Bukosa.
!\\
!\\
! !INTERFACE:
!
MODULE Isotope_Mod
!
! !USES:
!
  USE Error_Mod,     ONLY : Safe_Div
  USE Hco_Error_Mod, ONLY : HCO_SUCCESS, HCO_FAIL, HCO_WARNING, hp
  USE PhysConstants
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Emiss_Isotope
  PUBLIC :: Chem_Isotope
  PUBLIC :: Init_Isotope
  PUBLIC :: Cleanup_Isotope
!
! !REVISION HISTORY:
!  04 Apr 2022 - M.S. Long   - Initial version, based on work by B. Bukosa
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!

  ! Scalars
  INTEGER               :: id_CH4,     id_CH4_adv
  INTEGER               :: id_12CH4,   id_12CH4_adv
  INTEGER               :: id_13CH4,   id_13CH4_adv
  INTEGER               :: id_14CH4,   id_14CH4_adv
  INTEGER               :: id_12CH3D,  id_12CH3D_adv
  INTEGER               :: id_13CH3D,  id_13CH3D_adv
  INTEGER               :: id_14CO,    id_14CO_adv
  INTEGER               :: id_SF6,     id_SF6_adv
  INTEGER               :: id_OH,      id_O1D,       id_Cl

  ! Arrays
  REAL(fp), ALLOCATABLE :: sumOfCosSza(:,:)
!
! !DEFINED PARAMETERS:
!
  INTEGER,  PARAMETER   :: N_CH4_DIAGS = 16
  REAL(fp), PARAMETER   :: CM2perM2    = 1.0e+4_fp
  REAL(fp), PARAMETER   :: CM3perM3    = 1.0e+6_fp
  REAL(fp), PARAMETER   :: toMolecCm3  = ( AVO / AIRMW ) * 1.0e-3_fp

  ! 13C reference standard for Pee-Dee Belemnite [atom atom-1] [Craig, 1957]
  REAL(fp), PARAMETER   :: Rs_13C_VPDB = 0.01123720e0_fp

  ! D reference standard for Vienna Mean Standard Ocean Water [atom atom-1]
  REAL(fp), PARAMETER   :: Rs_D_VSMOW  = 1.5576e-4_fp

  ! Radioactive decay constant for 14C [s-1] = (8267 a)^-1
  REAL(fp), PARAMETER   :: lambda_14C  = 3.833082e-12_fp

  ! Absolute international standard activity defined for 1950 AD
  ! in 0.95 NBS oxalic acid [Bq/kgC] [Stuiver, 1980]
  REAL(fp), PARAMETER   :: A_abs       = 226_fp
  
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_isotope
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Emiss_Isotope( Input_Opt,  State_Chm, State_Diag,            &
                                 State_Grid, State_Met, RC                    )
!
! !USES:
!
    USE HCO_State_Mod,        ONLY : Hco_GetHcoId
    USE HCO_State_GC_Mod,     ONLY : HcoState
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_GetDiagn
    USE ErrCode_Mod
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Species_Mod,          ONLY : SpcConc
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    
    ! Strings
    CHARACTER(LEN=255)     :: thisLoc
    CHARACTER(LEN=512)     :: errMsg

    ! Arrays

    ! String arrays

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)
    REAL(f4),      POINTER :: Ptr2D(:,:)
!
! !DEFINED PARAMETERS:
!

    !========================================================================
    ! Emiss_Isotope begins here!
    !========================================================================

    ! Initialize
    RC       =  GC_SUCCESS
    Ptr2D    => NULL()
    Spc      => NULL()
    errMsg   =  ''
    thisLoc  =  &
     ' -> at Emiss_Isotope (in module GeosCore/isotope_mod.F90)'

    ! Exit with error if we can't find the HEMCO state object
    IF ( .NOT. ASSOCIATED( HcoState ) ) THEN
       errMsg = 'The HcoState object is undefined!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Free pointers for safety's sake
    Spc   => NULL()
    Ptr2D => NULL()

  END SUBROUTINE Emiss_Isotope
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_isotope
!
! !DESCRIPTION: Computes the chemical loss of carbon species (sources - sinks)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Chem_Isotope( Input_Opt,  State_Met,  State_Chm,            &
                                State_Grid, State_Diag, RC                   )
!
! !USES:
!
    USE carbon_Funcs
    USE gckpp_Global
    USE gckpp_Integrator,     ONLY : Integrate
    USE gckpp_Parameters
    USE gckpp_Precision
    USE gckpp_Rates,          ONLY : Update_Rconst
    USE ErrCode_Mod
    USE HCO_State_Mod,        ONLY : Hco_GetHcoId
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_HcoStateOK
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE Input_Opt_Mod,        ONLY : OptInput
    USE rateLawUtilFuncs,     ONLY : SafeDiv
    USE Species_Mod,          ONLY : SpcConc
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Met_Mod,        ONLY : MetState
    USE Time_Mod,             ONLY : Get_Ts_Chem, Get_Year
    USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
    !
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  CH4 SOURCES
!  ============================================================================
!  (1 ) Oxidation of methane, isoprene and monoterpenes (SRCO_fromHCs).
!  (2 ) Direct emissions of CO from fossil fuel combustion, biomass
!        burning and wood (for fuel) burning (SR SETEMIS).
!  (3 ) Emissions.
!                                                                             .
!  CH4 SINKS:
!  ============================================================================
!  (1 ) Removal of CO by OH (SR OHparam & CO_decay).
!  (2 ) CO uptake by soils (neglected).
!  (3 ) Transport of CO to stratosphere from troposphere
!        (in dynamical subroutines).
!  (4 ) Removal by OH (Clarissa's OH--climatol_OH.f and CO_decay.f)
!  (5 ) Transport of CH4 between troposphere and stratosphere, and
!        destruction in strat (CH4_strat.f).
!  (6 ) Removel by Cl
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL                :: first = .TRUE. 
   
    ! Scalars
    LOGICAL                :: failed,   found
    INTEGER                :: HcoID,    I,     SpcID,   KppID
    INTEGER                :: J,        L
    INTEGER                :: NA,       N
    INTEGER                :: IERR,     previous_units
    REAL(fp)               :: dtChem,   facDiurnal
    REAL(fp)               :: tsPerDay
    REAL(FP)               :: d13,    dD,   atomsH,   atomsD,    Rs, Rs_corr
    REAL(FP)               :: pMC,   A_S,     A_SN,  A_abs_c,  D14C
    REAL(FP)               :: gC

    ! Strings
    CHARACTER(LEN=63)      :: dgnName
    CHARACTER(LEN=512)     :: errMsg
    CHARACTER(LEN=255)     :: thisLoc

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    ! Arrays
    INTEGER                :: ICNTRL(20)
    INTEGER                :: ISTATUS(20)
    REAL(dp)               :: RCNTRL(20)
    REAL(dp)               :: RSTATE(20)
    REAL(fp)               :: OHdiurnalFac(State_Grid%NX, State_Grid%NY)

    ! Arrays for data read in via HEMCO
    REAL(fp) :: Global_OH(   State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: Global_Cl(   State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: Global_O1D(  State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: lsoil(       State_Grid%NX, State_Grid%NY               )
    REAL(fp) :: ksoil(       State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    
    !========================================================================
    ! Chem_Isotope begins here!
    !========================================================================

    ! Initialize
    RC       =  GC_SUCCESS
    dtChem   =  Get_Ts_Chem()
    tsPerDay =  86400.0_fp / dtChem
    Spc      => State_Chm%Species
    errMsg   = ''
    thisLoc  = &
     ' -> at Chem_Isotope (in module GeosCore/isotope_mod.F90)'

    ! First-time only safety checks
    IF ( first ) THEN

       ! Exit if the HEMCO state object has not yet been initialized
       IF ( .not. HCO_GC_HcoStateOK() ) THEN
          errMsg = 'HcoState object is not associated!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Reset first-time flag
       first = .FALSE.
    ENDIF

    !========================================================================
    ! Zero diagnostic archival arrays to make sure that we don't have any
    ! leftover values from the last timestep near the top of the chemgrid.
    !========================================================================
    IF ( State_Diag%Archive_D13CCH4          ) THEN
       State_Diag%D13CCH4 = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_D2HCH4           ) THEN
       State_Diag%D2HCH4 = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_PMCCH4           ) THEN
       State_Diag%PMCCH4 = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_D14CH4           ) THEN
       State_Diag%D14CH4 = 0.0_f4
    ENDIF    
    IF ( State_Diag%Archive_ClconcAfterChem  ) THEN
       State_Diag%ClconcAfterChem  = 0.0_f4
    ENDIF    
    IF ( State_Diag%Archive_OHconcAfterChem  ) THEN
       State_Diag%OHconcAfterChem  = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_O1DconcAfterChem ) THEN
       State_Diag%O1DconcAfterChem = 0.0_f4
    ENDIF

    !========================================================================
    ! Read chemical inputs (oxidant fields and soil uptake) via HEMCO
    !========================================================================

    !------------------------------------------------------------------------
    ! OH concentration:
    !------------------------------------------------------------------------
    DgnName = 'GLOBAL_OH'
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName,                     &
                         Global_OH, RC,         found=found                 )
    IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
       errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Convert orignal units [mol/mol dry air] to [molec/cm3] ! Already molec/cm3
    !Global_OH = ( Global_OH * State_Met%AirDen ) * toMolecCm3
    
    !------------------------------------------------------------------------
    ! Cl concentration:
    !------------------------------------------------------------------------
    DgnName = 'GLOBAL_Cl'
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName,                     &
                         Global_Cl, RC,         found=found                 )
    IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
       errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Convert orignal units [mol/mol dry air] to [molec/cm3] ! Already molec/cm3
    !Global_Cl = ( Global_Cl * State_Met%AirDen ) * toMolecCm3

    !------------------------------------------------------------------------
    ! O(1D) concentration:
    !------------------------------------------------------------------------
    DgnName = 'GLOBAL_O1D'
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName,                     &
                         Global_O1D, RC,        found=found                 )
    IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
       errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Convert orignal units [mol/mol dry air] to [molec/cm3] ! Already molec/cm3
    !Global_O1D = ( Global_O1D * State_Met%AirDen ) * toMolecCm3

    !------------------------------------------------------------------------
    ! Soil uptake
    !------------------------------------------------------------------------
    lsoil = 0d0
    ksoil = 0d0
    DgnName = 'CH4_SOILABSORB'
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'CH4_SOILABSORB',     &
                         lsoil,     RC,        found=found           )
    IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
       errMsg = 'Cannot get pointer to HEMCO field CH4_SOILABSORB'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Convert LSOIL [kg m-2 s-1] to KSOIL [s-1]
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       ! Convert from kg/m2/s to 1/s using total methane (still in units of kg here)
       kSOIL(I,J,1) = LSOIL(I,J) * State_Grid%Area_M2(I,J) / State_Chm%Species(id_CH4)%Conc(I,J,1)
    ENDDO
    ENDDO
    
    !========================================================================
    ! Compute OH diurnal cycle scaling factor
    ! (this scales OH by the position of the sun, and zeroes it at night)
    !========================================================================
    CALL Calc_Diurnal(                                                       &
         State_Grid   = State_Grid,                                          &
         State_Met    = State_Met,                                           &
         OHdiurnalFac = OHdiurnalFac                                        )

    !========================================================================
    ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
    !========================================================================

    IF ( State_Diag%Archive_OHconcAfterChem ) THEN
       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( I, J, L                                               )&
       !$OMP COLLAPSE( 3                                                    )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Archive OH if we are in the chemistry grid [molec/cm3]
          IF ( State_Met%InChemGrid(I,J,L) ) THEN
             IF ( State_Diag%Archive_OHconcAfterChem ) THEN
                State_Diag%OHconcAfterChem(I,J,L) = Global_OH(I,J,L)         &
                                                  * OHdiurnalFac(I,J)
             ENDIF
          ENDIF

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !========================================================================
    ! Convert species to [molec/cm3] (ewl, 8/16/16)
    !========================================================================

    ! Convert units of all species to molec/cm3 for KPP
    CALL Convert_Spc_Units(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         new_units      = MOLECULES_SPECIES_PER_CM3,                         &
         previous_units = previous_units,                                    &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'fullchem_mod.F90')
       RETURN
    ENDIF
    
    !========================================================================
    ! Main chemistry loop -- call KPP to integrate the mechanism forward
    !========================================================================

    ! KPP forward-Euler integrator settings
    ICNTRL     =  0
    ICNTRL(1)  =  1   ! Verbose error output
    ICNTRL(2)  =  0   ! Stop model on negative values
    ICNTRL(15) = -1   ! Do not call Update_SUN, Update_RCONST w/in integrator

    ! Set a flag to denote if the chemistry failed
    failed     = .FALSE.

    ! Loop over grid boxes
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J, L, N                                               )&
    !$OMP PRIVATE( SpcID, KppID                                             )&
    !$OMP COLLAPSE( 3                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 24                                             )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize PRIVATE and THREADPRIVATE loop variables
       C              = 0.0_dp                    ! Species conc. [molec/cm3]
       CFACTOR        = 1.0_dp                    ! Not used, set = 1
       NUMDEN         = State_Met%AIRNUMDEN(I,J,L)! Air density [molec/cm3]
       TEMP           = State_Met%T(I,J,L)        ! Temperature [K]
       INV_TEMP       = 1.0_dp / TEMP             ! 1/T  term for equations
       TEMP_OVER_K300 = TEMP / 300.0_dp           ! T/300 term for equations
       K300_OVER_TEMP = 300.0_dp / TEMP           ! 300/T term for equations
       SUNCOS         = State_Met%SUNCOSmid(I,J)  ! Cos(SZA) ) [1]

       !=====================================================================
       ! Initialize the KPP "C" vector of species concentrations [molec/cm3]
       !=====================================================================

       ! Isotopologues [molec/cm3]
       C(ind_CH4   ) = State_Chm%Species(id_CH4   )%Conc(I,J,L)
       C(ind_C12H4 ) = State_Chm%Species(id_12CH4 )%Conc(I,J,L)
       C(ind_C13H4 ) = State_Chm%Species(id_13CH4 )%Conc(I,J,L)
       C(ind_C14H4 ) = State_Chm%Species(id_14CH4 )%Conc(I,J,L)
       C(ind_C12H3D) = State_Chm%Species(id_12CH3D)%Conc(I,J,L)
       C(ind_C13H3D) = State_Chm%Species(id_13CH3D)%Conc(I,J,L)
       C(ind_C14O  ) = State_Chm%Species(id_14CO  )%Conc(I,J,L)
       
       ! Oxidant concentrations [molec/cm3]
       C(ind_FixedOH)  = Global_OH( I,J,L) * OHdiurnalFac(I,J)
       C(ind_FixedCl)  = Global_Cl( I,J,L)
       C(ind_FixedO1D) = Global_O1D(I,J,L)
       
       ! Update the array of rate constants for the KPP solver
       CALL Update_RCONST()

       !=====================================================================
       ! Call the KPP integrator
       !=====================================================================

       ! Integrate the mechanism forward in time
       CALL Integrate(                                                       &
            TIN      = 0.0_dp,                                               &
            TOUT     = dtChem,                                               &
            ICNTRL_U = ICNTRL,                                               &
            IERR_U   = IERR                                                 )

       ! Trap potential errors
       IF ( IERR /= 1 ) failed = .TRUE.

       ! Copy concentrations back into State_Chm
       State_Chm%Species(id_CH4   )%Conc(I,J,L) = C(ind_CH4   )
       State_Chm%Species(id_12CH4 )%Conc(I,J,L) = C(ind_C12H4 )
       State_Chm%Species(id_13CH4 )%Conc(I,J,L) = C(ind_C13H4 )
       State_Chm%Species(id_14CH4 )%Conc(I,J,L) = C(ind_C14H4 )
       State_Chm%Species(id_12CH3D)%Conc(I,J,L) = C(ind_C12H3D)
       State_Chm%Species(id_13CH3D)%Conc(I,J,L) = C(ind_C13H3D)
       State_Chm%Species(id_14CO  )%Conc(I,J,L) = C(ind_C14O  )

       !=====================================================================
       ! HISTORY: Archive KPP solver diagnostics
       !=====================================================================
       IF ( State_Diag%Archive_KppDiags ) THEN

          ! # of integrator calls
          IF ( State_Diag%Archive_KppIntCounts ) THEN
             State_Diag%KppIntCounts(I,J,L) = ISTATUS(1)
          ENDIF

          ! # of times Jacobian was constructed
          IF ( State_Diag%Archive_KppJacCounts ) THEN
             State_Diag%KppJacCounts(I,J,L) = ISTATUS(2)
          ENDIF

          ! # of internal timesteps
          IF ( State_Diag%Archive_KppTotSteps ) THEN
             State_Diag%KppTotSteps(I,J,L) = ISTATUS(3)
          ENDIF

          ! # of accepted internal timesteps
          IF ( State_Diag%Archive_KppTotSteps ) THEN
             State_Diag%KppAccSteps(I,J,L) = ISTATUS(4)
          ENDIF

          ! # of rejected internal timesteps
          IF ( State_Diag%Archive_KppTotSteps ) THEN
             State_Diag%KppRejSteps(I,J,L) = ISTATUS(5)
          ENDIF

          ! # of LU-decompositions
          IF ( State_Diag%Archive_KppLuDecomps ) THEN
             State_Diag%KppLuDecomps(I,J,L) = ISTATUS(6)
          ENDIF

          ! # of forward and backwards substitutions
          IF ( State_Diag%Archive_KppSubsts ) THEN
             State_Diag%KppSubsts(I,J,L) = ISTATUS(7)
          ENDIF

          ! # of singular-matrix decompositions
          IF ( State_Diag%Archive_KppSmDecomps ) THEN
             State_Diag%KppSmDecomps(I,J,L) = ISTATUS(8)
          ENDIF
       ENDIF
       
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    IF ( failed ) THEN
       errMsg = 'KPP integration failed!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

!$OMP PARALLEL DO                                                     &
!$OMP DEFAULT( SHARED                                                )&
!$OMP PRIVATE( I, J, L, Rs, Rs_corr, d13, dD, atomsH, atomsD         )&
!$OMP PRIVATE( A_S, A_SN, A_abs_c, D14C, pMC, gC                     )&
!$OMP COLLAPSE( 3                                                    )
      DO L = 1, State_Grid%NZ
      DO J = 1, State_Grid%NY
      DO I = 1, State_Grid%NX

         d13  = 0d0
         dD   = 0d0
         pMC  = 0d0
         D14C = 0d0
         
         !---------------------------------------------------------------
         ! d13C of Methane (permil)
         !---------------------------------------------------------------
         IF ( State_Diag%Archive_d13CCH4 ) THEN
       
            ! Calculate molecular ratio of 13C to 12C
            Rs = ( Spc(id_13CH4)%Conc( I, J, L ) + Spc(id_13CH3D)%Conc( I, J, L ) ) / &
                 ( Spc(id_12CH4)%Conc( I, J, L ) + Spc(id_12CH3D)%Conc( I, J, L ) )
            
            ! Compare to Pee-Dee Belemnite Standard (per mil)
            d13 = ( ( Rs / Rs_13C_VPDB ) - 1e0_fp ) * 1000e0_fp
            
            State_Diag%d13CCH4(I,J,L) = d13

         ENDIF

         !---------------------------------------------------------------
         ! dD of Methane (permil)
         !---------------------------------------------------------------
         IF ( State_Diag%Archive_d2HCH4 ) THEN
       
            ! Calculate atomic ratio of D to H
            atomsD =     Spc(id_12CH3D)%Conc( I, J, L ) + &
                         Spc(id_13CH3D)%Conc( I, J, L )
            atomsH = 3 * Spc(id_12CH3D)%Conc( I, J, L ) + &
                     3 * Spc(id_13CH3D)%Conc( I, J, L ) + &
                     4 * Spc(id_12CH4 )%Conc( I, J, L ) + &
                     4 * Spc(id_13CH4 )%Conc( I, J, L )
             
            Rs = atomsD / atomsH
             
            ! Compare to Vienna Standard Mean Ocean Water (per mil)
            dD = ( ( Rs / Rs_D_VSMOW ) - 1e0_fp ) * 1000e0_fp
             
            ! Archive
            State_Diag%d2HCH4(I,J,L) = dD

         ENDIF

         !---------------------------------------------------------------
         ! Percent Modern Carbon (%)
         !---------------------------------------------------------------
         IF ( State_Diag%Archive_pMCCH4 ) THEN
 
            ! Normalize to -25‰ (preindustrial d13) to account for
            ! atmospheric chemistry between 1950 standard and present-day
            IF ( d13 .eq. 0 ) THEN
               Rs  = ( Spc(id_13CH4)%Conc( I, J, L ) + Spc(id_13CH3D)%Conc( I, J, L ) ) / &
                     ( Spc(id_12CH4)%Conc( I, J, L ) + Spc(id_12CH3D)%Conc( I, J, L ) )
               d13 = ( ( Rs / Rs_13C_VPDB ) - 1e0_fp ) * 1000e0_fp
            ENDIF
 
            ! Calculate total g C of methane (using NIST values) (molec cm-3 -> g(C))
            gC = ( Spc(id_12CH4)%Conc(  I, J, L ) * 12.0d0           * State_Met%AIRVOL( I, J, L )*1d6 / AVO ) + &
                 ( Spc(id_13CH4)%Conc(  I, J, L ) * 13.00335483507d0 * State_Met%AIRVOL( I, J, L )*1d6 / AVO ) + &
                 ( Spc(id_12CH3D)%Conc( I, J, L ) * 12.0d0           * State_Met%AIRVOL( I, J, L )*1d6 / AVO ) + &
                 ( Spc(id_13CH3D)%Conc( I, J, L ) * 13.00335483507d0 * State_Met%AIRVOL( I, J, L )*1d6 / AVO ) + &
                 ( Spc(id_14CH4)%Conc(  I, J, L ) * 14.0032419884d0  * State_Met%AIRVOL( I, J, L )*1d6 / AVO ) 
            
            ! Calculate 14C activity (molec cm-3 -> Bq g(C)-1)
            A_S =   Spc(id_14CH4)%Conc( I, J, L ) * State_Met%AIRVOL( I, J, L ) * 1d6 * lambda_14C / gC
 
            ! Apply d13 normalization
            A_SN = A_S * ( 0.975d0 / ( 1d0 + (d13/1d3) ) )**2d0
 
            ! Calculate D14CH4 (permil deviation)
            D14C = ( ( A_SN / 0.2260d0 ) - 1d0 ) * 1000d0
            
            ! Archive
            State_Diag%D14CH4(I,J,L) = D14C
 
            ! pMC for simulation year
            pMC = 100d0 * ( 1d0 + (D14C/1d3) ) / ( EXP(lambda_14c*(1950-GET_YEAR()) ) )
 
            ! Archive
            State_Diag%pMCCH4(I,J,L) = pMC
            
         ENDIF
         
!         IF ( I .eq. 1 .and. J .eq. 1 .and. L .eq. 1 .and. Input_Opt%amIroot ) THEN
!100         FORMAT( 'ROC: CH4 :', F7.1, ' ppb | δ13-CH4: ', F7.1, ' ‰    |   δD-CH4: ', F7.1, ' ‰   |   pMC-CH4: ', F7.2, ' %  |  D14CH4: ',F7.2, '  ‰'  )
!            WRITE(6,100) 1e9*Spc(id_CH4)%Conc(I,J,L)/State_Met%AIRNUMDEN(I,J,L),  d13, dD, pMC, D14C
!         ENDIF
       
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
    
    !=======================================================================
    ! Convert species back to original units
    !=======================================================================

    ! Convert units of all species back to kg
    CALL Convert_Spc_Units(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         new_units  = previous_units,                                        &
         RC         = RC                                                    )
    
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'fullchem_mod.F90' )
       RETURN
    ENDIF
        
    ! Free pointers for safety's sake
    Spc => NULL()

  END SUBROUTINE Chem_Isotope
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_isotope
!
! !DESCRIPTION: Allocates and zeroes module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Isotope( Input_Opt,  State_Chm, State_Diag,             &
                                State_Grid, RC                                )
!
! !USES:
!
    USE gckpp_Global,   ONLY : MW, SR_MW, HENRY_CR, HENRY_K0
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState, Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine is called from GC_INIT_EXTRA (in GeosCore/input_mod.f)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: KppId,  N

    ! Strings
    CHARACTER(LEN=255) :: errMsg
    CHARACTER(LEN=255) :: thisLoc

    !========================================================================
    ! Initialize
    !========================================================================
    RC         = GC_SUCCESS
    errMsg     = ''
    thisLoc    = &
     ' -> at Init_Isotope (in module GeosCore/isotope_mod.F90)'

    !========================================================================
    ! Define GEOS-Chem species indices
    !
    ! Also denote which species are advected to facilitate single-tracer runs
    !========================================================================
    id_CH4        = Ind_( 'CH4'        )
    id_CH4_adv    = Ind_( 'CH4',    'A' )
    id_12CH4      = Ind_( 'C12H4'       )
    id_12CH4_adv  = Ind_( 'C12H4',  'A' )
    id_13CH4      = Ind_( 'C13H4'       )
    id_13CH4_adv  = Ind_( 'C13H4',  'A' )
    id_14CH4      = Ind_( 'C14H4'       )
    id_14CH4_adv  = Ind_( 'C14H4',  'A' )
    id_12CH3D     = Ind_( 'C12H3D'      )
    id_12CH3D_adv = Ind_( 'C12H3D', 'A' )    
    id_13CH3D     = Ind_( 'C13H3D'      )
    id_13CH3D_adv = Ind_( 'C13H3D', 'A' )
    id_14CO       = Ind_( 'C14O'        )
    id_14CO_adv   = Ind_( 'C14O',   'A' )
    id_SF6        = Ind_( 'SF6'         )
    id_SF6_adv    = Ind_( 'SF6',    'A' )
    id_OH         = Ind_( 'FixedOH'     )
    id_O1D        = Ind_( 'FixedO1D'    )
    id_Cl         = Ind_( 'FixedCl'     )
    
    !========================================================================
    ! Save physical parameters from the species_database.yml file into KPP
    ! arrays located in module gckpp_Global.F90.  These do not vary with
    ! (I,J,L) location, and so can be defined here in the init phase.
    !========================================================================
    DO KppId = 1, State_Chm%nKppSpc + State_Chm%nOmitted
       N                  = State_Chm%Map_KppSpc(KppId)
       IF ( N > 0 ) THEN
          MW(KppId)       = State_Chm%SpcData(N)%Info%MW_g
          SR_MW(KppId)    = SQRT( MW(KppId ) )
          HENRY_K0(KppId) = State_Chm%SpcData(N)%Info%Henry_K0
          HENRY_CR(KppId) = State_Chm%SpcData(N)%Info%Henry_CR
       ENDIF
    ENDDO

    !========================================================================
    ! Initialize variables OH
    !========================================================================
    ALLOCATE( sumOfCosSza( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'isotope_mod.F90:sumOfCosSza', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    sumOfCosSza = 0.0_fp
    
  END SUBROUTINE Init_Isotope
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_isotope
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_CH4 deallocates module arrays.
!  (bmy, 1/16/01)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Isotope( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! Cleanup_Isotope begins here!
    !=================================================================

    ! Initialize
    RC = GC_SUCCESS

    IF ( ALLOCATED( sumOfCosSza ) ) THEN
       DEALLOCATE( sumOfCosSza, STAT=RC )
       CALL GC_CheckVar( 'isotope_mod.F90:sumOfCosSza', 2, RC )
       RETURN
    ENDIF

  END SUBROUTINE Cleanup_Isotope
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_diurnal
!
! !DESCRIPTION: Subroutine CALC\_DIRUNAL computes the sume of the cosine
!  of the solar zenith angle over a 24 hour day as well as the total
!  length of daylight to scale the offline OH concentrations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_Diurnal( State_Grid, State_Met, OHdiurnalFac )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE TIME_MOD,       ONLY : ITS_A_NEW_DAY
    USE TIME_MOD,       ONLY : GET_MINUTE,    GET_SECOND,      GET_HOUR
    USE TIME_MOD,       ONLY : GET_TS_CHEM,   GET_DAY_OF_YEAR, GET_LOCALTIME
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid        ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met         ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT) :: OHdiurnalFac(                             &
                                    State_Grid%NX,                           &
                                    State_Grid%NY)   ! OH diurnal scaling [1]
!
! !REVISION HISTORY:
!  12 Mar 2014 - J. Fisher - Initial version, Copied from OHNO3TIME in
!                            carbon_mod and COSSZA in dao_mod
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Scalars
    INTEGER            :: I, J, N, NDYSTEP
    INTEGER            :: SECOND,  MINUTE, TS_SUN
    REAL*8             :: GMT_MID, TIMLOC, FACTOR
    REAL*8             :: R,       AHR,    DEC
    REAL*8             :: YMID_R,  SUNTMP_MID
    REAL*8             :: dtChem,  timestepsPerDay
!
! !DEFINED PARAMETERS:
!
    ! Coefficients for solar declination angle
    REAL*8,  PARAMETER :: A0 = 0.006918d0
    REAL*8,  PARAMETER :: A1 = 0.399912d0
    REAL*8,  PARAMETER :: A2 = 0.006758d0
    REAL*8,  PARAMETER :: A3 = 0.002697d0
    REAL*8,  PARAMETER :: B1 = 0.070257d0
    REAL*8,  PARAMETER :: B2 = 0.000907d0
    REAL*8,  PARAMETER :: B3 = 0.000148d0

    !=================================================================
    ! CALC_DIURNAL begins here!
    !=================================================================

    ! Only do at the start of a new day
    IF ( FIRST .or. ITS_A_NEW_DAY() ) THEN

       ! Zero array
       sumOfCosSza = 0d0

       ! Get time for central chemistry timestep
       TS_SUN = GET_TS_CHEM()                     ! Chemistry interval
       SECOND = GET_SECOND()                      ! Current seconds
       MINUTE = GET_MINUTE()                      ! Current minutes
       FACTOR = ( MINUTE * 60 + SECOND ) / TS_SUN ! Multiplying factor

       ! GMT at the midpoint of the chemistry time interval for first
       ! timestep of the day
       GMT_MID  = ( DBLE( GET_HOUR()        )        ) &
                + ( DBLE( TS_SUN * FACTOR ) / 3600d0 ) &
                + ( DBLE( TS_SUN / 2      ) / 3600d0 )

       ! Solar declination angle (low precision formula):
       ! Path length of earth's orbit traversed since Jan 1 [radians]
       R = ( 2d0 * PI / 365d0 ) * FLOAT( GET_DAY_OF_YEAR() - 1 )
       DEC = A0 - A1*COS(    R) + B1*SIN(    R) &
                - A2*COS(2d0*R) + B2*SIN(2d0*R) &
                - A3*COS(3d0*R) + B3*SIN(3d0*R)

       ! NDYSTEP is # of chemistry time steps
       NDYSTEP = INT( 24d0 * 3600d0 / GET_TS_CHEM() )

       ! Loop forward through NDYSTEP "fake" timesteps for this day
       DO N = 1, NDYSTEP

          ! Increment GMT (hours) to midpoint of next timestep
          IF ( N > 1 ) GMT_MID = GMT_MID + TS_SUN / 3600d0

          ! Loop over surface grid boxes
          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J, YMID_R, TIMLOC, AHR, SUNTMP_MID              )&
          !$OMP COLLAPSE( 2                                                 )
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Zero SUNTMP_MID
             SUNTMP_MID = 0d0

             ! Grid box latitude center [radians]
             YMID_R = State_Grid%YMid_R(I,J)

             ! Local time at box (I,J) [hours]
             TIMLOC = GET_LOCALTIME( I, J, 1, State_Grid, GMT=GMT_MID)

             ! Hour angle at box (I,J) [radians]
             AHR = ABS( TIMLOC - 12d0 ) * 15d00 * PI_180

             !===========================================================
             ! The cosine of the solar zenith angle (SZA) is given by:
             !
             !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)
             !
             ! where LAT = the latitude angle,
             !       DEC = the solar declination angle,
             !       AHR = the hour angle, all in radians.
             !
             ! If SUNCOS < 0, then the sun is below the horizon, and
             ! therefore does not contribute to any solar heating.
             !===========================================================

             ! Compute Cos(SZA)
             SUNTMP_MID = sin(YMID_R) * sin(DEC) + &
                          cos(YMID_R) * cos(DEC) * cos(AHR)

             ! sumOfCosSza is the sum of SUNTMP_MID at location (I,J)
             ! Do not include negative values of SUNTMP_MID
             sumOfCosSza(I,J) = sumOfCosSza(I,J) + MAX( SUNTMP_MID, 0d0 )

          ENDDO
          ENDDO
          !$OMP END PARALLEL DO
       ENDDO

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    !========================================================================
    ! Calculate the OH diurnal scaling factor
    !========================================================================

    ! Chemistry timestep [s] and timesteps per day
    dtChem          = GET_TS_CHEM()
    timestepsPerDay = 86400.0_fp / dtChem

    ! Loop over surface grid boxes
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J                                                     )&
    !$OMP COLLAPSE( 2                                                       )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize loop variables
       OHdiurnalFac(I,J) = 0.0_fp

       ! Scaling factor for OH diurnal cycles - zero at night
       IF ( State_Met%SUNCOSmid(I,J) > 0.0_fp  .and.                         &
            sumOfCosSza(I,J)         > 0.0_fp ) THEN
          OHdiurnalFac(I,J) = State_Met%SUNCOSmid(I,J)                       &
                            / sumOfCosSza(I,J)                               &
                            * timestepsPerDay
       ENDIF
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Calc_Diurnal
!EOC
END MODULE Isotope_Mod
