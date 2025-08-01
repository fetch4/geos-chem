!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fullchem_RateLawFuncs
!
! !DESCRIPTION: Provides rate-law functions used by the "fullchem" chemical
!  mechanism.  This will be referenced from within subroutine Update_RCONST.
!\\
!\\
! !INTERFACE:
!
MODULE isotope_RateLawFuncs
!
! !USES:
!
  USE gckpp_Global
  USE gckpp_Parameters
  USE gckpp_Precision
  USE rateLawUtilFuncs

  IMPLICIT NONE
  PUBLIC
!
! !DEFINED PARAMETERS:
!
! !REFERENCES:
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS

  !#########################################################################
  !#####          RATE-LAW FUNCTIONS FOR GAS-PHASE REACTIONS           #####
  !#####   Some common functions are defined in rateLawUtilFuncs.F90   #####
  !#########################################################################

  FUNCTION GCJPLPR_abab( a1, b1, a2, b2, fv ) RESULT( k )
    ! Third body effect for pressure dependence of rate coefficients.
    ! a1, b1 are the Arrhenius parameters for the lower-limit rate.
    ! a2, b2 are the Arrhenius parameters for the upper-limit rate.
    ! fv     is the falloff curve paramter, (see ATKINSON ET. AL (1992)
    !        J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
    !
    ! Used to compute the rate for these reactions:
    !    NO   + OH  {+M} = HNO2  {+M}
    !    HO2  + NO2 {+M} = HNO4
    !    NO2  + NO3 {+M} = N2O5
    !    ClO  + NO2 {+M} = ClNO3 {+M}
    !    MCO3 + NO2 {+M} = PAN
    !    RCO3 + NO2 {+M} = PPN
    !    PRPE + OH  {+M} = PO2
    !    MO2  + NO2 {+M} = MPN   {+M}
    !    BrO  + NO2 {+M} = BrNO3 {+M}
    !    NO2  + O   {+M} = NO3   {+M}
    !    H    + O2  {+M} = HO2   {+M}
    !    IO   + NO2 {+M} = IONO2 {+M}
    !
    ! For these reactions, these Arrhenius law terms evaluate to 1:
    !    EXP(c1/T)
    !    EXP(c2/T)
    ! because c1 = c2 = 0.  Therefore we can skip computing these
    ! terms.  Also, fct1 = fct2 = 0, so we will skip computing these
    ! terms as well.  This is more computationally efficient.
    ! (bmy, 1/25/20)
    !
    REAL(dp), INTENT(IN) :: a1,   b1,    a2,    b2,   fv
    REAL(dp)             :: rlow, rhigh, xyrat, blog, fexp, k
    !
    rlow  = a1 * ( K300_OVER_TEMP**b1 ) * NUMDEN
    rhigh = a2 * ( K300_OVER_TEMP**b2 )
    xyrat = rlow / rhigh
    blog  = LOG10( xyrat )
    fexp  = 1.0_dp / ( 1.0_dp + ( blog * blog ) )
    k     = rlow * ( fv**fexp ) / ( 1.0_dp + xyrat )
  END FUNCTION GCJPLPR_abab

  FUNCTION GCJPLAC_ababac( a1, b1, a2, b2, a3, c3, fv ) RESULT( k )
    ! Rate coefficient for activation reactions competing with a 
    ! termolecular association pathway
    ! a1, b1 are the Arrhenius parameters for the lower-limit rate.
    ! a2, b2 are the Arrhenius parameters for the upper-limit rate.
    ! a3, c3 are Arrhenius parameters for the activation path
    ! fv     is the falloff curve parameter, usually = 0.6.
    !
    ! Used to compute the rate for these reactions:
    !    NO2  + O  = O2 + NO
    !    HNO3 + OH = NO3 + H2O
    !    CO   + OH = HO2 + CO2
    !
    REAL(dp), INTENT(IN) :: a1, b1, a2, b2, a3, c3, fv
    REAL(dp)             :: rlow, rhigh, xyrat, blog
    REAL(dp)             :: fexp, k1, k2, k
    !
    rlow  = a1 * ( K300_OVER_TEMP**b1 ) * NUMDEN
    rhigh = a2 * ( K300_OVER_TEMP**b2 )
    xyrat = rlow / rhigh
    blog  = LOG10( xyrat )
    fexp  = 1.0_dp / ( 1.0_dp + ( blog * blog ) )
    k1    = rlow * ( fv**fexp ) / ( 1.0_dp + xyrat )
    k2    = a3 * EXP( c3 / TEMP )
    k     = k2 * (1.0_dp - (k1 / rhigh) )
  END FUNCTION GCJPLAC_ababac
  
END MODULE isotope_RateLawFuncs
!EOC
