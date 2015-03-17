SUBROUTINE forwardStep(coeffs,fluxQuad,flx,quadWeights,dLegVals,dxel,dt,nelem)
  ! ================================================================================
  ! Takes single forward Euler step applied to coefficient odes
  ! d a_kj / dt = forcingCoeffODE()
  ! Inputs:
  !         fluxQuad - flux function F(q) evaluated at quadrature nodes
  !         flx - numerical fluxes through interface
  !         quadWeghts - Gauss quadrature weights
  !         dLegVals - derivative of Legendre basis at quadrature nodes
  !         dxel - element spacing
  !         dt - time step size
  !         nelem - number of elements
  ! Outputs:
  !         coeffs - Legendre expansion coefficients
  ! ================================================================================
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nelem
  DOUBLE PRECISION, INTENT(IN) :: dt,dxel
  DOUBLE PRECISION, DIMENSION(0:nelem,1:meqn), INTENT(IN) :: flx
  DOUBLE PRECISION, DIMENSION(0:nQuad,1:nelem,1:meqn), INTENT(IN) :: fluxQuad
  DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadWeights
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: dLegVals
  ! Outputs
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem,1:meqn), INTENT(INOUT) :: coeffs
  ! Local variables
  INTEGER :: j,m
  DOUBLE PRECISION :: lam,fluxR,fluxL
  DOUBLE PRECISION, DIMENSION(0:nQuad) :: localFluxQuad

  INTERFACE
    FUNCTION forcingCoeffODE(fluxR,fluxL,fluxQuadVals,quadWeights,dLegVals)
      ! ==============================================================================
      ! Computes RHS forcing for all coefficients in current element and equation
      ! ==============================================================================
      USE commonTestParameters
      IMPLICIT NONE
      ! Inputs
      DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadWeights,fluxQuadVals
      DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: dLegVals
      DOUBLE PRECISION, INTENT(IN) :: fluxR,fluxL
      ! Outputs
      DOUBLE PRECISION, DIMENSION(0:maxPolyDegree) :: forcingCoeffODE
    END FUNCTION forcingCoeffODE
  END INTERFACE

  lam = dt/dxel

  DO m=1,meqn
    DO j=1,nelem
      fluxR = flx(j,m)
      fluxL = flx(j-1,m)
      localFluxQuad = fluxQuad(:,j,m)
      coeffs(:,j,m) = coeffs(:,j,m) + lam*forcingCoeffODE(fluxR,fluxL,localFluxQuad,quadWeights,dLegVals)
    ENDDO !j
  ENDDO !m

END SUBROUTINE forwardStep
