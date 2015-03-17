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
  ! Local variables
  INTEGER :: k

  forcingCoeffODE = 0D0
  DO k=0,maxPolyDegree
    forcingCoeffODE(k) = SUM(quadWeights(:)*dLegVals(k,:)*fluxQuadVals(:))
    forcingCoeffODE(k) = forcingCoeffODE(k)-fluxR+((-1D0)**k)*fluxL
    forcingCoeffODE(k) = (2D0*k+1D0)*forcingCoeffODE(k)
  ENDDO !k
END FUNCTION forcingCoeffODE
