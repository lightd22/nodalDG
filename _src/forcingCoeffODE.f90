FUNCTION forcingCoeffODE(fluxR,fluxL,fluxQuadVals,quadWeights,basisDeriv)
  ! ==============================================================================
  ! Computes RHS forcing for all coefficients in current element and equation
  ! ==============================================================================
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadWeights,fluxQuadVals
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: basisDeriv
  DOUBLE PRECISION, INTENT(IN) :: fluxR,fluxL
  ! Outputs
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree) :: forcingCoeffODE
  ! Local variables
  INTEGER :: k

  forcingCoeffODE = 0D0
  DO k=0,maxPolyDegree
    forcingCoeffODE(k) = SUM(quadWeights(:)*basisDeriv(k,:)*fluxQuadVals(:))
  ENDDO !k

  !write(*,*) 'b4-1',maxval(abs(forcingCoeffODE))

  forcingCoeffODE(0) = forcingCoeffODE(0)+fluxL

  !write(*,*) 'b4-2',maxval(abs(forcingCoeffODE))

  forcingCoeffODE(maxPolyDegree) = forcingCoeffODE(maxPolyDegree)-fluxR

  !write(*,*) 'b4-3',maxval(abs(forcingCoeffODE))

  forcingCoeffODE(:) = 2D0*forcingCoeffODE(:)/quadWeights(:)

  !write(*,*) 'after',maxval(abs(forcingCoeffODE))

END FUNCTION forcingCoeffODE
