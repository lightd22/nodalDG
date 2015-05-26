SUBROUTINE evaluateExpansion(coeffs,nelem,basisVals,qvals)
  ! ===========================================================================
  ! Evaluates polynomial expansion phi = \sum coeffs_k * P_k at local quad nodes
  ! Used for outputting solution values
  ! INPUTS:
  !         coeffs(0:maxPolyDegree,1:nelem,1:meqn)
  !         basisVals(0:maxPolyDegree,0:nQuad)
  !
  ! OUTPUTS: qvals
  ! ===========================================================================
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nelem
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem,1:meqn), INTENT(IN) :: coeffs
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: basisVals
  ! Outputs
  DOUBLE PRECISION, DIMENSION(0:nQuad,1:nelem,1:meqn), INTENT(OUT) :: qvals
  ! Local valriables
  INTEGER :: i,j,m

  DO m=1,meqn
    DO j=1,nelem
      DO i=0,nQuad
        qVals(i,j,m) = SUM(coeffs(:,j,m)*basisVals(:,i))
      ENDDO !i
    ENDDO!j
  ENDDO !m

END SUBROUTINE evaluateExpansion
