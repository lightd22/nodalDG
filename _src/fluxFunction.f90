SUBROUTINE fluxFunction(qvals,uvals,nx,nelem,fluxVals)
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nelem,nx
  DOUBLE PRECISION, DIMENSION(1:nx,1:nelem,1:meqn), INTENT(IN) :: qvals
  DOUBLE PRECISION, DIMENSION(1:nx,1:nelem), INTENT(IN) :: uvals
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nx,1:nelem,1:meqn), INTENT(OUT) :: fluxVals
  ! Local variables
  INTEGER :: m,j,i

  ! Advective flux: f(q) = u*q
  DO m=1,meqn
    DO j=1,nelem
      fluxVals(:,j,m) = uvals(:,j)*qvals(:,j,m)
    ENDDO !j
  ENDDO !m

END SUBROUTINE fluxFunction
