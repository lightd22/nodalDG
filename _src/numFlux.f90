SUBROUTINE numFlux(coeffs,uEdge,nelem,fluxes)
  ! ===========================================================================
  ! Returns upwind numerical fluxes
  ! ===========================================================================
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nelem
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem,1:meqn), INTENT(IN):: coeffs
  DOUBLE PRECISION, DIMENSION(0:nelem+1), INTENT(IN) :: uEdge
  ! Outputs
  DOUBLE PRECISION, DIMENSION(0:nelem,1:meqn), INTENT(OUT) :: fluxes
  ! Local Variables
  INTEGER :: i,m,j,whichSign,whichEl
  DOUBLE PRECISION, DIMENSION(0:1,0:nelem+1,1:meqn) :: qvals,fluxVals

  INTERFACE
    SUBROUTINE fluxFunction(qvals,uvals,nx,nelem,fluxVals)
      USE commonTestParameters
      IMPLICIT NONE
      ! Inputs
      INTEGER, INTENT(IN) :: nelem,nx
      DOUBLE PRECISION, DIMENSION(1:nx,1:nelem,1:meqn), INTENT(IN) :: qvals
      DOUBLE PRECISION, DIMENSION(1:nx,1:nelem), INTENT(IN) :: uvals
      ! Outputs
      DOUBLE PRECISION, DIMENSION(1:nx,1:nelem,1:meqn), INTENT(OUT) :: fluxVals
    END SUBROUTINE fluxFunction
  END INTERFACE

  DO m=1,meqn
    DO j=1,nelem
      qvals(1,j,m) = SUM(coeffs(:,j,m))
      qvals(0,j,m) = SUM(coeffs(:,j,m)*(/ ((-1D0)**i,i=0,maxPolyDegree) /))
    ENDDO !j
    qvals(:,0,m) = qvals(:,nelem,m)
    qvals(:,nelem+1,m) = qvals(:,1,m)
  ENDDO !m

  CALL fluxFunction(qvals,uEdge,2,nelem+2,fluxVals)

  ! NOTE: As written this is not as general as it should be.
  ! Only valid for fluxes of the form F(q) = u*g(q) for some function g
  DO j=0,nelem
    whichSign = 1-0.5D0*(1-(SIGN(1D0,uEdge(j))) )
    whichEl = j+0.5D0*(1-(SIGN(1D0,uEdge(j))) )
    fluxes(j,:) = fluxVals(whichSign,whichEl,:)
  ENDDO
END SUBROUTINE numFlux
