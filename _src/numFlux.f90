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
  DOUBLE PRECISION, DIMENSION(0:1,1:nelem,1:meqn) :: qvals,fluxVals
  DOUBLE PRECISION, DIMENSION(0:1,0:nelem+1,1:meqn) :: fluxValsPeriodic
  DOUBLE PRECISION, DIMENSION(0:1,1:nelem) :: uTmp

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
      qvals(1,j,m) = coeffs(maxPolyDegree,j,m)
      qvals(0,j,m) = coeffs(0,j,m)

      uTmp(1,j) = uEdge(j)
      uTmp(0,j) = uEdge(j-1)
    ENDDO !j
  ENDDO !m

  CALL fluxFunction(qvals,uTmp,2,nelem,fluxVals)

  ! Form periodic extension of flux values at element edges
  fluxValsPeriodic(:,1:nelem,:) = fluxVals
  fluxValsPeriodic(:,0,:) = fluxVals(:,nelem,:)
  fluxValsPeriodic(:,nelem+1,:) = fluxVals(:,1,:)

  ! NOTE: As written this is not as general as it should be.
  ! Only valid for fluxes of the form F(q) = u*g(q) for some function g
  DO j=0,nelem
    whichSign = 1-0.5D0*(1-(SIGN(1D0,uEdge(j))) )
    whichEl = j+0.5D0*(1-(SIGN(1D0,uEdge(j))) )
    fluxes(j,:) = fluxValsPeriodic(whichSign,whichEl,:)
  ENDDO
END SUBROUTINE numFlux
