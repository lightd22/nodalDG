SUBROUTINE updateSoln1d(q,u,uEdge,dt,dx,nelem,nx,quadWeights,avgOP,avgOP_LU,&
                        legVals,dlegVals,IPIV)
  ! ===========================================================================
  ! Takes full dt time step for one dimensional slice of subcell averages using SSPRK3
  ! integrator
  ! ===========================================================================
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nelem
  DOUBLE PRECISION, INTENT(IN) :: dxel,dt
  DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadWeights
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: legVals,&
    dlegVals
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:maxPolyDegree), INTENT(IN) :: avgOp,avgOp_LU
  INTEGER, DIMENSION(0:maxPolyDegree), INTENT(IN) :: IPIV
  DOUBLE PRECISION, DIMENSION(1:3,1:nx), INTENT(IN) :: u
  DOUBLE PRECISION, DIMENSION(1:3,1:nelem), INTENT(IN) :: uEdge
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nx,1:meqn), INTENT(INOUT) :: q
  ! Local variables
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem,1:meqn) :: coeffs,coeffsTmp,&
    qBar
  DOUBLE PRECISION, DIMENSION(1:3,0:nQuad,1:nelem) :: uTilde
  DOUBLE PRECISION, DIMENSION(1:3,0:nQuad,0:nelem+1) :: uEdgeTilde
  DOUBLE PRECISION, DIMENSION(0:nelem,1:meqn) :: numFlux
  DOUBLE PRECISION, DIMENSION(0:nQuad,1:meqn) :: localSolnQuad
  DOUBLE PRECISION, DIMENSION(0:nQuad) :: localVel
  INTEGER :: i,j,k,m

  ! Reshape incoming values
  DO j=1,nelem
    qBar(:,j,:) = q(1+(maxPolyDegree+1)*(j-1):(maxPolyDegree+1)*j,:)
    utild(1:3,:,j) = u(1:3,1+(maxPolyDegree+1)*(j-1) : (maxPolyDegree+1)*j)
  END DO
  ! Periodically extend edge velocities
  uedgeTilde(1:3,1:nelem) = uEdge(1:3,1:nelem)
  uedgeTilde(1:3,0) = uEdge(1:3,nelem)

  CALL projectAverages(coeffs,avgOP_LU,IPIV,qBar,nelem) ! Project incoming rhoq averages


END SUBROUTINE updateSoln1d
