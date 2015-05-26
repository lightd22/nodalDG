SUBROUTINE updateSoln1d(q,u,uEdge,dt,dxel,nelem,nx,quadWeights,&
                        basisVals,basisDeriv)
  ! ===========================================================================
  ! Takes full dt time step for one dimensional slice of subcell averages using SSPRK3
  ! integrator
  ! ===========================================================================
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nelem,nx
  DOUBLE PRECISION, INTENT(IN) :: dxel,dt
  DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadWeights
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: basisVals,&
    basisDeriv
  DOUBLE PRECISION, DIMENSION(1:3,1:nx), INTENT(IN) :: u
  DOUBLE PRECISION, DIMENSION(1:3,1:nelem), INTENT(IN) :: uEdge
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nx,1:meqn), INTENT(INOUT) :: q
  ! Local variables
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem,1:meqn) :: coeffs,coeffsTmp,&
    qBar
  DOUBLE PRECISION, DIMENSION(1:3,0:nQuad,1:nelem) :: uTilde
  DOUBLE PRECISION, DIMENSION(1:3,0:nelem+1) :: uEdgeTilde
  DOUBLE PRECISION, DIMENSION(0:nQuad,1:nelem) :: uQuadTmp
  DOUBLE PRECISION, DIMENSION(0:nQuad,1:nelem,1:meqn) :: quadVals,fluxQuad
  DOUBLE PRECISION, DIMENSION(0:nelem+1) :: uEdgeTmp
  DOUBLE PRECISION, DIMENSION(0:nelem,1:meqn) :: fluxes
  DOUBLE PRECISION, DIMENSION(0:nQuad,1:meqn) :: localSolnQuad
  DOUBLE PRECISION, DIMENSION(0:nQuad) :: localVel
  DOUBLE PRECISION, DIMENSION(1:nelem,1:meqn) :: elemAverages
  DOUBLE PRECISION :: cons0,consf
  INTEGER :: i,j,k,m,stage

  INTERFACE
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
    END SUBROUTINE evaluateExpansion

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
    END SUBROUTINE numFlux

    SUBROUTINE forwardStep(coeffs,fluxQuad,flx,quadWeights,basisDeriv,dxel,dt,nelem)
      ! ================================================================================
      ! Takes single forward Euler step applied to coefficient odes
      ! d a_kj / dt = forcingCoeffODE()
      ! Inputs:
      !         fluxQuad - flux function F(q) evaluated at quadrature nodes
      !         flx - numerical fluxes through interface
      !         quadWeghts - Gauss quadrature weights
      !         basisDeriv - derivative of Legendre basis at quadrature nodes
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
      DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: basisDeriv
      ! Outputs
      DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem,1:meqn), INTENT(INOUT) :: coeffs
    END SUBROUTINE forwardStep

    SUBROUTINE positivityLimiter(qBar,nelem,avgVals,quadWeights)
    	! Subroutine for mass filling within an element to remove negative cell averaged values
      USE commonTestParameters
    	IMPLICIT NONE
    	! Inputs
    	INTEGER, INTENT(IN) :: nelem
      DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadWeights
    	DOUBLE PRECISION, DIMENSION(1:nelem,1:meqn), INTENT(IN) :: avgVals
      ! Outputs
      DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem,1:meqn), INTENT(INOUT) :: qBar
    END SUBROUTINE positivityLimiter

    SUBROUTINE fluxCorrection(coeffs,flx,quadWeights,dxel,dt,nelem)
    	! Computes flux reductions factors to prevent total mass within each element from going negative
    	! Outputs fluxcf. fluxcf(j) is the reduction factor for the right face of element j,
    	!  with fluxcf(0) being the factor for the left domain interface
      USE commonTestParameters
    	IMPLICIT NONE
    	! -- Inputs
    	INTEGER, INTENT(IN) :: nelem
    	DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem,1:meqn), INTENT(IN) :: coeffs
      DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadWeights
    	DOUBLE PRECISION, INTENT(IN) :: dxel,dt
    	! -- Outputs
    	DOUBLE PRECISION, DIMENSION(0:nelem,1:meqn), INTENT(INOUT) :: flx
    END SUBROUTINE fluxCorrection
  END INTERFACE

  ! Reshape incoming values
  DO j=1,nelem
    coeffs(:,j,:) = q(1+(maxPolyDegree+1)*(j-1):(maxPolyDegree+1)*j,:)
    utilde(1:3,:,j) = u(1:3,1+(maxPolyDegree+1)*(j-1):(maxPolyDegree+1)*j)
  END DO
  ! Periodically extend edge velocities
  uedgeTilde(1:3,1:nelem) = uEdge(1:3,1:nelem)
  uedgeTilde(1:3,0) = uEdge(1:3,nelem)
  uedgeTilde(1:3,nelem+1) = uEdge(1:3,1)

  coeffsTmp = coeffs

!  ! Est. conservation
!  DO j=1,nelem
!    elemAverages(j,2) = 0.5D0*SUM(quadWeights(:)*coeffs(:,j,2))
!  ENDDO !j
!  cons0 = SUM(elemAverages(:,2))/DBLE(nelem)

  DO stage=1,3
    uQuadTmp = uTilde(stage,:,:)
    uEdgeTmp = uEdgeTilde(stage,:)

    ! Evaluate expansion and fluxes at quadrature nodes, compute numerical fluxes at interfaces
    CALL evaluateExpansion(coeffsTmp,nelem,basisVals,quadVals)
    CALL fluxFunction(quadVals,uQuadTmp,nQuad+1,nelem,fluxQuad)
    CALL numFlux(coeffsTmp,uEdgeTmp,nelem,fluxes)

    IF(doposlimit) CALL fluxCorrection(coeffsTmp,fluxes,quadWeights,dxel,dt,nelem)

    ! Take forward step
    CALL forwardStep(coeffsTmp,fluxQuad,fluxes,quadWeights,basisDeriv,dxel,dt,nelem)

    ! Update coefficients
    SELECT CASE(stage)
    CASE(2)
      coeffsTmp = 0.75D0*coeffs + 0.25D0*coeffsTmp
    CASE(3)
      coeffsTmp = coeffs/3d0 + 2D0*coeffsTmp/3D0
    END SELECT !stage
  ENDDO !stage

  ! Check conservation
!  DO m=1,meqn
!    DO j=1,nelem
!      elemAverages(j,m) = 0.5D0*SUM(quadWeights(:)*coeffsTmp(:,j,m))
!    ENDDO !j
!  ENDDO !m
!  consf = SUM(elemAverages(:,2))/DBLE(nelem)
!  write(*,*) 'cons estm =',consf-cons0


  IF(doposlimit) THEN
    DO m=1,meqn
      DO j=1,nelem
        elemAverages(j,m) = 0.5D0*SUM(quadWeights(:)*coeffsTmp(:,j,m))
      ENDDO !j
    ENDDO !m
    CALL positivityLimiter(coeffsTmp,nelem,elemAverages,quadWeights)
  ENDIF!doposlimit

  ! Reform original shaped arrays
  DO m=1,meqn
    DO j=1,nelem
      q(1+(maxPolyDegree+1)*(j-1):(maxPolyDegree+1)*j,m) = coeffsTmp(:,j,m)
    END DO !j
  ENDDO !m

END SUBROUTINE updateSoln1d
