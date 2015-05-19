SUBROUTINE init2d(q,u,v,uEdge,vEdge,xPlot,yPlot,quadNodes,dxel,dyel,&
                  dxPlot,dyPlot,elemCenterX,elemCenterY,reactiveCoeffs)
  ! ==============================================================================
  ! Computes initial conditions for q,u,v fields
  ! INPUTS: meqn - number of fields to evaluate
  !         nx,ny - number of points to evaluate q,u,v at
  !         quadNodes - local locations to evaluate velocities at
  !         elemCenterX, elemCenterY - element center locations
  !
  ! OUTPUTS: q(i,j,neq) - initial conditions evaluated for neqth field
  !          u(i,j),v(i,j) - velocities evaluated at quadrature locations
  !          uEdge,vEdge - velocities at edges of each element
  !          reactiveCoeffs - reaction coefficients at x(i),y(j)
  ! ==============================================================================
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  DOUBLE PRECISION, INTENT(IN) :: dxel, dyel,dxPlot,dyPlot
  DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadNodes
  DOUBLE PRECISION, DIMENSION(1:nxOut), INTENT(IN) :: xPlot
  DOUBLE PRECISION, DIMENSION(1:nyOut), INTENT(IN) :: yPlot
  DOUBLE PRECISION, DIMENSION(1:nex), INTENT(IN) :: elemCenterX
  DOUBLE PRECISION, DIMENSION(1:ney), INTENT(IN) :: elemCenterY
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn) :: q,reactiveCoeffs
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut) :: u,v
  DOUBLE PRECISION, DIMENSION(1:nex,1:nyOut), INTENT(OUT) :: uEdge
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:ney), INTENT(OUT) :: vEdge

  ! Local Variables
  INTEGER :: i,j,l
  DOUBLE PRECISION, DIMENSION(1:nxOut) :: DGx
  DOUBLE PRECISION, DIMENSION(1:nyOut) :: DGy
  DOUBLE PRECISION, DIMENSION(1:nxOut,0:1) :: xtilde
  DOUBLE PRECISION, DIMENSION(1:nyOut,0:1) :: ytilde
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,0:1) :: psiu,psiv
  DOUBLE PRECISION, DIMENSION(1:nex,1:nyOut,0:1) :: psiuEdge
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:ney,0:1) :: psivEdge

  INTERFACE
    SUBROUTINE qinit(xVals,yVals,nx,ny,q,reactiveCoeffs)
      USE commonTestParameters
      IMPLICIT NONE
      ! Inputs
      INTEGER,INTENT(IN) :: nx,ny
      DOUBLE PRECISION, DIMENSION(1:nx) :: xVals
      DOUBLE PRECISION, DIMENSION(1:ny) :: yVals
      ! Outputs
      DOUBLE PRECISION, DIMENSION(1:nx,1:ny,1:meqn) :: q,reactiveCoeffs
    END SUBROUTINE
  END INTERFACE

  ! Compute ICs on plotting grid
  CALL qinit(xPlot,yPlot,nxOut,nyOut,q,reactiveCoeffs)

  ! ====================================================
  ! Compute velocities at quad nodes from streamfunction
  ! ====================================================
  xtilde(:,0) = xPlot(:)-0.5D0*dxPlot
  xtilde(:,1) = xPlot(:)+0.5D0*dxPlot

  ytilde(:,0) = yPlot(:)-0.5D0*dyPlot
  ytilde(:,1) = yPlot(:)+0.5D0*dyPlot

  DO i=1,nex
    DGx(1+(i-1)*(maxPolyDegree+1):i*(maxPolyDegree+1)) = elemCenterX(i)+0.5D0*dxel*quadNodes(0:nQuad)
  ENDDO

  DO i=1,ney
    DGy(1+(i-1)*(maxPolyDegree+1):i*(maxPolyDegree+1)) = elemCenterY(i)+0.5D0*dyel*quadNodes(0:nQuad)
  ENDDO

  SELECT CASE(testID)
    CASE(0,1) ! uniform diagonal advection of a sine wave
      ! Evaluate stream function for horizontal velocities
      DO j=1,nyOut
        psiu(:,j,0) = -DGx(:) + ytilde(j,0)
        psiu(:,j,1) = -DGx(:) + ytilde(j,1)
        psiuEdge(:,j,0) = -(elemCenterX(:)+0.5D0*dxel) + ytilde(j,0)
        psiuEdge(:,j,1) = -(elemCenterX(:)+0.5D0*dxel) + ytilde(j,1)
      ENDDO!j

      ! Evaluate stream function for vertical velocities
      DO i=1,nxOut
        psiv(i,:,0) = -xtilde(i,0)+DGy(:)
        psiv(i,:,1) = -xtilde(i,1)+DGy(:)
        psivEdge(i,:,0) = -xtilde(i,0)+(elemCenterY(:)+0.5D0*dyel)
        psivEdge(i,:,1) = -xtilde(i,1)+(elemCenterY(:)+0.5D0*dyel)
      ENDDO!i
    CASE(99) ! no flow
      psiu = 0D0
      psiv = 0D0
      psiuEdge=0d0
      psivEdge=0d0
    CASE(2,5:7) ! LeVeque deformation flow
      ! Evaluate stream function for horizontal velocities (1/pi)*sin(pi*xf(i))**2 * sin(pi*yf(j))**2
      DO j=1,nyOut
        psiu(:,j,0) = (SIN(PI*DGx(:))**2 * SIN(PI*ytilde(j,0))**2 )/PI
        psiu(:,j,1) = (SIN(PI*DGx(:))**2 * SIN(PI*ytilde(j,1))**2 )/PI
        psiuEdge(:,j,0) = (SIN(PI*(elemCenterX(:)+0.5D0*dxel))**2 * SIN(PI*ytilde(j,0))**2)/PI
        psiuEdge(:,j,1) = (SIN(PI*(elemCenterX(:)+0.5D0*dxel))**2 * SIN(PI*ytilde(j,1))**2)/PI
      ENDDO!j

      ! Evaluate stream function for vertical velocities
      DO i=1,nxOut
        psiv(i,:,0) = (SIN(PI*xtilde(i,0))**2 * SIN(PI*DGy(:))**2)/PI
        psiv(i,:,1) = (SIN(PI*xtilde(i,1))**2 * SIN(PI*DGy(:))**2)/PI
        psivEdge(i,:,0) = (SIN(PI*xtilde(i,0))**2 * SIN(PI*(elemCenterY(:)+0.5D0*dyel))**2)/PI
        psivEdge(i,:,1) = (SIN(PI*xtilde(i,1))**2 * SIN(PI*(elemCenterY(:)+0.5D0*dyel))**2)/PI
      ENDDO!i
  END SELECT !testID

  ! Compute u velocities from stream function
  DO j=1,nyOut
    u(:,j) = (psiu(:,j,1)-psiu(:,j,0))/dyPlot
    uEdge(:,j) = (psiuEdge(:,j,1)-psiuEdge(:,j,0))/dyPlot
  ENDDO!j

  ! Compute v velocities from stream function
  DO i=1,nxOut
    v(i,:) = -1D0*(psiv(i,:,1)-psiv(i,:,0))/dxPlot
    vEdge(i,:) = -1D0*(psivEdge(i,:,1)-psivEdge(i,:,0))/dxPlot
  ENDDO!i
  
END SUBROUTINE init2d
