SUBROUTINE updateVelocities(uOut,vOut,uEdge,vEdge,xOut,elemEdgeX,DGx,yOut,elemEdgeY,DGy,time,dt)
  ! =========================================================
  ! Updates horizontal and vertical velocities
  ! at necessary grid points to time levels required by integrator
  ! =========================================================
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  DOUBLE PRECISION, INTENT(IN) :: time,dt
  DOUBLE PRECISION, DIMENSION(1:nxOut), INTENT(IN) :: xOut,DGx
  DOUBLE PRECISION, DIMENSION(1:nyOut), INTENT(IN) :: yOut,DGy
  DOUBLE PRECISION, DIMENSION(1:nex), INTENT(IN) :: elemEdgeX
  DOUBLE PRECISION, DIMENSION(1:ney), INTENT(IN) :: elemEdgeY
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:3,1:nxOut,1:nyOut), INTENT(INOUT) :: uOut,vOut
  DOUBLE PRECISION, DIMENSION(1:3,1:nex,1:nyOut), INTENT(INOUT) :: uEdge
  DOUBLE PRECISION, DIMENSION(1:3,1:nxOut,1:ney), INTENT(INOUT) :: vEdge

  ! Local Variables
  INTEGER :: stage
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut) :: uTmp,vTmp
  DOUBLE PRECISION, DIMENSION(1:nex,1:nyOut) :: uEdgeTmp
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:ney) :: vEdgeTmp
  DOUBLE PRECISION :: t_temp

  INTERFACE
    SUBROUTINE evalHorizVelocities(u,x,y,nx,ny,t)
      ! ======================================================
      ! User-specified subroutine which defines the non-divergent
      ! velocity field (u)
      ! ======================================================
      USE commonTestParameters
      IMPLICIT NONE
      ! Inputs
      INTEGER, INTENT(IN) :: nx,ny
      DOUBLE PRECISION, INTENT(IN) :: t
      DOUBLE PRECISION, DIMENSION(1:nx), INTENT(IN) :: x
      DOUBLE PRECISION, DIMENSION(1:ny), INTENT(IN) :: y

      ! Outputs
      DOUBLE PRECISION, DIMENSION(1:nx,1:ny), INTENT(OUT) :: u
    END SUBROUTINE evalHorizVelocities

    SUBROUTINE evalVertVelocities(v,x,y,nx,ny,t)
      ! ======================================================
      ! User-specified subroutine which defines the non-divergent
      ! velocity field (v)
      ! ======================================================
      USE commonTestParameters
      IMPLICIT NONE
      ! Inputs
      INTEGER, INTENT(IN) :: nx,ny
      DOUBLE PRECISION, INTENT(IN) :: t
      DOUBLE PRECISION, DIMENSION(1:nx), INTENT(IN) :: x
      DOUBLE PRECISION, DIMENSION(1:ny), INTENT(IN) :: y

      ! Outputs
      DOUBLE PRECISION, DIMENSION(1:nx,1:ny), INTENT(OUT) :: v
    END SUBROUTINE evalVertVelocities

  END INTERFACE

  DO stage = 1,3
    SELECT CASE(stage)
      CASE(1)
        t_temp = time
      CASE(2)
        t_temp = time + dt
      CASE(3)
        t_temp = time + 0.5D0*dt
    END SELECT
    ! Velocities at quadrature nodes
    CALL evalHorizVelocities(uTmp,DGx,yOut,nxOut,nyOut,t_temp)
    CALL evalVertVelocities(vTmp,xOut,DGy,nxOut,nyOut,t_temp)
    uOut(stage,:,:) = uTmp
    vOut(stage,:,:) = vTmp

    ! Velocities at element edges
    CALL evalHorizVelocities(uEdgeTmp,elemEdgeX,yOut,nex,nyOut,t_temp)
    CALL evalVertVelocities(vEdgeTmp,xOut,elemEdgeY,nxOut,ney,t_temp)
    uEdge(stage,:,:) = uEdgeTmp
    vEdge(stage,:,:) = vEdgeTmp

  ENDDO !stage
END SUBROUTINE updateVelocities
