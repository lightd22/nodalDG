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

  ! Local variables
  INTEGER :: j
  DOUBLE PRECISION :: timeFac

  INTERFACE
    function tfcn(t)
      USE commonTestParameters
      ! Inputs
      DOUBLE PRECISION, intent(in) :: t
      ! Outputs
      DOUBLE PRECISION :: tfcn
    END FUNCTION tfcn
  END INTERFACE

  timeFac = tfcn(t)

  SELECT CASE(testID)
    CASE(2,5:7) ! LeVeque deformation flow
      DO j=1,ny
        ! Horizontal velocities -- u
        u(:,j) = uMean + (SIN(PI*(x(:)-uMean*t))**2)*SIN(2D0*PI*(y(j)-vMean*t))&
                 *timeFac
      ENDDO !j
  END SELECT !testID
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

  ! Local variables
  INTEGER :: j
  DOUBLE PRECISION :: timeFac

  INTERFACE
    function tfcn(t)
      USE commonTestParameters
      ! Inputs
      DOUBLE PRECISION, intent(in) :: t
      ! Outputs
      DOUBLE PRECISION :: tfcn
    END FUNCTION tfcn
  END INTERFACE

  timeFac = tfcn(t)

  SELECT CASE(testID)
    CASE(2,5:7) ! LeVeque deformation flow
      DO j=1,ny
        ! Vertical velocities -- v
        v(:,j) = vMean - (SIN(PI*(y(j)-vMean*t))**2)*SIN(2D0*PI*(x(:)-uMean*t))&
                *timeFac
      ENDDO !j
  END SELECT !testID
END SUBROUTINE evalVertVelocities
