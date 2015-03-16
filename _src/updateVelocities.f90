SUBROUTINE updateVelocities(u,v,uEdge,vEdge,time)
  ! =========================================================
  ! Updates the edge and quadrature velocities to given time
  ! =========================================================
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  DOUBLE PRECISION, INTENT(IN) :: time
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut), INTENT(INOUT) :: u,v
  DOUBLE PRECISION, DIMENSION(1:nex,1:nyOut), INTENT(INOUT) :: uEdge
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:ney), INTENT(INOUT) :: vEdge

  ! Local Variables
  DOUBLE PRECISION :: timeFactor

  INTERFACE
    FUNCTION tfcn(t)
      REAL(KIND=8), INTENT(IN) :: t
    END FUNCTION tfcn
  END INTERFACE

  timeFactor = tfcn(time)
  
  u = u*timeFactor
  v = v*timeFactor
  uEdge = uEdge*timeFactor
  vEdge = vEdge*timeFactor

END SUBROUTINE updateVelocities
