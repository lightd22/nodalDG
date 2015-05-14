SUBROUTINE reactiveJacobian(jacobian,qVals,forcingCoeffs,nx,ny)
  ! ==============================================================================
  ! Computes right hand side forcing term for chemical reaction equation
  ! INPUTS:
  !         qVals(1:nx,1:ny,1:meqn) - solution values at given points
  !         forcingCoeffs(1:nx,1:ny,1:meqn) - forcing coefficients multiplying fields q1,..qmeqn
  !           evaluated at given grid
  ! OUTPUTS:
  !         jacobian(1:nx,1:ny,1:meqn,1:meqn) - jacobian matrix evaluated at grid points
  ! ==============================================================================

  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nx,ny
  DOUBLE PRECISION, DIMENSION(1:nx,1:ny,1:meqn), INTENT(IN) ::qVals,forcingCoeffs
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nx,1:ny,1:meqn,1:meqn), INTENT(OUT) :: jacobian
  ! Local variables
  INTEGER :: i,j

  ! Evaluate jacobian at grid points
  DO i=1,nx
    DO j=1,ny
!      jacobian(i,j,1,1) = -forcingCoeffs(i,j,1)
!      jacobian(i,j,1,2) = 2D0*forcingCoeffs(i,j,2)*qVals(i,j,2)
!      jacobian(i,j,2,1) = -2D0*jacobian(i,j,1,1)
!      jacobian(i,j,1,2) = -2D0*jacobian(i,j,1,2)
      jacobian(i,j,1,1) = -forcingCoeffs(i,j,1)*qVals(i,j,2)
      jacobian(i,j,1,2) = -forcingCoeffs(i,j,1)*qVals(i,j,1)
      jacobian(i,j,2,1) = forcingCoeffs(i,j,1)*qVals(i,j,2)
      jacobian(i,j,2,2) = forcingCoeffs(i,j,1)*qVals(i,j,1)
    ENDDO !j
  ENDDO !i

END SUBROUTINE reactiveJacobian
