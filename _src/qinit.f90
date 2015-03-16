SUBROUTINE qinit(xVals,yVals,nx,ny,q)
  ! ==============================================================================
  ! Computes initial conditions for q fields
  ! INPUTS: meqn - number of fields to evaluate
  !         nx,ny - number of points to evaluate q at
  !
  ! OUTPUTS: q(i,j,neq) - initial conditions evaluated at xvals(i) for neqth field
  ! ==============================================================================
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  INTEGER,INTENT(IN) :: nx,ny
  DOUBLE PRECISION, DIMENSION(1:nx) :: xVals
  DOUBLE PRECISION, DIMENSION(1:ny) :: yVals
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nx,1:ny,1:meqn) :: q
  ! Local variables
  INTEGER :: i,j
  DOUBLE PRECISION, DIMENSION(1:nx,1:ny) :: r

  SELECT CASE(testID)
    CASE(0,1) ! Uniform field
      q = 1D0
    CASE(5) ! Cosbell deformation from LeVeque
      DO j=1,ny
          r(:,j) = 4D0*SQRT( (xVals-0.25D0)**2 + (yVals(j)-0.25D0)**2 )
      ENDDO !j
      q = 0D0
      WHERE(r .lt. 1D0)
          q(:,:,1) = 0.25D0*(1D0+DCOS(PI*r))**2
      END WHERE
  END SELECT !testID
END SUBROUTINE qinit
