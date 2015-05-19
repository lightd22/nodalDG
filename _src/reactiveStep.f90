SUBROUTINE reactiveStep(q,dt,forcingCoeffs)
  ! ================================================================================
  ! Takes a single time step to update subcell average values to tn+1
  ! for "toy chemistry problem"
  ! dq1/dt = -k1 q1 + k2 q2^2
  ! dq2/dt = 2 k1 q1 - 2 k2 q2^2
  ! Currently uses 2-stage, 2nd order Rosenbock Runge-Kutta method with the following parameters:
  ! (see Durran "Numerical Methods for Fluid Dynamics")
  ! b1 = b2 = 0.5 ; alpha = 1+0.5*sqrt(2) ; a21 = 1 ; alpha21 = -2 alpha
  !
  ! INPUTS:
  ! OUTPUTS: q(i,j,m) - mth field subcell averages updated to new time
  ! ================================================================================
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  DOUBLE PRECISION, INTENT(IN) :: dt
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn), INTENT(IN) :: forcingCoeffs
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn), INTENT(INOUT) :: q
  ! Local variables
  INTEGER i,j,m,ierr
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn,1:meqn) :: jac
  DOUBLE PRECISION, DIMENSION(1:meqn) :: fRHS,localCoeffs,IPIV
  DOUBLE PRECISION, DIMENSION(1:meqn,1:meqn) :: A,eye
  DOUBLE PRECISION, DIMENSION(1:meqn) :: localQ,localQ1,localQ2
  DOUBLE PRECISION :: alpha,alpha21

  INTERFACE
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
    END SUBROUTINE reactiveJacobian

    SUBROUTINE reactiveForcing(forcing,qVals,forcingCoeffs)
      ! ==============================================================================
      ! Computes right hand side forcing term for chemical reaction equation
      ! INPUTS:
      !         qVals(1:meqn) - solution values at given points
      !         forcingCoeffs(1:meqn) - forcing coefficients multiplying fields q1,..qmeqn
      ! OUTPUTS: forcing(1:meqn) - RHS forcing function for fields q1,...,qmeqn
      !
      ! ==============================================================================
      USE commonTestParameters
      IMPLICIT NONE
      ! Inputs
      DOUBLE PRECISION, DIMENSION(1:meqn), INTENT(IN) :: qVals,forcingCoeffs
      ! Outputs
      DOUBLE PRECISION, DIMENSION(1:meqn), INTENT(OUT) :: forcing
    END SUBROUTINE reactiveForcing
  END INTERFACE

  !alpha = 1D0+0.5D0/sqrt(2D0)
  alpha = 1D0+0.5D0*sqrt(2D0)
  alpha21 = -2D0*alpha

  eye = 0D0
  DO m=1,meqn
    eye(m,m) = 1D0
  ENDDO !m

  ! Fill in jacobian for this time
  CALL reactiveJacobian(jac,q,forcingCoeffs,nxOut,nyOut)

  DO i=1,nxOut
    DO j=1,nyOut
      ! Form LHS A matrix
      A = eye - alpha*dt*jac(i,j,:,:)
      localCoeffs = forcingCoeffs(i,j,:)
      localQ = q(i,j,:)

      CALL reactiveForcing(fRHS,localQ,localCoeffs)
      ! Solve for first stage
      CALL DGESV(meqn,1,A,meqn,IPIV,fRHS,meqn,ierr)
      localQ1 = fRHS

      CALL reactiveForcing(fRHS,localQ+dt*localQ1,localCoeffs)
      fRHS = fRHS-2D0*localQ1
      ! Solve for second stage
      CALL DGESV(meqn,1,A,meqn,IPIV,fRHS,meqn,ierr)
      localQ2 = fRHS

      q(i,j,:) = q(i,j,:)+0.5D0*dt*(3D0*localQ1+localQ2)
    ENDDO !j
  ENDDO !i
END SUBROUTINE reactiveStep
