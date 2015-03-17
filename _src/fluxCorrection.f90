SUBROUTINE fluxCorrection(coeffs,flx,dxel,dt,nelem)
	! Computes flux reductions factors to prevent total mass within each element from going negative
	! Outputs fluxcf. fluxcf(j) is the reduction factor for the right face of element j,
	!  with fluxcf(0) being the factor for the left domain interface
  USE commonTestParameters
	IMPLICIT NONE
	! -- Inputs
	INTEGER, INTENT(IN) :: nelem
	DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem,1:meqn), INTENT(IN) :: coeffs
	DOUBLE PRECISION, INTENT(IN) :: dxel,dt
	! -- Outputs
	DOUBLE PRECISION, DIMENSION(0:nelem,1:meqn), INTENT(INOUT) :: flx
	! -- Local variables
	DOUBLE PRECISION :: Pj,Qj,eps
	DOUBLE PRECISION :: fluxcf
	DOUBLE PRECISION, DIMENSION(0:nelem+1) :: R ! Reduction ratio for outward fluxes so that element j has non-negative values (1D0 indicates no limiting needed)
	INTEGER :: j,m

	eps = 1D-6 ! Small parameter used to ensure no division by 0

  DO m=1,meqn
    DO j=1,nelem
      ! Compute maximum allowable flux out of element j
      Qj = (dxel/dt)*coeffs(0,j,m)

      Qj = MAX(Qj-epsilon(1D0),0D0)

      ! Compute actual flux out of element j
      Pj = MAX(0D0,flx(j,m)) - MIN(0D0,flx(j-1,m)) + eps

      ! Compute reduction ratio
      R(j) = MIN(1D0,Qj/Pj)
    END DO!j
    ! Periodicity
    R(0) = R(nelem)
    R(nelem+1) = R(1)
    ! Compute corrected factors
  	DO j=0,nelem
  		! If flux at right edge is negative, use limiting ratio in element to the right of current one
  		! (since that is where we are pulling mass from)
  		fluxcf = R(j) - 0.5D0*(1D0-INT(SIGN(1D0,flx(j,m))))*(R(j)-R(j+1))
  		flx(j,m) = flx(j,m)*fluxcf
  	END DO!j
  ENDDO !m
END SUBROUTINE fluxCorrection
