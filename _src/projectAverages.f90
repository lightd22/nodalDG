SUBROUTINE projectAverages(coeffs,avgOP_LU,IPIV,avgs,nelem)
  USE commonTestParameters
	IMPLICIT NONE
	! -- Inputs
	INTEGER, INTENT(IN) :: nelem
	DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:maxPolyDegree), INTENT(IN) :: avgOP_LU
	DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem,1:meqn), INTENT(IN) :: avgs
	INTEGER, DIMENSION(0:maxPolyDegree), INTENT(IN) :: IPIV
	! -- Outputs
	DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem,1:meqn) :: coeffs
	! -- Local variables
	INTEGER :: i,j,k,m
	DOUBLE PRECISION, DIMENSION(1:nelem) :: hold
	DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem) :: fooBAR
	DOUBLE PRECISION, DIMENSION(0:maxPolyDegree) :: FOO_y

  DO m=1,meqn
  	fooBAR = avgs(:,:,m)

  	DO i=0,maxPolyDegree ! Reorder RHS according to IPIV
      hold = fooBAR(i,1:nelem)
      fooBAR(i,:) = fooBAR(IPIV(i)-1,:)
      fooBAR(IPIV(i)-1,1:nelem) = hold
  	ENDDO

  	DO j=1,nelem
  		FOO_y = 0D0
  		! Solve Ly=RHS for y
  		FOO_y(0) = fooBAR(0,j)
  		DO k=1,maxPolyDegree
  			FOO_y(k) = fooBAR(k,j) - SUM(avgOP_LU(k,0:k-1)*FOO_y(0:k-1))
  		ENDDO
  		! Solve Ux=y for x
  		coeffs(maxPolyDegree,j,m) = (1D0/avgOP_LU(maxPolyDegree,maxPolyDegree))*FOO_y(maxPolyDegree)
  		DO k=maxPolyDegree-1,0,-1
  			coeffs(k,j,m) = (1D0/avgOP_LU(k,k))*(FOO_y(k) - SUM(avgOP_LU(k,k+1:maxPolyDegree)*coeffs(k+1:maxPolyDegree,j,m)))
  		ENDDO

  		IF(coeffs(0,j,m) .lt. 0D0) coeffs(0,j,m) = coeffs(0,j,m)+epsilon(1D0) ! To prevent numerical roundoff from causing negatives
  	ENDDO !j
  ENDDO !m

END SUBROUTINE projectAverages
