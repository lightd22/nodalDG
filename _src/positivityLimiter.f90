SUBROUTINE positivityLimiter(qBar,nelem,stat,avgVals)
	! Subroutine for mass filling within an element to remove negative cell averaged values
  USE commonTestParameters
	IMPLICIT NONE
	! Inputs
	INTEGER, INTENT(IN) :: nelem,stat
	DOUBLE PRECISION, DIMENSION(1:nelem,1:meqn), INTENT(IN) :: avgVals
  ! Outputs
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem,1:meqn), INTENT(INOUT) :: qBar
	! Local Variables
	INTEGER :: j,k,m
	DOUBLE PRECISION :: r,Mp,Mt,s,valMin,avg

	IF(stat .eq. 1) THEN
    DO m=1,meqn
      DO j=1,nelem
        Mp = 0D0
        Mt = 0D0

        DO k=0,maxPolyDegree
          Mt = Mt + qBar(k,j,m)
          qBar(k,j,m) = MAX(0D0,qBar(k,j,m)) ! Zero out negative masses
          Mp = Mp + qBar(k,j,m)
        ENDDO !k
        r = MAX(Mt,0D0)/MAX(Mp,TINY(1D0))
        qBar(:,j,m) = r*qBar(:,j,m) ! Reduce remaining positive masses by reduction factor
      ENDDO !j
    ENDDO !m
	ELSE
		! ===============================================================================================
		! ALTERNATIVE: Replace truncation mass redistribution with rescaling similar to Zhang and Shu (2010)
		! ===============================================================================================
    DO m=1,meqn
      DO j=1,nelem
        avg = avgVals(j,m) !SUM(qBar(:,j))/(N+1)
        valMin = MINVAL(qBar(:,j,m))-epsilon(1D0)
        r = MIN(avg/abs(valMin-avg),1D0)
        qBar(:,j,m) = r*(qBar(:,j,m)-avg)+avg
      ENDDO !j
    ENDDO !m
	ENDIF

END SUBROUTINE positivityLimiter
