SUBROUTINE positivityLimiter(qBar,nelem,avgVals)
	! Subroutine for mass filling within an element to remove negative cell averaged values
  USE commonTestParameters
	IMPLICIT NONE
	! Inputs
	INTEGER, INTENT(IN) :: nelem
	DOUBLE PRECISION, DIMENSION(1:nelem,1:meqn), INTENT(IN) :: avgVals
  ! Outputs
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:nelem,1:meqn), INTENT(INOUT) :: qBar
	! Local Variables
	INTEGER :: j,k,m
	DOUBLE PRECISION :: r,Mp,Mt,s,valMin,avg
  DOUBLE PRECISION, DIMENSION(1:meqn) :: rTmp

	IF(limitingType .eq. 1) THEN
    ! ===============================================================================================
    ! TYPE 1: TMAR (2015) Limiting
    ! ===============================================================================================
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
	ELSEIF(limitingType .eq. 2) THEN
		! ===============================================================================================
		! TYPE 2: Replace truncation mass redistribution with rescaling similar to Zhang and Shu (2010)
		! ===============================================================================================
    DO m=1,meqn
      DO j=1,nelem
        avg = avgVals(j,m) !SUM(qBar(:,j))/(N+1)
        valMin = MINVAL(qBar(:,j,m))-epsilon(1D0)
        r = MIN(avg/abs(valMin-avg),1D0)
        qBar(:,j,m) = r*(qBar(:,j,m)-avg)+avg
      ENDDO !j
    ENDDO !m
  ELSEIF(limitingType .eq. 3) THEN
    ! ===============================================================================================
    ! TYPE 3: 'Strict' ZS Rescaling
    ! ===============================================================================================
    DO j=1,nelem
      rTmp = 0D0
      ! Let each field determine its 'desired' limiting ratio
      DO m=1,meqn
        avg = avgVals(j,m) !SUM(qBar(:,j))/(N+1)
        valMin = MINVAL(qBar(:,j,m))-epsilon(1D0)
        rTmp(m) = MIN(avg/abs(valMin-avg),1D0)
      ENDDO!m
      ! Set ratio used to be most strict ratio
      r = MINVAL(rTmp)

      ! Apply worst ratio to all fields
      DO m=1,meqn
        avg = avgVals(j,m)
        qBar(:,j,m) = r*(qBar(:,j,m)-avg)+avg
      ENDDO !m
    ENDDO!j
	ENDIF

END SUBROUTINE positivityLimiter
