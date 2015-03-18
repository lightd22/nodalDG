SUBROUTINE qinit(xVals,yVals,nx,ny,q,reactiveCoeffs)
  ! ==============================================================================
  ! Computes initial conditions for q fields
  ! INPUTS: meqn - number of fields to evaluate
  !         nx,ny - number of points to evaluate q at
  !
  ! OUTPUTS: q(i,j,neq) - initial conditions evaluated at xvals(i) for neqth field
  !     reactiveCoeffs(i,j,neq) - reactive flow coefficient at xval(i),yval(j)
  ! ==============================================================================
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  INTEGER,INTENT(IN) :: nx,ny
  DOUBLE PRECISION, DIMENSION(1:nx), INTENT(IN) :: xVals
  DOUBLE PRECISION, DIMENSION(1:ny), INTENT(IN) :: yVals
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nx,1:ny,1:meqn), INTENT(OUT) :: q,reactiveCoeffs
  ! Local variables
  INTEGER :: i,j
  DOUBLE PRECISION, DIMENSION(1:nx,1:ny) :: r
  DOUBLE PRECISION :: x0,y0

  SELECT CASE(testID)
    CASE(0) ! Uniform field
      q = 1D0
    CASE(1) ! Sine wave
      DO j=1,ny
          q(:,j,1) = sin(2.d0*PI*xVals(:))*sin(2.d0*PI*yVals(j))
      ENDDO !j
    CASE(2) ! Two half plane initial conditions
      DO j=1,ny
        r(:,j) = ABS(xVals(:)-0.25D0)
      ENDDO !j
      q = 0D0
      WHERE(r .lt. 0.25D0)
        q(:,:,1) = 1D0
      END WHERE
      q(:,:,2) = 1D0-q(:,:,1)
      q(:,:,2) = 2D0*q(:,:,2)

    CASE(5) ! Cosbell deformation from LeVeque
      DO j=1,ny
          r(:,j) = 4D0*SQRT( (xVals-0.25D0)**2 + (yVals(j)-0.25D0)**2 )
      ENDDO !j
      q = 0D0
      WHERE(r .lt. 1D0)
          q(:,:,1) = 0.25D0*(1D0+DCOS(PI*r))**2
      END WHERE
      q(:,:,2) = q(:,:,1)
    CASE(7) ! Slotted cylinder in deformation flow
        x0 = 0.25D0
        y0 = 0.5D0
        DO j=1,ny
            r(:,j) = SQRT((xVals-x0)**2 + (yVals(j)-y0)**2)
        ENDDO !j
        q = 0D0
        WHERE(r .lt. .15D0)
            q(:,:,1) = 1D0
        END WHERE

        DO j=1,ny
          DO i=1,nx
            IF(ABS(xVals(i)-x0) .lt. 0.025D0 .AND. yVals(j) .gt.(y0-0.0625D0)) THEN
                q(i,j,1) = 0D0
            ENDIF
          ENDDO !i
        ENDDO !j
  END SELECT !testID

  reactiveCoeffs = 0D0
  IF(doreactive) THEN
    x0 = 0.25D0
    IF(meqn .ne. 2) then
      write(*,*) 'ERROR! In qinit().. reactive test only supports meqn = 2'
      STOP
    ENDIF

    reactiveCoeffs(:,:,2) = 1D0
    DO i=1,nx
      IF(abs(xVals(i)-x0) .le. 0.25D0) THEN
        reactiveCoeffs(i,:,1) = COS(2D0*PI*(xVals(i)-x0))
      ENDIF
    ENDDO !i
  ENDIF
END SUBROUTINE qinit
