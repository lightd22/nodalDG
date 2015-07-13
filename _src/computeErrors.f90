SUBROUTINE computeErrors(qOut,q0,quadWeights,e1,e2,ei,cons,qMax,qMin,tf,nRuns,nex0,ney0,nscale,stat)
  ! =============================================================================
  ! Prints error estimates and other useful information to screen
  ! INPUTS: qOut - current estimate solution
  !         q0   - initial conditions
  !         quadWeights - quadrature weights (used in conservation estimation)
  !         tf(p) - cput time for pth run
  !         stat - status integer
  !         ovrshoot(p,m) - maximum overshoot in appx soln at plotting times
  !         undrshoot(p,m) - maximum undershoot in appx soln at plotting times
  ! OUTPUTS: e1(p,m) - L1 error estimate
  !          e2(p,m) - L2 error estimate
  !          ei(p,m) - Linf error estimate
  !          cons(p,m) - conservation estimate
  ! =============================================================================
  USE commonTestParameters
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nRuns,stat,nscale,nex0,ney0
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn), INTENT(IN) :: q0,qOut
  DOUBLE PRECISION, DIMENSION(1:nRuns,1:meqn), INTENT(IN) :: qMax,qMin
  DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadWeights
  REAL(KIND=4), DIMENSION(1:nRuns),INTENT(IN) :: tf
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nRuns,1:meqn),INTENT(INOUT) :: e1,e2,ei,cons
  ! Local variables
  INTEGER :: p,m,startHoriz,startVert,endHoriz,endVert,i,j,l,currnex,currney
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:maxPolyDegree) :: coeffs,tmp
  CHARACTER(len=2) :: qname
  DOUBLE PRECISION :: cnvg1,cnvg2,cnvgi

  IF(stat == -1) THEN
    ! Write error output to screen
    DO m=1,meqn
      WRITE(qname,'(a,i1)') 'q',m
      WRITE(*,*) '===================='
      WRITE(*,'(a12)') qname
      WRITE(*,*) '===================='
    WRITE(*,'(A115)') &
'nex    ney    E1        E2         Einf        convergence rate  maximum   minimum       cons       cputime   tf'

      cnvg1 = 0D0
      cnvg2 = 0D0
      cnvgi = 0D0

      DO p=1,nRuns
        IF(p.gt.1) THEN
          cnvg1 = -log(e1(p,m)/e1(p-1,m))/log(dble(nscale))
          cnvg2 = -log(e2(p,m)/e2(p-1,m))/log(dble(nscale))
          cnvgi = -log(ei(p,m)/ei(p-1,m))/log(dble(nscale))
        ENDIF

        currnex = nex0*nscale**(p-1)
        currney = ney0*nscale**(p-1)

        WRITE(*,990) currnex, currney, e1(p,m), e2(p,m), ei(p,m), &
              cnvg1, cnvg2, cnvgi, &
              qMax(p,m), &
              qMin(p,m), &
              cons(p,m), tf(p),tfinal
      ENDDO!p
    ENDDO !m
    990    format(2i6,3e12.4,3f5.2,3e12.4,2f8.2)
  ELSE
    ! Compute error estimates for this run
    DO m=1,meqn
      ! Conservation estimate
      cons(stat,m) = 0D0
      e1(stat,m) = 0D0
      e2(stat,m) = 0D0
      DO i=1,nex
        DO j=1,ney
          startHoriz = 1+(maxPolyDegree+1)*(i-1)
          startVert = 1+(maxPolyDegree+1)*(j-1)
          endHoriz = startHoriz + maxPolyDegree
          endVert = startVert + maxPolyDegree

          coeffs(:,:) = qOut(startHoriz:endHoriz,startVert:endVert,m)-q0(startHoriz:endHoriz,startVert:endVert,m)
          DO l=0,maxPolyDegree
            tmp(:,l) = 0.25D0*coeffs(:,l)*quadWeights(:)*quadWeights(l)
          ENDDO !l
          cons(stat,m) = cons(stat,m) + SUM(tmp)

          DO l=0,maxPolyDegree
            tmp(:,l) = 0.25D0*abs(coeffs(:,l))*quadWeights(:)*quadWeights(l)
          ENDDO !l
          e1(stat,m) = e1(stat,m) + SUM(tmp)

          DO l=0,maxPolyDegree
            tmp(:,l) = 0.25D0*quadWeights(:)*quadWeights(l)*coeffs(:,l)**2
          ENDDO !l
          e2(stat,m) = e2(stat,m) + SUM(tmp)
        ENDDO !j
      ENDDO!i
      cons(stat,m) = cons(stat,m)/DBLE(nex*ney)
      e1(stat,m) = e1(stat,m)/DBLE(nex*ney)
      e2(stat,m) = SQRT(e2(stat,m)/DBLE(nex*ney))
      ei(stat,m) = MAXVAL(ABS( qOut(:,:,m)-q0(:,:,m) ))
    ENDDO !m

  ENDIF

END SUBROUTINE computeErrors
