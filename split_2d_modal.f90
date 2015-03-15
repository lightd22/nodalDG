! =====================================
! Split 2D modal algorithm for tracer transport
! Strang splitting and Legendre basis to simulate 2D tracer transport equations with variable windspeeds
!
! Dependencies:
! 	netCDF
!	LAPACK
!	mDGmod.f90 ; mDGsweep.f90 ; tfcn.f90 ; inputs.nl
! By: Devin Light ; Sept. 2014
! =====================================

PROGRAM EXECUTE
    USE commonTestParam
    USE netCDF

    IMPLICIT NONE
    INTEGER :: startRes,noutput,nRuns,nScale
    LOGICAL :: debug
    REAL(KIND=8) :: muMAX,cflCoeff
    INTEGER, PARAMETER :: inUnit=20,outUnit=21

    INTERFACE
      SUBROUTINE DRIVER(testID,nex0,ney0,nscale,nruns,noutput,maxCFL)
        INTEGER, INTENT(IN) :: testID,nex0,ney0,nscale,nruns,noutput
        REAL(KIND=8), INTENT(IN) :: maxCFL
      END SUBROUTINE DRIVER
    END INTERFACE

    NAMELIST /inputs/ startRes,nRuns,nScale,maxPolyDegree,cflCoeff,noutput,meqn, &
                      testID,TRANSIENT,DEBUG

    OPEN(unit=inUnit,file="inputs.nl",action="read")
    READ(inUnit,NML=inputs)

    doposlimit = .FALSE.

    write(*,*) meqn
    muMAX  = determineCFL(maxPolyDegree,cflCoeff)

    write(*,*) '======================================================'
    write(*,*) '             BEGINNING RUN OF MODAL TESTS             '
    write(*,'(A27,F7.4)') 'muMAX=',muMAX
    write(*,'(A10,L5)') 'TRANSIENT =',transient
    write(*,'(A10,I5)') 'meqn = ',meqn
    write(*,*) '======================================================'

    write(*,*) '======'
    SELECT CASE(testID)
        CASE(0)
        	write(*,*) 'TEST 0: Consistency test'
        CASE(1)
        	write(*,*) 'TEST 1: Uniform advection (u=v=1)'
        CASE(2)
        	write(*,*) 'TEST 2: Solid body rotation of cylinder'
        CASE(5)
          write(*,*) 'TEST 5: LeVeque Cosbell Deformation Test'
        CASE(6)
        	write(*,*) 'TEST 6: LeVeque Smoother Cosbell Deformation Test'
        CASE(7)
        	write(*,*) 'TEST 7: Slotted Cylinder Deformation Test'
    END SELECT
	write(*,*) '======'
	CALL test2d_modal(whichTest,startRes,startRes,nScale,nRuns,noutput,muMAX) !1D0/(2D0*4D0-1D0) !0.3D0/sqrt(2d0)
  CLOSE(inUnit)
  WRITE(*,*) 'PROGRAM COMPLETE!'

CONTAINS
  FUNCTION determineCFL(maxPolyDegree,cflCoeff)
    ! ===============================================================
    ! Determines max CFL for SSPRK3 timestepping as a function of
    ! maximum reconstructing polynomial degree
    ! ===============================================================
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: maxPolyDegree
    DOUBLE PRECISION :: cflCoeff
    ! Outputs
    DOUBLE PRECISION determineCFL
    ! Local variables

    SELECT CASE(maxPolyDegree)
      CASE(2)
        determineCFL = 0.209D0
      CASE(3)
        determineCFL = 0.130D0
      CASE(4)
        determineCFL = 0.089D0
      CASE(5)
        determineCFL = 0.066D0
      CASE(6)
        determineCFL = 0.051D0
      CASE(7)
        determineCFL = 0.04D0
      CASE(8)
        determineCFL = 0.033D0
      CASE(9)
        determineCFL = 0.026D0
    END SELECT
    determineCFL = determineCFL*cflCoeff

  END FUNCTION determineCFL
    SUBROUTINE test2d_modal(ntest,nex0,ney0,nscale,nlevel,noutput,maxcfl)
      IMPLICIT NONE
      ! Inputs
      INTEGER, INTENT(IN) :: ntest,nex0,ney0,nscale,nlevel,noutput
      REAL(KIND=8), INTENT(IN) :: maxcfl
      ! Outputs
      ! Local variables
      INTEGER, DIMENSION(10) :: tmp_method
      REAL(KIND=8), DIMENSION(nlevel) :: e1, e2, ei
	    REAL(KIND=8) :: cnvg1, cnvg2, cnvgi, cons
      INTEGER :: nmethod,nmethod_final,imethod,ierr,nmethodx,nmethody
      INTEGER :: i,j,l,n,nOrder,p,nex,ney,nstep,nx,ny,nout
      INTEGER, DIMENSION(1:2) :: horizBnd,vertBnd
      LOGICAL :: oddstep,dogllGrid

	    CHARACTER(len=40) :: cdf_out
      CHARACTER(len=9) :: outdir

	    REAL(KIND=4) :: t0,tf,t1,t2,t3
      REAL(KIND=8) :: maxTime, minTime,totTime

	    REAL(KIND=8) :: PI,dxel,dyel,tfinal,tmp_umax,tmp_vmax,tmp_qmax,tmp_qmin,dxm,dym,dt,time,calculatedMu,dxPlot,dyPlot
      REAL(KIND=8), DIMENSION(1:2) :: xEdge,yEdge,domainCenter
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: q,q0,u0,v0,uEdge0,vEdge0,avgXferOp,avgXferOpLU,legendreVal,legendreDeriv
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: elemCenterX,elemCenterY,quadNodes,quadWeights,DGx,DGy,xPlot,yPlot,FOO,&
                                                   gllNodes,gllWeights,posWeight
	    INTEGER, ALLOCATABLE, DIMENSION(:) :: IPIV

        if(nlevel.lt.1) STOP 'nlev should be at least 1 in calls to test2d_modal'

        PI = DACOS(-1D0)

        ! Set domain size
        xEdge(1) = 0
        xEdge(2) = 1
        yEdge = xEdge
        domainCenter(1) = (xEdge(2) + xEdge(1))/2D0
        domainCenter(2) = (yEdge(2) + yEdge(1))/2D0

        write(*,*) 'Domain is: [',xEdge(1),',',xEdge(2),'].'
        write(*,*) 'Warning: Not all tests have been implemented for non-unit square domain!'
        IF(doTimeTest) THEN
            write(*,*) 'Warning: doing timing tests overwrites number of outputs!'
        ENDIF

        nmethod_final = 1
        tmp_method = 0
        tmp_method(1) = 1 ! Split modal DG, no limiting
        tmp_method(1) = 2 ! Split modal DG, mass redistribution limiting for positivity
        tmp_method(3) = 3 ! Split modal DG, evaluate at GLL points for comparison
        tmp_method(4) = 4 ! Split modal DG, evaluate at GLL points positivity redistribution

        DO nmethod = 1,nmethod_final
            imethod = tmp_method(nmethod)

            SELECT CASE(imethod)
    		    CASE(1)
        		  write(*,*) 'DG, averages, no limiting'
		          write(*,*) 'WARNING: Should only be used with periodic BCs'
        		  doposlimit = .false.
              dogllGrid = .false.
		          outdir = '_dgnolim/'
              nmethodx = 99
              nmethody = 99
              nOrder = polyDegree
              write(*,FMT='(A5,i1)') ' N = ',nOrder
            CASE(2)
		          write(*,*) 'DG, averages, element mass redist'
		          write(*,*) 'WARNING: Should only be used with periodic BCs'
		          doposlimit = .true.
              dogllGrid = .false.
		          outdir = '_dgmfill/'
		          nmethodx = 99
		          nmethody = 99
              nOrder = polyDegree
		          write(*,FMT='(A5,i1)') ' N = ',nOrder
            CASE(3)
              write(*,*) 'DG, solution at nodal values, no limiting'
              write(*,*) 'WARNING: Should only be used with periodic BCs'
              doposlimit = .false.
              dogllGrid = .true.
              outdir = '_dgnodes/'
              nmethodx = 99
              nmethody = 99
              nOrder = polyDegree
              write(*,FMT='(A5,i1)') ' N = ',nOrder
            CASE(4)
              write(*,*) 'DG, solution at nodal values, positivity limiting'
              write(*,*) 'WARNING: Should only be used with periodic BCs'
              doposlimit = .true.
              dogllGrid = .true.
              outdir = '_posnods/'
              nmethodx = 99
              nmethody = 99
              nOrder = polyDegree
              write(*,FMT='(A5,i1)') ' N = ',nOrder
            END SELECT !imethod

            ! Fill in quadrature nodes and weights
            ALLOCATE(quadNodes(0:nOrder),quadWeights(0:nOrder),legendreVal(0:nOrder,0:nOrder),&
                     legendreDeriv(0:nOrder,0:nOrder),posWeight(0:nOrder),STAT=ierr)

            CALL quad_nodes(nOrder+1,quadNodes)
            CALL quad_weights(nOrder+1,quadNodes,quadWeights)

        		! Fill array of Legendre polynomials evaluated at quad nodes + Leg. derivative at quad nodes
        		DO i=0,nOrder
        			DO j=0,nOrder
        				legendreVal(i,j) = legendre(quadNodes(j),i)
        				legendreDeriv(i,j) = dlegendre(quadNodes(j),i)
        			ENDDO !j
        		ENDDO !i

            posWeight = 1D0
            DO p=1,nlevel
                CALL cpu_time(t0)

                nex = nex0*nscale**(p-1) ! Number of x elements
                ney = ney0*nscale**(p-1)
                nx = nex*(nOrder+1) ! Number of x subcells (for plotting final solution)
                ny = ney*(nOrder+1)

                dxel = (xEdge(2)-xEdge(1))/DBLE(nex)
	              dyel = (yEdge(2)-yEdge(1))/DBLE(ney)
                dxPlot = (xEdge(2)-xEdge(1))/DBLE(nx)
                dyPlot = (yEdge(2)-yEdge(1))/DBLE(ny)

                ALLOCATE(xPlot(1:nx),yPlot(1:ny),DGx(1:nx),DGy(1:ny),elemCenterX(1:nex),elemCenterY(1:ney),&
                        u0(1:nx,1:ny),v0(1:nx,1:ny),uEdge0(1:nex,1:ny),vEdge0(1:nx,1:ney),&
                        STAT=ierr)

                ALLOCATE(q(1:nx,1:ny),q0(1:nx,1:ny),avgXferOp(0:nOrder,0:nOrder),avgXferOpLU(0:nOrder,0:nOrder),&
                         IPIV(0:nOrder),FOO(0:nOrder),&
                         STAT=ierr)

                ! =====================================================================================================
                ! Set up grids
                ! =====================================================================================================
                ! Set up DG grids for ICs
                elemCenterX(1) = xEdge(1)+0.5D0*dxel
                DO i=2,nex
                    elemCenterX(i) = elemCenterX(i-1)+dxel
                ENDDO !i

                elemCenterY(1) = yEdge(1)+0.5D0*dyel
                DO i=2,ney
                    elemCenterY(i) = elemCenterY(i-1)+dyel
                ENDDO !i

                DO i=1,nex
                    DGx(1+(i-1)*(nOrder+1):i*(nOrder+1)) = elemCenterX(i)+0.5D0*dxel*quadNodes(0:nOrder)
                ENDDO !i

                DO i=1,ney
                    DGy(1+(i-1)*(nOrder+1):i*(nOrder+1)) = elemCenterY(i)+0.5D0*dyel*quadNodes(0:nOrder)
                ENDDO !i

                ! Set up plotting grids
                IF(dogllGrid) THEN
                ! Output solution on tensor product of GLL nodes
                    ALLOCATE(gllNodes(0:nOrder),gllWeights(0:nOrder),stat=ierr)
                    CALL gllquad_nodes(nOrder,gllNodes)
                    CALL gllquad_weights(nOrder,gllNodes,gllWeights)
                    posWeight = gllWeights

                    DO j=1,nex
                        xPlot(1+(nOrder+1)*(j-1):(nOrder+1)*j) = 0.5D0*dxel*gllNodes(:)+elemCenterX(j)
                    ENDDO !j
                    DO j=1,nex
                        yPlot(1+(nOrder+1)*(j-1):(nOrder+1)*j) = 0.5D0*dyel*gllNodes(:)+elemCenterY(j)
                    ENDDO !j

                ELSE
                ! Use equally spaced subcells
                    xPlot(1) = xEdge(1)+0.5D0*dxPlot
                    DO j=2,nx
                        xPlot(j) = xPlot(j-1)+dxPlot
                    ENDDO!j
                    yPlot(1) = yEdge(1)+0.5D0*dyPlot
                    DO j=2,ny
                        yPlot(j) = yPlot(j-1)+dyPlot
                    ENDDO!j
                ENDIF !dogllGrid


                ! =====================================================================================================
                ! Fill in operator used to transfer between subcell averages on plotting grid and DG modal coeffs
                ! =====================================================================================================
                IF(dogllGrid) THEN
                    CALL Cmat_FILL(nOrder,gllNodes,gllWeights,dxPlot,dxel,avgXferOp, 'gll')   ! Assumes number of GLL points matches DOFs. Also has repeated points along bdry
                ELSE
                    CALL Cmat_FILL(nOrder,quadNodes,quadWeights,dxPlot,dxel,avgXferOp,'avgs') ! Assumes an evenly spaced sub-grid
                ENDIF

            	  ! Compute LU decomposition of avgXferOp, stored in avgXferOpLU
            	  avgXferOpLU = avgXferOp
            	  FOO = 0D0
            	  CALL DGESV(nOrder+1,1,avgXferOpLU,nOrder+1,IPIV,FOO,nOrder+1,ierr)

                ! =====================================================================================================
                ! Compute ICs for q,u, and v
                ! =====================================================================================================
                CALL init2d(ntest,nx,ny,xPlot,yPlot,nex,ney,DGx,DGy,elemCenterX,elemCenterY,dxel,dyel,&
                            q0,u0,v0,uEdge0,vEdge0,tfinal,cdf_out)
                q = q0

                ! =====================================================================================================
                ! Set up time step
                ! =====================================================================================================
                time = 0D0

                dxm = dxel
                dym = dyel

                tmp_umax = MAXVAL(u0)
                tmp_vmax = MAXVAL(v0)

                IF(noutput .eq. -1) THEN
                    nstep = CEILING( tfinal*MAX(tmp_umax/dxel,tmp_vmax/dyel)/maxcfl )
                    nout = nstep
                ELSE
                    nstep = noutput*CEILING( tfinal*MAX(tmp_umax/dxel,tmp_vmax/dyel)/maxcfl/DBLE(noutput) )
                    nout = noutput
                ENDIF

                IF(doTimeTest) THEN
                    nout = 1
                    nstep = 800*nscale**(p-1)
                ENDIF

                dt = tfinal/DBLE(nstep)
                calculatedMu = MAX(tmp_umax/dxm,tmp_vmax/dym)*dt
                write(*,*) 'Mu used=',calculatedMu

                IF(p.eq.1) THEN ! Set up netCDF output file
                    write(*,*) 'Maximum velocity: |u| = ',maxval(sqrt(u0**2+v0**2))
                    ! Make output file to capture intputs and screen output
                    cdf_out = outdir//TRIM(cdf_out)
                    OPEN(unit=outUnit,file=TRIM(cdf_out)//'.out',status="REPLACE",action="WRITE")
                    WRITE(outUnit,NML=inputs)
                    WRITE(outUnit,FMT='(A)') ''
                    WRITE(outUnit,FMT='(A)') '======================'

                    cdf_out = TRIM(cdf_out)//'.nc'
                    CALL output2d(q,xPlot,yPlot,quadWeights,quadNodes,nex,ney,nOrder,nx,ny,tfinal,calculatedMu,&
                                  cdf_out,nout,-1)
                ENDIF
                ! Set up variables for this value of p ; Write x, y, and initial conditions
                CALL output2d(q,xPlot,yPlot,quadWeights,quadNodes,nex,ney,nOrder,nx,ny,tfinal,calculatedMu,&
                              cdf_out,p,0)

                ! =====================================================================================================
                ! Time integration
                ! =====================================================================================================
                tmp_qmax = MAXVAL(q)
                tmp_qmin = MINVAL(q)

                maxTime = 0D0
                minTime = 1D0
                totTime = 0D0

                oddstep = .TRUE.
                DO n=1,nstep
                    CALL CPU_TIME(t1)
                    CALL strangSplitUpdate(q,u0,v0,uEdge0,vEdge0,quadNodes,quadWeights,time,&
                                 legendreVal,legendreDeriv,avgXferOp,avgXferOpLU,IPIV,&
                                 dt,dxel,dyel,nOrder,nx,ny,nex,ney,oddstep,doposlimit,maxTime,minTime,totTime,posWeight)
                    CALL CPU_TIME(t2)
                    !write(*,*) 'Post-step time:',t2-t1

                    time = time + dt

		              IF((MOD(n,nstep/nout).eq.0).OR.(n.eq.nstep)) THEN ! Write output variables
			              CALL output2d(q,xPlot,yPlot,quadWeights,quadNodes,nex,ney,nOrder,nx,ny,tfinal,calculatedMu,&
                                  cdf_out,p,2)
                  ENDIF

                    tmp_qmin = MIN(tmp_qmin,MINVAL(q(1:nx,1:ny)))
                    tmp_qmax = MAX(tmp_qmax,MAXVAL(q(1:nx,1:ny)))

                    oddstep = .not. oddstep

                ENDDO !n
                CALL cpu_time(tf)
                tf = tf - t0

                IF(dogllGrid) THEN
                  CALL computeErrors(q,q0,dxel,dyel,gllWeights,e1,e2,ei,p,nlevel,nx,ny,nex,ney,nOrder,'gll')
                ELSE
                  CALL computeErrors(q,q0,dxel,dyel,gllWeights,e1,e2,ei,p,nlevel,nx,ny,nex,ney,nOrder,'avgs')
                ENDIF

                IF(p.eq.1) THEN
                    write(UNIT=6,FMT='(A119)') &
' nex  ney       E1          E2         Einf   convergence rate overshoot   undershoot    cons     cputime    step    tf '

                    write(unit=outUnit,FMT='(A)') &
' nex  ney       E1          E2         Einf   convergence rate overshoot   undershoot    cons     cputime    step    tf '
                    cnvg1 = 0.d0
                    cnvg2 = 0.d0
                    cnvgi = 0.d0
                ELSE
                    cnvg1 = -log(e1(p)/e1(p-1))/log(dble(nscale))
                    cnvg2 = -log(e2(p)/e2(p-1))/log(dble(nscale))
                    cnvgi = -log(ei(p)/ei(p-1))/log(dble(nscale))
                ENDIF
                IF(dogllGrid) THEN
                    ! Estimate conservation using GLL quadrature
                    cons = 0D0
                    vertBnd(2) = -1
                    DO i=1,nex
                        horizBnd(1) = 1+(nOrder+1)*(i-1)
                        horizBnd(2) = (nOrder+1)*i
                        DO j=1,ney
                            DO l=0,nOrder
                                vertBnd(1) = 1+(nOrder+1)*(j-1)+l
                                cons = cons + 0.25D0*dxel*dyel*gllWeights(l)*&
                                SUM( gllWeights(:)*(q(horizBnd(1):horizBnd(2),vertBnd(1))-q0(horizBnd(1):horizBnd(2),vertBnd(1))) )
                            ENDDO
                        ENDDO
                    ENDDO
                    cons = cons/DBLE(nx*ny)
                ELSE
                    ! Estimate conservation using subcell averages
                    cons = SUM(q(1:nx,1:ny)-q0)/DBLE(nx*ny)
                ENDIF
                write(*,FMT=990) nex, ney, e1(p), e2(p), ei(p), &
                             cnvg1, cnvg2, cnvgi, &
                             tmp_qmax-MAXVAL(q0), &
                             MINVAL(q0)-tmp_qmin, &
                             cons, tf, nstep,tfinal
                write(unit=outUnit,FMT=990) nex, ney, e1(p), e2(p), ei(p), &
                             cnvg1, cnvg2, cnvgi, &
                             tmp_qmax-MAXVAL(q0), &
                             MINVAL(q0)-tmp_qmin, &
                             cons, tf, nstep,tfinal
!write(*,'(3f8.2)') minTime, maxTime, totTime/tf

				IF(p .eq. nlevel) THEN
                    CALL output2d(q,xPlot,yPlot,quadWeights,quadNodes,nex,ney,nOrder,nx,ny,tfinal,calculatedMu,&
                              cdf_out,p,1) ! Close netCDF files
				ENDIF
                DEALLOCATE(xPlot,yPlot,DGx,DGy,elemCenterX,elemCenterY,u0,v0,uEdge0,vEdge0,q,q0,avgXferOp,avgXferOpLU,&
                           IPIV,FOO,&
                           STAT=ierr)

            ENDDO !p
            DEALLOCATE(quadNodes,quadWeights,legendreVal,legendreDeriv,posWeight)
            IF(dogllGrid) DEALLOCATE(gllNodes,gllWeights,STAT=ierr)
            CLOSE(outUnit)
        ENDDO ! nmethod

990    format(2i5,3e12.4,3f5.2,3e12.4,f8.2,i8,f8.2)

    END SUBROUTINE test2d_modal

    SUBROUTINE computeErrors(soln,ics,dxel,dyel,quadWeight,e1,e2,ei,p,nlvl,nx,ny,nex,ney,N,stat)
    ! ===========================================================================================
    ! computeErrors() computes the error at the final time of the integration in one of two ways:
    ! soln and ics arrays are assumed to contain approximations to the final and initial solutions
    ! the nature of these approximations depends on stat string:
    !   1. stat = 'gll'
    !     approximations are nodal values on a GLL quadrature grid with N+1 nodes along each coordinate per element.
    !     L1 and L2 errors are estimated using GLL quadrature.
    !
    !   2. stat != 'gll'
    !     approximations are cell averages using a uniform cartesian mesh with N+1 cells along each coordinate per element.
    !     L1 and L2 errors are estimated using grid norms.
    ! ===========================================================================================
      IMPLICIT NONE
      ! Inputs
  		INTEGER, INTENT(IN) :: p,nx,ny,nlvl,N,nex,ney
  		REAL(KIND=8), INTENT(IN) :: dxel,dyel
  		REAL(KIND=8), DIMENSION(0:N), INTENT(IN) :: quadWeight
  		REAL(KIND=8), DIMENSION(1:nx,1:ny), INTENT(IN) :: soln,ics
  		CHARACTER(len=*), INTENT(IN) :: stat
  		! Outputs
  		REAL(KIND=8), DIMENSION(1:nlvl), INTENT(OUT) :: e1,e2,ei
  		! Local Variables
  		INTEGER :: i,j,k
      REAL(KIND=8) :: errijk
  		INTEGER, DIMENSION(1:2) :: hBnd,vBnd

  		IF(TRIM(stat) .eq. 'gll') THEN
  		! Compute errors using gll quadrature estimation
  			e1(p) = 0D0
  			e2(p) = 0D0
  			ei(p) = MAXVAL(ABS(soln(:,:)-ics(:,:)))

        vBnd(2) = -1
        errijk = 0D0
  			DO i=1,nex
          hBnd(1) = 1+(i-1)*(N+1) ! Lower bound of element index
          hBnd(2) = i*(N+1)				! Upper bound of element index
          DO j=1,ney
            DO k=0,N
              vBnd(1) = 1+k+(j-1)*(N+1)

              errijk = SUM(quadWeight*quadWeight(k)*ABS(soln(hBnd(1):hBnd(2),vBnd(1))-ics(hBnd(1):hBnd(2),vBnd(1))))
              e1(p) = e1(p) + errijk

              errijk = SUM(quadWeight*quadWeight(k)*ABS(soln(hBnd(1):hBnd(2),vBnd(1))-ics(hBnd(1):hBnd(2),vBnd(1)))**2)
              e2(p) = e2(p) + errijk

            ENDDO
          ENDDO
  			ENDDO
  			e1(p) = 0.25D0*dxel*dyel*e1(p)
  			e2(p) = SQRT(0.25D0*dxel*dyel*e2(p))

  		ELSE
  		! Compute errors using equal area cells estimation
      e1(p) = SUM(ABS(soln(1:nx,1:ny)-ics))/DBLE(nx)/DBLE(ny)
      e2(p) = SQRT(SUM((soln(1:nx,1:ny)-ics)**2)/DBLE(nx)/DBLE(ny))
      ei(p) = MAXVAL(ABS(soln(1:nx,1:ny)-ics))

  		ENDIF
    END SUBROUTINE computeErrors

    SUBROUTINE strangSplitUpdate(q,u0,v0,uEdge0,vEdge0,quadNodes,quadWeights,time,&
                                 legendreVal,legendreDeriv,avgXferOp,avgXferOpLU,IPIV,&
                                 dt,dxel,dyel,nOrder,nx,ny,nex,ney,oddstep,doposlimit,&
                                 maxTime,minTime,totTime,posWeight)
    ! =====================================================================================================
    ! strangSplitUpdate is responsible for selecting which slice of subcell volumes is sent to mDGsweep for update to time
    ! level tn+1 following a Strang splitting.
    ! For Strang splitting:
    !   - Each slice is updated
    !   - Odd steps: x-slices are updated first (horizontal advection) then y-slices are updated (vertical advection)
    !   - Even steps: y-slices are updated first then x-slices are updated (vertical advection)
    ! =====================================================================================================

        IMPLICIT NONE
        ! Inputs
        REAL(KIND=8), EXTERNAL :: tfcn
        INTEGER, INTENT(IN) :: nOrder,nx,ny,nex,ney
        REAL(KIND=8), INTENT(IN) :: dt,dxel,dyel,time
        REAL(KIND=8), DIMENSION(1:nx,1:ny), INTENT(IN) :: u0,v0
        REAL(KIND=8), DIMENSION(1:nex,1:ny), INTENT(IN) :: uEdge0
        REAL(KIND=8), DIMENSION(1:nx,1:ney), INTENT(IN) :: vEdge0
        REAL(KIND=8), DIMENSION(0:nOrder), INTENT(IN) :: quadNodes,quadWeights,posWeight
        REAL(KIND=8), DIMENSION(0:nOrder,0:nOrder), INTENT(IN) :: legendreVal,legendreDeriv,avgXferOp,avgXferOpLU
        INTEGER, DIMENSION(0:nOrder), INTENT(IN) :: IPIV
        LOGICAL, INTENT(IN) :: oddstep,doposlimit
        ! Outputs
        REAL(KIND=8), DIMENSION(1:nx,1:ny), INTENT(INOUT) :: q
        REAL(KIND=8), INTENT(INOUT) :: maxTime, minTime, totTime
        ! Local variables
        INTEGER :: i,j,k
        REAL(KIND=8) :: t_temp,time_factor
        REAL(KIND=8), DIMENSION(1:3,1:nx,1:ny) :: u,v
        REAL(KIND=8), DIMENSION(1:3,1:nex,1:ny) :: uEdge
        REAL(KIND=8), DIMENSION(1:3,1:nx,1:ney) :: vEdge
        REAL(KIND=8), DIMENSION(1:nx) :: q1dx
        REAL(KIND=8), DIMENSION(1:ny) :: q1dy
        REAL(KIND=8), DIMENSION(1:3,1:nx) :: u1dx
        REAL(KIND=8), DIMENSION(1:3,1:ny) :: v1dy
        REAL(KIND=8), DIMENSION(1:3,1:nex) :: uEdge1dx
        REAL(KIND=8), DIMENSION(1:3,1:ney) :: vEdge1dy

        ! Update velocities at times required for ssprk3 update
        u(1,:,:) = u0
        u(2,:,:) = u0
        u(3,:,:) = u0

        uEdge(1,:,:) = uEdge0
        uEdge(2,:,:) = uEdge0
        uEdge(3,:,:) = uEdge0

        v(1,:,:) = v0
        v(2,:,:) = v0
        v(3,:,:) = v0

        vEdge(1,:,:) = vEdge0
        vEdge(2,:,:) = vEdge0
        vEdge(3,:,:) = vEdge0

        IF(transient) THEN
            DO i=1,3
                SELECT CASE(i)
                CASE(1)
                    t_temp = time
                CASE(2)
                    t_temp = time+dt
                CASE(3)
                    t_temp = time+0.5D0*dt
                END SELECT
                time_factor = tfcn(t_temp)
                u(i,:,:) = u(i,:,:)*time_factor
                v(i,:,:) = v(i,:,:)*time_factor
                uEdge(i,:,:) = uEdge(i,:,:)*time_factor
                vEdge(i,:,:) = vEdge(i,:,:)*time_factor
            ENDDO !i
        ENDIF !transient
        IF(oddstep) THEN
            ! ===================================
            ! Perform sweeps in x-direction first
            ! ===================================
            DO j=1,ny
                q1dx = q(:,j)
                u1dx(1:3,:) = u(1:3,:,j)
                uEdge1dx(1:3,:) = uEdge(1:3,:,j)
                CALL mDGsweep(q1dx,u1dx,uEdge1dx,dxel,nex,nOrder,quadWeights,avgXferOp,avgXferOpLU, &
                              legendreVal,legendreDeriv,IPIV,dt,doposlimit,posWeight,maxTime,minTime,&
                              totTime)
                ! Update solution
                q(:,j) = q1dx
            ENDDO!j

            DO i=1,nx
                q1dy = q(i,:)
                v1dy(1:3,:) = v(1:3,i,:)
                vEdge1dy(1:3,:) = vEdge(1:3,i,:)

                CALL mDGsweep(q1dy,v1dy,vEdge1dy,dyel,ney,nOrder,quadWeights,avgXferOp,avgXferOPLU,&
                              legendreVal,legendreDeriv,IPIV,dt,doposlimit,posWeight,maxTime,minTime,&
                              totTime)
                ! Update solution
                q(i,:) = q1dy
            ENDDO !i

        ELSE
            ! ===================================
            ! Perform sweeps in y-direction first
            ! ===================================
            DO i=1,nx
                q1dy = q(i,:)
                v1dy(1:3,:) = v(1:3,i,:)
                vEdge1dy(1:3,:) = vEdge(1:3,i,:)

                CALL mDGsweep(q1dy,v1dy,vEdge1dy,dyel,ney,nOrder,quadWeights,avgXferOp,avgXferOPLU,&
                              legendreVal,legendreDeriv,IPIV,dt,doposlimit,posWeight,maxTime,minTime,&
                              totTime)
                ! Update solution
                q(i,:) = q1dy
            ENDDO !i

            DO j=1,ny
                q1dx = q(:,j)
                u1dx(1:3,:) = u(1:3,:,j)
                uEdge1dx(1:3,:) = uEdge(1:3,:,j)

                CALL mDGsweep(q1dx,u1dx,uEdge1dx,dxel,nex,nOrder,quadWeights,avgXferOp,avgXferOpLU, &
                              legendreVal,legendreDeriv,IPIV,dt,doposlimit,posWeight,maxTime,minTime,&
                              totTime)
                ! Update solution
                q(:,j) = q1dx
            ENDDO!j
        ENDIF !oddstep
    END SUBROUTINE strangSplitUpdate

    SUBROUTINE init2d(ntest,nx,ny,xPlot,yPlot,nex,ney,DGx,DGy,elemCenterX,elemCenterY,dxel,dyel,&
                      q0,u0,v0,uEdge0,vEdge0,tfinal,cdf_out)
        IMPLICIT NONE
        ! Inputs
        INTEGER, INTENT(IN) :: ntest,nx,ny,nex,ney
        REAL(KIND=8), INTENT(IN) :: dxel,dyel
        REAL(KIND=8), DIMENSION(1:nx), INTENT(IN) :: xPlot,DGx
        REAL(KIND=8), DIMENSION(1:ny), INTENT(IN) :: yPlot,DGy
        REAL(KIND=8), DIMENSION(1:nex), INTENT(IN) :: elemCenterX
        REAL(KIND=8), DIMENSION(1:ney), INTENT(IN) :: elemCenterY

        ! Outputs
        REAL(KIND=8), INTENT(OUT) :: tfinal
        REAL(KIND=8), DIMENSION(1:nx,1:ny), INTENT(OUT) :: q0,u0,v0
        REAL(KIND=8), DIMENSION(1:nex,1:ny), INTENT(OUT) :: uEdge0
        REAL(KIND=8), DIMENSION(1:nx,1:ney), INTENT(OUT) :: vEdge0
        CHARACTER(LEN=40), INTENT(OUT) :: cdf_out
        ! Local Variables
        INTEGER :: i,j,l
        REAL(KIND=8) :: PI,waveSpd,dxPlot,dyPlot,x0,y0
        REAL(KIND=8), DIMENSION(1:nx,0:1) :: xtilde
        REAL(KIND=8), DIMENSION(1:ny,0:1) :: ytilde
        REAL(KIND=8), DIMENSION(1:nx,1:ny,0:1) :: psiu,psiv
        REAL(KIND=8), DIMENSION(1:nex,1:ny,0:1) :: psiuEdge
        REAL(KIND=8), DIMENSION(1:nx,1:ney,0:1) :: psivEdge

        REAL(KIND=8), DIMENSION(1:nx,1:ny) :: r


        PI = DACOS(-1D0)
        waveSpd = 1D0

        dxPlot = xPlot(2)-xPlot(1)
        xtilde(:,0) = xPlot(:)-0.5D0*dxPlot
        xtilde(:,1) = xPlot(:)+0.5D0*dxPlot

        dyPlot = yPlot(2)-yPlot(1)
        ytilde(:,0) = yPlot(:)-0.5D0*dyPlot
        ytilde(:,1) = yPlot(:)+0.5D0*dyPlot

        SELECT CASE(ntest)
            CASE(0) ! Consistency test: uniform field in deformation flow
                cdf_out =  'spltMod2d_consistency'
                waveSpd = 1D0
                tfinal = 5D0
                q0 = 1D0

            CASE(1) ! Constant sine adv
                cdf_out =  'spltMod2d_adv_sine'
                waveSpd = 1D0
                tfinal = 1D0
                DO j=1,ny
                    q0(:,j) = sin(2.d0*PI*xPlot(:))*sin(2.d0*PI*yPlot(j))
                ENDDO !j
            CASE(2) ! Solid body rotation of cylinder
                cdf_out = 'spltMod2d_rot_cylinder'
                waveSpd = 2D0*PI
                tfinal = 1D0
                DO j=1,ny
                    r(:,j) = SQRT((xPlot-0.3D0)**2 + (yPlot(j)-0.3D0)**2)
                ENDDO !j
                q0 = 0D0
                WHERE(r .lt. 0.125D0)
                    q0 = 1D0
                END WHERE

            CASE(5) ! Cosbell deformation from LeVeque
                cdf_out = 'spltMod2d_def_cosinebell'
                waveSpd = 1D0
                tfinal = 5D0

                DO j=1,ny
                    r(:,j) = 4D0*SQRT( (xPlot-0.25D0)**2 + (yPlot(j)-0.25D0)**2 )
                ENDDO !j
                q0 = 0D0
                WHERE(r .lt. 1D0)
                    q0 = 0.25D0*(1D0+DCOS(PI*r))**2
!                    q0 = 0.5D0*(1D0+DCOS(PI*r))
                END WHERE

            CASE(6) ! Smoother cosbell deformation flow (LeVeque)
                cdf_out = 'spltMod2d_def_smth_cosbell'
                waveSpd = 1D0
                tfinal = 5D0
                DO j=1,ny
                    r(:,j) = 4D0*SQRT( (xPlot-0.25D0)**2 + (yPlot(j)-0.25D0)**2 )
                ENDDO !j
                q0 = 0D0
                WHERE(r .lt. 1D0)
                    q0 = (0.5D0*(1D0+DCOS(PI*r)))**3
                END WHERE

            CASE(7) ! Slotted cylinder in deformation flow
                cdf_out = 'spltMod2d_def_cyl'
                waveSpd = 1D0
                tfinal = 5D0

                x0 = 0.25D0
                y0 = 0.5D0
                DO j=1,ny
                    r(:,j) = SQRT((xPlot-x0)**2 + (yPlot(j)-y0)**2)
                ENDDO !j
                q0 = 0D0
                WHERE(r .lt. .15D0)
                    q0 = 1D0
                END WHERE

                DO j=1,ny
                    DO i=1,nx
                        IF(ABS(xPlot(i)-x0) .lt. 0.025D0 .AND. yPlot(j) .gt.(y0-0.0625D0)) THEN
                            q0(i,j) = 0D0
                        ENDIF
                    ENDDO !i
                ENDDO !j
        END SELECT

        SELECT CASE(ntest)
            CASE(1) ! uniform diagonal advection of a sine wave
              ! Evaluate stream function for horizontal velocities
              DO j=1,ny
                psiu(:,j,0) = -DGx(:) + ytilde(j,0)
                psiu(:,j,1) = -DGx(:) + ytilde(j,1)
                psiuEdge(:,j,0) = -(elemCenterX(:)+0.5D0*dxel) + ytilde(j,0)
                psiuEdge(:,j,1) = -(elemCenterX(:)+0.5D0*dxel) + ytilde(j,1)
              ENDDO!j

          		! Evaluate stream function for vertical velocities
              DO i=1,nx
                psiv(i,:,0) = -xtilde(i,0)+DGy(:)
                psiv(i,:,1) = -xtilde(i,1)+DGy(:)
                psivEdge(i,:,0) = -xtilde(i,0)+(elemCenterY(:)+0.5D0*dyel)
                psivEdge(i,:,1) = -xtilde(i,1)+(elemCenterY(:)+0.5D0*dyel)
              ENDDO!i

            CASE(2) ! Solid body rotation
              ! Evaluate stream function for horizontal velocities
              DO j=1,ny
                psiu(:,j,0) = 0.5D0*waveSpd*( (DGx(:)-0.5D0)**2+(ytilde(j,0)-0.5D0)**2)
                psiu(:,j,1) = 0.5D0*waveSpd*( (DGx(:)-0.5D0)**2+(ytilde(j,1)-0.5D0)**2)
                psiuEdge(:,j,0) = 0.5D0*waveSpd*( (elemCenterX(:)+0.5D0*dxel-0.5D0)**2 + (ytilde(j,0)-0.5D0)**2)
                psiuEdge(:,j,1) = 0.5D0*waveSpd*( (elemCenterX(:)+0.5D0*dxel-0.5D0)**2 + (ytilde(j,1)-0.5D0)**2)
              ENDDO!j

              ! Evaluate stream function for vertical velocities
              DO i=1,nx
                psiv(i,:,0) = 0.5D0*waveSpd*((xtilde(i,0)-0.5D0)**2+(DGy(:)-0.5D0)**2)
                psiv(i,:,1) = 0.5D0*waveSpd*((xtilde(i,1)-0.5D0)**2+(DGy(:)-0.5D0)**2)
                psivEdge(i,:,0) = 0.5D0*waveSpd*((xtilde(i,0)-0.5D0)**2+(elemCenterY(:)+0.5D0*dyel-0.5D0)**2)
                psivEdge(i,:,0) = 0.5D0*waveSpd*((xtilde(i,1)-0.5D0)**2+(elemCenterY(:)+0.5D0*dyel-0.5D0)**2)
              ENDDO!i

            CASE(0,5:7) ! LeVeque deformation flow

              ! Evaluate stream function for horizontal velocities (1/pi)*sin(pi*xf(i))**2 * sin(pi*yf(j))**2
              DO j=1,ny
                psiu(:,j,0) = (SIN(PI*DGx(:))**2 * SIN(PI*ytilde(j,0))**2 )/PI
                psiu(:,j,1) = (SIN(PI*DGx(:))**2 * SIN(PI*ytilde(j,1))**2 )/PI
                psiuEdge(:,j,0) = (SIN(PI*(elemCenterX(:)+0.5D0*dxel))**2 * SIN(PI*ytilde(j,0))**2)/PI
                psiuEdge(:,j,1) = (SIN(PI*(elemCenterX(:)+0.5D0*dxel))**2 * SIN(PI*ytilde(j,1))**2)/PI
              ENDDO!j

              ! Evaluate stream function for vertical velocities
              DO i=1,nx
                psiv(i,:,0) = (SIN(PI*xtilde(i,0))**2 * SIN(PI*DGy(:))**2)/PI
                psiv(i,:,1) = (SIN(PI*xtilde(i,1))**2 * SIN(PI*DGy(:))**2)/PI
                psivEdge(i,:,0) = (SIN(PI*xtilde(i,0))**2 * SIN(PI*(elemCenterY(:)+0.5D0*dyel))**2)/PI
                psivEdge(i,:,1) = (SIN(PI*xtilde(i,1))**2 * SIN(PI*(elemCenterY(:)+0.5D0*dyel))**2)/PI
              ENDDO!i

        END SELECT !ntest

        ! Compute u velocities from stream function
        DO j=1,ny
          u0(:,j) = (psiu(:,j,1)-psiu(:,j,0))/dyPlot
          uEdge0(:,j) = (psiuEdge(:,j,1)-psiuEdge(:,j,0))/dyPlot
        ENDDO!j

        ! Compute v velocities from stream function
        DO i=1,nx
          v0(i,:) = -1D0*(psiv(i,:,1)-psiv(i,:,0))/dxPlot
          vEdge0(i,:) = -1D0*(psivEdge(i,:,1)-psivEdge(i,:,0))/dxPlot
        ENDDO!i

    END SUBROUTINE init2d

	SUBROUTINE output2d(q,xPlot,yPlot,quadWeights,quadNodes,nex,ney,nOrder,nx,ny,tval_in,mu,cdf_out,ilvl,stat)
		IMPLICIT NONE

		! Inputs
		INTEGER, INTENT(IN) :: nOrder,nex,ney,nx,ny,stat,ilvl
		CHARACTER(len=40), INTENT(IN) :: cdf_out
		REAL(KIND=8), INTENT(IN) :: tval_in,mu
		REAL(KIND=8), DIMENSION(1:nx), INTENT(IN) :: xPlot
		REAL(KIND=8), DIMENSION(1:ny), INTENT(IN) :: yPlot
		REAL(KIND=8), DIMENSION(1:nx,1:ny), INTENT(IN) :: q
        REAL(KIND=8), DIMENSION(0:nOrder), INTENT(IN) :: quadWeights,quadNodes

		! Outputs

		! Local variables
		INTEGER :: cdfid ! ID for netCDF file
		INTEGER, PARAMETER :: NDIMS = 3
		INTEGER :: ierr
	    INTEGER :: idq,idt,idx,idy,dimids(NDIMS),idweight,idnode,idmu
	    INTEGER :: x_dimid, y_dimid, t_dimid,node_dimid
		INTEGER, DIMENSION(1:NDIMS) :: start, count
		CHARACTER(len=8) :: nxname,xname,nyname,yname,qname,muname

		REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: tmp
        REAL(KIND=8), DIMENSION(1:nx) :: qSlice
		INTEGER :: i,j,l,m,ylvl

	    SAVE cdfid, idq, t_dimid, start, count
		IF(stat .eq. -1) THEN
			! Create netCDF file and time variables
			ierr = NF90_CREATE(TRIM(cdf_out),NF90_CLOBBER,cdfid)

			ierr = NF90_REDEF(cdfid)
			ierr = NF90_DEF_DIM(cdfid, "nt", ilvl+1, t_dimid)
            ierr = NF90_DEF_DIM(cdfid, "nnodes", nOrder+1, node_dimid)

			ierr = NF90_DEF_VAR(cdfid, "qweights",NF90_FLOAT, node_dimid, idweight)
			ierr = NF90_DEF_VAR(cdfid, "qnodes",NF90_FLOAT, node_dimid, idnode)
			ierr = NF90_DEF_VAR(cdfid, "time", NF90_FLOAT, t_dimid,idt)

			ierr = NF90_ENDDEF(cdfid)

			! Calculate time at output levels (note ilvl=noutput)
			ALLOCATE(tmp(1:ilvl+1), STAT=ierr)
			DO i=0,ilvl
				tmp(i+1) = DBLE(i)*tval_in/DBLE(ilvl)
			ENDDO

			! Write t values
			ierr = NF90_PUT_VAR(cdfid,idt,tmp)
			ierr = NF90_PUT_VAR(cdfid,idweight,quadWeights)
            ierr = NF90_PUT_VAR(cdfid,idnode,quadNodes)

			DEALLOCATE(tmp, STAT=ierr)

			RETURN

		ELSEIF(stat .eq. 0) THEN
			! Create dimensions and variables for this level of runs (ilvl = p)
			start = 1
			count = 1

			! Define names of variables
			WRITE(nxname,'(a2,i1)') 'nx',ilvl
			WRITE(nyname,'(a2,i1)') 'ny',ilvl
			WRITE(xname, '(a1,i1)') 'x',ilvl
			WRITE(yname, '(a1,i1)') 'y',ilvl
			WRITE(qname, '(a1,i1)') 'Q',ilvl
            WRITE(muname, '(a2,i1)') 'mu',ilvl

			ierr = NF90_REDEF(cdfid)

			ierr = NF90_DEF_DIM(cdfid, TRIM(nxname), nx, x_dimid)
			ierr = NF90_DEF_DIM(cdfid, TRIM(nyname), ny, y_dimid)

			dimids(1) = x_dimid
			dimids(2) = y_dimid
			dimids(3) = t_dimid

			ierr = NF90_DEF_VAR(cdfid, TRIM(qname),NF90_FLOAT,dimids,idq)
			ierr = NF90_DEF_VAR(cdfid, TRIM(xname),NF90_FLOAT,x_dimid,idx)
			ierr = NF90_DEF_VAR(cdfid, TRIM(yname),NF90_FLOAT,y_dimid,idy)
            ierr = NF90_DEF_VAR(cdfid, TRIM(muname),NF90_FLOAT,idmu)

			ierr = NF90_enddef(cdfid)

			! Write x and y values
			ierr = NF90_PUT_VAR(cdfid, idx, xPlot)
			ierr = NF90_PUT_VAR(cdfid, idy, yPlot)
            ierr = NF90_PUT_VAR(cdfid,idmu,mu)

			start(3) = 1

		ELSEIF(stat .eq. 1) THEN
			ierr = NF90_CLOSE(cdfid)
			RETURN
		ENDIF

		! Write out concentration field
		count(1) = nx
        DO ylvl=1,ny
            start(2) = ylvl
            qSlice = q(:,ylvl)
            ierr = NF90_PUT_VAR(cdfid,idq,qSlice,start,count)
        ENDDO !ylvl

		! Increment t level
		start(3) = start(3) + 1

	END SUBROUTINE output2d
END PROGRAM EXECUTE
