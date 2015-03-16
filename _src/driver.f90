SUBROUTINE DRIVER(nex0,ney0,nscale,nruns,noutput,maxCFL)
  ! ===============================================================
  ! Main driver subroutine for DG simulations
  ! Inputs:
  !   testID    : which test is being run
  !   nex0,ney0 : number of initial spatial cells
  !   nscale    : multiple of nex0 used for subsequent runs (nex = nex0*nscale**p)
  !   nruns     : number of total runs to make
  !   maxCFL    : maximal CFL number to use throughout integration
  ! ===============================================================
  USE commonTestParameters
  USE mDGmod
  USE netCDF

  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nex0,ney0,nscale,nruns,noutput
  REAL(KIND=8), INTENT(IN) :: maxCFL
  ! Outputs
  ! Local variables
  CHARACTER(len=40) :: cdfOut
  INTEGER, DIMENSION(10) :: tmp_method
  INTEGER :: nmethod,nmethod_final,imethod,ierr,i,j,p,n,m
  INTEGER :: nstep,nout
  REAL(KIND=4) :: t0,tf
  LOGICAL :: oddstep
  DOUBLE PRECISION :: dxel,dyel,dxPlot,dyPlot,tmp_umax,tmp_vmax,calculatedMu,dt,time
  DOUBLE PRECISION, DIMENSION(1:meqn) :: tmp_qmax,tmp_qmin
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: quadNodes,quadWeights,elemCenterX,elemCenterY,&
      xPlot,yPlot,DGx,DGy,FOO
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: u0,v0,uEdge0,vEdge0,&
                                  avgXferOp,avgXferOpLU,legendreVal,legendreDeriv
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: q,q0
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IPIV

  ! =================================================================================
  ! Subroutine Interfaces
  ! =================================================================================

  INTERFACE
    SUBROUTINE init2d(q,u,v,uEdge,vEdge,xPlot,yPlot,quadNodes,dxel,dyel,&
                      dxPlot,dyPlot,elemCenterX,elemCenterY)
      ! ==============================================================================
      ! Computes initial conditions for q,u,v fields
      ! INPUTS: meqn - number of fields to evaluate
      !         nx,ny - number of points to evaluate q,u,v at
      !         quadNodes - local locations to evaluate velocities at
      !         elemCenterX, elemCenterY - element center locations
      !
      ! OUTPUTS: q(i,j,neq) - initial conditions evaluated for neqth field
      !          u(i,j),v(i,j) - velocities evaluated at quadrature locations
      !          uEdge,vEdge - velocities at edges of each element
      ! ==============================================================================
      USE commonTestParameters
      IMPLICIT NONE
      ! Inputs
      DOUBLE PRECISION, INTENT(IN) :: dxel,dyel,dxPlot,dyPlot
      DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadNodes
      DOUBLE PRECISION, DIMENSION(1:nxOut), INTENT(IN) :: xPlot
      DOUBLE PRECISION, DIMENSION(1:nyOut), INTENT(IN) :: yPlot
      DOUBLE PRECISION, DIMENSION(1:nex), INTENT(IN) :: elemCenterX
      DOUBLE PRECISION, DIMENSION(1:ney), INTENT(IN) :: elemCenterY
      ! Outputs
      DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn) :: q
      DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut) :: u,v
      DOUBLE PRECISION, DIMENSION(1:nex,1:nyOut), INTENT(OUT) :: uEdge
      DOUBLE PRECISION, DIMENSION(1:nxOut,1:ney), INTENT(OUT) :: vEdge
    END SUBROUTINE init2d

    SUBROUTINE output2d(q,xOut,yOut,timeOut,muOut,cdfOut,ilvl,stat)
      ! ============================================================================
      ! output2d - Creates netCDF output files and writes out different output fields
      ! INPUTS: q(nx,ny,meqn)
      !         xOut(nx),yOut(ny)
      !         timeOut,muOut
      ! OUTPUTS: -None-
      ! ============================================================================
      USE commonTestParameters
      USE netCDF
      IMPLICIT NONE
      ! Inputs
      INTEGER, INTENT(IN) :: ilvl,stat
      CHARACTER(len=40), INTENT(IN) :: cdfOut
      DOUBLE PRECISION, INTENT(IN) :: muOut,timeOut
      DOUBLE PRECISION, DIMENSION(1:nxOut), INTENT(IN) :: xOut
      DOUBLE PRECISION, DIMENSION(1:nyOut), INTENT(IN) :: yOut
      DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn), INTENT(IN) :: q
      ! Outputs
    END SUBROUTINE output2d

    SUBROUTINE strangSplit(q,u0,v0,uEdge0,vEdge0,quadNodes,quadWeights,time,&
                           legendreVal,legendreDeriv,avgXferOp,avgXferOpLU,IPIV,&
                           dt,dxel,dyel,oddstep)
    ! =====================================================================================================
    ! strangSplitUpdate is responsible for selecting which slice of subcell volumes is sent to mDGsweep for update to time
    ! level tn+1 following a Strang splitting.
    ! For Strang splitting:
    !   - Each slice is updated
    !   - Odd steps: x-slices are updated first (horizontal advection) then y-slices are updated (vertical advection)
    !   - Even steps: y-slices are updated first then x-slices are updated (vertical advection)
    ! =====================================================================================================
        USE commonTestParameters
        IMPLICIT NONE
        ! Inputs
        REAL(KIND=8), INTENT(IN) :: dt,dxel,dyel,time
        REAL(KIND=8), DIMENSION(1:nxOut,1:nyOut), INTENT(IN) :: u0,v0
        REAL(KIND=8), DIMENSION(1:nex,1:nyOut), INTENT(IN) :: uEdge0
        REAL(KIND=8), DIMENSION(1:nxOut,1:ney), INTENT(IN) :: vEdge0
        REAL(KIND=8), DIMENSION(0:nQuad), INTENT(IN) :: quadNodes,quadWeights
        REAL(KIND=8), DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: legendreVal,legendreDeriv,avgXferOp,avgXferOpLU
        INTEGER, DIMENSION(0:maxPolyDegree), INTENT(IN) :: IPIV
        LOGICAL, INTENT(IN) :: oddstep
        ! Outputs
        REAL(KIND=8), DIMENSION(1:nxOut,1:nyOut,1:meqn), INTENT(INOUT) :: q
      END SUBROUTINE strangSplit
  END INTERFACE
  ! =================================================================================
  ! END SUBROUTINE INTERFACES
  ! =================================================================================

  if(nRuns.lt.1) STOP 'nRuns should be at least 1 in DRIVER()'
  PI = DACOS(-1D0)

  ! Set spatial domain
  xDomain(1) = 0D0
  xDomain(2) = 1D0
  yDomain = xDomain

  nmethod_final = 1
  tmp_method = 0
  tmp_method(1) = 1 ! Split modal DG, no limiting
  tmp_method(1) = 2 ! Split modal DG, mass redistribution limiting for positivity

  SELECT CASE(testID)
    CASE(0)
      cdfOut = 'spltMod2d_consistency'
    CASE(1)
      cdfOut = 'spltMod2d_uniform'
    CASE(5)
      cdfOut = 'spltMod2d_def_cosinebell'
    CASE(7)
      cdfOut = 'spltMod2d_def_cyl'
  END SELECT !testID

  DO nmethod = 1,nmethod_final
    imethod = tmp_method(nmethod)

    SELECT CASE(imethod)
      CASE(1)
        write(*,*) 'DG, averages, no limiting'
        doposlimit = .false.
        outdir = '_modal/'
      CASE(2)
        write(*,*) 'DG, averages, element mass redist'
        doposlimit = .true.
        outdir = '_pdModal/'
    END SELECT !imethod

    write(*,FMT='(A5,i1)') ' N = ',maxPolyDegree

    nQuad = maxPolyDegree

    ALLOCATE(quadNodes(0:nQuad),quadWeights(0:nQuad),legendreVal(0:maxPolyDegree,0:nQuad),&
            legendreDeriv(0:maxPolyDegree,0:nQuad),STAT=ierr)

    CALL quad_nodes(nQuad+1,quadNodes)
    CALL quad_weights(nQuad+1,quadNodes,quadWeights)

    ! Fill array of Legendre polynomials evaluated at quad nodes + Leg. derivative at quad nodes
    DO i=0,maxPolyDegree
      DO j=0,nQuad
        legendreVal(i,j) = legendre(quadNodes(j),i)
        legendreDeriv(i,j) = dlegendre(quadNodes(j),i)
      ENDDO !j
    ENDDO !i

    DO p=1,nRuns
      CALL cpu_time(t0)

      nex = nex0*nscale**(p-1) ! Number of x elements
      ney = ney0*nscale**(p-1)
      nxOut = nex*(nQuad+1) ! Number of local subcells for plotting final solution
      nyOut = ney*(nQuad+1)

      ALLOCATE(elemCenterX(1:nex),elemCenterY(1:ney),xPlot(1:nxOut),yPlot(1:nyOut),&
              DGx(1:nxOut),DGy(1:nyOut),q(1:nxOut,1:nyOut,1:meqn),q0(1:nxOut,1:nyOut,1:meqn),&
              u0(1:nxOut,1:nyOut),v0(1:nxOut,1:nyOut),uEdge0(1:nex,1:nyOut),vEdge0(1:nxOut,1:ney),STAT=ierr)

      ALLOCATE(avgXferOp(0:maxPolyDegree,0:maxPolyDegree),avgXferOpLU(0:maxPolyDegree,0:maxPolyDegree),&
               IPIV(0:maxPolyDegree),FOO(0:maxPolyDegree),STAT=ierr)

      ! Create plotting grids
      CALL makeGrid(dxel,dyel,elemCenterX,elemCenterY,dxPlot,dyPlot,xPlot,yPlot)

      ! =====================================================================================================
      ! Fill in operator used to transfer between subcell averages on plotting grid and DG modal coeffs
      ! =====================================================================================================
      CALL Cmat_FILL(maxPolyDegree,quadNodes,quadWeights,dxPlot,dxel,avgXferOp,'avgs') ! Assumes an evenly spaced sub-grid
      ! Compute LU decomposition of avgXferOp, stored in avgXferOpLU
      avgXferOpLU = avgXferOp
      FOO = 0D0
      CALL DGESV(maxPolyDegree+1,1,avgXferOpLU,maxPolyDegree+1,IPIV,FOO,maxPolyDegree+1,ierr)

      ! =====================================================================================================
      ! Initialize q, u, and v arrays.
      ! =====================================================================================================
      CALL init2d(q0,u0,v0,uEdge0,vEdge0,xPlot,yPlot,quadNodes,dxel,dyel,&
                  dxPlot,dyPlot,elemCenterX,elemCenterY)
      q = q0

      ! =====================================================================================================
      ! Set up time step size
      ! =====================================================================================================
      time = 0D0

      tmp_umax = MAXVAL(u0)
      tmp_vmax = MAXVAL(v0)

      IF(noutput .eq. -1) THEN
          nstep = CEILING( tfinal*MAX(tmp_umax/dxel,tmp_vmax/dyel)/maxcfl )
          nout = nstep
      ELSE
          nstep = noutput*CEILING( tfinal*MAX(tmp_umax/dxel,tmp_vmax/dyel)/maxcfl/DBLE(noutput) )
          nout = noutput
      ENDIF

      dt = tfinal/DBLE(nstep)
      calculatedMu = MAX(tmp_umax/dxel,tmp_vmax/dyel)*dt
      write(*,*) 'Mu used=',calculatedMu

      IF(p==1) THEN ! Set up netCDF output file
        cdfOut = TRIM(outdir)//TRIM(cdfOut)//'.nc'
        CALL output2d(q,xPlot,yPlot,tfinal,calculatedMu,cdfOut,nout,-1)
      ENDIF
      CALL output2d(q,xPlot,yPlot,tfinal,calculatedMu,cdfOut,p,0)

      ! =====================================================================================================
      ! Time integration
      ! =====================================================================================================
      DO m=1,meqn
        tmp_qmax(m) = MAXVAL(q(:,:,m))
        tmp_qmin(m) = MINVAL(q(:,:,m))
      ENDDO !m

      oddstep = .TRUE.
      DO n=1,nstep
        ! Call update
        CALL strangSplit(q,u0,v0,uEdge0,vEdge0,quadNodes,quadWeights,time,&
                               legendreVal,legendreDeriv,avgXferOp,avgXferOpLU,IPIV,&
                               dt,dxel,dyel,oddstep)

        time = time + dt
        ! Check if this is output time
        IF((MOD(n,nstep/nout).eq.0).OR.(n.eq.nstep)) THEN ! Write output variables
          CALL output2d(q,xPlot,yPlot,tfinal,calculatedMu,cdfOut,p,2)
        ENDIF
        DO m=1,meqn
          tmp_qmin(m) = MIN(tmp_qmin(m),MINVAL(q(:,:,m)))
          tmp_qmax(m) = MAX(tmp_qmax(m),MAXVAL(q(:,:,m)))
        ENDDO !m

        oddstep = .NOT. oddstep
      ENDDO !n

      CALL cpu_time(tf)
      tf = tf - t0
      write(*,*) 'tf=', tf

      IF(p == nRuns) THEN
        ! Close output files
        CALL output2d(q,xPlot,yPlot,tfinal,calculatedMu,cdfOut,nout,1)
      ENDIF

      DEALLOCATE(elemCenterX,elemCenterY,xPlot,yPlot,DGx,DGy,avgXferOp,avgXferOpLU,IPIV,FOO,&
                q,q0,u0,v0,uEdge0,vEdge0,STAT=ierr)
    ENDDO !p

    DEALLOCATE(quadNodes,quadWeights,legendreVal,legendreDeriv,STAT=ierr)

    ENDDO !nmethod

CONTAINS
  SUBROUTINE makeGrid(dxel,dyel,xCenter,yCenter,dxPlot,dyPlot,xPlot,yPlot)
    ! =============================================================================
    ! Computes cell width and initializes cell centers and quadrature grid
    ! INPUTS:   nex,ney - number of elements
    !           nQuad - number of quadrature nodes
    !           quadNodes(0:nQuad) - Gauss-Legendre quadrature nodes
    ! OUTPUTS:  dxel,dyel - width of elements
    !           xCenter(j),yCenter(j) - location of jth element center
    ! =============================================================================
    USE commonTestParameters, ONLY: xDomain,yDomain,nxOut,nyOut
    IMPLICIT NONE
    ! Inputs
    ! Outputs
    DOUBLE PRECISION, INTENT(OUT) :: dxel,dyel,dxPlot,dyPlot
    DOUBLE PRECISION, DIMENSION(1:nex), INTENT(OUT) :: xCenter
    DOUBLE PRECISION, DIMENSION(1:ney), INTENT(OUT) :: yCenter
    DOUBLE PRECISION, DIMENSION(1:nxOut), INTENT(OUT) :: xPlot
    DOUBLE PRECISION, DIMENSION(1:nyOut), INTENT(OUT) :: yPlot

    ! Local variables
    DOUBLE PRECISION :: domainWidth
    INTEGER :: k,j

    domainWidth = xDomain(2)-xDomain(1)
    dxel = domainWidth/DBLE(nex)
    dxPlot = domainWidth/DBLE(nxOut)

    domainWidth = yDomain(2)-yDomain(1)
    dyel = domainWidth/DBLE(ney)
    dyPlot = domainWidth/DBLE(nyOut)

    xCenter(1) = xDomain(1)+0.5D0*dxel
    DO j=2,nex
      xCenter(j) = xCenter(j-1)+dxel
    ENDDO!j

    xPlot(1) = xDomain(1)+0.5D0*dxPlot
    DO j=2,nxOut
      xPlot(j) = xPlot(j-1)+dxPlot
    ENDDO!j

    yCenter(1) = yDomain(1)+0.5D0*dyel
    DO j=2,ney
      yCenter(j) = yCenter(j-1)+dyel
    ENDDO!j

    yPlot(1) = xDomain(1)+0.5D0*dyPlot
    DO j=2,nyOut
      yPlot(j) = yPlot(j-1)+dyPlot
    ENDDO!j


  END SUBROUTINE makeGrid
END SUBROUTINE DRIVER
