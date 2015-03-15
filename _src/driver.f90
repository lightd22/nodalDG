SUBROUTINE DRIVER(testID,nex0,ney0,nscale,nruns,noutput,maxCFL)
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
  INTEGER, INTENT(IN) :: testID,nex0,ney0,nscale,nruns,noutput
  REAL(KIND=8), INTENT(IN) :: maxCFL
  ! Outputs
  ! Local variables
  DOUBLE PRECISION :: xLeft,xRight,yLeft,yRight
  INTEGER, DIMENSION(10) :: tmp_method
  INTEGER :: nmethod,nmethod_final,imethod,ierr,i,j
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: q,q0,u0,v0,uEdge0,vEdge0,&
                                  avgXferOp,avgXferOpLU,legendreVal,legendreDeriv

  if(nRuns.lt.1) STOP 'nRuns should be at least 1 in DRIVER()'
  PI = DACOS(-1D0)

  ! Set spatial domain
  xLeft = 0D0
  xRight = 1D0
  yLeft = xLeft
  yRight = xRight

  nmethod_final = 1
  tmp_method = 0
  tmp_method(1) = 1 ! Split modal DG, no limiting
  tmp_method(1) = 2 ! Split modal DG, mass redistribution limiting for positivity

  DO nmethod = 1,nmethod_final
    imethod = tmp_method(nmethod)

    SELECT CASE(imethod)
      CASE(1)
        write(*,*) 'DG, averages, no limiting'
        doposlimit = .false.
        dogllGrid = .false.
        outdir = '_modal/'
      CASE(2)
        write(*,*) 'DG, averages, element mass redist'
        doposlimit = .true.
        dogllGrid = .false.
        outdir = '_pdModal/'
    END SELECT !imethod

    write(*,FMT='(A5,i1)') ' N = ',nOrder
    write(*,*) 'WARNING: Only periodic BCs are implemented'

    nQuad = maxPolyDegree

    ALLOCATE(quadNodes(0:nQuad),quadWeights(0:nQuad),legendreVal(0:maxPolyDegree,0:nQuad),&
            legendreDeriv(0:maxPolyDegree,0:nQuad),STAT=ierr)

    CALL quad_nodes(nOrder+1,quadNodes)
    CALL quad_weights(nOrder+1,quadNodes,quadWeights)

    ! Fill array of Legendre polynomials evaluated at quad nodes + Leg. derivative at quad nodes
    DO i=0,nOrder
      DO j=0,nOrder
        legendreVal(i,j) = legendre(quadNodes(j),i)
        legendreDeriv(i,j) = dlegendre(quadNodes(j),i)
      ENDDO !j
    ENDDO !i

    DEALLOCATE(quadNodes,quadWeights,legendreVal,legendreDeriv,STAT=ierr)

    ENDDO !nmethod

END SUBROUTINE DRIVER
