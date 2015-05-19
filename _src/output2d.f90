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
  CHARACTER(len=60), INTENT(IN) :: cdfOut
  DOUBLE PRECISION, INTENT(IN) :: muOut,timeOut
  DOUBLE PRECISION, DIMENSION(1:nxOut), INTENT(IN) :: xOut
  DOUBLE PRECISION, DIMENSION(1:nyOut), INTENT(IN) :: yOut
  DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn), INTENT(IN) :: q
  ! Outputs
  ! Local variables
  INTEGER, PARAMETER :: NDIMS = 4
  INTEGER :: idt, meqn_dimid, x_dimid, y_dimid, t_dimid,idx, idy, idmu, idmpd,idq
  CHARACTER(len=8) :: nxname,xname,nyname,yname,qname,muname,meqnName
  INTEGER, DIMENSION(1:NDIMS) :: start, count,dimids(NDIMS)
  INTEGER :: i,ierr,m,j
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp
  DOUBLE PRECISION, DIMENSION(1:nxOut) :: qSlice

  SAVE start,count,idq,t_dimid,meqn_dimid

  IF(stat == -1) THEN
    ! Create netCDF file and time variables
    ierr = NF90_CREATE(TRIM(cdfOut),NF90_CLOBBER,cdfID)

    ierr = NF90_REDEF(cdfID)
    ierr = NF90_DEF_DIM(cdfID, "nt", ilvl+1, t_dimid)
    ierr = NF90_DEF_DIM(cdfID, "meqn", meqn, meqn_dimid)
    ierr = NF90_DEF_VAR(cdfID, "time", NF90_FLOAT, t_dimid,idt)
    ierr = NF90_DEF_VAR(cdfID, "maxPoly",NF90_INT,idmpd)

    ierr = NF90_ENDDEF(cdfID)

    ! Calculate time at output levels (note ilvl=noutput)
    ALLOCATE(tmp(1:ilvl+1), STAT=ierr)
    DO i=0,ilvl
      tmp(i+1) = DBLE(i)*timeOut/DBLE(ilvl)
    ENDDO

    ! Write t values
    ierr = NF90_PUT_VAR(cdfID,idt,tmp)
    ierr = NF90_PUT_VAR(cdfID,idmpd,maxPolyDegree)

    DEALLOCATE(tmp, STAT=ierr)
    RETURN
  ELSEIF(stat == 0) THEN
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

    ierr = NF90_REDEF(cdfID)

    ierr = NF90_DEF_DIM(cdfID, TRIM(nxname), nxOut, x_dimid)
    ierr = NF90_DEF_DIM(cdfID, TRIM(nyname), nyOut, y_dimid)

    dimids(1) = x_dimid
    dimids(2) = y_dimid
    dimids(3) = t_dimid
    dimids(4) = meqn_dimid

    ierr = NF90_DEF_VAR(cdfid, TRIM(qname),NF90_FLOAT,dimids,idq)
    ierr = NF90_DEF_VAR(cdfid, TRIM(xname),NF90_FLOAT,x_dimid,idx)
    ierr = NF90_DEF_VAR(cdfid, TRIM(yname),NF90_FLOAT,y_dimid,idy)
    ierr = NF90_DEF_VAR(cdfid, TRIM(muname),NF90_FLOAT,idmu)

    ierr = NF90_ENDDEF(cdfid)

    ! Write x and y values
    ierr = NF90_PUT_VAR(cdfid, idx, xOut)
    ierr = NF90_PUT_VAR(cdfid, idy, yOut)
    ierr = NF90_PUT_VAR(cdfid,idmu,muOut)

    start(3) = 1
  ELSEIF(stat == 1) THEN
    ! Close netCDF output file
    ierr = NF90_CLOSE(cdfid)
    RETURN
  ENDIF ! stat
  ! Write out concentration field
  count(1) = nxOut
  DO m=1,meqn
    start(4) = m
    DO j=1,nyOut
        start(2) = j
        qSlice = q(:,j,m)
        ierr = NF90_PUT_VAR(cdfid,idq,qSlice,start,count)
    ENDDO !ylvl
  ENDDO !meqn

  ! Increment t level
  start(3) = start(3) + 1


END SUBROUTINE output2d
