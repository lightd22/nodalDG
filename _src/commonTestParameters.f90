MODULE commonTestParameters
  ! ========================================================================
  ! This module will contain information about the current physical domain and test parameters.
  ! Used to simplify passing of this information throughout modal subroutines and functions
  ! ========================================================================
  IMPLICIT NONE
  INTEGER :: nex,ney,nQuad,nxOut,nyOut,meqn,testID,maxPolyDegree,limitingType
  INTEGER :: cdfID
  INTEGER :: inUnit
  DOUBLE PRECISION, DIMENSION(1:2) :: xDomain,yDomain
  DOUBLE PRECISION :: PI,tfinal
  LOGICAL :: transient,doposlimit,doreactive
  CHARACTER(len=24) :: outdir
  SAVE

END MODULE commonTestParameters
