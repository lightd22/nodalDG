MODULE commonTestParameters
  ! ========================================================================
  ! This module will contain information about the current physical domain and test parameters.
  ! Used to simplify passing of this information throughout modal subroutines and functions
  ! ========================================================================
  IMPLICIT NONE
  INTEGER :: nex,nQuad,nxOut,meqn,testID,maxPolyDegree
  DOUBLE PRECISION, DIMENSION(1:2) :: xDomain
  DOUBLE PRECISION :: PI,tfinal,u0,K0,rho0
  LOGICAL :: transient,doposlimit
  CHARACTER(len=24) :: outdir
  SAVE

END MODULE commonTestParameters
