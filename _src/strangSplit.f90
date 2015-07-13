SUBROUTINE strangSplit(q,u,v,uEdge,vEdge,quadNodes,quadWeights,time,&
                       basisPolyVal,basisPolyDeriv,avgOP,avgOP_LU,IPIV,&
                       dt,dxel,dyel,reactiveCoeffs,oddstep)
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
    DOUBLE PRECISION, INTENT(IN) :: dt,dxel,dyel,time
    DOUBLE PRECISION, DIMENSION(1:3,1:nxOut,1:nyOut), INTENT(IN) :: u,v
    DOUBLE PRECISION, DIMENSION(1:3,1:nex,1:nyOut), INTENT(IN) :: uEdge
    DOUBLE PRECISION, DIMENSION(1:3,1:nxOut,1:ney), INTENT(IN) :: vEdge
    DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadNodes,quadWeights
    DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: basisPolyVal,basisPolyDeriv
    DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:maxPolyDegree),INTENT(IN) :: avgOP,avgOp_LU
    DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn), INTENT(IN) :: reactiveCoeffs
    INTEGER, DIMENSION(0:maxPolyDegree), INTENT(IN) :: IPIV
    LOGICAL, INTENT(IN) :: oddstep
    ! Outputs
    DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn), INTENT(INOUT) :: q
    ! Local variables
    INTEGER :: i,j,k
    DOUBLE PRECISION, DIMENSION(1:nxOut,1:meqn) :: q1dx
    DOUBLE PRECISION, DIMENSION(1:nyOut,1:meqn) :: q1dy
    DOUBLE PRECISION, DIMENSION(1:3,1:nxOut) :: u1dx
    DOUBLE PRECISION, DIMENSION(1:3,1:nyOut) :: v1dy
    DOUBLE PRECISION, DIMENSION(1:3,1:nex) :: uEdge1dx
    DOUBLE PRECISION, DIMENSION(1:3,1:ney) :: vEdge1dy

    INTERFACE
      SUBROUTINE updateSoln1d(q,u,uEdge,dt,dxel,nelem,nx,quadWeights,&
                              basisVals,basisDeriv)
        ! ===========================================================================
        ! Takes full dt time step for one dimensional slice of subcell averages using SSPRK3
        ! integrator
        ! ===========================================================================
        USE commonTestParameters
        IMPLICIT NONE
        ! Inputs
        INTEGER, INTENT(IN) :: nelem,nx
        DOUBLE PRECISION, INTENT(IN) :: dxel,dt
        DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadWeights
        DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: basisVals,&
          basisDeriv
        DOUBLE PRECISION, DIMENSION(1:3,1:nx), INTENT(IN) :: u
        DOUBLE PRECISION, DIMENSION(1:3,1:nelem), INTENT(IN) :: uEdge
        ! Outputs
        DOUBLE PRECISION, DIMENSION(1:nx,1:meqn), INTENT(INOUT) :: q
      END SUBROUTINE updateSoln1d

      SUBROUTINE reactiveStep(q,dt,reactiveCoeffs)
        ! ================================================================================
        ! Takes a single time step to update subcell average values to tn+1
        ! for "toy chemistry problem"
        ! dq1/dt = -k1 q1 + k2 q2^2
        ! dq2/dt = 2 k1 q1 - 2 k2 q2^2
        ! Currently uses 2-stage, 2nd order Rosenbock Runge-Kutta method with the following parameters:
        ! (see Durran "Numerical Methods for Fluid Dynamics")
        ! b1 = b2 = 0.5 ; alpha = 1+1/(2 sqrt(2)) ; a21 = 1 ; alpha21 = -2 alpha
        !
        ! INPUTS:
        ! OUTPUTS: q(i,j,m) - mth field subcell averages updated to new time
        ! ================================================================================
        USE commonTestParameters
        IMPLICIT NONE
        ! Inputs
        DOUBLE PRECISION, INTENT(IN) :: dt
        DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn), INTENT(IN) :: reactiveCoeffs
        ! Outputs
        DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn), INTENT(INOUT) :: q
      END SUBROUTINE reactiveStep
    END INTERFACE

    IF(oddstep) THEN
        ! ===================================
        ! Perform sweeps in x-direction first
        ! ===================================
        DO j=1,nyOut
          q1dx = q(:,j,:)
          u1dx(1:3,:) = u(1:3,:,j)
          uEdge1dx(1:3,:) = uEdge(1:3,:,j)
          CALL updateSoln1d(q1dx,u1dx,uEdge1dx,dt,dxel,nex,nxOut,quadWeights,&
                            basisPolyVal,basisPolyDeriv)
          ! Update solution
          q(:,j,:) = q1dx
        ENDDO!j

        DO i=1,nxOut
          q1dy = q(i,:,:)
          v1dy(1:3,:) = v(1:3,i,:)
          vEdge1dy(1:3,:) = vEdge(1:3,i,:)
          CALL updateSoln1d(q1dy,v1dy,vEdge1dy,dt,dyel,ney,nyOut,quadWeights,&
                            basisPolyVal,basisPolyDeriv)
          ! Update solution
          q(i,:,:) = q1dy
        ENDDO !i

        IF(doreactive) CALL reactiveStep(q,dt,reactiveCoeffs)

    ELSE
        ! ===================================
        ! Perform sweeps in y-direction first
        ! ===================================
        IF(doreactive) CALL reactiveStep(q,dt,reactiveCoeffs)

        DO i=1,nxOut
          q1dy = q(i,:,:)
          v1dy(1:3,:) = v(1:3,i,:)
          vEdge1dy(1:3,:) = vEdge(1:3,i,:)
          CALL updateSoln1d(q1dy,v1dy,vEdge1dy,dt,dyel,ney,nyOut,quadWeights,&
                            basisPolyVal,basisPolyDeriv)
          ! Update solution
          q(i,:,:) = q1dy
        ENDDO !i

        DO j=1,nyOut
          q1dx = q(:,j,:)
          u1dx(1:3,:) = u(1:3,:,j)
          uEdge1dx(1:3,:) = uEdge(1:3,:,j)
          CALL updateSoln1d(q1dx,u1dx,uEdge1dx,dt,dxel,nex,nxOut,quadWeights,&
                            basisPolyVal,basisPolyDeriv)
          ! Update solution
          q(:,j,:) = q1dx
        ENDDO!j
    ENDIF !oddstep
END SUBROUTINE strangSplit
