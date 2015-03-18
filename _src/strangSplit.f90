SUBROUTINE strangSplit(q,u0,v0,uEdge0,vEdge0,quadNodes,quadWeights,time,&
                       legendreVal,legendreDeriv,avgOP,avgOP_LU,IPIV,&
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
    DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut), INTENT(IN) :: u0,v0
    DOUBLE PRECISION, DIMENSION(1:nex,1:nyOut), INTENT(IN) :: uEdge0
    DOUBLE PRECISION, DIMENSION(1:nxOut,1:ney), INTENT(IN) :: vEdge0
    DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadNodes,quadWeights
    DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: legendreVal,legendreDeriv
    DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:maxPolyDegree),INTENT(IN) :: avgOP,avgOp_LU
    DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn), INTENT(IN) :: reactiveCoeffs
    INTEGER, DIMENSION(0:maxPolyDegree), INTENT(IN) :: IPIV
    LOGICAL, INTENT(IN) :: oddstep
    ! Outputs
    DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut,1:meqn), INTENT(INOUT) :: q
    ! Local variables
    INTEGER :: i,j,k
    DOUBLE PRECISION :: t_temp
    DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut) :: utmp,vtmp
    DOUBLE PRECISION, DIMENSION(1:nex,1:nyOut) :: uEdgetmp
    DOUBLE PRECISION, DIMENSION(1:nxOut,1:ney) :: vEdgetmp
    DOUBLE PRECISION, DIMENSION(1:3,1:nxOut,1:nyOut) :: u,v
    DOUBLE PRECISION, DIMENSION(1:3,1:nex,1:nyOut) :: uEdge
    DOUBLE PRECISION, DIMENSION(1:3,1:nxOut,1:ney) :: vEdge
    DOUBLE PRECISION, DIMENSION(1:nxOut,1:meqn) :: q1dx
    DOUBLE PRECISION, DIMENSION(1:nyOut,1:meqn) :: q1dy
    DOUBLE PRECISION, DIMENSION(1:3,1:nxOut) :: u1dx
    DOUBLE PRECISION, DIMENSION(1:3,1:nyOut) :: v1dy
    DOUBLE PRECISION, DIMENSION(1:3,1:nex) :: uEdge1dx
    DOUBLE PRECISION, DIMENSION(1:3,1:ney) :: vEdge1dy

    INTERFACE
      SUBROUTINE updateVelocities(u,v,uEdge,vEdge,time)
        ! =========================================================
        ! Updates the edge and quadrature velocities to given time
        ! =========================================================
        USE commonTestParameters
        IMPLICIT NONE
        ! Inputs
        DOUBLE PRECISION, INTENT(IN) :: time
        ! Outputs
        DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut), INTENT(INOUT) :: u,v
        DOUBLE PRECISION, DIMENSION(1:nex,1:nyOut), INTENT(INOUT) :: uEdge
        DOUBLE PRECISION, DIMENSION(1:nxOut,1:ney), INTENT(INOUT) :: vEdge
      END SUBROUTINE updateVelocities

      SUBROUTINE updateSoln1d(q,u,uEdge,dt,dxel,nelem,nx,quadWeights,avgOP,avgOP_LU,&
                              legVals,dlegVals,IPIV)
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
        DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: legVals,&
          dlegVals
        DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:maxPolyDegree), INTENT(IN) :: avgOp,avgOp_LU
        INTEGER, DIMENSION(0:maxPolyDegree), INTENT(IN) :: IPIV
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

    ! Update velocities at times required for ssprk3 update
    IF(transient) THEN
      DO i=1,3
        utmp = u0
        vtmp = v0
        uEdgetmp = uEdge0
        vEdgetmp = vEdge0
        SELECT CASE(i)
          CASE(1)
            t_temp = time
          CASE(2)
            t_temp = time+dt
          CASE(3)
            t_temp = time+0.5D0*dt
        END SELECT

        CALL updateVelocities(utmp,vtmp,uEdgetmp,vEdgetmp,t_temp)

        u(i,:,:) = utmp
        v(i,:,:) = vtmp
        uEdge(i,:,:) = uEdgetmp
        vEdge(i,:,:) = vEdgetmp
      ENDDO !i
    ELSE
      DO i=1,3
        u(i,:,:) = u0
        v(i,:,:) = v0
        uEdge(i,:,:) = uEdge0
        vEdge(i,:,:) = vEdge0
      ENDDO !i
    ENDIF !transient
    IF(oddstep) THEN
        ! ===================================
        ! Perform sweeps in x-direction first
        ! ===================================
        DO j=1,nyOut
          q1dx = q(:,j,:)
          u1dx(1:3,:) = u(1:3,:,j)
          uEdge1dx(1:3,:) = uEdge(1:3,:,j)
          CALL updateSoln1d(q1dx,u1dx,uEdge1dx,dt,dxel,nex,nxOut,quadWeights,&
                            avgOP,avgOP_LU,legendreVal,legendreDeriv,IPIV)
          ! Update solution
          q(:,j,:) = q1dx
        ENDDO!j

        DO i=1,nxOut
          q1dy = q(i,:,:)
          v1dy(1:3,:) = v(1:3,i,:)
          vEdge1dy(1:3,:) = vEdge(1:3,i,:)
          CALL updateSoln1d(q1dy,v1dy,vEdge1dy,dt,dyel,ney,nyOut,quadWeights,&
                            avgOP,avgOP_LU,legendreVal,legendreDeriv,IPIV)
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
                            avgOP,avgOP_LU,legendreVal,legendreDeriv,IPIV)
          ! Update solution
          q(i,:,:) = q1dy
        ENDDO !i

        DO j=1,nyOut
          q1dx = q(:,j,:)
          u1dx(1:3,:) = u(1:3,:,j)
          uEdge1dx(1:3,:) = uEdge(1:3,:,j)
          CALL updateSoln1d(q1dx,u1dx,uEdge1dx,dt,dxel,nex,nxOut,quadWeights,&
                            avgOP,avgOP_LU,legendreVal,legendreDeriv,IPIV)
          ! Update solution
          q(:,j,:) = q1dx
        ENDDO!j
    ENDIF !oddstep
END SUBROUTINE strangSplit
