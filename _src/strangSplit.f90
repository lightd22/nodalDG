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
    DOUBLE PRECISION, INTENT(IN) :: dt,dxel,dyel,time
    DOUBLE PRECISION, DIMENSION(1:nxOut,1:nyOut), INTENT(IN) :: u0,v0
    DOUBLE PRECISION, DIMENSION(1:nex,1:nyOut), INTENT(IN) :: uEdge0
    DOUBLE PRECISION, DIMENSION(1:nxOut,1:ney), INTENT(IN) :: vEdge0
    DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadNodes,quadWeights
    DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: legendreVal,legendreDeriv,avgXferOp,avgXferOpLU
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
    ENDIF !transient
    IF(oddstep) THEN
        ! ===================================
        ! Perform sweeps in x-direction first
        ! ===================================
        DO j=1,nyOut
            q1dx = q(:,j,:)
            u1dx(1:3,:) = u(1:3,:,j)
            uEdge1dx(1:3,:) = uEdge(1:3,:,j)
!            CALL mDGsweep(q1dx,u1dx,uEdge1dx,dxel,nex,nOrder,quadWeights,avgXferOp,avgXferOpLU, &
!                          legendreVal,legendreDeriv,IPIV,dt,doposlimit,posWeight,maxTime,minTime,&
!                          totTime)
            ! Update solution
            q(:,j,:) = q1dx
        ENDDO!j

        DO i=1,nxOut
            q1dy = q(i,:,:)
            v1dy(1:3,:) = v(1:3,i,:)
            vEdge1dy(1:3,:) = vEdge(1:3,i,:)

!            CALL mDGsweep(q1dy,v1dy,vEdge1dy,dyel,ney,nOrder,quadWeights,avgXferOp,avgXferOPLU,&
!                          legendreVal,legendreDeriv,IPIV,dt,doposlimit,posWeight,maxTime,minTime,&
!                          totTime)
            ! Update solution
            q(i,:,:) = q1dy
        ENDDO !i

    ELSE
        ! ===================================
        ! Perform sweeps in y-direction first
        ! ===================================
        DO i=1,nxOut
            q1dy = q(i,:,:)
            v1dy(1:3,:) = v(1:3,i,:)
            vEdge1dy(1:3,:) = vEdge(1:3,i,:)

!            CALL mDGsweep(q1dy,v1dy,vEdge1dy,dyel,ney,nOrder,quadWeights,avgXferOp,avgXferOPLU,&
!                          legendreVal,legendreDeriv,IPIV,dt,doposlimit,posWeight,maxTime,minTime,&
!                          totTime)
            ! Update solution
            q(i,:,:) = q1dy
        ENDDO !i

        DO j=1,nyOut
            q1dx = q(:,j,:)
            u1dx(1:3,:) = u(1:3,:,j)
            uEdge1dx(1:3,:) = uEdge(1:3,:,j)

!            CALL mDGsweep(q1dx,u1dx,uEdge1dx,dxel,nex,nOrder,quadWeights,avgXferOp,avgXferOpLU, &
!                          legendreVal,legendreDeriv,IPIV,dt,doposlimit,posWeight,maxTime,minTime,&
!                          totTime)
            ! Update solution
            q(:,j,:) = q1dx
        ENDDO!j
    ENDIF !oddstep
END SUBROUTINE strangSplit
