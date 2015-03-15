! ##############################################################
! Module containing Modal Discontinuous Galerkin basis functions, quadrature data
! and transformation matrix for use in 1D simulations.
! By : Devin Light 04.18.2013
! ##############################################################

MODULE mDGmod
	IMPLICIT NONE
	INTEGER, PARAMETER, PRIVATE :: DOUBLE=KIND(1D0)

	CONTAINS

		! ########################################################################
    ! N-choose-k Function
    ! ########################################################################
       REAL(KIND=DOUBLE) FUNCTION choose(alpha,k)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: k
            REAL(KIND=DOUBLE), INTENT(IN) :: alpha
            INTEGER :: i
            REAL(KIND=DOUBLE) :: HOLDER

            HOLDER = 1D0

            DO i = 1,k
                HOLDER = HOLDER*((alpha-DBLE(k-i))/(DBLE(i)))
            END DO
            choose = HOLDER
        END FUNCTION choose

    ! ########################################################################
    ! Legendre Polynomial function of degree N
    ! ########################################################################
        REAL(KIND=DOUBLE) FUNCTION legendre(x,N)
            IMPLICIT NONE
            REAL(KIND=DOUBLE), INTENT(IN) :: x
            REAL(KIND=DOUBLE) :: HOLDER
            INTEGER, INTENT(IN) :: N
            INTEGER :: k

            HOLDER = 0.D0
            DO k = 0,N
                HOLDER = HOLDER + choose(DBLE(N),k)*choose((N+k-1)/2D0,N)*x**k
            END DO

            legendre = HOLDER*(2**N)

        END FUNCTION legendre

    ! ########################################################################
    ! Derivative of Legendre Polyomial of degree N
    ! ########################################################################
        REAL(KIND=DOUBLE) FUNCTION dlegendre(x,N)
            IMPLICIT NONE
            REAL(KIND=DOUBLE),INTENT(IN) :: x
            REAL(KIND=DOUBLE) :: HOLDER
            INTEGER, INTENT(IN) :: N ! Order of legendre polynomial
            INTEGER :: k

            HOLDER = 0.D0
            DO k = 1,N
                HOLDER = HOLDER + k*choose(DBLE(N),k)*choose((N+k-1)/2D0,N)*x**(k-1)
            END DO

            dlegendre = HOLDER*(2**N)

        END FUNCTION dlegendre

    ! ###########################################################################################################
    ! Subroutine for computing Gaussian quadrature nodes based on the derivative of Mth Order Legendre Polynomial
		! For Modal DG, we require M=N+1 nodes, where N is the highest order of Legendre polynomial being used
    ! ###########################################################################################################
        SUBROUTINE quad_nodes(M,nodes)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M
            REAL(KIND=DOUBLE), DIMENSION(0:M-1), INTENT(OUT) :: nodes
            REAL(KIND=DOUBLE) :: xnew,xold,error,tol, PI
            INTEGER :: k

            PI = DACOS(-1D0)

            tol  = 10D-10

            DO k = 0,M-1
                error = 1D0
                xold = -1D0*DCOS(((2*k+1)/(2D0*M))*PI)

                DO WHILE (error>tol)
                    xnew = xold - (legendre(xold,M))/(1D0*dlegendre(xold,M))
                    error = DABS(xnew-xold)
                    xold = xnew
                END DO
                nodes(k) = xold
            END DO
        END SUBROUTINE quad_nodes

		! #######################################################################################################
    ! Computing weights associated with N+1 nodes for quadratures on [-1,1]
		! For Modal DG, we require M=N+1 weights, where N is the highest order of Legendre polynomial being used
    ! #######################################################################################################
        SUBROUTINE quad_weights(M,nodes,wghts)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M
            REAL(KIND=DOUBLE), DIMENSION(0:M-1), INTENT(IN) :: nodes
            REAL(KIND=DOUBLE), DIMENSION(0:M-1), INTENT(OUT) :: wghts
            INTEGER :: k

            DO k = 0,M-1
                wghts(k) = 2D0*(1-nodes(k)**2)/((M*legendre(nodes(k),M-1))**2)
				!wghts(k) = 2D0/( (1-nodes(k)**2)*(dlegendre(nodes(k),M))**2 )
            END DO

        END SUBROUTINE quad_weights

				! #######################################################################
					! Subroutine for filling in the C-matrix. Used for interchanging solution
					! between DG and PPM, computed using Gaussian quadrature
				! #######################################################################
					SUBROUTINE Cmat_FILL(N,nodes,wghts,dx,dxelem,output,stat)
							IMPLICIT NONE
							INTEGER, INTENT(IN) :: N
							REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: nodes, wghts
							REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(OUT) :: output
							CHARACTER(len=*), INTENT(IN) :: stat
							REAL(KIND=DOUBLE), INTENT(IN) :: dx,dxelem
							INTEGER :: i,k,m
							REAL(KIND=DOUBLE) :: dz,zi
							REAL(KIND=DOUBLE), DIMENSION(0:N) :: foo

							dz = dx/dxelem

							foo = 0D0

							IF(TRIM(stat) .eq. 'gll') THEN
								DO k=0,N
									DO m=0,N
										output(k,m) = legendre(nodes(k),m)
									ENDDO
								ENDDO

							ELSE
								DO k=0,N
									DO m=0,N
										foo = 0D0
									DO i=0,N
										zi = -1D0 + dz*(2*k+nodes(i)+1D0)
										foo(i) = wghts(i)*legendre(zi,m)
									ENDDO
										output(k,m) = 0.5D0*SUM(foo)
									ENDDO
								ENDDO
							ENDIF

					END SUBROUTINE Cmat_FILL
					
	! ###############################################################################################################
	! phitld(xi,j,A,N,nelem) computes the complete series expansion form of solution based on given coefficents and element
	! ###############################################################################################################

		REAL(KIND=DOUBLE) FUNCTION phitld(xi,j,Ain,N,nelem)
			IMPLICIT NONE
			REAL(KIND=DOUBLE), INTENT(IN) :: xi
			INTEGER, INTENT(IN) :: j, N, nelem
			REAL(KIND=DOUBLE), DIMENSION(0:N,0:nelem+1), INTENT(IN) :: Ain
			INTEGER :: i ! Looping variable
			REAL(KIND=DOUBLE), DIMENSION(0:N) :: foo

			DO i = 0,N
				foo(i) = Ain(i,j)*legendre(xi,i)
			END DO

			phitld = SUM(foo)

		END FUNCTION phitld



END MODULE mDGmod
