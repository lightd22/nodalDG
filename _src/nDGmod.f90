! ##############################################################
! Module containing 4th and 5th order Legendre polynomials
! and their derivatives used in ppmwrap.f90 as part of test_advection_slskam_2d.f90
! By : Devin Light 11.29.2012
! ##############################################################

MODULE nDGmod
	IMPLICIT NONE
    INTEGER, PARAMETER :: DOUBLE = KIND(1D0)


    PRIVATE :: DOUBLE

    ! #######################################################################
    ! node4 is vector of zeros of Jacobi polynomial
    ! associated with the problem, and are also the location of the GLL nodes
    ! #######################################################################

    REAL(KIND = DOUBLE), DIMENSION(0:4) :: node4 = (/ &
            -1D0, &
            -0.654653670707978D0, &
            0D0, &
            0.654653670707977D0, &
            1D0 /)


    ! ##########################################################
    ! w4 is the weights used in the GLL quadrature
    ! with w4(i) being the weight of the ith term in the sum
    ! ##########################################################

    REAL(KIND = DOUBLE), DIMENSION(0:4) :: w4 = (/ &
            1D0/10D0, &
            0.5444444444444456D0, &
            32D0/45D0,&
            0.5444444444444456D0,&
            1D0/10D0 /)

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
            INTEGER, INTENT(IN) :: N
            INTEGER :: k

            HOLDER = 0.D0
            DO k = 1,N
                HOLDER = HOLDER + k*choose(DBLE(N),k)*choose((N+k-1)/2D0,N)*x**(k-1)
            END DO

            dlegendre = HOLDER*(2**N)

        END FUNCTION dlegendre

    ! ########################################################################
    ! 2nd Derivative of Legendre Polyomial of degree N
    ! ########################################################################
        REAL(KIND=DOUBLE) FUNCTION ddlegendre(x,N)
            IMPLICIT NONE
            REAL(KIND=DOUBLE), INTENT(IN) :: x
            REAL(KIND=DOUBLE) :: HOLDER
            INTEGER, INTENT(IN) :: N
            INTEGER :: k

            HOLDER = 0.D0
            DO k = 2,N
                HOLDER = HOLDER + k*(k-1)*choose(DBLE(N),k)*choose((N+k-1)/2D0,N)*x**(k-2)
            END DO

            ddlegendre = HOLDER*(2**N)

        END FUNCTION ddlegendre

    ! ########################################################################
    ! Subroutine for computing GLL nodes based on the derivative of N'th Order Legendre Polynomial
    ! ########################################################################
        SUBROUTINE gllquad_nodes(N,nodes)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: N
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(OUT) :: nodes
            REAL(KIND=DOUBLE) :: xnew,xold,error,tol, PI
            INTEGER :: k

            PI = DACOS(-1D0)

            tol  = 10.D0**(-8)

            nodes(0) = -1D0
            nodes(N) = 1D0

            DO k = 1,N-1
                error = 1D0
                xold = -1D0*DCOS( ((2*k-1)/(2D0*(N-1)))*PI)

                DO WHILE (error>tol)
                    xnew = xold - (dlegendre(xold,N))/(1D0*ddlegendre(xold,N))
                    error = DABS(xnew-xold)
                    xold = xnew
                END DO
                nodes(k) = xold
            END DO
        END SUBROUTINE gllquad_nodes

    ! ###########################################################################################################
    ! Subroutine for computing Gaussian quadrature nodes based on the derivative of Mth Order Legendre Polynomial
	! For Modal DG, we require M=N+1 nodes, where N is the highest order of Legendre polynomial being used
    ! ###########################################################################################################
        SUBROUTINE gaussquad_nodes(M,nodes)
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
        END SUBROUTINE gaussquad_nodes

    ! ########################################################################
    ! Computing weights associated with N+1 nodes for quadratures on [-1,1]
    ! ########################################################################
        SUBROUTINE gllquad_weights(N,nodes,wghts)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: N
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: nodes
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(OUT) :: wghts
            INTEGER :: k

            DO k = 0,N
                wghts(k) = 2D0/(N*(N+1)*(legendre(nodes(k),N))**2)
            END DO

        END SUBROUTINE gllquad_weights

	! #######################################################################################################
    ! Computing weights associated with N+1 nodes for quadratures on [-1,1]
	! For Modal DG, we require M=N+1 weights, where N is the highest order of Legendre polynomial being used
    ! #######################################################################################################
        SUBROUTINE gaussquad_weights(M,nodes,wghts)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M
            REAL(KIND=DOUBLE), DIMENSION(0:M-1), INTENT(IN) :: nodes
            REAL(KIND=DOUBLE), DIMENSION(0:M-1), INTENT(OUT) :: wghts
            INTEGER :: k

            DO k = 0,M-1
                wghts(k) = 2D0*(1-nodes(k)**2)/((M*legendre(nodes(k),M-1))**2)
				!wghts(k) = 2D0/( (1-nodes(k)**2)*(dlegendre(nodes(k),M))**2 )
            END DO

        END SUBROUTINE gaussquad_weights


    ! ###########################################################
    ! baryWeights computes the set of barycentric weights for the Lagrange interpolating polynomial,
    ! used to evaluate the basis functions
    ! ###########################################################

        SUBROUTINE fillBaryWeights(lambda,nodes,N)
            IMPLICIT NONE
            ! Inputs
            INTEGER, INTENT(IN) :: N
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: nodes
            ! Outputs
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(OUT) :: lambda ! lambda(k) is the kth barycentric weight
            ! Local variables
            INTEGER :: k
            LOGICAL, DIMENSION(0:N) :: MASK

            DO k=0,N
                MASK = .TRUE.
                MASK(k) = .FALSE.

                lambda(k) = 1D0/PRODUCT(nodes(k)-nodes,MASK)
            ENDDO !k

        END SUBROUTINE fillBaryWeights

    ! ###########################################################
    ! phi computes the k'th basis function, a Lagrange interpolating
    ! polynomial, for a given set of nodes
    ! ###########################################################
        REAL(KIND=DOUBLE) FUNCTION lagrange(xi,k,N,nodes,lambda)
            IMPLICIT NONE
            REAL(KIND=DOUBLE), INTENT(IN) :: xi
            INTEGER, INTENT(IN) :: k,N
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: nodes,lambda
            LOGICAL, DIMENSION(0:N) :: MASK
            REAL(KIND=DOUBLE), DIMENSION(0:N) :: l
            INTEGER :: i = 0

            MASK = .TRUE.
            MASK(k) = .FALSE.

            l = xi - nodes
            lagrange = PRODUCT(l,MASK)*lambda(k)

       END FUNCTION lagrange


    ! #############################################################################
    ! Subroutine D(N,nodes,output) computes the matrix of values of diff(phi(k,x),x)
    ! evaluated at x = nodes(n); Used in calc of the Galerkin step
    ! #############################################################################
        SUBROUTINE Dmat(N,nodes,output)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: N
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: nodes
            REAL(KIND=DOUBLE), DIMENSION(0:N, 0:N), INTENT(OUT) :: output
            INTEGER :: i = 0, j = 0

            DO i = 0, N
                DO j = 0, N
                    IF (i /= j) THEN
                      output(i,j) = (legendre(nodes(j),N))/(legendre(nodes(i),N)*(nodes(j)-nodes(i)))
                    ELSE
                      output(i,j) = 0.D0
                    END IF
                END DO
            END DO

            output(0,0) = -(N)*(N+1)/4D0
            output(N,N) = (N)*(N+1)/4D0

        END SUBROUTINE Dmat

        SUBROUTINE GEPP_INV (M,N,Minv)
        !
        ! Subroutine to perform the partial-pivoting Gaussian elimination.
        ! A(N,N) is the original matrix in the input and transformed matrix
        ! plus the pivoting element ratios below the diagonal in the output.
        ! INDX(N) records the pivoting order.
        ! Ainv is found by performing the same operations on I(N,N) and then
        ! solving the related systems for each column of Ainv. [A|I] -> [A|Y]
        !
          IMPLICIT NONE
          INTEGER, INTENT (IN) :: N
          INTEGER :: I,J,K,ITMP
          INTEGER, DIMENSION (N) :: INDX
          REAL :: C1,PI,PI1,PJ
          REAL(KIND=DOUBLE), INTENT (IN), DIMENSION (N,N) :: M
          REAL(KIND=DOUBLE), INTENT(OUT), DIMENSION(N,N) :: Minv
          REAL(KIND=DOUBLE), DIMENSION(N,N) :: Y,Tmp1,Tmp2,A
          REAL, DIMENSION (N) :: C

        ! Initialize A and Minv
        A(:,:) = M(:,:)
        Minv(:,:) = 0.D0

        ! Initialize Y as I(N,N)
          Y(:,:) = 0D0
          DO I = 1,N
            Y(I,I) = 1D0
          END DO

        !
        ! Initialize the index
        !
          DO I = 1, N
            INDX(I) = I
          END DO
        !
        ! Select largest absval element, one from each row
        !
!          DO I = 1, N
!            C1= 0.0
!            DO J = 1, N
!              C1 = DMAX1(C1,ABS(A(I,J)))
!            END DO
!            C(I) = C1
!          END DO


          DO J = 1, N-1

            ! Select pivoting (largest) element from each column
            PI1 = 0.0
            DO I = J, N
              PI = DABS(A(INDX(I),J)) !/C(INDX(I))
              IF (PI.GT.PI1) THEN
                PI1 = PI
                K   = I
              END IF
            END DO
        !
        ! Interchange the rows via INDX(N) to record pivoting order
        !
            ITMP    = INDX(J)
            INDX(J) = INDX(K)
            INDX(K) = ITMP
            DO I = J+1, N
              PJ  = A(INDX(I),J)/A(INDX(J),J)
        !
        ! Record pivoting ratios below the diagonal
        !
              A(INDX(I),J) = 0D0!PJ
        !
        ! Modify other elements accordingly
        !
              DO K = J+1, N
                A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
              END DO
              DO K = 1,N
                Y(INDX(I),K) = Y(INDX(I),K)-PJ*Y(INDX(J),K)
              END DO
            END DO
          END DO

        ! Swap rows to get it back to the correct form
        Tmp1 = A
        Tmp2 = Y
        DO I=1,N
            A(I,:) = Tmp1(INDX(I),:)
            Y(I,:) = Tmp2(INDX(I),:)
        END DO

        ! To find Minv (for the PIVOTED matrix), solve n-systems using back substitution
        ! Anew*Ainv(:,k) = Y(:,k) k=1..n
          DO K = 1,N
            Minv(N,K) = Y(N,K)/A(N,N)
            DO I = N-1,1,-1
                Minv(I,K) = Y(I,K)
                DO J = I+1,N
                    Minv(I,K) = Minv(I,K) - A(I,J)*Minv(J,K)
                END DO
            Minv(I,K) = Minv(I,K)/A(I,I)
            END DO
          END DO

        END SUBROUTINE GEPP_INV

END MODULE nDGmod
