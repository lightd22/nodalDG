&inputs
    ! Spatial element parameters
    startRes = 12,      ! Which resolution is run first
    nRuns = 3,          ! How many runs are done
    nScale = 2,         ! Ratio between number of elements in successive runs
    maxPolyDegree = 4,  ! Degree of reconstructing polynomial

    ! Time stepping paramteters
    cflCoeff = 0.9D0    ! Ratio of used CFL number to maximum stable CFL

    ! Outputting parameters
    noutput = 30         ! Number of times to output output, including final time (must be >= 1) (automatically includes ICs)

    ! Testing parmeters
    meqn = 2            ! Number of tracers being simulated

    testID = 5          ! 0 = Consistency test
                        ! 1 = Uniform diagonal advection
                        ! 2 = Reactive def. ICs
                        ! 5 = LeVeque deformation of C^3 cosinebell
                        ! 6 = LeVeque deformation of C^5 cosinebell
                        ! 7 = LeVeque deformation of slotted cylinder

    tfinal = 5D0        ! Final time of integration

    TRANSIENT = .TRUE.  ! Time-dependent flow
    DOREACTIVE = .TRUE. ! Reactive flow

    ! Misc parameters
    DEBUG = .FALSE.

/
