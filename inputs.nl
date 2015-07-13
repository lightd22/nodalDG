&inputs
    ! Spatial element parameters
    startRes = 20,      ! Which resolution is run first
    nRuns = 2,          ! How many runs are done
    nScale = 2,         ! Ratio between number of elements in successive runs
    maxPolyDegree = 5,  ! Degree of reconstructing polynomial

    ! Time stepping paramteters
    cflCoeff = 0.45D0    ! Ratio of used CFL number to maximum stable CFL

    ! Outputting parameters
    noutput = 30         ! Number of times to output output, including final time (must be >= 1) (automatically includes ICs)

    ! Testing parmeters
    meqn = 1            ! Number of tracers being simulated

    testID = 5          ! 0 = Consistency test
                        ! 1 = Uniform diagonal advection
                        ! 2 = Reactive def. ICs
                        ! 5 = LeVeque deformation of C^3 cosinebell
                        ! 6 = LeVeque deformation of C^5 cosinebell
                        ! 7 = LeVeque deformation of slotted cylinder
                        ! 8 = Solid body rotation of cylinder

    tfinal = 5D0        ! Final time of integration

    TRANSIENT = .TRUE.  ! Time-dependent flow
    DOREACTIVE = .FALSE. ! Reactive flow

    ! Misc parameters
    DEBUG = .FALSE.
    uMean = 0D0
    vMean = 0D0

/
