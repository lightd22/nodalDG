  function tfcn(t)
    USE commonTestParameters
    ! Inputs
    DOUBLE PRECISION, intent(in) :: t
    ! Outputs
    DOUBLE PRECISION :: tfcn
    ! Local variables
    DOUBLE PRECISION, parameter :: t_period = 5.d0

    tfcn = dcos(pi*t/t_period)

  end function tfcn
