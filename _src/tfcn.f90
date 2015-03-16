  function tfcn(t)
    real(kind=8), intent(in) :: t
    real(kind=8) :: pi, tfcn
    real(kind=8), parameter :: t_period = 5.d0

    pi = atan2(0.d0,-1.d0)
    tfcn = dcos(pi*t/t_period)

  end function tfcn
