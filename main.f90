!     solve senkin problem using VODE

program test_main
      use chemkin_params, only : initialize_chemkin_workarray

      !   ------- start of user input data ---------

      integer,parameter :: num_spec = 53
      
      ! phisycal values from CFD calculation
      real(8) :: pressure = 1.01325d5 ! Pa
      real(8) :: temperature = 1000d0  ! K
      real(8)    y(num_spec)          ! Mass fractions
      real(8) :: time_end = 1.0d0     ! s
      real(8) :: tolerances(2)        ! Tolerances

      ! output transport data
      real(8) :: D_mix(num_spec) ! mixture diffusion coefficient [CM**2/S]
      real(8) :: Lambda_mix ! mixture thermal conductivity [ERG/CM*K*S]
      real(8) c_p !  mean specific heat at constant pressure [ergs/(gm*K)]
      
      ! Assurme Y has a same secuence as species in ckout
      data y          /0.00d+00,0.00d+00,0.00d+00,2.20d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                       0.00d+00,0.00d+00,0.00d+00,0.00d+00,5.51d-02,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                       0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                       0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                       0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                       0.00d+00,0.00d+00,7.24d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00/
      
      data tolerances /1.d-5, 1.d-5/ ! absolute and relative error tolerances
      
      !   ------- end of user input data ---------

      call initialize_chemkin_workarray(pressure)

      call get_tranport_data(temperature, pressure, y, num_spec, &
                             D_mix, Lambda_mix, c_p)

      ! write(6, *) D_mix
      ! write(6, *) Lambda_mix
      ! write(6, *) c_p

      call solve_senkin(temperature, y, time_end, tolerances, num_spec)

      write(6, *) 'write from test_main'
      write(6, *) temperature, y

end program test_main

!   ----------------------------------------

subroutine get_tranport_data(t_cfd, p_cfd, y_cfd, num_spec, &
                             D_mix, Lambda_mix, c_p)
      use chemkin_params
      
      ! input values
      real(8), intent(in) :: t_cfd    ! K
      real(8), intent(in) :: p_cfd    ! Pa
      real(8) :: y_cfd(num_spec) ! Mass fractions
      integer, intent(in) :: num_spec

      ! output transport data
      ! mixture diffusion coefficient [CM**2/S]
      real(8), intent(out) :: D_mix(num_spec) 
      ! mixture thermal conductivity [ERG/CM*K*S]
      real(8), intent(out) :: Lambda_mix
      !  mean specific heat at constant pressure [ergs/(gm*K)]
      real(8), intent(out) :: c_p

      ! variables for calculations
      real(8) p_calc ! dyne/cm**2
      real(8) :: x_calc(num_spec) ! Mole fractions
      p_calc = p_cfd*10.0d0       ! Pa to dyne/cm**2
      call ckytx(y_cfd, int_ckwk, real_ckwk, x_calc)

      call mcadif(p_calc, t_cfd, x_calc, real_tpwk, D_mix) 
      call mcacon (t_cfd, x_calc, real_tpwk, Lambda_mix)
      call ckcpbs(t_cfd, y_cfd, int_ckwk, real_ckwk, c_p)
end subroutine get_tranport_data

!   ----------------------------------------

subroutine solve_senkin(temperature, y, time_end, tolerances, num_spec)
      use dvode_f90_m

      real(8), intent(inout) :: temperature
      real(8), intent(inout) :: y(num_spec)
      real(8), intent(in)    :: time_end
      real(8), intent(in)    :: tolerances(2)

      real(8) :: z(num_spec+1) ! vector of variables for ODE
      real(8) :: zdot(num_spec+1) ! vector of variables for ODE
      real(8) :: time = 0.0d0  ! s
      real(8) :: dt = 1.0d-3   ! s
      real(8) :: time_out      ! s

      integer :: itask = 2
      integer :: istate = 1

      type(vode_opts) :: options

      external rconp_fex

      ! vode options
      options = set_normal_opts(dense_j=.true., relerr=tolerances(1), &
                                abserr=tolerances(2))

      ! put variables to vector
      z(1) = temperature
      z(2:num_spec+1) = y
      
      time_out = time + dt
      
      do while (time < time_end)
            
            call dvode_f90(f=rconp_fex, neq=num_spec+1, y=z, t=time,  &
                          tout=time_out, itask=itask, istate=istate, & 
                          opts=options)
            
            write(6, *) time, itask, istate

      enddo

      ! call rconp_fex(num_spec+1, time, z, zdot)


      ! return variables from vector
      temperature = z(1)
      y = z(2:num_spec+1)

end subroutine solve_senkin
   
!   ----------------------------------------

subroutine rconp_fex(neq, time, z, zdot)
      use chemkin_params

      real(8), intent(in) :: z(neq)
      real(8), intent(out) :: zdot(neq)

      real(8) rho             ! gm/cm**3
      real(8) c_pb            ! ergs/(gm*K) 
      real(8) volume          ! cm**3/gm
      real(8) :: wdot(kk)     ! moles/(cm**3*sec)
      real(8) :: enthalpy(kk) ! ergs/gm
      integer i

      ! ---------- prepare phisical values ----------
      ! pertubation factor
      ! do i = 1, ii
      !       call ckrdex(-i, real_ckwk, real_ckwk(iprd+i-1))
      ! enddo
      
      ! get rho, c_pb, wdot, enthalpy, wt
      call ckrhoy(pressure, z(1), z(2:), int_ckwk, real_ckwk, rho) 
      call ckcpbs(z(1), z(2), int_ckwk, real_ckwk, c_pb)
      call ckwyp(pressure, z(1), z(2), int_ckwk, real_ckwk, wdot)
      call ckhms(z(1), int_ckwk, real_ckwk, enthalpy)
      
      if (rho < 0.0) then
            write (mout,*) 'Stop, zero density.'
      endif
      volume = 1.0d0/rho

      ! ---------- energey equation ----------
      zdot(1) = zdot(1) + volume*sum(enthalpy*wdot*wt)/c_pb
      
      ! ---------- species equations ----------
      zdot(2:neq) = zdot(2:neq) - wdot*wt*volume

end subroutine rconp_fex