!     test code for obtaining chemkin transport data

!     general parameter indempendent from each cells 
module chemkin_params
      ! unit nuber for input
      integer, parameter :: unit_tplink = 20
      integer, parameter :: unit_cklink = 21

      ! chemkin work array
      integer,      allocatable :: int_ckwk(:)
      real(8),      allocatable :: real_ckwk(:)
      character(6), allocatable :: char_ckwk(:)*16
      
      ! transport work array
      integer,      allocatable :: int_tpwk(:)
      real(8),      allocatable :: real_tpwk(:)
end module chemkin_params

!   ----------------------------------------

program test_main

      !   ------- start of user input data ---------

      integer,parameter :: num_spec = 53
      
      ! phisycal values from CFD calculation
      real(8) :: t_cfd = 298d0       ! K
      real(8) :: p_cfd = 1.01325d5   ! Pa
      real(8) y_cfd(num_spec)        ! Mass fractions

      ! output transport data
      ! mixture diffusion coefficient [CM**2/S]
      real(8) :: D_mix(num_spec) 
      ! mixture thermal conductivity [ERG/CM*K*S]
      real(8) :: Lambda_mix
      !  mean specific heat at constant pressure [ergs/(gm*K)]
      real(8) c_p
      
      ! Assurme Y has a same secuence as species in ckout
      data y_cfd /0.00d+00,0.00d+00,0.00d+00,2.20d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,5.51d-02,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,7.24d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00/

      !   ------- end of user input data ---------

      call initialize_chemkin_workarray

      call get_tranport_data(t_cfd, p_cfd, y_cfd, num_spec, &
                             D_mix, Lambda_mix, c_p)

      write(6, *) D_mix
      write(6, *) Lambda_mix
      write(6, *) c_p

end program test_main

!   ----------------------------------------

subroutine initialize_chemkin_workarray()
      use chemkin_params

      integer, parameter :: unit_stdout = 6

      ! length of work array
      integer len_logi_ckwk
      integer len_int_ckwk
      integer len_real_ckwk
      integer len_char_ckwk

      ! open input unit
      open(unit_cklink, form='unformatted', file='./link/cklink')
      open(unit_tplink, form='unformatted', file='./link/tplink')

      !   ------- initialize chemkin work array ---------

      CALL cklen(unit_cklink, unit_stdout, len_int_ckwk, len_real_ckwk, len_char_ckwk)

      allocate(int_ckwk(len_int_ckwk))
      allocate(real_ckwk(len_real_ckwk))
      allocate(char_ckwk(len_char_ckwk))

      call ckinit(len_int_ckwk, len_real_ckwk, len_char_ckwk, unit_cklink, &
                  unit_stdout, int_ckwk, real_ckwk, char_ckwk)

      !   ------- initialize transport work array ---------

      call mclen(unit_tplink, unit_stdout, len_int_tpwk, len_real_tpwk)
      
      allocate(int_tpwk(len_int_tpwk))
      allocate(real_tpwk(len_real_tpwk))

      call mcinit(unit_tplink, unit_stdout, len_int_tpwk, len_real_tpwk, &
                  int_tpwk, real_tpwk)

end subroutine initialize_chemkin_workarray

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