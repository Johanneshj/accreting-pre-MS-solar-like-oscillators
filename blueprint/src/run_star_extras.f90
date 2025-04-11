! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff
      
      implicit none

      logical, save :: save_model
      real(dp), save :: next_mass, Lacc, Ladd
      real(dp):: accretion_factor


      logical :: table_loaded = .false.
      logical :: reached_c12 = .false.
      logical :: need_to_save
      real(dp) :: Xc_save, Xc_save_step
      logical :: cahnged_mesh_options = .false.


      type acc_table_t
      character(len=255) :: filename
      integer            :: n
      real(kind=8)       :: min_time, max_time
      real(kind=8)       :: total_mass
      real(kind=8)     :: final_accretion
      real(kind=8), allocatable, dimension(:) :: time, rate
      end type
      

      type(acc_table_t), target :: acc_table

! My own stuff (Johannes)
      real(dp) :: Age_ZAMS 
      logical :: zams_age_saved = .false.
      logical :: final_mass_reached = .false.

      real(dp) :: init_M_Johannes
      real(dp) :: star_mass_Johannes 
      real(dp) :: star_mass_Johannes_old

      ! these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         s% other_energy => energy_routine
         s% other_adjust_mdot => other_adjust_mdot

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

         s% job% initial_h2 = s% x_ctrl(1) * 1d-6
         s% job% initial_he3 = s% x_ctrl(2) * 1d-6
         s% job% initial_he4 = s% initial_y
         s% job% initial_h1 = 1 - s% initial_y - s% job% initial_h2 - s% job% initial_he3 - s% initial_z
      
         s% accretion_h2 = s% x_ctrl(1) * 1d-6
         s% accretion_he3 = s% x_ctrl(2) * 1d-6
         s% accretion_he4 = s% initial_y
         s% accretion_h1 = 1 - s% initial_y - s% job% initial_h2 - s% job% initial_he3 - s% initial_z

         s% kap_rq% Zbase = s% initial_z

         Xc_save = s% job% initial_h1 - 0.001
         Xc_save_step = 0.01

      end subroutine extras_controls



      
      subroutine other_adjust_mdot(id, ierr)
            use star_def
            implicit none
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            integer :: step, i
            type(star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr) ! retrieve star id


            if (s% star_age > 1d-10) then
               s% mstar_dot =  acc_rate(s% time_old, s% time_old + s% dt)*Msun/secyer
            else
               s% mstar_dot = acc_table% rate(1)*Msun/secyer
            end if    

            DO i=1, s% nz
                  IF (s% star_age < acc_table% time(i) ) THEN
                     step = i-1
                     exit
                  END IF
            END DO



            accretion_factor = get_accretion_factor(s% time_old + s% dt, s% mstar / Msun, s% x_ctrl(3))

            write(*,*) "current Accretion rate: " ,  s% mstar_dot
            s% mstar_dot = s% mstar_dot*accretion_factor

            write(*,*) "current Accretion Factor: " , accretion_factor
            write(*,*) "scaled Accretion rate: " ,  s% mstar_dot    
            
            !if (s% mstar_dot/Msun*secyer > 1d-4) then
            !   s% use_dPrad_dm_form_of_T_gradient_eqn = .True.
            !   s% eps_mdot_factor = 0.0
	         !   s% eps_mdot_leak_frac_factor = 0.0
            !else if (s% mstar_dot/Msun*secyer > 5d-5) then
            !   s% use_dPrad_dm_form_of_T_gradient_eqn = .false.
            !   s% eps_mdot_factor = 0.5
            !   s% eps_mdot_leak_frac_factor = 0.5
            !else 
            !   s% use_dPrad_dm_form_of_T_gradient_eqn = .false.
            !   s% eps_mdot_factor = 1d0
	         !   s% eps_mdot_leak_frac_factor = 1d0
            !end if
            
            if (s% star_mass .ge. s% x_ctrl(3)) then
               s% use_dPrad_dm_form_of_T_gradient_eqn = .false.
	            s% use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn = .false.
               s% eps_mdot_factor = 1.0
	            s% eps_mdot_leak_frac_factor = 1.0
            else
               s% use_dPrad_dm_form_of_T_gradient_eqn = .true.
	            s% use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn = .true.
               s% eps_mdot_factor = 0.0
	            s% eps_mdot_leak_frac_factor = 0.0
            end if
            
            !if (s% mstar_dot/Msun*secyer > 5d-5) then
            !   s% delta_mdot_limit = 0.1
            !   s% delta_mdot_hard_limit = 0.25 
            if (s% mstar_dot/Msun*secyer > 5d-6) then
               s% delta_mdot_limit = 0.2
               s% delta_mdot_hard_limit = 0.5 
            else
               s% delta_mdot_limit = -1
               s% delta_mdot_hard_limit = -1
            end if
    
            if (MODULO(s% model_number, s% terminal_interval) == 0) then  
                  write(*,*) "Adjust_mdot: ", s% mstar_dot/Msun*secyer
            end if

        end subroutine other_adjust_mdot

      subroutine energy_routine(id, ierr)
            use const_def
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s

            real(dp) :: mass, mdot, mdot_msun, alpha, radius, mcumul
            real(dp) :: xposk, dif, xpos, thingie, thingie2, Mstar, mass_dif
            integer :: k, xposk_orig

            real(dp), parameter :: mdot_break = 3.5e-5
            real(dp), parameter :: alpha_high = 0.2
            real(dp), parameter :: alpha_low = 0.005
            real(dp), parameter :: delta_break = 2.0/3.0*1e-6
            real(dp) :: numerator, denominator

            ierr = 0
            call star_ptr(id, s, ierr)

               mass = s% mstar
               mdot = s% mstar_dot
               mdot_msun = s% mstar_dot / Msun* secyer
               numerator = alpha_low * exp(mdot_break/delta_break) + alpha_high * exp(mdot_msun/delta_break)
               denominator = exp(mdot_break/delta_break) + exp(mdot_msun/delta_break)

               alpha = numerator / denominator
  
            write(*,*) "Alpha is set to", alpha, denominator, numerator, mdot_msun
            write(*,*) 'accretion rate in Msun units:', mdot_msun

            mass = s% mstar
            mdot = s% mstar_dot

            if (mdot == 0.) return

            radius = s% r(1)

            mcumul = 0.
            k = 0
            Lacc = 0.5 * (standard_cgrav * (s% mstar_dot) * (mass)) / radius

            Ladd = alpha/2*standard_cgrav * (s% star_mass *Msun) * (s% mstar_dot) / (s% r(1)) ![erg/s]
            Lacc = safe_log10((1-alpha)/2*standard_cgrav * (s% star_mass *Msun) * (s% mstar_dot) / (s% r(1)) / Lsun) ![erg/s]
            
            Mstar=(s% star_mass *Msun) ![g]


   
            !xpos = s% x_ctrl(5)
            xpos = ( 2.5d-3 * exp(-6/0.5) + 2.0d-1 * exp(safe_log10(mdot_msun)/0.5) ) / ( exp(-6/0.5) + exp(safe_log10(mdot_msun)/0.5) )
            write(*,*) 'xpos:', xpos
            do k = 1, s% nz
               if (s% star_mdot > 0.0d0 .and. s% m(k)>=(1.0d0-xpos)*(s% star_mass)*Msun ) then
                     s% extra_heat(k) = 2.0d0*Ladd/(Mstar*xpos**2)*(s% m(k)/Mstar-(1.0d0-xpos))
               end if
            end do

  
      end subroutine energy_routine


      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         integer :: nlines, io, i
         character(2) :: history_num
         character(len = 15) :: intermediate_step
         type (star_info), pointer :: s
         real(kind=8) :: offset
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         if (s% x_integer_ctrl(1) < 10) then
            Write( history_num, '(i1)')  s% x_integer_ctrl(1)
            history_num = '0' // history_num 
         else
            Write( history_num, '(i2)')  s% x_integer_ctrl(1)
         end if

         intermediate_step = '../ext_50myr/' // history_num
         acc_table% filename = intermediate_step // '_ex50.dat'

         write(*,*) acc_table % filename


         nlines = 0
         open(unit = 10, file = acc_table% filename)
         DO 
            READ(10,*,iostat=io)
            IF (io/=0) EXIT
            nlines = nlines + 1 
         END DO 
         close(10)


         acc_table%n = nlines

         open(unit = 10, file = acc_table% filename)

         allocate(acc_table % time(acc_table % n))
         allocate(acc_table % rate(acc_table % n))
         DO i = 1, acc_table %n
            READ(10,*) acc_table % time(i), acc_table % rate(i)
         END DO

         offset = acc_table % time(1)
         DO i = 1, acc_table %n
            acc_table% time(i) = acc_table% time(i) - offset
         END DO

         acc_table% time = acc_table% time * secyer ! Time is in Myrs

         acc_table% min_time = acc_table% time(1)
         acc_table% max_time = acc_table% time(acc_table% n)

         acc_table% total_mass = 0.01
         DO i = 1, acc_table %n-1
            acc_table% total_mass = acc_table% total_mass + (acc_table % rate(i)+acc_table % rate(i+1))* &
            (acc_table% time(i+1) - acc_table% time(i))/(2* secyer)
         END DO

         i = nlines
         DO 
         if (acc_table % rate (i) > 0) then
            acc_table% final_accretion = acc_table % time(i)/secyer
            exit
         end if
         i = i-1
         END DO
      end subroutine extras_startup


      function get_timescale(age, mass) result(max_timestep)

         implicit none
         integer :: k, idn, i4, i5, i6
         real(dp) :: max_timestep
         real(dp), intent(in) :: age, mass


         max_timestep = 20*secyer
         if (mass/acc_table% total_mass < 0.05) return
         if (mass/acc_table% total_mass > 0.99) then
            max_timestep = 1000*secyer
            return
         end if  

         idn = 0
         i6 = acc_table%n 
         do k=1, acc_table %n
            if (acc_table % time(k) > age) then
               idn = k 
               exit
            end if 
         end do
         
         if (idn == 0) return
         do k=idn, acc_table %n
            if (acc_table % rate(k) > 1d-5) then
               i6 = k
               exit
            end if
         end do

         write(*,*) 'i6 weird thiingie:', i6
         if (i6 == 0) return
         max_timestep = min(max(20.0*secyer, (acc_table %time(i6) - acc_table % time(idn))/500), 10.0*secyer)

      end function get_timescale
      


     function acc_rate(tstart, tend)
         use star_def
         implicit none
         real(kind=8) :: tstart, tend, acc_rate, lb, mp, ub
         real(kind=8), pointer, dimension(:) :: t, acc
         real(kind=8) :: epsilon, mdot_max, mdot_min, epsilon_limit
         real(kind=8) :: lower_side, upper_side
         integer, save:: ilast=1, ioff=1
         integer      :: i, step, n
         integer :: ilast_old, ioff_old

         if (tend > acc_table% max_time ) then
         acc_rate = 0.0
         return
         elseif (tend < acc_table% min_time) then
         acc_rate = acc_table% rate(1)
         return
         endif

         n   =  acc_table% n
         t   => acc_table% time
         acc => acc_table% rate

         DO i=1, n
            IF (tend < t(i)) THEN
                  step = i-1
                  exit
            END IF
         END DO

         acc_rate = acc(step) + (tend - t(step))/(t(step+1)-t(step))*(acc(step+1)-acc(step))

      end function acc_rate


      function get_accretion_factor(time, star_mass, end_mass)  result(accretion_factor)
         use star_def
         implicit none
         real(kind=8) :: end_mass, time, accretion_factor, star_mass, to_accrete, mass_needed
         real(kind=8), pointer, dimension(:) :: t, acc
         real(kind=8) :: epsilon, mdot_max, mdot_min, epsilon_limit
         real(kind=8) :: lower_side, upper_side
         integer      :: i, step, n
         integer :: ilast_old, ioff_old

         n   =  acc_table% n
         t   => acc_table% time
         acc => acc_table% rate

      
         DO i=1, n
            IF (time < t(i)) THEN
                  step = i-1
                  exit
            END IF
         END DO

         to_accrete = 0
         DO i = step, acc_table %n-1
         to_accrete = to_accrete + (acc(i)+acc(i+1)) * (t(i+1) - t(i)) / (2 * secyer)
         
         END DO

         write(*,*) 'current_expected total mass (not scaled):', star_mass + to_accrete

         mass_needed = end_mass - star_mass

         write(*,*) 'current mass:', star_mass, 'out of', end_mass
         write(*,*) 'mass needed:', mass_needed

         if (mass_needed < 0) then
            accretion_factor = 0
         else
            accretion_factor = mass_needed / to_accrete
         end if

         write(*,*) 'current_expected total mass with factor (scaled!):', star_mass + accretion_factor * to_accrete

      end function get_accretion_factor


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0

         if (s% star_age < acc_table % final_accretion*0.99) then        
            s% delta_lgL_limit = 1d-3 
            s% delta_lgL_hard_limit = 0.0025
            s% max_timestep = get_timescale(s% star_age *secyer ,s% star_mass)
            if (acc_rate(s% time_old, s% time_old + s% dt) > s% x_ctrl(6)) then
                s% max_timestep = min(s% max_timestep, 2.5*secyer)
                s% delta_lgL_limit = 5d-4 
                s% delta_lgL_hard_limit = 0.00125
                write(*,*) 'max_timestep:', s% max_timestep/secyer
            end if
         else
            s% max_timestep = 0
            if (.not. cahnged_mesh_options) then
              s % prune_bad_cz_min_Hp_height = 0
              s % prune_bad_cz_min_log_eps_nuc = -99
              s % min_convective_gap = -1
              s % max_resid_jump_limit = 1d6
              s % corr_coeff_limit = 0.05d0

              s % mesh_logX_species(1) = 'h1'
              s % mesh_logX_min_for_extra(1) = -10
              s % mesh_dlogX_dlogP_extra(1) = 0.15             ! resol coeff for chemical gradients
              s % mesh_dlogX_dlogP_full_on(1) = 1d-12           ! additional resol on for gradient larger than this
              s % mesh_dlogX_dlogP_full_off(1) = 1d-12         ! additional resol off for gradient smaller than this             
              cahnged_mesh_options = .true.
            end if
         end if

      if (s% star_mass < 0.01) then
               s% atm_option = 'table'
               s% atm_table = 'tau_100'
      else if (s% star_mass < 0.02) then
               s% atm_option = 'T_tau'
               s% atm_T_tau_relation = 'Eddington'
               s% atm_T_tau_opacity = 'iterated' ! 'varying'
               s% atm_T_tau_errtol = 1d-7
               s% atm_T_tau_max_steps = 50000  
      else if (s% star_mass < 0.08) then
               s% atm_option = 'T_tau'
               s% atm_T_tau_relation = 'Eddington'
               s% atm_T_tau_opacity = 'iterated'
               s% atm_T_tau_errtol = 1d-9
      else
               s% atm_option = 'T_tau'
               s% atm_T_tau_relation = 'Eddington'
               s% atm_T_tau_opacity = 'iterated'
               s% atm_T_tau_errtol = 1d-9
      end if

      end function extras_start_step

      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr

         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         

         s% xtra(3) = s% center_h1

         if (.not. zams_age_saved) then
            if ( s% xtra(3) .lt. 0.999*(s% job% initial_h1) ) then
               Age_ZAMS = s% star_age
               zams_age_saved = .true.
               s% need_to_save_profiles_now = .true.
            end if
         end if

         s% xtra(1) = s% star_mass
         init_M_Johannes = s% x_ctrl(3)
         star_mass_Johannes = s% xtra(1) / init_M_Johannes
         star_mass_Johannes_old = s% xtra_old(1) / init_M_Johannes

         if (.not. zams_age_saved) then
            if (star_mass_Johannes .lt. init_M_Johannes) then
               if ( ABS( floor(star_mass_Johannes*100) - floor(star_mass_Johannes_old*100) ) .ne. 0 ) then
                  s% need_to_save_profiles_now = .true.
               end if
            else if (star_mass_Johannes .ge. init_M_Johannes) then
               final_mass_reached = .true.
            end if
         end if

         s% xtra(4) = 10**(s% log_surface_radius)

         if (.not. zams_age_saved) then
            if (final_mass_reached) then
               if ( ABS( floor(s% xtra(4)*100) - floor(s% xtra_old(4)*100) ) .ne. 0 ) then
                  write(*,*) "Radius changed by 1% on the ZAMS --> Saving profile"
                  s% need_to_save_profiles_now = .true.
               end if 
            end if
         end if

         s% xtra(2) = s% star_age
         
         if (.not. zams_age_saved) then
            if ( ABS( floor( (s% xtra(2))/1e6 ) - floor( (s% xtra_old(2))/1e6 ) ) .ne. 0 ) then
               write(*,*) "Age changed by 1 Myr on the ZAMS --> Saving profile"
               s% need_to_save_profiles_now = .true.
            end if 
         end if
       
         if (zams_age_saved) then
            if ( ABS( floor( (s% xtra(3))*10 ) - floor( (s% xtra_old(3))*10 ) ) .ne. 0 ) then
               write(*,*) "Central h1 dropped by 10% --> saving profile"
               s% need_to_save_profiles_now = .true.
            end if 
         end if
         
         if (zams_age_saved) then
            if ( s% xtra(3) .le. 0.99*(s% job% initial_h1)) then
               s% need_to_save_profiles_now = .true.
               extras_check_model = terminate
               write(*, *) 'pre-ms has ended'
            end if
         end if


         if (MODULO(s% model_number, s% terminal_interval) == 0) then  
               write(*,*) "Expected total mass:", acc_table % total_mass, "current ratio:", s% star_mass/acc_table % total_mass
               write(*,*) "Star age:", s% star_age, "accretion time ratio:", s% star_age /acc_table % final_accretion
         end if

         if (extras_check_model == terminate) s% termination_code = t_extras_check_model

      end function extras_check_model



      subroutine save_at_Xc_c12 (id, Xc_save, make_save, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: Xc_save
         logical, intent(out) :: make_save
         integer, intent(out) :: ierr

         type(star_info), pointer :: s
         real(dp) :: Xc_hi, Xc_lo, Xc
         real(dp) :: dXc

         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'Error: save_at_Xc: star_ptr failed'
            return
         endif

         if (Xc_save .le. 0) then
            make_save = .false.
            return
         endif

         Xc  = s% center_c12
         if (ABS(Xc - Xc_save) > 5d-4) then
            dXc = 5d-7
         else
            dXc = 5d-9
         endif

         Xc_hi = Xc_save+5d-4
         Xc_lo = Xc_save-2d-5
         call set_grid_dXc(id, Xc,  Xc_hi, Xc_lo, dXc, Xc_save, ierr)
         if (ABS(Xc - Xc_save) < 1d-9) then
            make_save = .true.
            return
         endif
         make_save = .false.
         return
      end subroutine save_at_Xc_c12

      subroutine save_at_Xc (id, Xc_save, make_save, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: Xc_save
         logical, intent(out) :: make_save
         integer, intent(out) :: ierr

         type(star_info), pointer :: s
         real(dp) :: Xc_hi, Xc_lo, Xc
         real(dp) :: dXc

         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'Error: save_at_Xc: star_ptr failed'
            return
         endif

         if (Xc_save .le. 0) then
            make_save = .false.
            return
         endif

         Xc  = s% center_h1
         if (ABS(Xc - Xc_save) > 0.04) then
            dXc = 0.005
         else
            dXc = 0.000005
         endif

         Xc_hi = Xc_save+0.005
         Xc_lo = Xc_save-0.002
         call set_grid_dXc(id, Xc,  Xc_hi, Xc_lo, dXc, Xc_save, ierr)
         if (ABS(Xc - Xc_save) < 1d-6) then
            make_save = .true.
            return
         endif
         make_save = .false.
         return
      end subroutine save_at_Xc


    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! if Xc outside the range specified by Xc_hi and Xc_lo, then we
    ! resort back to the values in mesa/star/defaults/controls.defaults
      subroutine set_grid_dXc(id, Xc,  Xc_hi, Xc_lo, dXc, Xc_save, ierr)
         integer, intent(in) :: id
         real(dp), optional, intent(in) :: Xc_hi, Xc_lo
         real(dp), intent(in) :: Xc, dXc, Xc_save
         integer, intent(out) :: ierr

         type(star_info), pointer :: s
         real(dp) :: Xc_next, max_Xc_chb, delta_XH_cntr_limit, delta_XC_cntr_limit
         real(dp), parameter :: min_Xc_chb = 1d-12
         real(dp), parameter :: lowest_delta_XH_cntr_limit = 5d-7
         real(dp), parameter :: lowest_delta_XC_cntr_limit = 5d-10
         logical :: is_chb

         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'Error: set_grid_dXc: star_ptr failed'
            return
         endif

         s% delta_XH_cntr_limit  = 1d-3  
         s% delta_XC_cntr_limit  = 1d-4  

         max_Xc_chb = (s% job% initial_h1 + s% job% initial_h2)

         is_chb = ((Xc .ge. min_Xc_chb) .and. (Xc <= max_Xc_chb))

         Xc_next = Xc - dXc

         if (Xc_save > 1d-3) then
            delta_XH_cntr_limit     = (Xc - Xc_save)/20
            delta_XH_cntr_limit     = max(delta_XH_cntr_limit, lowest_delta_XH_cntr_limit)
            s% delta_lg_XH_cntr_limit  = delta_XH_cntr_limit
            write(*,*) s% delta_XH_cntr_limit, Xc, Xc_next,  dXc
         else
            delta_XC_cntr_limit     = (Xc - Xc_save)/20
            delta_XC_cntr_limit     = max(delta_XC_cntr_limit, lowest_delta_XC_cntr_limit)
            s% delta_XC_cntr_limit  = delta_XC_cntr_limit
            write(*,*) s% delta_XC_cntr_limit, Xc, Xc_next,  dXc
         end if

      end subroutine set_grid_dXc


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 6
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: max_temp 
         ! integral
         integer :: max_pos, k
          !i, list_size
         !double precision, allocatable, dimension(:) :: dcdr, product, steps
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         

         names(1) = 'Ladd'
         vals(1) = Ladd/Lsun

         names(2) = 'inversion_h2'
         names(3) = 'inversion_li7'
         names(4) = 'inversion_c12'
         names(5) = 'inversion_h1'
         names(6) = 'inversion_position'

         max_temp = -99 
         do k=1, s% nz
            if (s % T(k) > max_temp) then
              max_temp = s% T(k)
              max_pos = k
            end if
         end do

         vals(2) = s% xa(2, max_pos)
         vals(3) = s% xa(5, max_pos)
         vals(4) = s% xa(8, max_pos)
         vals(4) = s% xa(1, max_pos)
         vals(6) = s% m(max_pos) / s% m(1)

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 3
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         names(1) = 'log_opacity'
         names(2) = 'grad_temperature'
         names(3) = 'grad_L'
         do k = 1, s% nz
            vals(k,1) = safe_log10(s% opacity(k))
            vals(k,2) = safe_log10(s% grad_temperature(k))
            vals(k,3) = safe_log10(s% gradL(k))
         end do
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return
      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve


      end module run_star_extras
      
