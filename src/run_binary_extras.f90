! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
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

      module run_binary_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use math_lib
      
      implicit none

      ! debugging flag
      logical, parameter :: dbg = .false.

      ! b% xtra(x_inclination) contains the inclination value (changing at each timestep)
      integer, parameter :: x_inclination = 1
      ! b% xtra(x_idot) contains the derivative of the inclination
      integer, parameter :: x_idot = 2
      ! b% xtra(x_omega_eq) contains the value of the equilibrium spin
      integer, parameter :: x_omega_eq = 3

      ! b% lxtra(lx_convective_envelope) is true when envelope of star is convective
      integer, parameter :: lx_convective_envelope = 1
      ! b% lxtra(lx_use_qin_E2) is true when using fits of E2 in dynamical tides from Qin+ 2018
      integer, parameter :: lx_use_qin_E2 = 2

      ! inclination value
      real(dp) :: initial_inclination, idot_step, di_step, i_evolve
      real(dp), parameter :: min_inclination = 0d0
      real(dp), parameter :: max_abs_di = 1d0
      real(dp), parameter :: fi_hard = 0.1
      real(dp), parameter :: fi_limit = 1d-3

      ! minimum value for the fraction of convective envelope
      real(dp), parameter :: min_convective_fraction = 0.1d0
      
      contains
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! apply tides in non-coplanar case
         b% other_sync_spin_to_orbit => sync_non_coplanar
         b% other_edot_tidal => edot_non_coplanar

         ! Set these function pointers to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
         b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup => extras_binary_startup
         b% extras_binary_start_step => extras_binary_start_step
         b% extras_binary_check_model => extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve => extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message, 
          b% warn_binary_extra =.false.
         
      end subroutine extras_binary_controls


      subroutine sync_non_coplanar(id, nz, osep, qratio, rl, dt_next, Ftid, sync_type, sync_mode, ierr)
         integer, intent(in) :: id
         integer, intent(in) :: nz
         real(dp), intent(in) :: osep  ! orbital separation (cm)
         real(dp), intent(in) :: qratio  ! mass_other_star/mass_this_star
         real(dp), intent(in) :: rl  ! roche lobe radius (cm)
         real(dp), intent(in) :: dt_next  ! next timestep
         real(dp), intent(in) :: Ftid  ! efficiency of tidal synchronization. (time scale/Ftid)
         character (len = strlen), intent(in) :: sync_type  ! synchronization timescale
         character (len = strlen), intent(in) :: sync_mode  ! where to put/take angular momentum
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         integer :: k
         real(dp) :: moment_of_inertia, rGyr_squared, m, r_phot, porb
         real(dp) :: omega_orb
         real(dp) :: a1, a2, f_sync, t_sync, eff_factor
         real(dp), dimension(nz) :: j_sync, delta_j
         real(dp) :: i_step
         logical :: convective_envelope

         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_ptr'
            return
         end if

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         s% extra_jdot(:) = 0d0

         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if
         porb = b% period
         
         omega_orb = 2d0*pi/b% period
         do k = 1, nz
            j_sync(k) = omega_orb*s% i_rot(k)
         end do
         
         i_step = b% xtra(x_inclination)
         
         ! compute gyration radius
         moment_of_inertia = dot_product(s% i_rot(:s% nz), s% dm_bar(:s% nz))
         rGyr_squared = (moment_of_inertia/(m*r_phot*r_phot))

         ! compute sync timescale based on different energy transport on the star envelope
         if (sync_type == 'Hut_conv') then
            t_sync = 3d0*k_div_T(b, s, .true.) * (qratio*qratio/rGyr_squared) * pow6(r_phot/osep)
         else if (sync_type == 'Hut_rad') then
            t_sync = 3d0*k_div_T(b, s, .false.) * (qratio*qratio/rGyr_squared) * pow6(r_phot/osep)
         else if (sync_type == 'Orb_period') then
            t_sync = b% period
         else
            ierr = -1
            write(*,*) 'unrecognized sync_type', sync_type
            return
         end if
         t_sync = (1d0/t_sync) / Ftid

         ! set change in angular momentum on donor star as MESA does.
         ! one thing to note is that it follows the work of Wellstein (2001) where the
         ! ang mom (J) is syncronized assuming constant j_orb & t_sync during a timestep
         ! thus delta_j is a decaiying 
         !! a1 = f2(b% eccentricity) * cos(i_step)
         !! a2 = pow(1-pow2(b% eccentricity), 1.5d0) * f5(b% eccentricity) * &
         !!    (3d0+cos(2*i_step)) / 4d0
         a1 = f2(b% eccentricity)
         a2 = pow(1-pow2(b% eccentricity), 1.5d0) * f5(b% eccentricity)
         eff_factor = 1 / a2
         do k = 1, nz
            delta_j(k) = (1d0-exp(-dt_next/(t_sync*eff_factor)))*(s% j_rot(k) - a1/a2*j_sync(k))
         end do

         if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
            b% t_sync_1 = t_sync
            if (b% model_twins_flag) b% t_sync_2 = t_sync
         else
            b% t_sync_2 = t_sync
         end if

         if (.not. b% doing_first_model_of_run) then
            do k = 1, nz
               s% extra_jdot(k) = s% extra_jdot(k) - delta_j(k)/dt_next
            end do
         end if

      end subroutine sync_non_coplanar


      subroutine edot_non_coplanar(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b
         
         include 'formats.inc'

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         b% edot_tidal = 0d0

         if (b% point_mass_i /= 1) then
            if (b% circ_type_1 == 'Hut_conv') then
               b% edot_tidal = edot_tidal_non_coplanar(b, b% s1, .true.)
            else if (b% circ_type_1 == 'Hut_rad') then
               b% edot_tidal = edot_tidal_non_coplanar(b, b% s1, .false.)
            else
               write(*,*) 'unrecognized circ_type_1', b% circ_type_1
            end if
         end if
         if (b% point_mass_i /= 2) then
            if (b% circ_type_2 == 'Hut_conv') then
               b% edot_tidal = b% edot_tidal + edot_tidal_non_coplanar(b, b% s2, .true.)
            else if (b% circ_type_2 == 'Hut_rad') then
               b% edot_tidal = b% edot_tidal + edot_tidal_non_coplanar(b, b% s2, .false.)
            else
               write(*,*) 'unrecognized circ_type_2', b% circ_type_2
            end if
         end if

         if (b% model_twins_flag) then
            b% edot_tidal = b% edot_tidal + b% edot_tidal
         end if  
         
         if (dbg) write(*,1) 'edot', b% edot_tidal

      end subroutine edot_non_coplanar


      real(dp) function edot_tidal_non_coplanar(b, s, convective_envelope) result(edot)
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         logical, intent(in) :: convective_envelope
         integer :: k
         real(dp) :: m, porb, r_phot, osep, qratio, omega_s, omega_sync
         real(dp) :: i_step

         include 'formats.inc'
         
         i_step = b% xtra(x_inclination)

         edot = 0d0

         porb = b% period
         omega_sync = 2d0*pi/b% period
         omega_s = s% omega_avg_surf
         osep = b% separation
         
         qratio = b% m(b% a_i) / b% m(b% d_i)
         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            qratio = 1.0d0/qratio
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if
         
         ! eq. (10) of Hut, P. 1981, A&A, 99, 126
         edot = -27.0d0 * qratio * (1+qratio) * pow8(r_phot/osep) &
            * b% eccentricity * pow(1-pow2(b% eccentricity), -6.5d0) * b% Ftid_1
         ! add multiplication by (k/T), eq. (29) of Hurley et al. 2002
         edot = edot * k_div_T(b, s, convective_envelope)
         ! add terms dependant on omega & inclination
         edot = edot * (f3(b% eccentricity) - &
            11d0/18d0 * omega_s/omega_sync * f4(b% eccentricity) * &
            pow(1-pow2(b% eccentricity), 1.5d0) * cos(i_step))

      end function edot_tidal_non_coplanar


      real(dp) function idot_non_coplanar(b, s, coplanar_type) result(idot)
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         character (len = strlen), intent(in) :: coplanar_type
         real(dp) :: i_step, par
         real(dp) :: porb, omega_sync, omega_s, osep, qratio, m, r_phot
         real(dp) :: rGyr_squared, moment_of_inertia
         logical :: convective_envelope
         
         include 'formats.inc'

         if (b% doing_first_model_of_run .or. b% xtra(x_inclination) <= min_inclination) then
            idot = 0d0
            return
         end if

         i_step = b% xtra(x_inclination)
         
         porb = b% period
         omega_sync = 2d0*pi/b% period
         omega_s = s% omega_avg_surf
         osep = b% separation
         
         qratio = b% m(b% a_i) / b% m(b% d_i)
         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            qratio = 1.0d0/qratio
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if
         
         if (coplanar_type == 'Hut_conv') then
            convective_envelope = .true.
         else if (coplanar_type == 'Hut_rad') then
            convective_envelope = .false.
         else
            write(*,*) 'unrecognized coplanar_type', coplanar_type
         end if

         ! calculate the gyration radius squared
         moment_of_inertia = dot_product(s% i_rot(:s% nz), s% dm_bar(:s% nz))
         rGyr_squared = (moment_of_inertia/(m*r_phot*r_phot))

         idot = -3d0 * k_div_T(b, s, convective_envelope) * &
            (qratio*qratio/rGyr_squared) * pow6(r_phot/osep) * (omega_sync/omega_s) * &
            pow(1-pow2(b% eccentricity), -6d0) * sin(i_step)

         par = f2(b% eccentricity) - 0.5d0 * f5(b% eccentricity) * &
            ((omega_s/omega_sync) * cos(i_step) * pow(1-pow2(b% eccentricity), 1.5d0) + &
            r_phot*r_phot*osep*omega_s*omega_s*rGyr_squared * (1-pow2(b% eccentricity)) / &
            (b% m(b% a_i) * standard_cgrav))

         idot = par*idot

         if (dbg) write(*,1) 'idot', idot

      end function idot_non_coplanar


      logical function is_convective(s)
         type(star_info), pointer :: s
         real(dp) :: m_conv, f_conv
         integer :: k

         include 'formats.inc'

         is_convective = .false.

         ! two different ways to check for convective envelope. One for MS evolution, 
         ! the other for post-MS evolution
         if (s% he_core_mass > 0d0 .and. s% center_h1 < 1d-6) then
            if (s% he_core_k == 1 .or. (s% star_mass-s% he_core_mass) < 1d-3) then
               m_conv = 0d0
               f_conv = 0d0
            else
               m_conv = 0d0
               do k = 1, s% he_core_k
                  if (s% mixing_type(k) == convective_mixing) m_conv = m_conv+s% dm(k)
               end do
               f_conv = (m_conv/Msun) / (s% star_mass-s% he_core_mass)
            end if
         else
            m_conv = 0d0
            do k = 1, s% nz
               if (s% m(k) / s% mstar < 0.9) exit
               if (s% mixing_type(k) == convective_mixing) m_conv = m_conv+s% dm(k)
            end do
            f_conv = (m_conv/Msun) / s% star_mass
         end if
         
         if (f_conv > min_convective_fraction) is_convective = .true.

         if (dbg) then
            write(*,1) 'f_conv', f_conv
            write(*,14) 'is_convective', is_convective
         end if

      end function is_convective


      real(dp) function k_div_T(b, s, has_convective_envelope)
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         logical, intent(in) :: has_convective_envelope
         integer :: k
         real(dp) osep, qratio, m, r_phot, porb, m_env, r_env, tau_conv, P_tid, f_conv
         real(dp) :: E2

         ! k/T computed as in Hurley, J., Tout, C., Pols, O. 2002, MNRAS, 329, 897
         ! Kudos to Francesca Valsecchi for help implementing and testing this

         k_div_T = 0d0

         osep = b% separation
         qratio = b% m(b% a_i) / b% m(b% d_i)
         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            qratio = 1d0/qratio
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if
         porb = b% period

         if (has_convective_envelope) then
            m_env = 0d0
            r_env = 0d0
            do k = 1, s% nz
               if (s% mixing_type(k) /= convective_mixing .and. &
                   s% rho(k) > 1d5*s% rho(1)) then
                  r_env = (r_phot-s% r(k))/Rsun
                  m_env = (s% m(1) - s% m(k))/Msun
                  exit
               end if
            end do
            tau_conv = 0.431d0*pow(m_env*r_env* &
               (r_phot/Rsun-r_env/2d0)/3d0/s% L_phot, one_third) * secyer
            if (1d0/porb /= s% omega_avg_surf/(2d0*pi)) then
               P_tid = 1d0/abs(1d0/porb-s% omega_avg_surf/(2d0*pi))
               f_conv = min(1.0d0, pow(P_tid/(2d0*tau_conv),b% tidal_reduction))
            else
               f_conv = 1d0
            end if
            ! Belczynski+2008
            k_div_T = 2d0/21d0*f_conv/tau_conv*m_env/(m/Msun)
         else
            if (b% lxtra(lx_use_qin_E2)) then
               do k = s% nz, 1, -1
                  if (s% brunt_N2(k) >= 0) exit
               enddo
               E2 = 10**(-0.42) * (s% r(k)/r_phot)**(7.5)
               k_div_T = sqrt(standard_cgrav*m*r_phot**2/pow5(osep)/(Msun/pow3(Rsun)))
               k_div_T = k_div_T*pow(1d0+qratio,5d0/6d0)
               k_div_T = k_div_T * E2
            else
               !NOTE:There is a typo in eq. (42) of Hurley+2002, 
               !correct expression is given in footnote 3 of
               !Sepinsky+2007
               k_div_T = 1.9782d4*sqrt(m*r_phot*r_phot/pow5(osep)/(Msun/pow3(Rsun)))
               k_div_T = k_div_T*pow(1d0+qratio, 5d0/6d0)
               !!k_div_T = k_div_T*1.592d-9*pow(m/Msun, 2.84d0)/secyer
               E2 = 1.592d-9*pow(m/Msun, 2.84d0)/secyer
               k_div_T = k_div_T*E2
            end if
         end if

      end function k_div_T


      real(dp) function f2(e)
         real(dp), intent(in) :: e

         f2 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f2 after eq. 11
         if (e > 0d0) then
             f2 = 1d0+15d0/2d0*pow2(e) + 45d0/8d0*pow4(e) + 5d0/16d0*pow6(e)
         end if

      end function f2


      real(dp) function f3(e)
         real(dp), intent(in) :: e

         f3 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f3 after eq. 11
         if (e > 0d0) then
             f3 = 1d0+15d0/4d0*pow2(e) + 15d0/8d0*pow4(e) + 5d0/64d0*pow6(e)
         end if

      end function f3


      real(dp) function f4(e)
         real(dp), intent(in) :: e

         f4 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f4 after eq. 11
         if (e > 0d0) then
             f4 = 1d0+3d0/2d0*pow2(e) + 1d0/8d0*pow4(e)
         end if

      end function f4


      real(dp) function f5(e)
         real(dp), intent(in) :: e

         f5 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f5 after eq. 11
         if (e > 0d0) then
             f5 = 1d0+3d0*pow2(e) + 3d0/8d0*pow4(e)
         end if

      end function f5


      integer function how_many_extra_binary_history_header_items(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_header_items = 0
      end function how_many_extra_binary_history_header_items


      subroutine data_for_extra_binary_history_header_items( &
           binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id, n
         character (len = maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
      end subroutine data_for_extra_binary_history_header_items


      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id

         how_many_extra_binary_history_columns = 5

      end function how_many_extra_binary_history_columns


      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len = maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: mu, I1, I2

         real(dp) :: t_circ_1, t_circ_2

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         names(1) = 'inclination'
         names(2) = 'idot'

         vals(1) = b% xtra(x_inclination) * rad2a
         vals(2) = b% xtra(x_idot) * rad2a
         

         ! compute circ timescale based on different energy transport on the star envelope
         names(3) =  'lg_t_circ_1'
         names(4) =  'lg_t_circ_2'
        
         t_circ_1 = 0d0
         if (b% have_star_1) then
            t_circ_1 = get_t_circ(b, b% s1, b% circ_type_1, b% Ftid_1, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in get_t_circ for star 1'
               return
            end if
         end if
         vals(3) = safe_log10(abs(t_circ_1)/secyer)

         t_circ_2 = 0d0
         if (b% have_star_2) then
            t_circ_2 = get_t_circ(b, b% s2, b% circ_type_2, b% Ftid_2, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in get_t_circ for star 2'
               return
            end if
         end if
         vals(4) = safe_log10(abs(t_circ_2)/secyer)

         ! Darwin unstable separation
         names(5) = 'a_Darwin'

         mu = b% m(1) * b% m(2) / (b% m(1) + b% m(2))
         I1 = 0d0
         I2 = 0d0
         if (b% point_mass_i /= 1) &
            I1 = I1 + dot_product(b% s1% dm_bar(1:b% s1% nz), b% s1% i_rot(1:b% s1% nz))
         if (b% point_mass_i /= 2) &
            I2 = I2 + dot_product(b% s1% dm_bar(1:b% s1% nz), b% s1% i_rot(1:b% s1% nz))

         vals(5) = sqrt(3 * (I1 + I2) / mu) / Rsun

         
      end subroutine data_for_extra_binary_history_columns

      
      real(dp) function get_t_circ(b, s, circ_type, Ftid, ierr) result(t_circ)
         character (len=strlen) :: circ_type
         real(dp), intent(in) :: Ftid
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         type (star_info), pointer :: s

         real(dp) :: osep, qratio, m, r_phot

         ierr = 0

         osep = b% separation
         qratio = b% m(b% a_i) / b% m(b% d_i)
         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            qratio = 1.0d0/qratio
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if

         if (circ_type == 'Hut_conv') then
            t_circ = (21d0/2d0)*k_div_T(b, s, .true.)*qratio*(1+qratio)*pow8(r_phot/osep)
         else if (circ_type == 'Hut_rad') then
            t_circ = (21d0/2d0)*k_div_T(b, s, .false.)*qratio*(1+qratio)*pow8(r_phot/osep)
         else if (circ_type == 'Orb_period') then
            t_circ = b% period
         else
            ierr = -1
            write(*,*) 'unrecognized circ_type', circ_type
            return
         end if

         t_circ = (1d0/t_circ) / Ftid
         
      end function get_t_circ 
      
      
      integer function extras_binary_startup(binary_id, restart, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then  ! failure in  binary_ptr
            return
         end if
         
         ! for now, just use it as input
         initial_inclination = b% s1% x_ctrl(1)

         if (.not. restart) then
            b% lxtra(lx_convective_envelope) = .false.
            b% lxtra(lx_use_qin_E2) = b% s1% x_logical_ctrl(1)
            b% xtra(x_inclination) = initial_inclination*a2rad
            b% xtra_old(x_inclination) = initial_inclination*a2rad
            b% xtra(x_idot) = 0d0
            b% xtra(x_omega_eq) = 0d0
         end if

         extras_binary_startup = keep_going

      end function  extras_binary_startup
     

      integer function extras_binary_start_step(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         real(dp) :: i_step
         character (len = strlen) :: coplanar_type

         include 'formats.inc'

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then  ! failure in  binary_ptr
            return
         end if
         
         extras_binary_start_step = keep_going
        
         i_step = b% xtra(x_inclination)

         ! to update inclination based on Repetto & Nelemans (2014)
         ! first we compute di/dt
         idot_step = 0d0
         if (b% point_mass_i /= 1) then
            ! convective envelope check for star 1
            if (is_convective(b% s_donor)) then
               b% circ_type_1 = 'Hut_conv'
               b% sync_type_1 = 'Hut_conv'
               b% Ftid_1 = 50
            else
               b% circ_type_1 = 'Hut_rad'
               b% sync_type_1 = 'Hut_rad'
               b% Ftid_1 = 1
            end if

            if (b% circ_type_1 == 'Hut_conv') then
               coplanar_type = 'Hut_conv'
            else if (b% circ_type_1 == 'Hut_rad') then
               coplanar_type = 'Hut_rad'
            else
               write(*,*) 'unrecognized coplanar_type'
               ierr = -1
               return
            end if
            idot_step = idot_non_coplanar(b, b% s1, coplanar_type)
         end if
         if (b% point_mass_i /= 2) then
            ! convective envelope check for star 2
            if (is_convective(b% s_donor)) then
               b% circ_type_2 = 'Hut_conv'
               b% sync_type_2 = 'Hut_conv'
               b% Ftid_2 = 50
            else
               b% circ_type_2 = 'Hut_rad'
               b% sync_type_2 = 'Hut_rad'
               b% Ftid_2 = 1
            end if

            if (b% circ_type_2 == 'Hut_conv') then
               coplanar_type = 'Hut_conv'
            else if (b% circ_type_2 == 'Hut_rad') then
               coplanar_type = 'Hut_rad'
            else
               write(*,*) 'unrecognized coplanar_type'
               ierr = -1
               return
            end if
            idot_step = idot_step + idot_non_coplanar(b, b% s2, coplanar_type)
         end if
         
         ! compute di in radians. use to check for retry later on
         di_step = idot_step * b% time_step * secyer

         ! update inclination in the timestep. not saved in xtra array - yet
         i_step = i_step + di_step
         i_evolve =  i_step
         
      end function extras_binary_start_step
      

      !Return either keep_going, retry or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         real(dp) :: rel_i_change

         include 'formats.inc'

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then  ! failure in  binary_ptr
            return
         end if

         extras_binary_check_model = keep_going
        
         rel_i_change = abs(di_step / max(b% xtra(x_inclination), fi_limit))
         if (rel_i_change > fi_hard) then
            write(*,1) 'retry because of fi_hard limit', rel_i_change
            extras_binary_check_model = retry
            return
         end if

         ! if (abs(di_step) * rad2a > max_abs_di)  then
         !    write(*,1) 'retry due to di > max_abs_di', abs(di_step) * rad2a
         !    extras_binary_check_model = retry
         !    return
         ! end if
         
         ! after checking inclination value is OK, store it in xtra array
         b% xtra(x_idot) = idot_step
         
         ! be sure that we do not go beyond minimum inclination value
         if (i_evolve < min_inclination) then
            b% xtra(x_inclination) = min_inclination
         else
            b% xtra(x_inclination) = i_evolve
         end if
         
      end function extras_binary_check_model
      
      
      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         real(dp) :: omega_orb, a1, a2
         real(dp) :: mu, I1, I2, a_Darwin

         include 'formats.inc'

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then  ! failure in  binary_ptr
            return
         end if

         extras_binary_finish_step = keep_going

         ! compute the equilibrium value of spin for star
         omega_orb = 2d0*pi/b% period
         a1 = f2(b% eccentricity)
         a2 = pow(1-pow2(b% eccentricity), 1.5d0) * f5(b% eccentricity)
         b% xtra(x_omega_eq) = abs(a1/a2*omega_orb)

         write(*,'(a)')
         write(*,1) 'Ftid_1', b% Ftid_1
         write(*,1) 'idot', b% xtra(x_idot) * rad2a
         write(*,1) 'inclination', b% xtra(x_inclination) * rad2a
         write(*,1) 'edot_tidal', b% edot_tidal
         write(*,1) 'eccentricity', b% eccentricity
         write(*,1) 'Omega', b% s_donor% omega_avg_surf
         write(*,1) 'Omega_eq', b% xtra(x_omega_eq)
         write(*,1) 'Omega_div_Omega_eq', b% s_donor% omega_avg_surf/b% xtra(x_omega_eq)
         write(*,1) 'Jrot_div_Jorb', b% s_donor% total_angular_momentum/b% angular_momentum_j
         write(*,'(a)')

         ! check for Darwin unstable binary
         mu = b% m(1) * b% m(2) / (b% m(1) + b% m(2))
         I1 = 0d0
         I2 = 0d0
         if (b% point_mass_i /= 1) &
            I1 = I1 + dot_product(b% s1% dm_bar(1:b% s1% nz), b% s1% i_rot(1:b% s1% nz))
         if (b% point_mass_i /= 2) &
            I2 = I2 + dot_product(b% s1% dm_bar(1:b% s1% nz), b% s1% i_rot(1:b% s1% nz))
         a_Darwin = sqrt(3 * (I1 + I2) / mu)

         if (b% doing_first_model_of_run .and. b% separation < a_Darwin) then
            write(*,'(a)') 'system is Darwin unstable'
            b% s_donor% termination_code = t_xtra1
            termination_code_str(t_xtra1) = 'Darwin unstable'
            extras_binary_finish_step = terminate
            return
         end if

         ! terminate when equilibrium is reached
         if (.false. .and. abs(b% s_donor% omega_avg_surf/b% xtra(x_omega_eq) - 1) < 1d-2 .and. &
            abs(b% eccentricity - 0d0) < 1d-6 .and. abs(b% xtra(x_inclination) - 0d0) < 1d-6) then
            write(*,'(a)') 'reach equilibrium conditions'
            extras_binary_finish_step = terminate
            return
         end if

         ! end at beginning of MT
         if(.false. .and. b% rl_relative_gap(b% d_i) >= 0d0) then
            write(*,'(a)') 'reach MT phase'
            extras_binary_finish_step = terminate
            return
         end if

         ! stop at He depletion
         if (b% s1% center_h1 < 1d-5 .and. b% s1% center_he4 < 1d-5) then
            write(*,'(a)') 'reach he depletion'
            extras_binary_finish_step =  terminate
            return
         end if
         
      end function extras_binary_finish_step
     

      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then  ! failure in  binary_ptr
            return
         end if 
         
      end subroutine extras_binary_after_evolve     
      

      end module run_binary_extras
