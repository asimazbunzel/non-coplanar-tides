! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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
      use chem_def
      use num_lib
      use binary_def
      
      implicit none

      ! these routines are called by the standard run_star check_model
      contains
      
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% other_wind => brott_wind

      end subroutine extras_controls

      
      subroutine brott_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr

         integer :: h1, he4
         real(dp) :: Xs, Ys, Z_div_Z_solar, Teff_jump, alfa, L1, M1, R1, T1, &
            vink_wind, nieu_wind, hamann_wind
         type (star_info), pointer :: s
         type (binary_info), pointer :: b
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) return

         L1 = Lsurf
         M1 = Msurf
         R1 = Rsurf
         T1 = Tsurf

         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         Xs = s% xa(h1,1)
         Ys = s% xa(he4,1)
         ! Z=0.017 is Z from Grevesse et al. 1996
         Z_div_Z_solar = s% kap_rq% Zbase / 0.017d0
         ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
         Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z_div_Z_solar)))

         vink_wind = 0d0
         nieu_wind = 0d0
         hamann_wind = 0d0
         w = 0

         call eval_Vink_wind(vink_wind)
         call eval_Nieuwenhuijzen_wind(nieu_wind)
         call eval_Hamann_wind(hamann_wind)

         ! use 1/10 hamann
         hamann_wind = hamann_wind/10d0

         if (T1 < Teff_jump) then
            ! low T wind
            w = max(vink_wind, nieu_wind)
         else
            ! high T wind
            alfa = 0d0
            if (Xs > 0.7d0) then
               alfa = 1d0
            else if (Xs > 0.4d0 .and. Xs < 0.7d0) then
               alfa = (Xs - 0.4d0)/0.3d0
            end if
            w = alfa * vink_wind + (1d0-alfa) * hamann_wind
         end if

         ierr = 0

         contains

         subroutine eval_Vink_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc

            ! alfa = 1 for hot side, = 0 for cool side
            if (T1 > 27500d0) then
               alfa = 1
            else if (T1 < 22500d0) then
               alfa = 0
            else
               dT = 100d0
               if (T1 > Teff_jump + dT) then
                  alfa = 1
               else if (T1 < Teff_jump - dT) then
                  alfa = 0
               else
                  alfa = (T1 - (Teff_jump - dT)) / (2*dT)
               end if
            end if
            
            if (alfa > 0) then ! eval hot side wind (eqn 24)
               vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10(L1/Lsun/1d5) &
                  - 1.313d0*log10(M1/Msun/30) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.933d0*log10(T1/4d4) &
                  - 10.92d0*pow2(log10(T1/4d4)) &
                  + 0.85d0*log10(Z_div_Z_solar)
               w1 = exp10(logMdot)
            else
               w1 = 0
            end if
            
            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.688d0 &
                  + 2.210d0*log10(L1/Lsun/1d5) &
                  - 1.339d0*log10(M1/Msun/30) &
                  - 1.601d0*log10(vinf_div_vesc/2d0) &
                  + 1.07d0*log10(T1/2d4) &
                  + 0.85d0*log10(Z_div_Z_solar)
               w2 = exp10(logMdot)
            else
               w2 = 0
            end if
            
            w = alfa*w1 + (1 - alfa)*w2
            
         end subroutine eval_Vink_wind

         subroutine eval_Nieuwenhuijzen_wind(w)
            ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -14.02d0 &
                     +1.24d0*log10(L1/Lsun) &
                     +0.16d0*log10(M1/Msun) &
                     +0.81d0*log10(R1/Rsun) &
                     +0.85d0*log10(Z_div_Z_solar)
            w = exp10(log10w)
         end subroutine eval_Nieuwenhuijzen_wind

         subroutine eval_Hamann_wind(w)
            ! Hamann, W.-R.; Koesterke, L.; Wessolowski, U. 1995, A&A, 299, 151
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -11.95d0 &
                     +1.5d0*log10(L1/Lsun) &
                     -2.85d0*Xs &
                     + 0.85d0*log10(Z_div_Z_solar)
            w = exp10(log10w)
         end subroutine eval_Hamann_wind

      end subroutine brott_wind
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         extras_start_step = 0

         if (s% doing_relax) then
            s% job% pgstar_flag = .false.
         else
            s% job% pgstar_flag = .true.
         end if

      end function extras_start_step


      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         include 'formats'

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         extras_check_model = keep_going

      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id

         how_many_extra_history_columns = 15

      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: dt
         type(star_info), pointer :: s
         type(binary_info), pointer :: b
         real(dp) :: m_conv, Ebind
         integer :: k

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) return

         names(1) = 'Mconv_env'
         names(2) = 'fconv_env'

         if (s% he_core_mass > 0d0 .and. s% center_h1 < 1d-6) then
            m_conv = 0d0
            if (s% he_core_k == 1 .or. (s% star_mass - s% he_core_mass) < 1d-3) then
               vals(1) = 0d0
               vals(2) = 0d0
            else
               do k = 1, s% he_core_k
                  if (s% mixing_type(k) == convective_mixing) m_conv = m_conv + s% dm(k)
               end do
               vals(1) = m_conv / Msun
               vals(2) = (m_conv / Msun) / (s% star_mass - s% he_core_mass)
            end if
         else
            vals(1) = 0d0
            vals(2) = 0d0
         end if

         names(3) = 'Ebind'

         Ebind = 0d0
         do k=1, s% nz
            Ebind = Ebind + s% dm(k) * (-standard_cgrav * s% m(k) / s% r(k) + s% energy(k))
         end do
         vals(3) = Ebind

         names(4) = 'spin_magnitude'
         vals(4) = clight * s% total_angular_momentum / (standard_cgrav * s% m(1)**2)

         names(5) = 'log_mtransfer_rate'
         vals(5) = safe_log10(abs(b% step_mtransfer_rate)/Msun*secyer)

         names(6) = 'eff_xfer_fraction'
         if (b% component_mdot(b% d_i) == 0d0) then
            vals(6) = 1d0
         else
            vals(6) = (-b% component_mdot(b% a_i))/(b% component_mdot(b% d_i))
         end if

         names(7) = 'inclination'
         vals(7) = b% xtra(1) * rad2a

         names(8) = 'eccentricity'
         vals(8) = b% eccentricity

         names(9) = 'Prot'
         vals(9) = (2d0 * pi / s% omega_avg_surf) / (60d0 * 60d0 * 24d0)

         names(10) = 'Porb'
         vals(10) = b% period / (60d0 * 60d0 * 24d0)

         names(11) = 'Prot_div_Porb'
         vals(11) = vals(9) / vals(10)

         names(12) = 'Omega'
         vals(12) = s% omega_avg_surf

         names(13) = 'Omega_div_Omega_crit'
         vals(13) = s% w_div_w_crit_avg_surf

         names(14) = 'Omega_eq'
         vals(14) = b% xtra(6)

         names(15) = 'Omega_div_Omega_eq'
         vals(15) = s% omega_avg_surf / b% xtra(6)

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id

         how_many_extra_profile_columns = 1

      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         type(star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = 'Ebind'
         vals(:,1) = 0d0
         vals(1,1) = s% dm(1) * (-standard_cgrav * s% m(1) / s% r(1) + s% energy(1))
         do k=2, s% nz
            vals(k,1) = vals(k-1,1) + s% dm(k) * (-standard_cgrav * s% m(k) / s% r(k) + s% energy(k))
         end do

      end subroutine data_for_extra_profile_columns
      

      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: center_h1, center_he4
         real(dp) :: center_h1_old, center_he4_old

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         extras_finish_step = keep_going
         
         ! saved profiles at specific stages
         center_h1 = s% xa(s% net_iso(ih1), s% nz)
         center_h1_old = s% xa_old(s% net_iso(ih1), s% nz_old)
         center_he4 = s% xa(s% net_iso(ihe4), s% nz)
         center_he4_old = s% xa_old(s% net_iso(ihe4), s% nz_old)
         if (center_h1 < 0.5 .and. center_h1_old > 0.5) then
            call star_write_profile_info(id, 'LOGS/prof_h50.data', ierr)
            if (ierr /= 0) return
         else if (center_h1 < 1d-6 .and. center_h1_old > 1d-6) then
            call star_write_profile_info(id, 'LOGS/prof_h00.data', ierr)
            if (ierr /= 0) return
         else if (center_h1 < 1d-6) then
            if (center_he4 < 0.5 .and. center_he4_old > 0.5) then
               call star_write_profile_info(id, 'LOGS/prof_he50.data', ierr)
               if (ierr /= 0) return
            else if (center_he4 < 1d-5 .and. center_he4_old > 1d-5) then
               call star_write_profile_info(id, 'LOGS/prof_he00.data', ierr)
               if (ierr /= 0) return
            end if
         end if

         ! stop at C depletion
         if (s% center_c12 < 1d-6 .and. s% center_h1 < 1d-6) then
            call star_write_profile_info(id, 'LOGS/prof_final.data', ierr)
            if (ierr /= 0) return
            s% termination_code = t_xtra1
            termination_code_str(t_xtra1) = 'reach C depletion'
            extras_finish_step = terminate
            return
         end if

      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

      end subroutine extras_after_evolve

      end module run_star_extras
      
