subroutine diagnostics(diagnostic_mode)

!
! pcc - 19/10/2020
!
! This code prints out things that the user might like to know, such as information 
! on timesteps, energy, momentum, (V)COM, masses. We've got two modes:
! 	diagnostic_mode = 0 --> initial information from startup / ICs
!	diagnostic_mode = 0 --> information from inside the time integration loop
!

use global_variables
use global_parameters

implicit none

integer :: diagnostic_mode

call get_energies_and_momentum

if ( diagnostic_mode.eq.0 ) then
   !
   ! set the counter to zero
   !
   diagnostic_counter = 0

   !
   ! print out the initial information. Can add to this as more
   ! useful quantities become apparent!
   !
   print *, '------------------------------------------- '
   print *, ' >>>>> the N Body Dynamics (NBD) code <<<<<' 
   print *, '------------------------------------------- '
   print *, 'Using PEC-n scheme with n = ', n_order
   if ( int_mode .eq. 1 ) then
      print *, 'Using global (but variable) timesteps'
   else if (int_mode .eq. 2 ) then
      print *, 'Using block (i.e. variable and individual) timesteps' 
   else
      print *, 'Using a constant timestep for all particles'
   end if 
   print *, ' '
   print *, '>>> Information on INITIAL conditions <<<'
   print *, 'Number of bodies (N): ', N
   print *, 'mass units (g):      ', units_mass_g
   print *, 'distance units (cm): ', units_length_cm
   print *, 'time units (s)     : ', units_time_s
   print *, 'Gravtiational constant: ', G_int
   print *, 'centre of mass:      ', xcom, ycom, zcom
   print *, 'centre of velocity   ', vxcom, vycom, vzcom
   print *, 'total spec angular momentum: ', spec_ang_mom_tot
   print *, 'gravitational and kinetic energies: ', e_grav, e_kin
   print *, 'total energy:                       ', e_grav + e_kin
else if (diagnostic_mode.eq.1) then
   !
   ! update the counter 
   !   
   diagnostic_counter = diagnostic_counter + 1
   !
   ! is it time to write out information?
   !
   if (diagnostic_counter.lt.diagnostic_freq) then
      return
   else
      diagnostic_counter = 0
   end if

   !
   ! print out the diagnostics from the loop
   !
   print *, ' '
   print *, 'INFO AT STEP ', integration_loop_counter
   print *, 'time and dt', global_time, dt
   print *, 'centre of mass:      ', xcom, ycom, zcom
   print *, 'centre of velocity   ', vxcom, vycom, vzcom
   print *, 'total spec angular momentum: ', spec_ang_mom_tot
   print *, 'gravitational and kinetic energies: ', e_grav, e_kin
   print *, 'total energy:                       ', e_grav + e_kin      
else
   !
   ! doesn't exist!
   !
   print *, "wrong diagnostic mode!"
   stop
end if

end subroutine diagnostics
