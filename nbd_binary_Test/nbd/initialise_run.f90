subroutine initialise_run

!
! pcc 18/10/2020
!
! sets lots of parameters before we start
!

use global_variables
use global_parameters

implicit none

integer :: i, nn
integer :: ic_filename_length
character*4 :: filename_tail
integer :: t_power
double precision :: dt_initial(1:N)

!
! set the snapshot time stamps and filecounter
!
time_at_last_snapshot = global_time
time_at_next_snapshot = global_time + time_between_snapshots
if (restart_flag .eq. 1) then
   !
   ! restarted from a snapshot, so need to get the filestamp from the tail
   ! of the ic filename
   !
   ic_filename_length = len_trim( initial_conditions )
   filename_tail = initial_conditions(ic_filename_length-4:ic_filename_length)
   read(filename_tail, *) snap_file_counter
   print *, "Started from a snapshot ending in number ",  snap_file_counter
else
   !
   ! started from a true IC file, so can start the snapshot numbering from 0
   !
   snap_file_counter = 0
end if

!
! Important one... Set G! :) 
!
G_int = gg * units_mass_g * units_time_s**2 / units_length_cm**3

!
! set the time arrays
!
t(:) = global_time

!
! get the hermite block time steps scheme up and running, if it's being used.
!
if (int_mode .eq. 2) then

   !
   ! set all particles to be active XXX: this could be removed!
   !
   N_act = N
   do i = 1, N
      i_active(i) = i
      isok_change_dt(i) = 1
   end do

   !
   ! get the accels at the current time
   !
   call get_accel_for_active(N, N_act, i_active, x, y, z, vx, vy, vz, m, soft, ax0, ay0, az0, axdot0, aydot0, azdot0)
   a = sqrt(ax0*ax0 + ay0*ay0 + az0*az0)
   adot = sqrt(axdot0*axdot0 + aydot0*aydot0 + azdot0*azdot0)
   dt_initial = min(dt_max, eta_dt * a/adot)

   !
   ! prevent going past output time...
   !
   do i = 1, N 
      if ( (global_time + dt_initial(i)) .gt. time_at_next_snapshot ) then
         dt_initial(i) = time_at_next_snapshot - global_time
      end if
   end do

   !
   ! define log2_dtmax
   !
   log2_dtmax = log(dt_max)/log2

   !
   ! Now convert all dts into a power of 2... 
   !
   t_power_bin = int(log2_dtmax - log(dt_initial)/log2) + 1
   delta_t = dt_max / 2**t_power_bin
   
   ! make the time between snapshots an integer multiple of dt_max
   ! XXX put a print statement here to notify user and check that this works!  
   !t_power = int(log2_dtmax - log(time_between_snapshots)/log2)
   !print *, "INITIALIZE: modifying snap cadence from ", time_between_snapshots
   !print *, 'to ', dt_max / 2**t_power
   !time_between_snapshots = dt_max / 2**t_power
end if

end subroutine initialise_run
