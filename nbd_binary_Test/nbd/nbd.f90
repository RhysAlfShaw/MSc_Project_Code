program nbd

!
! pcc - 18/10/2020
!
! N-Body Dynamics (NBD) code with at least 4th order Hermite.
!
! Takes an initial condition file and parameter file as input
!

use global_parameters
use global_variables

implicit none

!
! load parameters file
!
call read_param_file

!
! load the initial conditions from file
!
call get_ics

!
! initialise the run
!
call initialise_run

!
! initial diagnostics
!
call diagnostics(0)

!
! The main integration loop
!
integration_loop_counter = 0
do while (global_time .lt. end_time)
   !print *, 'in loop ', integration_loop_counter, global_time, dt 
   !
   ! Is it time for a snapshot? Do energy checks here
   !
   call write_a_snapshot

   !
   ! Do the integration step here. This subroutine takes care of timesteps and
   ! performs calls to get_accel_info.
   !
   if ( int_mode .eq. 1 ) then
      call hermite_4thorder_PECn_global_dt
   else if ( int_mode .eq. 2 ) then
      call hermite_4thorder_PECn_block_dt
   else
      call hermite_4thorder_PECn_fixed_dt
   end if

   !
   ! diagnostics
   !
   call diagnostics(1)

   !
   ! update the integration loop counter
   !
   integration_loop_counter = integration_loop_counter + 1
   !if ( integration_loop_counter.eq. 3) stop
end do

end program nbd
