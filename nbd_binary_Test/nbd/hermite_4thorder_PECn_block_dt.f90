subroutine hermite_4thorder_PECn_block_dt

use global_variables
use global_parameters

implicit none

double precision :: t_diff(1:N), t_diff2(1:N), t_diff3(1:N)
integer :: n_iter
integer :: i, nn
integer :: dt_bin_try

!
! set the new global_time (i.e. the time at the END of the upcoming step) 
! and various 'dt' arrays
!

global_time = minval(t + delta_t)
t_diff = global_time - t
t_diff2 = t_diff**2
t_diff3 = t_diff**3

!
! Now get a list of active particles (i.e. those that have t = global_time)
!
N_act = 0
do i = 1, N
   if ( abs(t(i) + delta_t(i) - global_time) < 2*tiny(global_time) ) then
      N_act = N_act + 1
      i_active(N_act) = i
   end if
end do
if( N_act .eq. 0) then
   print *, 'PROBLEM: no active particles! Aborting! :-/'
   print *, 'time loop counter', integration_loop_counter
   stop
end if

!
! using the current position and velocity, predict (p) new point
!
xp = x + t_diff*vx + 0.5d0*t_diff2*ax0 + sixth*t_diff3*axdot0
yp = y + t_diff*vy + 0.5d0*t_diff2*ay0 + sixth*t_diff3*aydot0
zp = z + t_diff*vz + 0.5d0*t_diff2*az0 + sixth*t_diff3*azdot0

vxp = vx + t_diff*ax0 + 0.5d0*t_diff2*axdot0
vyp = vy + t_diff*ay0 + 0.5d0*t_diff2*aydot0
vzp = vz + t_diff*az0 + 0.5d0*t_diff2*azdot0

!
! apply the iterative corrector loop as discussed in Kokubo et al. 1998.
! NOTE: this is only really worth it when using a constant global timestep
! but we've kept it here for testing purposes
!
do n_iter = 1, n_order
   !
   ! get the accels based on the predicted or corrected values:
   !
   if (n_iter .eq. 1 ) then
      call get_accel_for_active(N, N_act, i_active, xp, yp, zp, vxp, vyp, vzp, m, soft, ax1, ay1, az1, axdot1, aydot1, azdot1)
   else
      call get_accel_for_active(N, N_act, i_active, xc, yc, zc, vxc, vyc, vzc, m, soft, ax1, ay1, az1, axdot1, aydot1, azdot1)
   end if

   !
   ! now work out the 2nd and 3rd order "corrected" accels, which are based on a0 and a1:
   ! NOTE1: we don't divide by t here, to prevent rounding error in correction)
   !        but rather do it at the end of the step, before working out new dts 
   !
   do nn = 1, N_act
      i = i_active(nn)
      axc2(i) = -6.0d0*(ax0(i) - ax1(i)) - t_diff(i)*(4.0d0*axdot0(i) + 2.0d0*axdot1(i))
      ayc2(i) = -6.0d0*(ay0(i) - ay1(i)) - t_diff(i)*(4.0d0*aydot0(i) + 2.0d0*aydot1(i))
      azc2(i) = -6.0d0*(az0(i) - az1(i)) - t_diff(i)*(4.0d0*azdot0(i) + 2.0d0*azdot1(i))

      axc3(i) = 12.0d0*(ax0(i) - ax1(i)) + 6.0d0*t_diff(i)*(axdot0(i) + axdot1(i))
      ayc3(i) = 12.0d0*(ay0(i) - ay1(i)) + 6.0d0*t_diff(i)*(aydot0(i) + aydot1(i))
      azc3(i) = 12.0d0*(az0(i) - az1(i)) + 6.0d0*t_diff(i)*(azdot0(i) + azdot1(i))
   end do

   !
   ! finally, correct the predicted step using these corrected accels
   !
   xc = xp
   yc = yp
   zc = zp
   vxc = vxp
   vyc = vyp
   vzc = vzp
   do nn = 1, N_act
      i = i_active(nn)
      xc(i) = xp(i)  + t_diff2(i)*axc2(i)/24.0d0 + t_diff2(i)*axc3(i)/120.0d0
      yc(i) = yp(i)  + t_diff2(i)*ayc2(i)/24.0d0 + t_diff2(i)*ayc3(i)/120.0d0
      zc(i) = zp(i)  + t_diff2(i)*azc2(i)/24.0d0 + t_diff2(i)*azc3(i)/120.0d0

      vxc(i) = vxp(i)  + t_diff(i)*axc2(i)/6.0d0 + t_diff(i)*axc3(i)/24.0d0
      vyc(i) = vyp(i)  + t_diff(i)*ayc2(i)/6.0d0 + t_diff(i)*ayc3(i)/24.0d0
      vzc(i) = vzp(i)  + t_diff(i)*azc2(i)/6.0d0 + t_diff(i)*azc3(i)/24.0d0
   end do

   !
   ! divide ac2 and ac3 by t for use in the timestep criterion
   ! 
   axc2 = axc2 / t_diff2
   ayc2 = ayc2 / t_diff2
   azc2 = azc2 / t_diff2
   axc3 = axc3 / t_diff3
   ayc3 = ayc3 / t_diff3
   azc3 = azc3 / t_diff3
end do

!
! update pos / vel for *active* particles with new iterative solution and set the 'predicted' values
! to these too, for use in the accel function
!
do nn = 1, N_act
   i = i_active(nn)
   x(i) = xc(i)
   y(i) = yc(i)
   z(i) = zc(i)
   vx(i) = vxc(i)
   vy(i) = vyc(i)
   vz(i) = vzc(i)
   xp(i) = xc(i)
   yp(i) = yc(i)
   zp(i) = zc(i)
   vxp(i) = vxc(i)
   vyp(i) = vyc(i)
   vzp(i) = vzc(i)
end do

!
! get new accel info at this position / velocity for the active particles
! Note we use the pos / vel values stored in the 'predicted' arrays, since
! these are 'drifted' in space between active times.
!

call get_accel_for_active(N, N_act, i_active, xp, yp, zp, vxp, vyp, vzp, m, soft, ax0, ay0, az0, axdot0, aydot0, azdot0)

!
! update time and  get new timesteps for ACTIVE particles. Note that we
! don't prevent particles evolving past the output time here. Instead, 
! we recalculate the positions for output separately.
! 
do nn = 1, N_act
   i = i_active(nn)
   !
   ! set time of particle to global_time + current dt
   !
   t(i) = t(i) + delta_t(i)

   ! 
   ! flip the isok_change_dt flag between pos/neg
   !
   isok_change_dt(i) = isok_change_dt(i) * (-1) 

   !
   ! use the 'mystical' (according to W. Denon) Aarseth criterion
   ! to set a new dt, if allowed! 
   !
   if ( isok_change_dt(i).gt.0 ) then
      a1 = sqrt( ax1(i)*ax1(i) + ay1(i)*ay1(i) + az1(i)*az1(i) )
      a1dot = sqrt( axdot1(i)*axdot1(i) + aydot1(i)*aydot1(i) + azdot1(i)*azdot1(i) ) 
      ac3 = sqrt( axc3(i)*axc3(i) + ayc3(i)*ayc3(i) + azc3(i)*azc3(i) )
      acx2_ext = axc2(i) + t_diff(i)*axc3(i)
      acy2_ext = ayc2(i) + t_diff(i)*ayc3(i)
      acz2_ext = azc2(i) + t_diff(i)*azc3(i)
      ac2_ext = sqrt(acx2_ext*acx2_ext + acy2_ext*acy2_ext + acz2_ext*acz2_ext)
      dt = min(dt_max, sqrt( eta_dt * (a1*ac2_ext + a1dot*a1dot) / (a1dot*ac3 + ac2_ext*ac2_ext) ))
      !
      ! Now convert this into a power of 2...
      !
      dt_bin_try = int(log2_dtmax - log(dt)/log2) + 1 
      if (dt_bin_try .gt. t_power_bin(i)) then
         t_power_bin(i) = t_power_bin(i) + 1
      else if ( dt_bin_try .gt. t_power_bin(i) ) then
         t_power_bin(i) = max(0, t_power_bin(i) - 1)
      end if
      
      delta_t(i) = dt_max / 2**t_power_bin(i)
   end if
end do


end subroutine hermite_4thorder_PECn_block_dt
