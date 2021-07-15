subroutine hermite_4thorder_PECn_global_dt

use global_variables
use global_parameters

implicit none

double precision :: dt_try(1:N)
integer :: n_iter
integer :: i

!
! get the accels at the current time
!
!print *, 'getting accels'
call get_accel_info(N, x, y, z, vx, vy, vz, m, soft, ax0, ay0, az0, axdot0, aydot0, azdot0)
if ( integration_loop_counter.eq. 0 ) then
   print *, 'Initiallising dts for global time-stepping algortithm' 
   do i = 1, N
      a = sqrt(ax0(i)*ax0(i) + ay0(i)*ay0(i) + az0(i)*az0(i))
      adot = sqrt(axdot0(i)*axdot0(i) + aydot0(i)*aydot0(i) + azdot0(i)*azdot0(i))
      dt_try(i) = eta_dt * a(i)/adot(i)
   end do
   dt = minval(dt_try)
   dt = min(dt, dt_max)
   if ( (global_time + dt) .gt. time_at_next_snapshot ) then
      dt = time_at_next_snapshot - global_time
      print *, "HERMITE: dt limited by snapshot time", global_time, dt
   end if   
end if

!
! using the current position and velocity, predict (p) new point
!
xp = x + dt*vx + 0.5d0*dt*dt*ax0 + sixth*dt*dt*dt*axdot0
yp = y + dt*vy + 0.5d0*dt*dt*ay0 + sixth*dt*dt*dt*aydot0
zp = z + dt*vz + 0.5d0*dt*dt*az0 + sixth*dt*dt*dt*azdot0

vxp = vx + dt*ax0 + 0.5d0*dt*dt*axdot0
vyp = vy + dt*ay0 + 0.5d0*dt*dt*aydot0
vzp = vz + dt*az0 + 0.5d0*dt*dt*azdot0

!
! apply the iterative corrector loop as discussed in Kokubo et al. 1998
!
do n_iter = 1, n_order
   !
   ! get the accels based on the predicted or corrected values:
   !
   if (n_iter .eq. 1 ) then
       call get_accel_info(N, xp, yp, zp, vxp, vyp, vzp, m, soft, ax1, ay1, az1, axdot1, aydot1, azdot1)
   else
      call get_accel_info(N, xc, yc, zc, vxc, vyc, vzc, m, soft, ax1, ay1, az1, axdot1, aydot1, azdot1)
   end if

   !
   ! now work out the 2nd and 3rd order "corrected" accels, which are based on a0 and a1:
   ! (note that we don't divide by t here, to prevent rounding error in correction)
   !§
   axc2 = -6.0d0*(ax0 - ax1) - dt*(4.0d0*axdot0 + 2.0d0*axdot1)
   ayc2 = -6.0d0*(ay0 - ay1) - dt*(4.0d0*aydot0 + 2.0d0*aydot1)
   azc2 = -6.0d0*(az0 - az1) - dt*(4.0d0*azdot0 + 2.0d0*azdot1)

   axc3 = 12.0d0*(ax0 - ax1) + 6.0d0*dt*(axdot0 + axdot1)
   ayc3 = 12.0d0*(ay0 - ay1) + 6.0d0*dt*(aydot0 + aydot1)
   azc3 = 12.0d0*(az0 - az1) + 6.0d0*dt*(azdot0 + azdot1)

   !
   ! finally, correct the predicted step using these corrected accels
   !
   xc = xp + (dt**2)*axc2/24.0d0 + (dt**2)*axc3/120.0d0
   yc = yp + (dt**2)*ayc2/24.0d0 + (dt**2)*ayc3/120.0d0
   zc = zp + (dt**2)*azc2/24.0d0 + (dt**2)*azc3/120.0d0

   vxc = vxp + dt*axc2/6.0d0 + dt*axc3/24.0d0
   vyc = vyp + dt*ayc2/6.0d0 + dt*ayc3/24.0d0
   vzc = vzp + dt*azc2/6.0d0 + dt*azc3/24.0d0


   !
   ! divide ac2 and ac3 by t for use in the timestep criterion
   ! 
   axc2 = axc2 / (dt**2)
   ayc2 = ayc2 / (dt**2)
   azc2 = azc2 / (dt**2)
   axc3 = axc3 / (dt**3)
   ayc3 = ayc3 / (dt**3)
   azc3 = azc3 / (dt**3)

end do

!
! update pos / vel with new iterative solution
!
x = xc
y = yc
z = zc
vx = vxc
vy = vyc
vz = vzc

!
! update global time
!
global_time = global_time + dt
t(:) = t(:) + dt

!
! get timesteps
!
dt_old = dt
do i = 1, N
   a1 = sqrt( ax1(i)*ax1(i) + ay1(i)*ay1(i) + az1(i)*az1(i) )
   a1dot = sqrt( axdot1(i)*axdot1(i) + aydot1(i)*aydot1(i) + azdot1(i)*azdot1(i) )
   ac3 = sqrt( axc3(i)*axc3(i) + ayc3(i)*ayc3(i) + azc3(i)*azc3(i) )
   acx2_ext = axc2(i) + dt_old*axc3(i)
   acy2_ext = ayc2(i) + dt_old*ayc3(i)
   acz2_ext = azc2(i) + dt_old*azc3(i)
   ac2_ext = sqrt(acx2_ext*acx2_ext + acy2_ext*acy2_ext + acz2_ext*acz2_ext)
   dt_try(i) = sqrt( eta_dt * (a1*ac2_ext + a1dot*a1dot) / (a1dot*ac3 + ac2_ext*ac2_ext) )
end do
dt = minval(dt_try)
! apply the Aarseth stability criterion
dt = min(dt, dt_old*dt_stab_factor)
! prevent going past output time...
if ( (global_time + dt) .gt. time_at_next_snapshot ) then
   dt = time_at_next_snapshot - global_time
   print *, "HERMITE: dt limited by snapshot time", global_time, dt
end if

end subroutine hermite_4thorder_PECn_global_dt
