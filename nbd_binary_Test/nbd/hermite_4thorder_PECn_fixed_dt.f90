subroutine hermite_4thorder_PECn_fixed_dt

use global_variables
use global_parameters

implicit none

double precision :: dt_try
integer :: n_iter

!
! get the accels at the current time
!
call get_accel_info(N, x, y, z, vx, vy, vz, m, soft, ax0, ay0, az0, axdot0, aydot0, azdot0)

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

end subroutine hermite_4thorder_PECn_fixed_dt
