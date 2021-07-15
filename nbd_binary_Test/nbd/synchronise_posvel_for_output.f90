subroutine synchronise_posvel_for_output

use global_variables
use global_parameters

implicit none

double precision :: t_diff(1:N), t_diff2(1:N), t_diff3(1:N)
integer :: i

print *, 'Creating time-synchonous pos/vel for output!'

!
! set the various 'dt' arrays based on the snapshot output time
!

print *, 't', t
print *, 'time_at_next_snapshot', time_at_next_snapshot
t_diff = time_at_next_snapshot - t
t_diff2 = t_diff**2
t_diff3 = t_diff**3

!
! Sanity check: make sure there are no negative diffs! 
! no particle should ever get the chance to evolve beyond the snapshot
! time before this is called, but we should check! 
!
do i = 1, N
   if ( t_diff(i) < 0 ) then
      print *, 'Negative t-diff in snapshot sync.Must stop! i = ', i
      print *, 't_diff:', t_diff
      stop
   end if
end do

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
! get the accels based on the predicted values
!
call get_accel_info(N, xp, yp, zp, vxp, vyp, vzp, m, soft, ax1, ay1, az1, axdot1, aydot1, azdot1)

!
! now work out the 2nd and 3rd order "corrected" accels, which are based on a0 and a1:
! NOTE1: we don't divide by t here, to prevent rounding error in correction)
!        but rather do it at the end of the step, before working out new dts 
!
axc2 = -6.0d0*(ax0 - ax1) - t_diff*(4.0d0*axdot0 + 2.0d0*axdot1)
ayc2 = -6.0d0*(ay0 - ay1) - t_diff*(4.0d0*aydot0 + 2.0d0*aydot1)
azc2 = -6.0d0*(az0 - az1) - t_diff*(4.0d0*azdot0 + 2.0d0*azdot1)

axc3 = 12.0d0*(ax0 - ax1) + 6.0d0*t_diff*(axdot0 + axdot1)
ayc3 = 12.0d0*(ay0 - ay1) + 6.0d0*t_diff*(aydot0 + aydot1)
azc3 = 12.0d0*(az0 - az1) + 6.0d0*t_diff*(azdot0 + azdot1)

!
! finally, correct the predicted step using these corrected accels
!
xo = xp  + t_diff2*axc2/24.0d0 + t_diff2*axc3/120.0d0
yo = yp  + t_diff2*ayc2/24.0d0 + t_diff2*ayc3/120.0d0
zo = zp  + t_diff2*azc2/24.0d0 + t_diff2*azc3/120.0d0

vxo = vxp  + t_diff*axc2/6.0d0 + t_diff*axc3/24.0d0
vyo = vyp  + t_diff*ayc2/6.0d0 + t_diff*ayc3/24.0d0
vzo = vzp  + t_diff*azc2/6.0d0 + t_diff*azc3/24.0d0

end subroutine synchronise_posvel_for_output 
