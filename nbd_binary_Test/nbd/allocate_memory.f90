subroutine allocate_memory

use global_variables
use global_parameters

!
! this is called after N is read, to set all the arrays in the code
!

! the basic pos, vel, mass
allocate ( x(1:N) )
allocate ( y(1:N) )
allocate ( z(1:N) )
allocate ( vx(1:N) )
allocate ( vy(1:N) )
allocate ( vz(1:N) )
allocate ( m(1:N) )
allocate ( t(1:N) )
if (int_mode .eq. 2 ) then
   !
   ! extra arrays with block timesteps.
   !
   allocate ( delta_t(1:N) )
   allocate ( t_power_bin(1:N) )
   allocate ( i_active(1:N) )
   allocate ( isok_change_dt(1:N) )
end if

! the predicted and corrected pos / vel
allocate ( xp(1:N) )
allocate ( yp(1:N) )
allocate ( zp(1:N) )
allocate ( vxp(1:N) )
allocate ( vyp(1:N) )
allocate ( vzp(1:N) )
allocate ( xc(1:N) )
allocate ( yc(1:N) )
allocate ( zc(1:N) )
allocate ( vxc(1:N) )
allocate ( vyc(1:N) )
allocate ( vzc(1:N) )


! the acceleration arrays
allocate( ax0(1:N) )
allocate( ay0(1:N) )
allocate( az0(1:N) )
allocate( axdot0(1:N) )
allocate( aydot0(1:N) )
allocate( azdot0(1:N) )
allocate( ax1(1:N) )
allocate( ay1(1:N) )
allocate( az1(1:N) )
allocate( axdot1(1:N) )
allocate( aydot1(1:N) )
allocate( azdot1(1:N) )

allocate( a(1:N) )
allocate( adot(1:N) )

allocate( axc2(1:N) )
allocate( ayc2(1:N) )
allocate( azc2(1:N) )
allocate( axc3(1:N) )
allocate( ayc3(1:N) )
allocate( azc3(1:N) )

! the snapshot arrays
allocate ( xo(1:N) )
allocate ( yo(1:N) )
allocate ( zo(1:N) )
allocate ( vxo(1:N) )
allocate ( vyo(1:N) )
allocate ( vzo(1:N) )


end subroutine allocate_memory
