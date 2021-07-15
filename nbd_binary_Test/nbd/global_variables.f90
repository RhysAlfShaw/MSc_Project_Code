module global_variables

!
! body data
!
integer :: N, N_act
double precision, allocatable :: x(:), y(:), z(:)
double precision, allocatable :: vx(:), vy(:), vz(:)
double precision, allocatable :: m(:)
double precision, allocatable :: t(:), delta_t(:)
integer, allocatable :: t_power_bin(:)
double precision :: soft ! only global softening for just now... need to change this
integer, allocatable :: i_active(:)
integer, allocatable :: isok_change_dt(:)
!
! temporary pos vel, acc, arrays used by hermite scheme
!
double precision, allocatable :: xp(:), yp(:), zp(:), vxp(:), vyp(:), vzp(:)
double precision, allocatable :: xc(:), yc(:), zc(:), vxc(:), vyc(:), vzc(:)
double precision, allocatable :: ax0(:), ay0(:), az0(:), axdot0(:), aydot0(:), azdot0(:)
double precision, allocatable :: ax1(:), ay1(:), az1(:), axdot1(:), aydot1(:), azdot1(:)
double precision, allocatable :: axc2(:), ayc2(:), azc2(:)
double precision, allocatable :: axc3(:), ayc3(:), azc3(:)
double precision, allocatable :: a(:), adot(:)

!
! which integration mode?
!
integer :: int_mode ! 0 = fixed dt; 1 = global dt; 2 = block dt

!
! accel info for timestep control
!
double precision :: acx2_ext, acy2_ext, acz2_ext, ac2_ext
double precision :: a1, a1dot, ac3

!
! time / loop stuff
!
double precision :: global_time
double precision :: dt, dt_old
double precision :: end_time
integer :: integration_loop_counter
double precision :: dt_max, log2_dtmax

!
! filenames
!
character*100 :: initial_conditions
character*100 :: snapshot_base, snapshot_filename


!
! snapshots
!
double precision :: time_at_last_snapshot
double precision :: time_at_next_snapshot
double precision :: time_between_snapshots
double precision, allocatable :: xo(:), yo(:), zo(:), vxo(:), vyo(:), vzo(:)
integer :: restart_flag
integer :: snap_file_counter

!
! diagnostic control (controls how frequently the code spews infomation
! via the standard output).
!
integer :: diagnostic_freq  ! how many passes through main loop before info
integer :: diagnostic_counter
double precision :: mean_dt, min_dt, max_dt
double precision :: xcom, ycom, zcom, vxcom, vycom, vzcom, total_mass
double precision :: spec_ang_mom_tot
double precision :: e_kin, e_grav

!
! units
!
double precision :: units_mass_g, units_length_cm, units_time_s

end module global_variables
