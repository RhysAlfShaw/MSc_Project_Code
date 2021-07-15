module global_parameters

!
! Integration Scheme
!
integer, parameter :: n_order = 1

!
! Time-stepping accuracy parameters
!
double precision :: eta_dt = 0.01d0
double precision :: dt_stab_factor = 1.2d0

!
! Physical constants
!
double precision, parameter :: gg = 6.672d-08
double precision, parameter :: year = 3600.0d0*24.0d0*365.2425d0
double precision, parameter :: solar_mass = 1.98847d33

!
! Internal version of G (kept here for convenience)
!
double precision :: G_int

!
! numerical constants
!
double precision, parameter :: sixth = 1.0d0 / 6.0d0
double precision, parameter :: log2 = log(2.0d0)

end module global_parameters
