subroutine get_ics

use global_variables
use global_parameters

!
! Open the file and read in the data. The snapshots use the same format as the ics,
! ics, excepts that they hold additional information in a footer
!
open(20, file = initial_conditions, form='unformatted')

!
! First read in the number of bodies and other header info
!
read(20) N
read(20) global_time
read(20) units_mass_g, units_length_cm, units_time_s

!
! now allocate all the arrays for the program
!
call allocate_memory

!
! read in the arrays
!
read(20) x(1:N)
read(20) y(1:N)
read(20) z(1:N)
read(20) vx(1:N)
read(20) vy(1:N)
read(20) vz(1:N)
read(20) m(1:N)

!
! close the file
!
close(20)

print *, 'pos 1', x(1), y(1), z(1)
print *, 'pos 2', x(2), y(2), z(2)
print *, 'vel 1', vx(1), vy(1), vz(1)
print *, 'vel 2', vx(2), vy(2), vz(2)

end subroutine get_ics
