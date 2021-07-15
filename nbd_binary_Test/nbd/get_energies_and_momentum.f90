subroutine get_energies_and_momentum

!
! pcc - 18/10/2020
!

use global_variables
use global_parameters

implicit none

integer :: i, j
double precision :: spec_ang_mom_x(1:N), spec_ang_mom_y(1:N), spec_ang_mom_z(1:N)

!
! get COM / VCOM
!
total_mass = sum(m)
xcom = sum(x*m) / total_mass
ycom = sum(y*m) / total_mass
zcom = sum(z*m) / total_mass
vxcom = sum(vx*m) / total_mass
vycom = sum(vy*m) / total_mass
vzcom = sum(vz*m) / total_mass

!
! get angular momentum
!
spec_ang_mom_x = (y - ycom)*(vz - vzcom) - (z - zcom)*(vy - vycom)
spec_ang_mom_y = - (x - xcom)*(vz - vzcom) + (z - zcom)*(vx - vxcom)
spec_ang_mom_z = (x - xcom)*(vy - vycom) - (y - ycom)*(vx - vxcom)
spec_ang_mom_tot = sqrt((sum(spec_ang_mom_x))**2 + (sum(spec_ang_mom_y))**2 + (sum(spec_ang_mom_z))**2)

!
! get kinetic energy
!
e_kin = 0.5*sum(((vx-vxcom)**2 + (vy-vycom)**2 + (vz-vzcom)**2)*m )

!
! get grav energy
!
e_grav = 0
do i = 1, N-1
   do j = i+1, N
      e_grav = e_grav - G_int*m(i)*m(j) / sqrt( (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2 )
   end do
end do



end subroutine get_energies_and_momentum
