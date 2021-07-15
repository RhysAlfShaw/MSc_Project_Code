subroutine write_a_snapshot
!
! pcc = 18/10/2010
!
! Writes the snapshots in the same format as the ic file, but can have extra
! information stored in the header (means that snapshots can be used as ics).
! Files are written in standard fortran binary.
!

use global_variables
use global_parameters

implicit none
character*4 :: filename_tail
integer :: snapfile_len
integer :: i
double precision :: dx, dy, dz
double precision :: dvx, dvy, dvz
double precision :: v_rel2, r_rel, rdotv
double precision :: e_bin, a_bin, p_bin, energy_bin
double precision :: m_sys, m_red
double precision :: t_predict

!
! find the time that our next integration step will take us to 
! 
if (int_mode .eq. 2) then
   t_predict = minval(t + delta_t)
else
   ! minval is just to get a variable, rather than an array
   ! they should all be the same!
   t_predict = minval(t + dt)
end if

!
! decide whether or not we need to write a snapshot
!
if ( t_predict .ge. time_at_next_snapshot ) then
   print *, "WRITESNAP, T-PREDICT:", t_predict, time_at_next_snapshot
   !
   ! set the filename and open
   !
1000 format(I4.4)
   write(filename_tail, 1000) snap_file_counter
   snapfile_len = len_trim(snapshot_base)
   snapshot_filename = snapshot_base(1:snapfile_len) // filename_tail
   open(30, file=snapshot_filename, form='unformatted')

   !
   ! We need to synchronise all the particles at the output time
   ! Note that this is done into new variables, xo, vxo, etc, so the orginal
   ! pos/vels are not disturbed.
   !
   print *, 'WRITE x before sync', x
   call synchronise_posvel_for_output
   print *, 'WRITE x after sync', x
   print *, 'WRITE xo', xo
   print *, 'Should be using xo, vxo, etc for output'

   !
   ! notify user
   !
   print *, 'WRITING SNAPSHOT to file ', snapshot_filename, ' at time ', time_at_next_snapshot

   !
   ! write the data and close. Note that we use the *predicted* value here. In the
   ! block time-step code, we need to synchonise the particles for output. The
   ! predicted value arrays already store this, so we make use of them here.
   ! when block timesteps are not in use, we first pass the standard x, vx, etc arrays
   ! to the predicted value arrays.
   !
   ! now write...
   write(30) N
   write(30) time_at_next_snapshot
   write(30) units_mass_g, units_length_cm, units_time_s
   write(30) xo(1:N)
   write(30) yo(1:N)
   write(30) zo(1:N)
   write(30) vxo(1:N)
   write(30) vyo(1:N)
   write(30) vzo(1:N)
   write(30) m(1:N) 
   close(30)   

   !
   ! write an ascii too if N is small
   !
   if (N.gt.2  .and.  N.lt.30) then
      if ( snap_file_counter.eq.0 ) then
         open(40, file='body_evolution.dat')
      end if
      write(40, *) N
      write(40, *) global_time
      do i = 1, N
         write(40, 4000) xo(i),yo(i),zo(i),vxo(i),vyo(i),vzo(i),m(i)
      end do
   end if
   if (N.le.10) then
      !
      ! write out pos, vel, m
      !
      if ( snap_file_counter.eq.0 ) then
         open(50, file='orbit_evolution.dat')
         open(60, file='binary_parameters.dat')
      end if
      write(50, 4000) global_time,xo(2),yo(2),xo(3),yo(3),xo(4),yo(4),xo(5),yo(5),xo(10),yo(10)
      !
      ! binary properties to a separate file!
      !
      dx = xo(2) - xo(1)
      dy = yo(2) - yo(1)
      dz = zo(2) - zo(1)
      dvx = vxo(2) - vxo(1)
      dvy = vyo(2) - vyo(1)
      dvz = vzo(2) - vzo(1)
      r_rel = sqrt(dx*dx + dy*dy + dz*dz)
      v_rel2 = dvx*dvx + dvy*dvy + dvz*dvz
      m_sys = m(1) + m(2)
      m_red = m(1) * m(2) / m_sys
      energy_bin = 0.5d0 * m_red * v_rel2 - G_int * m(1) * m(2) / r_rel
      a_bin = - G_int * m(1) * m(2) / 2.0d0 / energy_bin
      p_bin = sqrt( a_bin**3 / m_sys )
      rdotv = dx*dvx + dy*dvy + dz*dvz
      e_bin = sqrt((1.0d0 - r_rel/a_bin)**2 + (rdotv**2)/a_bin/G_int/m_sys)
      write(60, 4000) global_time, dt, energy_bin, a_bin, p_bin, e_bin, r_rel, sqrt(v_rel2)
   end if
 4000 format(16(1X,E16.10))

   !
   ! update the time stamps and snapcounter
   !
   time_at_last_snapshot = time_at_next_snapshot
   time_at_next_snapshot = time_at_last_snapshot + time_between_snapshots
   print *, 'Next snapshot will be written at time ', time_at_next_snapshot

   !
   ! update the snapshot file counter
   !
   snap_file_counter = snap_file_counter + 1
end if

end subroutine write_a_snapshot
