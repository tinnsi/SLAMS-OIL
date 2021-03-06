SUBROUTINE interpolate(t1,t2,z1,z2,temp,z)
   
   DOUBLE PRECISION, INTENT(IN) :: t1, t2, z1, z2, z
   DOUBLE PRECISION, INTENT(OUT) :: temp
   DOUBLE PRECISION :: a, b, temp_0

   IF (t1 .eq. t2) THEN
      temp = t1
! 	  write(*,*) 'temp', temp
   ELSE
      a = (t1-t2)/(z1-z2)
      b = t2 - a*z2
      temp = a*z + b
! 	  write(*,*) 'temp', temp
	  temp_0 = (t2*(z-z1) + t1*(z2-z))/(z2-z1)
! 	  write(*,*) 'temp0', temp_0
   ENDIF

END SUBROUTINE interpolate
!============================================================
SUBROUTINE check_index(i, max_size, x)
   USE the_info
   USE  data_sim

! arguments to this subroutine are the index that is being used
! and the size of the dimension in use: SIZE(array, dim).  x is
! 0 if all is well and 1 if it is being overwritten.

   INTEGER, INTENT(IN) :: i, max_size
   INTEGER, INTENT(OUT) :: x

   IF (i .le. max_size) THEN
      x = 0
   ELSE
      x = 1
   ENDIF

END SUBROUTINE check_index
