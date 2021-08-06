program main
  implicit none

  real(8), dimension(:,:), allocatable :: U, Uold
  real(8) :: err, errmax, tolerance
  integer :: k, i, j, N

  N = 200
  allocate(U(N,N), Uold(N,N))
  tolerance = 1.e-4
  err = 2. * tolerance
  U = 0.
  k = 0

  U(N/2-10 : N/2+10 , N/2-10 : N/2+10) = 1.0


  do while (err > tolerance) 

    Uold = U
    errmax = 0.
    err = 0.
    
    do j = 2 , N-1
       do i = 2 , N-1

          if (i > N/2-10 .and. i < N/2+10 .and. j > N/2-10 .and. j < N/2+10) cycle

          U(i,j) = (Uold(i-1, j) + Uold(i, j-1) + Uold(i+1, j) + Uold(i, j+1))/4.

          err = abs(Uold(i,j) - U(i,j))

          if (err > errmax) then
             errmax = err

          end if

       end do
    end do

    err = errmax

    print*, k, err
    
    k = k + 1
    

  end do


  open(101, file='sol.dat')
  do j = 1, N
     do i = 1, N
        write(101, *) U(i,j)
     end do
     write(101,*)
  end do
  close(101)































end program main
