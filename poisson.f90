program poisson
  use precision
  implicit none

  real(dp), dimension(:,:), allocatable :: U, Uold
  real(dp) :: tolerance, a, err, maxerr, w
  integer :: i,j,k, N, error
  character(1) :: M
  character(20) :: arg

  if (iargc()<3) then
     write(*,*) 'poisson Method N tolerance w'
     stop
  endif
  
  call getarg(1,arg)
  read(arg,*) M

  call getarg(2,arg)
  read(arg,*) N

  call getarg(3,arg)
  read(arg,*) tolerance

  call getarg(4,arg)
  read(arg,*) w

  allocate(U(N,N), stat=error)
  allocate(Uold(N,N), stat=error)
  if (error /= 0) then
     write(*,*) 'allocatio error'
     stop
  endif

  U=0.0_dp
  
  U(N/2-10:N/2+10, N/2-10:N/2+10) = 1.0_dp
  
  select case (M)
  case("J")
    write(*,*) "Jacobi iteration"
 
    err = 2.0_dp * tolerance
    k = 1
    do while (err > tolerance)
      Uold = U 
      maxerr=0.0_dp        !lo devo riniziaizzare perche' altrimenti mi porto lo strascico dell'errore vecchio (che e' sempre piu' grosso)
 
      do j = 2, N-1         !j parte da 2 perche' sul bordo resta fissa al valore 0
        do i = 2, N-1
         
          if (i >= N/2-10 .and. i<=N/2+10 .and. j>=N/2-10 .and. j<=N/2+10) cycle  ! cycle vuol dire salta quell'iterazione del ciclo
          
          U(i,j) = (Uold(i-1,j) + Uold(i+1,j) + Uold(i,j-1) + Uold(i,j+1))/4.0_dp  !c'e' U old, e' un'approssimazione
          
          if (abs(Uold(i,j)-U(i,j)) > maxerr ) then
             maxerr = abs(Uold(i,j) - U(i,j))        !e' la massima differenza tra gli elementi delle 2 matrici, ha senso perche' voglio che converga componente per componente
          end if
     
        end do
      end do
 
      write(*,*) 'iter: ',k, maxerr
      err = maxerr
      k = k + 1
    end do   
 
  case("G") 
    write(*,*) "Gauss-Siedel iteration"
 
    err = 2.0_dp * tolerance
    k = 1
    do while (err > tolerance)
      Uold = U 
      maxerr=0.0_dp
 
      do j = 2, N-1
        do i = 2, N-1
         
          if (i >= N/2-10 .and. i<=N/2+10 .and. j>=N/2-10 .and. j<=N/2+10) cycle  
          
          U(i,j) = (U(i-1,j) + Uold(i+1,j) + U(i,j-1) + Uold(i,j+1))/4.0_dp
          
          if (abs(Uold(i,j)-U(i,j)) > maxerr ) then
             maxerr = abs(Uold(i,j) - U(i,j))
          end if
          
          U(i,j) = (1-w)*Uold(i,j) + w*U(i,j)
     
        end do
      end do
 
      write(*,*) 'iter: ',k, maxerr
      err = maxerr
      k = k + 1
    end do  

  case("N")
  
    U=0.0_dp
    U(N/2-10:N/2+10,1) = -1.0_dp
    U(N/2-10:N/2+10,N) = +1.0_dp

    err = 2.0_dp * tolerance
    k = 1
    do while (err > tolerance)

      Uold = U 
      maxerr=0.0_dp
 
      do j = 2, N-1
        do i = 2, N-1
         
          U(i,j) = (U(i-1,j) + Uold(i+1,j) + U(i,j-1) + Uold(i,j+1))/4.0_dp
          
          if (abs(Uold(i,j)-U(i,j)) > maxerr ) then
             maxerr = abs(Uold(i,j) - U(i,j))
          end if
          
          U(i,j) = (1-w)*Uold(i,j) + w*U(i,j)
     
        end do
        !Neumann di O(a^2) lungo x==1 e x==N
        U(1,j) = 4.0_dp/3.0_dp * U(2,j) - 1.0_dp/3.0_dp * U(3,j)
        U(N,j) = 4.0_dp/3.0_dp * U(N-1,j) - 1.0_dp/3.0_dp * U(N-2,j) 
        U(1,1) = 4.0_dp/3.0_dp * U(2,1) - 1.0_dp/3.0_dp * U(3,1)
        U(N,N) = 4.0_dp/3.0_dp * U(N-1,N) - 1.0_dp/3.0_dp * U(N-2,N) 
      end do
      !Neumann di O(a^2) lungo y==1 e y==N
      U(1:N/2-11,1) = 4.0_dp/3.0_dp * U(1:N/2-11,2) - 1.0_dp/3.0_dp * U(1:N/2-11,3)
      U(N/2+11:N,1) = 4.0_dp/3.0_dp * U(N/2+11:N,2) - 1.0_dp/3.0_dp * U(N/2+11:N,3)
      U(1:N/2-11,N) = 4.0_dp/3.0_dp * U(1:N/2-11,N-1) - 1.0_dp/3.0_dp * U(1:N/2-11,N-2)
      U(N/2+11:N,N) = 4.0_dp/3.0_dp * U(N/2+11:N,N-1) - 1.0_dp/3.0_dp * U(N/2+11:N,N-2)
 
      write(*,*) 'iter: ',k, maxerr
      err = maxerr
      k = k + 1
    end do  

 end select


 open(101, file='sol.dat')
 do j = 1, N
    do i = 1, N
       write(101, *) U(i,j)
    end do
    write(101,*)
 end do
 close(101)
 
  

end program poisson
