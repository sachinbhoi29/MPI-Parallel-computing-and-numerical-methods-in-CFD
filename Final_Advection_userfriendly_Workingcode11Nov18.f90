
 program advection

 implicit none
 integer :: xmax,xmin,nmax, i, j, k, itr, method, fun             
 real :: dx, dt, cfl, length, c, t_current, t_final, start, finish         
 real, dimension(:), allocatable :: xcord, f0, f1, Uan, Unew, Ucur,Ucurmac, L_1, L_infinity   ! Unew (n+1), Ucur (n) !double precision
 print*, "Enter number of grid points: "
 read(*,*) nmax
 print*, "Enter CFL number: "
 read(*,*) cfl
 print*, "Enter final time: "
 read(*,*) t_final
 xmin=-40
 xmax=40 
 length=xmax-xmin
 c=1.5
 dx=length/(nmax-1)       			      ! 80/100
 t_current = 0.0 
 print*, "1-FTCS, 2-Upwind, 3-Lax-Friedrichs, 4-Lax-Wendroff, 5-MacCormack," 
 print*, "Enter the method: "
 read(*,*) method    ! 1 - FTCS, 2 - Upwind, 3 - Lax-Friedrichs, 4 - Lax-Wendroff, 5 - MacCormack
 print*, "1 - sign, 2 - exponential"
 print*, "Enter the function: "
 read(*,*) fun


!allocate or declare the size of the array
allocate (xcord(nmax-1))
allocate (f0(nmax-1))
allocate (f1(nmax-1))
allocate (Unew(nmax-1))
allocate (Ucur(nmax-1))
allocate (Uan (nmax-1))
allocate (L_1 (nmax-1))
allocate (L_infinity (nmax-1))
allocate (Ucurmac (nmax-1))

!run the loop

call cpu_time(start)

! Spatial discretisation (grid generation) ! it is an array of 100 numbers 
do i=0 , nmax-1
  xcord(i)=xmin+i*dx

!	print*, 'xcord',xcord
  !write(*,*) xcord	
end do

! Set intial conditions or initial condition array 
if (fun.eq.1) then
	do i=0,nmax-1
		f0(i)=0.5*(sign(1.0,xcord(i))+1.0)
	end do
else if (fun.eq.2) then 
	do i=0,nmax-1
		f1(i)=0.5*exp(-(xcord(i))**2.0)
	end do
else 
	write(*,*) "Invalid function choice."
	stop
end if												! 



! Comnpute time-step
dt = dx*cfl/c

! Initialise numerical array

if (fun.eq.1) then
		Ucur = f0
else if (fun.eq.2) then 
		Ucur = f1
else 
	print*, "Invalid choice"
end if	

																	

! Set boundary conditions
if (fun.eq.1) then
		Ucur(0)      = 0.0
		Ucur(nmax-1) = 1.0
else if (fun.eq.2) then 
		Ucur(0)   =0.0
		Ucur(nmax-1) = 0.0
else 
	print*, "Invalid choice"
end if	


! Begin time-loop
do while (t_current .le. t_final)
	
	if (fun.eq.1) then
		Ucur(0)      = 0.0
		Ucur(nmax-1) = 1.0
	else if (fun.eq.2) then 
		Ucur(0)   =0.0
		Ucur(nmax-1) = 0.0
	else 
	print*, "Invalid choice"
	end if
	! Spatial loop
	do i=1,nmax-2
		if (method .eq. 1) then
			! FTCS
			Unew(i) = Ucur(i) - (c*dt/(2.0*dx))*(Ucur(i+1) - Ucur(i-1))
			!print*, 'UNEW: ', Unew(i)
		else if (method .eq. 2) then
			! Upwind scheme
			Unew(i) = Ucur(i) - (c*dt/dx)*(Ucur(i) - Ucur(i-1))
			!print*, 'UNEW: ', Unew(i)
		else if (method .eq. 3) then
			! Lax-Fredrichs
			Unew(i) = 0.5*(Ucur(i+1)+Ucur(i-1))- (c*dt/(2*dx))*(Ucur(i+1) - Ucur(i-1))
			!print*, 'UNEW: ', Unew(i)
		else if (method .eq. 4) then
			! Lax-Wendroff
			Unew(i) = Ucur(i)- c*dt*((Ucur(i+1)-Ucur(i-1))/(2*dx)) + 0.5*c*c*dt*dt*((Ucur(i+1)-2*Ucur(i)+Ucur(i-1))/dx*dx)
			!print*, 'UNEW: ', Unew(i)
		else if (method .eq. 5) then
			! MacCormack
			Ucurmac(i)=Ucur(i)-(c*dt/dx)*(Ucur(i+1)-Ucur(i))
			Unew(i) = 0.5*((Ucur(i)+Ucurmac(i))-(c*dt/dx)*(Ucurmac(i)-Ucurmac(i-1)))
			!print*, 'UNEW: ', Unew(i)
		else 
			write(*,*) "Invalid method choice."
			stop
		end if
	end do
	
	! Update numerical values
	do i=1,nmax-2
		Ucur(i) = Unew(i)
	end do

	! Re-enforce boundary conditions
!	Ucur(0)      = 0.0														!Apparently not needed to update here
!	Ucur(nmax-1) = 1.0
!	Ucur(0)      = 0.0														!Re-enforce set boundary condotions for the exponential function
!       Ucur(nmax-1) = 0.0
															
	t_current = t_current + dt	! update time for each successive loop
end do

do i=0,nmax-1

		if (fun.eq.1) then
		Uan(i) = 0.5*(sign(1.0,xcord(i)-1.5*t_current)+1.0)
		else if (fun.eq.2) then 
		Uan(i) = 0.5*exp(-(xcord(i)-1.5*t_final)**2)
		else 
		print*, "Invalid choice"
		end if	

end do


!Computing the errora by L_1 and L_infinity norms
!do i=1, nmax-1
!	L_1(i)=(Uan(i)-Ucur(i))
	
!print*,'L_1: ', L_1(i)
!end do

call cpu_time(finish)

!printing results
do i=0,nmax-1
   write(*,*), i,xcord(i),Ucur(i),Uan(i)
end do

print*, "        Time to run the code is ", finish-start, " seconds."	

open(1,file = 'OutputProfileNS1D.dat')

	do i=1,Nmax-1
	  write (1,*) i,xcord(i),Ucur(i),Uan(i)
 	end do
	close (1)
 
end program advection
