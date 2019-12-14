                                                       !Parallel code using HPC

PROGRAM Parallel_code_advection 

! Send & Receive Data using MPI

IMPLICIT NONE
include 'mpif.h'

! Declarations
integer myid, numprocs, i, j, k, xmax, xmin, nmax,  numproc,Snmax,First_cell,remainder,SSnmax, fun, method
real :: dx, dt, cfl, length, c, t_current, t_final, start, finish         
real, dimension(:), allocatable :: xcord, f0, f1, Uan, Unew, Ucur,L_1, L_infinity, SUcur,SUnew,SUan, finaldata, Sxcord,SSxcord,SUcurmac    ! Unew (n+1), Ucur (n) !double precision

! Declarations and initialisation
real :: Leftboundary 
real :: Rightboundary   
integer :: IERROR
integer :: istatus(mpi_status_size) 

! Initialise MPI libraries
call MPI_INIT(IERROR)
call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, IERROR)
call MPI_COMM_RANK( MPI_COMM_WORLD, myid, IERROR)  

 nmax=101 				         !Grid size 
 length=80.0
 c=1.5
 dx=length/(nmax-1)     		    ! 80/100
 cfl=0.5                                ! numerical VISCOSITY
 t_current = 0.0 
 t_final = 5.0                  !question 2 output time equal 5 and 10, try 
! print*,"numprocs",numprocs,nmax-1,nmax/numprocs 
 xmin=-40
 xmax=40
 !print*, "Enter the function: "
 !read(*,*) fun
 fun =2                   ! define the function   1 sign 2 exponential
!print*, "1-FTCS, 2-Upwind, 3-Lax-Friedrichs, 4-Lax-Wendroff, 5-MacCormack," 
! print*, "Enter the method: "
! read(*,*) 
method = 1   ! 1 - FTCS, 2 - Upwind, 3 - Lax-Friedrichs, 4 - Lax-Wendroff, 5 - MacCormack

!allocate or declare the size of the array
allocate(FinalData(numprocs))
allocate (xcord(nmax))
allocate (f0(nmax))
allocate (f1(nmax))
allocate (Unew(nmax))
allocate (Ucur(nmax))
allocate (Uan (nmax))

		dx=length/(nmax-1)     		        ! 80/100
		Snmax=(nmax/numprocs) 			!smallcontainersize
		remainder = mod(nmax, numprocs)
		SSnmax= Snmax+remainder
		allocate (SSxcord(SSnmax))
		allocate (SUcur(Snmax))
		allocate (SUnew(Snmax))
		allocate (SUan(Snmax))
		allocate (Sxcord(Snmax))
		allocate (SUcurmac(Snmax))









!if (remainder.eq.0) then
   if(myid.eq.0) then

	!print*,snmax

		do i=1, SSnmax-1
        		Sxcord(i)=xmin+(dx*i)-dx
		!	print*,i,myid,xcord(i)
 		end do
!print*, numprocs-1
		do i=1,numprocs-1
			!print*, i*Snmax+1
			First_cell=xmin+i*Snmax*dx
			!print*, i, myid, 'check'
			call MPI_SEND(First_cell, 1, MPI_REAL, i, i, MPI_COMM_WORLD, IERROR)
		end do
   end if
!end if

call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	do i=1, numprocs-1	
		if (myid .eq. i) then
!print*, 'hello', myid
			call MPI_RECV(First_cell, 1,MPI_REAL,0, i, MPI_COMM_WORLD,istatus, IERROR)
			print*,First_cell
			
			do j=1,Snmax
				Sxcord(j)=First_cell+dx*j
			end do

			do j=1, Snmax
!				print*,"myid ", myid,"Sxcord", Sxcord(j)
			end do
		end if   
	end do



!--------------------------------------calculating the function-----------------------------------------------------------------------
if (fun.eq.1) then

  do i=0, numprocs-1
	if (myid.eq.i) then
		do j=1, Snmax
			f0(j)=0.5*(sign(1.0,Sxcord(j))+1.0)
			!print*,"myid ", myid,"f0", f0(j)
		end do
	end if
  end do

end if



if (fun.eq.2) then


do i=0, numprocs-1
	if (myid.eq.i) then
		do j=1, Snmax
			f1(j)=0.5*exp(-(Sxcord(j))**2.0)
		!	print*,"myid ", myid,"f1", f0(j)
		end do
	end if
end do

end if

!--------------------------------------------------------Copying in Ucurrent------------------------------------------------------------------

if (fun.eq.1) then

do i=0, numprocs-1
	if (myid.eq.i) then
		do j=1, Snmax
			Ucur(j)=f0(j)
			!print*,"myid ", myid,"f0", f0(j)
		end do
	end if
end do
	
end if



if (fun.eq.1) then

do i=0, numprocs-1
	if (myid.eq.i) then
		do j=1, Snmax
			Ucur(j)=f1(j)
			!print*,"myid ", myid,"f0", f0(j)
		end do
	end if
end do
	
end if

!--------------------------------------------------------------------------------------------------------------------------------------------

! Set boundary conditions
if (fun.eq.1) then



	do i=0, 0
		if (myid.eq.i) then
				Ucur(1)=0.0	
		end if
	end do


	do i=numprocs-1,numprocs-1 
		if (myid.eq.i) then
				Ucur(Snmax)=1.0	
		end if
	end do	



else if (fun.eq.2) then


	do i=0, 0
		if (myid.eq.i) then
				Ucur(1)=0.0	
		end if
	end do


	do i=numprocs-1,numprocs-1 
		if (myid.eq.i) then
			Ucur(Snmax)=0.0	
		end if
	end do
end if

!print*,"Check"
!---------------------------------------------------------Main function caluclation-----------------------------------------------------------



! Begin time-loop
do i=0,numprocs-1
	  do while (t_current .le. t_final)
	
		! Set boundary conditions
		if (fun.eq.1) then

!print*,i,j

			do k=0, 0
				if (myid.eq.i) then
					Ucur(1)=0.0	
				end if
			end do


			do k=numprocs-1,numprocs-1 
				if (myid.eq.i) then
					Ucur(Snmax)=1.0	
				end if
			end do	


		else if (fun.eq.2) then


			do k=0, 0
				if (myid.eq.i) then
					Ucur(1)=0.0	
				end if
			end do


			do k=numprocs-1,numprocs-1 
				if (myid.eq.i) then
					Ucur(Snmax)=0.0	
				end if
			end do
		end if
print*,"Check"
	! Spatial loop
	do j=2,Snmax-2
		if (method .eq. 1) then
			! FTCS
			SUnew(j) = SUcur(j) - (c*dt/(2.0*dx))*(SUcur(j+1) - SUcur(j-1))
			!print*, 'UNEW: ', Unew(i)
		else if (method .eq. 2) then
			! Upwind scheme
			SUnew(j) = SUcur(j) - (c*dt/dx)*(SUcur(j) - SUcur(j-1))
			!print*, 'UNEW: ', Unew(i)
		else if (method .eq. 3) then
			! Lax-Fredrichs
			SUnew(j) = 0.5*(SUcur(j+1)+SUcur(j-1))- (c*dt/(2*dx))*(SUcur(j+1) - SUcur(j-1))
			!print*, 'UNEW: ', Unew(i)
		else if (method .eq. 4) then
			! Lax-Wendroff
			SUnew(j) = SUcur(j)- c*dt*((SUcur(j+1)-SUcur(j-1))/(2*dx)) + 0.5*c*c*dt*dt*((SUcur(j+1)-2*SUcur(j)+SUcur(j-1))/dx*dx)
			!print*, 'UNEW: ', Unew(i)
		else if (method .eq. 5) then
			! MacCormack
			SUcurmac(j)=SUcur(j)-(c*dt/dx)*(SUcur(j+1)-SUcur(j))
			SUnew(j) = 0.5*((SUcur(j)+SUcurmac(j))-(c*dt/dx)*(SUcurmac(j)-SUcurmac(j-1)))
			!print*, 'UNEW: ', Unew(j)
		else 
			write(*,*) "Invalid method choice."
			stop
		end if
	end do
	
	! Update numerical values
	do j=1,Snmax
		SUcur(i) = SUnew(i)
	end do

	! Re-enforce boundary conditions
!	Ucur(0)      = 0.0														!Apparently not needed to update here
!	Ucur(nmax-1) = 1.0
!	Ucur(0)      = 0.0														!Re-enforce set boundary condotions for the exponential function
!       Ucur(nmax-1) = 0.0
															
	t_current = t_current + dt	! update time for each successive loop
    end do
   do j=2,Snmax-1

		if (fun.eq.1) then
		SUan(j) = 0.5*(sign(1.0,xcord(j)-1.5*t_current)+1.0)
		else if (fun.eq.2) then 
		SUan(j) = 0.5*exp(-(xcord(j)-1.5*t_final)**2)
		else 
		print*, "Invalid choice"
		end if	

   end do
end do

do i=0, numprocs-1
	do j=1, Snmax
		print*, SUcur(j)
	end do
end do
	


!--------------------------------------------------------------------------------------------------------------------------------------

!if (remainder.ne.0) 	then		    !special small xcord to accommodate more nodes  
	
!	if (myid .eq. 0) then
!			do j=1,Snmax
!				Sxcord(j)=Xmin+dx*j
!			end do
!		do i=1, numprocs-1
		!do j=1, Snmax
			!print*, i*Snmax+1
!			First_cell=i*Snmax+1
!			call MPI_SEND(First_cell, 1, MPI_REAL, i, i*456, MPI_COMM_WORLD, IERROR)
			!end do
!		end do
 !  	end if
!call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

 !   do i=1, numproc-2
!	if (myid .eq.i) then
!			call MPI_RECV(First_cell, 1,MPI_REAL,0, i*456, MPI_COMM_WORLD, IERROR)
!			do j=1,Snmax
!				Sxcord(j)=First_cell+dx*j
!			end do
		!do j=1, Snmax
			!print*,"myid ", myid,"Sxcord", Sxcord(j)
		!end do
			
!	end if	
 !   end do
  !      if (myid.eq.numproc-1) then
!		call MPI_RECV(First_cell, 1,MPI_REAL,0, i*456, MPI_COMM_WORLD, IERROR)
!			do j=1,Snmax
!				Sxcord(j)=First_cell+dx*j
!			end do
!	end if
!end if
		    

 	


!do j=1, Snmax
!			print*,"myid ", myid,"Sxcord", Sxcord(j)
!end do	
	
!print*, "check"

! Set intial conditions or initial condition array 
!if (myid.eq.0) then
!  do i=1,nmax
!	f0(i)=0.5*(sign(1.0,xcord(i))+1.0)
	!f1(i)=0.5*(exp(-xSUcurcord(i)**2))						! set initial boundary conditions for the exponential function
!   write(*,*), i,xcord(i),f0(i)
!  end do
!do i=1,nmax
!	print*,i,myid,xcord(i), f0(i)
!end do
!end if
!print*,"Snmax", Snmax

!call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!till here program is working

!Initial condition, sending the data of f0 to small Ucur, SUcur
!do i=0,numprocs-1
!	if (myid.eq.i) then
		!do j=1, Snmax			
			!SUcur(j)=f0(i*Snmax+j)							!Initial Condition for sign
			!SUcur(j)=f1(i*Snmax+j)							!Initial condition for exp
			!call MPI_BARRIER(MPI_COMM_WORLD,IERROR)print*,i,j
			!print*,i,j,myid,xcord(i),SUcur(j)
			
!			do j=1, Snmax
!				SUcur(j) = f0(i*Snmax+j)
!				print*,"myid",myid,SUcur(j)
!			end do
		!end do
		
!	end if
	
!end do
		 





!boundary condition initilisation  

!do i=0, Numprocs-1

	! Boundary conditions at the first and last cell 
!	if (myid.eq.0) then
!		SUnew(1)=0.0									!BC for sign
		!SUnew(Snmax)=0.0								!BC for exp
!	else if (myid.eq.numprocs-1) then
!		SUnew(nmax)=1.0									!BC for sign
		!SUnew(Snmax)=0.0								!BC for exp
!	end if

!end do






! Seinding the Ucurrent boundary data by mpi to leftboundary and rightboundary for the calculation of the Unew on the boundary 

!do i=0, Numprocs-2 !numprocs-2 because number of processors start from 0, so -1, and last processor won't send data so one more -1 =>-1-1=>-2
! send the data from Ucur array, from previous processor to the rightbob, last processor won't send the data cuz no next processor 
!	if (myid.eq.i) then     
			 			!     for 100 grid cells and 4 procs
!	        		call MPI_SEND(SUcur(Snmax-1), 1, MPI_REAL, i+1, i+1, MPI_COMM_WORLD, IERROR)	!sending the last cell data from previous proc for the calculation of next proc first cell
									!send tonext processor is i+1
!	end if
!end do

!do i=1, Numprocs-1 !from second processor to last processor, send the first cell value
!	if(myid.eq.i) then					!   
!				call MPI_SEND(SUcur(1), 1, MPI_REAL, i-1, 200+i, MPI_COMM_WORLD, IERROR)
								  !send to Previous processor i-1    		
!	end if
!end do
	
!call MPI_BARRIER(MPI_COMM_WORLD,IERROR)    ! let every processor send the data, if you put recv the data without putting barrier, cuz of dif proc speedit might happen that proc will try to recv beforing sending

!call MPI_BARRIER(MPI_COMM_WORLD,IERROR)





! receiving the data sent by the previous processors   LFboundary============================
!Leftboundary
!do i=1, numprocs-1 ! recieving the data sent by the previous processor so it'll start from the second processor (i=1 and not i=0) till the end processor
!	if (myid.eq.i) then
!		call MPI_RECV(Leftboundary, 1, MPI_INT, i, i+1, MPI_COMM_WORLD, istatus, IERROR)
!	end if		!adrs to store data, no of data, integer,source proc, tag
!end do

! Receiving the data sent next processor to the previous processor
!Rightboundary
!do i=0, numprocs-2 ! receicall MPI_BARRIER(MPI_COMM_WORLD,IERROR)ve the data sent by the next processor and storing it as rightboundary of the previous processor      =============================RTboundary
!	if (myid.eq.i)	then
!		call MPI_RECV(Rightboundary, 1, MPI_INT, i+1, 20000+i, MPI_COMM_WORLD, istatus, IERROR)
!	end if
!end do







! we have the data for the boubdaries, now we need to calculate the SUnew on the boundaries only  LFboundary +==================================+  RTboundary

!using rghtboundary to calculate the array last cell, processor will start from 0,1,2
!do i=0, numprocs-2
!	if (myid.eq.i)	then
		!time loop
!		do while (t_current.le.t_final)
			!spatial loop
!			do j=1, Snmax-1    ! j=1 to snmax-1 because it is already defined on bc defined on the final boundary condition
!				SUnew((j+1)*Snmax) = SUcur(j) - (c*dt/(2.0*dx))*(Rightboundary - SUcur(j-1))	
!			end do
!		end do
!	end if
!end do
		
!Using right boundary to calculate the array first cell, processor will start from 1,2,3
!do i=1, numprocs-1
!	if (myid.eq.i)	 then
		!time loop
!		do while (t_current.le.t_final)
			!scall MPI_BARRIER(MPI_COMM_WORLD,IERROR)patial loop
!			do j=1, Snmax-1    ! j=1 to snmax-1 because it is already defined on bc defined on the final boundary condition
!				SUnew((j+1)*Snmax) = SUcur(j) - (c*dt/(2.0*dx))*(SUcur(j+1) - Leftboundary)	
!			end do
!		end do
!	end if
!end do







! calculate all the data between the boundaries
!do i=1, numprocs-1
 !do the calculation for the unew except on the boundaries
!	if (myid.eq.i) then
	!time loop 
!		do while (t_current.le.t_final)
		!spatial loop
!			do j=1, Snmax-1    ! j=1 to snmax-1 because it is already defined on bc defined on the final boundary condition
!				SUnew(j) = SUcur(j) - (c*dt/(2.0*dx))*(SUcur(j+1) - SUcur(j-1))	
!			end do
!		end do
!	end if
!end do







! Blocking future tasks to ensure the previous task is complete.
! This is good practice, especially in parallel programming.
! Data cannot be received before it is sent!
!call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

!call MPI_ALLREDUCE(NewData, global_sum, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, IERROR)

!call MPI_ALLGATHER(SUnew, Snmax, MPI_REAL, UCur, nmax-1, MPI_REAL, MPI_COMM_WORLD, IERROR)

 
!printing results
!do i=1,nmax-1
!   print*, i,xcord(i),Ucur(i),Uan(i)
!end do
call MPI_FINALIZE(IERROR)

END PROGRAM Parallel_code_advection

   

