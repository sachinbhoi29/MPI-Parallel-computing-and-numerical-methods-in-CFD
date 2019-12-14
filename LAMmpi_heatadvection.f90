program mpi_heatadvection
    implicit none
    include 'mpif.h'

    integer myid, numprocs
    integer :: IERROR
    integer :: istatus(mpi_status_size)

    ! declarations of variable
    real*8 :: u, leftdomain, rightdomain ! u: speed; leftdomain & rightdomain: left & right end of domain
    real*8 :: leftbound, rightbound ! boundary condition at left & right end of domain
    integer :: i, j, k ! counters
    integer :: quotient, remainder ! quotient & remainder, used for allocate indexs to processors
    integer :: start_index, numnodes ! start index and number of nodes for each process
    integer :: problemset, solver, node ! problemset: number of problem set; solver: number of solver; node: number of node in domain
    real*8 :: t, cfl, dx, dt, timer ! t: target time; cfl: courant number; dt: time step; dx: spacing between node; timer: simulation time
    real*8, dimension(:), allocatable :: x, f0, f1, fs ! x: x coordinate; f0, f1 & fs: numerical steps
    real*8 :: leftvalue, rightvalue ! leftvalue & rightvalue: values at ends of domain in each processors
    real*8 :: lefthalo, righthalo ! lefthalo & righthalo: values from adjacent left & right processors
    real :: t_start, t_end, t_comm, t_comp ! t_start & t_end: real time process starts & ends; t_comm & t_comp: time for communication & computation
    character(len=100) :: filename, timelog, cproblemset, csolver, cnode, ct, ccfl ! for writing file

    ! solver
    ! 1: FTCS
    ! 2: Upwind
    ! 3: Lax-Friedrichs
    ! 4: Lax-Wendroff
    ! 5: MacCormack

    problemset = 0
    solver = 5
    node = 10000
    t = 10
    cfl = 1

    ! static parameters
    u = 1.5
    leftdomain = -40
    rightdomain = 40

    ! allocate problem set
    if (problemset .eq. 0) then
        leftbound = 0
        rightbound = 1
    else if (problemset .eq. 1) then
        leftbound = 0
        rightbound = 0
    end if

    ! calculate basic parameters
    dx = (rightdomain - leftdomain) / (node-1)
    dt = cfl * dx / u

    ! set timer for allocation, communication & computation time
    t_comm = 0
    t_comp = 0

    ! strings for writing numerical solution
    write (cproblemset, *) problemset
    write (csolver, *) solver
    write (cnode, *) node
    write (ct, '(f5.1)') t
    write (ccfl, '(f5.1)') cfl
    filename = "mpi_set"//trim(adjustl(cproblemset))//"_sol"//trim(adjustl(csolver))//"_n"//trim(adjustl(cnode))//"_t"//trim(adjustl(ct))//"_cfl"//trim(adjustl(ccfl))//".dat"
    timelog = "mpitime_set"//trim(adjustl(cproblemset))//"_sol"//trim(adjustl(csolver))//"_n"//trim(adjustl(cnode))//"_t"//trim(adjustl(ct))//"_cfl"//trim(adjustl(ccfl))//".dat"

    open(1, file = timelog)

!---------------------------------------------------------------------------------------------------------------------------------------
    call MPI_INIT(IERROR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, IERROR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, IERROR)

    ! start counting communication time
    t_start = MPI_WTIME()

    ! allocate indexs to different processors
    ! for master
    if (myid .eq. 0) then
        quotient = node / numprocs
        remainder = mod(node, numprocs)
        do i=1, numprocs-1
            if (remainder .eq. 0) then
                start_index = i * quotient + 1
                numnodes = quotient
            else
                if (i .le. quotient) then
                    start_index = i * quotient + 1 + i
                    numnodes = quotient + 1
                else
                    start_index = i * quotient + 1 + quotient
                    numnodes = quotient
                end if
            end if
            
            call MPI_SEND(start_index, 1, MPI_INT, i, 0, MPI_COMM_WORLD, IERROR)
            call MPI_SEND(numnodes, 1, MPI_INT, i, 1, MPI_COMM_WORLD, IERROR)
        end do
        
        ! set start_index & numnodes for master
        start_index = 1
        if (remainder .eq. 0) then
            numnodes = quotient
        else
            numnodes = quotient + 1
        end if
    end if

    ! wait for sending out information
    call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

    ! for workers
    do i=1, numprocs-1
        if (myid .eq. i) then
            call MPI_RECV(start_index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, istatus, IERROR)
            call MPI_RECV(numnodes, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, istatus, IERROR)
        end if
    end do

    ! wait for receiving information
    call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

    ! compute communication time
    t_end = MPI_WTIME()
    t_comm = t_comm + (t_end - t_start)
    write (1,*) "comm: send & receive index", t_end - t_start
    t_start = MPI_WTIME()

    ! allocate arrays for all processors, initialization & setup boundary conditions 
    do i=0, numprocs-1
        ! allocate x coordinate array
        if (myid .eq. i) then
            allocate(x(numnodes))
            allocate(f0(numnodes))
            allocate(f1(numnodes))
            if (solver .eq. 5) then
                allocate(fs(numnodes))
            end if

            x(1:numnodes) = (/ (leftdomain + (start_index + i - 2) * dx, i=1, numnodes) /)
            if (problemset .eq. 0) then
                f0(1:numnodes) = (/ (0.5 * (sign(1.0, x(i)) + 1), i=1, numnodes) /)
            else if (problemset .eq. 1) then
                f0(1:numnodes) = (/ (0.5 * exp(-x(i) ** 2), i=1, numnodes) /)
            end if
        end if

        ! setup boundary conditions
        if (myid .eq. 0) then
            f0(1) = leftbound
        else if (myid .eq. numprocs-1) then
            f0(numnodes) = rightbound
        end if
    end do

    ! wait for initialization finished
    call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

    ! compute computation time
    t_end = MPI_WTIME()
    t_comp = t_comp + (t_end - t_start)
    write (1,*) "comp: allocate array, assign x coord, set up boundary conditions", t_end - t_start
    t_start = MPI_WTIME()

!---------------------------------------------------------------------------------------------------------------------------------------
    ! numerical calculation
    timer = 0
    do while (timer .le. t)
        ! pass information to adjcent processors
        t_start = MPI_WTIME()
        do i=0, numprocs-1
            if ((solver .eq. 1) .or. (solver .eq. 5)) then
                ! solver 1 & 5: pass leftmost data to processor on the left
                if (myid .eq. 0) then
                else if (myid .eq. i) then
                    leftvalue = f0(1)
                    call MPI_SEND(leftvalue, 1, MPI_REAL, i-1, 0, MPI_COMM_WORLD, IERROR)
                end if
            else if (solver .eq. 2) then
                ! solver 2: pass rightmost data to processor on the right
                if (myid .eq. numprocs-1) then
                else if (myid .eq. i) then
                    rightvalue = f0(numnodes)
                    call MPI_SEND(rightvalue, 1, MPI_REAL, i+1, 1, MPI_COMM_WORLD, IERROR)
                end if
            else if ((solver .eq. 3) .or. (solver .eq. 4)) then
                ! solver 3 & 4: pass leftmost & rightmost data to processor on both sides
                if (myid .eq. 0) then
                    rightvalue = f0(numnodes)
                    call MPI_SEND(rightvalue, 1, MPI_REAL, i+1, 1, MPI_COMM_WORLD, IERROR)
                else if (myid .eq. numprocs-1) then
                    leftvalue = f0(1)
                    call MPI_SEND(leftvalue, 1, MPI_REAL, i-1, 0, MPI_COMM_WORLD, IERROR)
                else if (myid .eq. i) then
                    leftvalue = f0(1)
                    rightvalue = f0(numnodes)
                    call MPI_SEND(leftvalue, 1, MPI_REAL, i-1, 0, MPI_COMM_WORLD, IERROR)
                    call MPI_SEND(rightvalue, 1, MPI_REAL, i+1, 1, MPI_COMM_WORLD, IERROR)
                end if
            end if
        end do

        ! wait for sending out information
        call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

        ! receive information from adjcent processors
        do i=0, numprocs-1
            if ((solver .eq. 1) .or. (solver .eq. 5)) then
                ! solver 1 & 5: receive leftmost data from processor on the right
                if (myid .eq. numprocs-1) then
                else if (myid .eq. i) then
                    call MPI_RECV(righthalo, 1, MPI_REAL, i+1, 0, MPI_COMM_WORLD, IERROR)
                end if
            else if (solver .eq. 2) then
                ! solver 2: receive rightmost data from processor on the left
                if (myid .eq. 0) then
                else if (myid .eq. i) then
                    call MPI_RECV(lefthalo, 1, MPI_REAL, i-1, 1, MPI_COMM_WORLD, IERROR)
                end if
            else if ((solver .eq. 3) .or. (solver .eq. 4)) then
                ! solver 3 & 4: receive leftmost & rightmost data from processor from both sides
                if (myid .eq. 0) then
                    call MPI_RECV(righthalo, 1, MPI_REAL, i+1, 0, MPI_COMM_WORLD, IERROR)
                else if (myid .eq. numprocs-1) then
                    call MPI_RECV(lefthalo, 1, MPI_REAL, i-1, 1, MPI_COMM_WORLD, IERROR)
                else if (myid .eq. i) then
                    call MPI_RECV(righthalo, 1, MPI_REAL, i+1, 0, MPI_COMM_WORLD, IERROR)
                    call MPI_RECV(lefthalo, 1, MPI_REAL, i-1, 1, MPI_COMM_WORLD, IERROR)
                end if
            end if
        end do

        ! wait for receiving information
        call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

        ! compute communication time
        t_end = MPI_WTIME()
        t_comm = t_comm + (t_end - t_start)
        write (1,*) "comm: pass and receive data from neighbours", timer, t_end - t_start
        t_start = MPI_WTIME()

        !*******************************************************************************************************************************
        ! calculation
        do i=0, numprocs-1
            if (myid .eq. 0) then
                ! left end boundary
                f1(1) = leftbound
                if (solver .eq. 1) then
                    f1(2:numnodes-1) = (/ (f0(i) - cfl * (f0(i+1) - f0(i)), i=2, numnodes-1) /)
                    f1(numnodes) = f0(numnodes) - cfl * (righthalo - f0(numnodes))
                else if (solver .eq. 2) then
                    f1(2:numnodes) = (/ (f0(i) - cfl * (f0(i) - f0(i-1)), i=2, numnodes) /)
                else if (solver .eq. 3) then
                    f1(2:numnodes-1) = (/ (0.5 * (f0(i+1) + f0(i-1)) - 0.5 * cfl * (f0(i+1) - f0(i-1)), i=2, numnodes-1) /)
                    f1(numnodes) = 0.5 * (righthalo + f0(i-1)) - 0.5 * cfl * (righthalo - f0(i-1))
                else if (solver .eq. 4) then
                    f1(2:numnodes-1) = (/ (f0(i) - 0.5 * cfl * (f0(i+1) - f0(i-1)) + 0.5 * cfl ** 2 * (f0(i+1) - 2 * f0(i) + f0(i-1)), i=2, numnodes-1) /)
                    f1(numnodes) = f0(i) - 0.5 * cfl * (righthalo - f0(i-1)) + 0.5 * cfl ** 2 * (righthalo - 2 * f0(i) + f0(i-1))
                else if (solver .eq. 5) then
                    fs(1) = leftbound ! apply boundary condition to left end of domain
                    fs(2:numnodes-1) = (/ (f0(i) - cfl * (f0(i+1) - f0(i)), i=2, numnodes-1) /)
                    fs(numnodes) = f0(i) - cfl * (righthalo - f0(i))
                end if
            else if (myid .eq. numprocs-1) then
                ! right end boundary
                f1(numnodes) = rightbound
                if (solver .eq. 1) then
                    f1(1:numnodes-1) = (/ (f0(i) - cfl * (f0(i+1) - f0(i)), i=2, numnodes-1) /)
                else if (solver .eq. 2) then
                    f1(1) = f0(i) - cfl * (f0(i) - lefthalo)
                    f1(2:numnodes) = (/ (f0(i) - cfl * (f0(i) - f0(i-1)), i=2, numnodes) /)
                else if (solver .eq. 3) then
                    f1(1) = 0.5 * (f0(i+1) + lefthalo) - 0.5 * cfl * (f0(i+1) - lefthalo)
                    f1(2:numnodes-1) = (/ (0.5 * (f0(i+1) + f0(i-1)) - 0.5 * cfl * (f0(i+1) - f0(i-1)), i=2, numnodes-1) /)
                else if (solver .eq. 4) then
                    f1(1) = f0(i) - 0.5 * cfl * (f0(i+1) - lefthalo) + 0.5 * cfl ** 2 * (f0(i+1) - 2 * f0(i) + lefthalo)
                    f1(2:numnodes-1) = (/ (f0(i) - 0.5 * cfl * (f0(i+1) - f0(i-1)) + 0.5 * cfl ** 2 * (f0(i+1) - 2 * f0(i) + f0(i-1)), i=2, numnodes-1) /)
                else if (solver .eq. 5) then
                    fs(numnodes) = rightbound ! apply boundary condition to right end of domain
                    fs(1:numnodes-1) = (/ (f0(i) - cfl * (f0(i+1) - f0(i)), i=2, numnodes-1) /)
                end if
            else if (myid .eq. i) then
                ! middle domains
                if (solver .eq. 1) then
                    f1(1:numnodes-1) = (/ (f0(i) - cfl * (f0(i+1) - f0(i)), i=2, numnodes-1) /)
                    f1(numnodes) = f0(numnodes) - cfl * (righthalo - f0(numnodes))
                else if (solver .eq. 2) then
                    f1(1) = f0(i) - cfl * (f0(i) - lefthalo)
                    f1(2:numnodes) = (/ (f0(i) - cfl * (f0(i) - f0(i-1)), i=2, numnodes) /)
                else if (solver .eq. 3) then
                    f1(1) = 0.5 * (f0(i+1) + lefthalo) - 0.5 * cfl * (f0(i+1) - lefthalo)
                    f1(2:numnodes-1) = (/ (0.5 * (f0(i+1) + f0(i-1)) - 0.5 * cfl * (f0(i+1) - f0(i-1)), i=2, numnodes-1) /)
                    f1(numnodes) = 0.5 * (righthalo + f0(i-1)) - 0.5 * cfl * (righthalo - f0(i-1))
                else if (solver .eq. 4) then
                    f1(1) = f0(i) - 0.5 * cfl * (f0(i+1) - lefthalo) + 0.5 * cfl ** 2 * (f0(i+1) - 2 * f0(i) + lefthalo)
                    f1(2:numnodes-1) = (/ (f0(i) - 0.5 * cfl * (f0(i+1) - f0(i-1)) + 0.5 * cfl ** 2 * (f0(i+1) - 2 * f0(i) + f0(i-1)), i=2, numnodes-1) /)
                    f1(numnodes) = f0(i) - 0.5 * cfl * (righthalo - f0(i-1)) + 0.5 * cfl ** 2 * (righthalo - 2 * f0(i) + f0(i-1))
                else if (solver .eq. 5) then
                    fs(1:numnodes-1) = (/ (f0(i) - cfl * (f0(i+1) - f0(i)), i=2, numnodes-1) /)
                    fs(numnodes) = f0(i) - cfl * (righthalo - f0(i))
                end if
            end if
        end do
        
        ! compute computation time
        t_end = MPI_WTIME()
        t_comp = t_comp + (t_end - t_start)
        write (1,*) "comm: calulation", timer, t_end - t_start
        t_start = MPI_WTIME()

        !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ! communication for additional time steps for solver 5
        if (solver .eq. 5) then
            ! pass rightmost data to processor on the right
            do i=0, numprocs-1
                if (myid .eq. numprocs-1) then
                else if (myid .eq. i) then
                    rightvalue = fs(numnodes)
                    call MPI_SEND(rightvalue, 1, MPI_REAL, i+1, 0, MPI_COMM_WORLD, IERROR)
                end if
            end do

            ! wait for sending out information
            call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

            ! receive rightmost data from processor on the left
            do i=0, numprocs-1
                if (myid .eq. 0) then
                else if (myid .eq. i) then
                    call MPI_RECV(lefthalo, 1, MPI_REAL, i-1, 0, MPI_COMM_WORLD, IERROR)
                end if
            end do

            ! wait for receiving information
            call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        
            ! compute communication time
            t_end = MPI_WTIME()
            t_comm = t_comm + (t_end - t_start)
            write (1,*) "comm: pass and receive data from neighbours for extra time step, solver 5", timer, t_end - t_start
            t_start = MPI_WTIME()

            ! compute for additional time steps for solver 5
            if (solver .eq. 5) then
                do i=0, numprocs-1
                    if (myid .eq. 0) then
                        f1(1) = leftbound
                        f1(2:numnodes) = (/ (0.5 * ((f0(i) + fs(i)) - cfl * (fs(i) - fs(i-1))), i=2, numnodes) /)
                    else if (myid .eq. i) then
                        f1(1) = 0.5 * ((f0(i) + fs(i)) - cfl * (fs(i) - lefthalo))
                        f1(2:numnodes) = (/ (0.5 * ((f0(i) + fs(i)) - cfl * (fs(i) - fs(i-1))), i=2, numnodes) /)
                    end if
                end do
            end if

            ! compute computation time
            t_end = MPI_WTIME()
            t_comp = t_comp + (t_end - t_start)
            write (1,*) "comm: calulation for extra time step, solver 5", timer, t_end - t_start
            t_start = MPI_WTIME()
        end if
        !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        do i=0, numprocs-1
            if (myid .eq. i) then
                f0 = f1
            end if
        end do

        timer = timer + dt

        ! compute computation time
        t_end = MPI_WTIME()
        t_comp = t_comp + (t_end - t_start)
        write (1,*) "comm: calulation, arrage arrays & timer", timer, t_end - t_start
        t_start = MPI_WTIME()
    end do

    write (1,*) "total time"
    write (1,*) "comm:", t_comm
    write (1,*) "comp:", t_comp
    close (1)

!---------------------------------------------------------------------------------------------------------------------------------------
    ! write numerical solution
    open(2, file = filename)

    do i=0, numprocs-1
        if (myid .eq. i) then
            do j=1, numnodes
                write (1, *) x(i), f0(i)
            end do
        end if
    end do

    close (2)

!---------------------------------------------------------------------------------------------------------------------------------------
    ! deallocate arrays and finalize mpi
    do i=0, numprocs-1
        if (myid .eq. i) then
            deallocate(x, f0, f1)
            if (solver .eq. 5) then
                deallocate(fs)
            end if
        end if
    end do

    call MPI_FINALIZE(IERROR)

    print *, "finish"
end program mpi_heatadvection
