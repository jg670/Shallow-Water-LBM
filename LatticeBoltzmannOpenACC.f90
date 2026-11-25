module lattice_boltzmann_OpenACC
    implicit none

    ! type used for velocity vectors !
    type :: velocity
        real(8) :: x
        real(8) :: y
    end type velocity

    ! Lattice dimensions, will be overwritten by command line arguments !
    integer :: global_width = 15, global_height = 15, instance_height = 17, instance_width = 17
    
    ! Simulation parameters !
    integer, parameter :: directions = 9
    real(8) :: omega = 1
    real(8), parameter :: shift_directions_x(directions) = [0.0,  1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0]
    real(8), parameter :: shift_directions_y(directions) = [0.0,  0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0]
    real(8), parameter :: weights(directions) = [4.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/9.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8, 1.0_8/36.0_8]

    ! Persistent MPI communication buffers !
    real(8), allocatable :: send_w_buf(:,:), recv_w_buf(:,:)
    real(8), allocatable :: send_e_buf(:,:), recv_e_buf(:,:)
    real(8), allocatable :: send_n_buf(:,:), recv_n_buf(:,:)
    real(8), allocatable :: send_s_buf(:,:), recv_s_buf(:,:)
    
    ! Wall velocity, used when simulation has a moving wall !
    type(velocity) :: wall_speed = velocity(0.1_8, 0.0_8)

    ! Width of coarray grid !
    integer :: coarray_dimensions = 1

    contains

        ! Perform one time step of the LBM simulation !
        subroutine perform_one_time_step(lattice_in, lattice_out, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
            
            ! Input variable declarations !
            real(8), intent(inout) :: lattice_in(instance_width, instance_height, directions)
            real(8), intent(inout) :: lattice_out(instance_width, instance_height, directions)
            integer, intent(in) :: comm_cart, n_n, n_s, n_e, n_w
            integer, intent(in) :: coords(2), dims(2)

            ! Update all ghost nodes on input array !
            call sync_ghost_nodes(lattice_in, comm_cart, n_n, n_s, n_e, n_w)

            ! Apply boundary conditions to ghost nodes !
            call sliding_lid_boundary(lattice_in, coords, dims)

            ! Combined streaming and collision (pulls from lattice_in, computes, writes to lattice_out) !
            call stream_and_collide(lattice_in, lattice_out)

        end subroutine perform_one_time_step

        ! Synchronize ghost nodes across all images !
        subroutine sync_ghost_nodes(lattice, comm, n_n, n_s, n_e, n_w)
            use mpi

            ! Input variable declarations !
            real(8), intent(inout) :: lattice(instance_width, instance_height, directions)
            integer, intent(in) :: comm, n_n, n_s, n_e, n_w
    
            ! MPI request variables !
            integer :: req(8), req_count
            integer :: stat(MPI_STATUS_SIZE, 8), ierr

            ! Iterators !
            integer :: d, j, k
    
            req_count = 0

            ! Pack buffers on GPU and send to CPU!
            if (n_w /= MPI_PROC_NULL) then
                !$acc parallel loop collapse(2) present(lattice, send_w_buf)
                do k = 1, instance_height
                    do d = 1, directions
                        send_w_buf(k, d) = lattice(2, k, d)
                    end do
                end do
                !$acc update host(send_w_buf)
            end if

            if (n_e /= MPI_PROC_NULL) then
                !$acc parallel loop collapse(2) present(lattice, send_e_buf)
                do k = 1, instance_height
                    do d = 1, directions
                        send_e_buf(k, d) = lattice(instance_width - 1, k, d)
                    end do
                end do
                !$acc update host(send_e_buf)
            end if

            if (n_n /= MPI_PROC_NULL) then
                !$acc parallel loop collapse(2) present(lattice, send_n_buf)
                do j = 1, instance_width
                    do d = 1, directions
                        send_n_buf(j, d) = lattice(j, 2, d)
                    end do
                end do
                !$acc update host(send_n_buf)
            end if

            if (n_s /= MPI_PROC_NULL) then
                !$acc parallel loop collapse(2) present(lattice, send_s_buf)
                do j = 1, instance_width
                    do d = 1, directions
                        send_s_buf(j, d) = lattice(j, instance_height - 1, d)
                    end do
                end do
                !$acc update host(send_s_buf)
            end if

            ! X-direction sync !
            if (n_w /= MPI_PROC_NULL) then
                req_count = req_count + 1
                call MPI_Isend(send_w_buf, directions*instance_height, MPI_DOUBLE_PRECISION, n_w, 0, comm, req(req_count), ierr)
                req_count = req_count + 1
                call MPI_Irecv(recv_w_buf, directions*instance_height, MPI_DOUBLE_PRECISION, n_w, 1, comm, req(req_count), ierr)
            end if

            if (n_e /= MPI_PROC_NULL) then
                req_count = req_count + 1
                call MPI_Isend(send_e_buf, directions*instance_height, MPI_DOUBLE_PRECISION, n_e, 1, comm, req(req_count), ierr)
                req_count = req_count + 1
                call MPI_Irecv(recv_e_buf, directions*instance_height, MPI_DOUBLE_PRECISION, n_e, 0, comm, req(req_count), ierr)
            end if

            ! Y-direction sync !
            if (n_n /= MPI_PROC_NULL) then
                req_count = req_count + 1
                call MPI_Isend(send_n_buf, directions*instance_width, MPI_DOUBLE_PRECISION, n_n, 2, comm, req(req_count), ierr)
                req_count = req_count + 1
                call MPI_Irecv(recv_n_buf, directions*instance_width, MPI_DOUBLE_PRECISION, n_n, 3, comm, req(req_count), ierr)
            end if

            if (n_s /= MPI_PROC_NULL) then
                req_count = req_count + 1
                call MPI_Isend(send_s_buf, directions*instance_width, MPI_DOUBLE_PRECISION, n_s, 3, comm, req(req_count), ierr)
                req_count = req_count + 1
                call MPI_Irecv(recv_s_buf, directions*instance_width, MPI_DOUBLE_PRECISION, n_s, 2, comm, req(req_count), ierr)
            end if

            ! Wait for all communication to stop !
            if (req_count > 0) call MPI_Waitall(req_count, req(1:req_count), stat, ierr)

            ! Send buffers to GPU and unpack to lattice !
            if (n_w /= MPI_PROC_NULL) then
                !$acc update device(recv_w_buf)
                !$acc parallel loop collapse(2) present(lattice, recv_w_buf)
                do k = 1, instance_height
                    do d = 1, directions
                        lattice(1, k, d) = recv_w_buf(k, d)
                    end do
                end do
            end if

            if (n_e /= MPI_PROC_NULL) then
                !$acc update device(recv_e_buf)
                !$acc parallel loop collapse(2) present(lattice, recv_e_buf)
                do k = 1, instance_height
                    do d = 1, directions
                        lattice(instance_width, k, d) = recv_e_buf(k, d)
                    end do
                end do
            end if

            if (n_n /= MPI_PROC_NULL) then
                !$acc update device(recv_n_buf)
                !$acc parallel loop collapse(2) present(lattice, recv_n_buf)
                do j = 1, instance_width
                    do d = 1, directions
                        lattice(j, 1, d) = recv_n_buf(j, d)
                    end do
                end do
            end if

            if (n_s /= MPI_PROC_NULL) then
                !$acc update device(recv_s_buf)
                !$acc parallel loop collapse(2) present(lattice, recv_s_buf)
                do j = 1, instance_width
                    do d = 1, directions
                        lattice(j, instance_height, d) = recv_s_buf(j, d)
                    end do
                end do
            end if

        end subroutine sync_ghost_nodes

        ! Combined streaming and collision step !
        subroutine stream_and_collide(lattice_in, lattice_out)

            ! Input variable declarations !
            real(8), intent(in) :: lattice_in(instance_width, instance_height, directions)
            real(8), intent(inout) :: lattice_out(instance_width, instance_height, directions)
            
            ! Iterators and data variables !
            integer :: j, k
            real(8) :: f0, f1, f2, f3, f4, f5, f6, f7, f8
            real(8) :: density, vel_x, vel_y, u_sq, c_dot_u

            ! Loop through lattice on GPU !
            !$acc parallel loop collapse(2) present(lattice_in, lattice_out) private(f0, f1, f2, f3, f4, f5, f6, f7, f8, density, vel_x, vel_y, u_sq, c_dot_u)
            do k = 2, instance_height - 1
                do j = 2, instance_width - 1
                    
                    ! Take data that will "stream in" from neighbors in input lattice !
                    f0 = lattice_in(j, k, 1)
                    f1 = lattice_in(j - 1, k, 2)     ! x: 1, y: 0   -> pulled from left node !
                    f2 = lattice_in(j, k + 1, 3)     ! x: 0, y: 1   -> pulled from bottom node !
                    f3 = lattice_in(j + 1, k, 4)     ! x:-1, y: 0   -> pulled from right node !
                    f4 = lattice_in(j, k - 1, 5)     ! x: 0, y:-1   -> pulled from top node !
                    f5 = lattice_in(j - 1, k + 1, 6) ! x: 1, y: 1   -> pulled from bottom left node !
                    f6 = lattice_in(j + 1, k + 1, 7) ! x:-1, y: 1   -> pulled from bottom rght node !
                    f7 = lattice_in(j + 1, k - 1, 8) ! x:-1, y:-1   -> pulled from top right node !
                    f8 = lattice_in(j - 1, k - 1, 9) ! x: 1, y:-1   -> pulled from top left node !
                    
                    ! Values used for collision calculations !
                    density = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8
                    vel_x = (f1 - f3 + f5 - f6 - f7 + f8) / density
                    vel_y = (f2 - f4 + f5 + f6 - f7 - f8) / density
                    u_sq = vel_x**2 + vel_y**2

                    ! Do collision calculations and save to output lattice !
                    lattice_out(j, k, 1) = f0 - omega * (f0 - weights(1) * density * (1.0_8 - 1.5_8 * u_sq))

                    c_dot_u = vel_x
                    lattice_out(j, k, 2) = f1 - omega * (f1 - weights(2) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = vel_y
                    lattice_out(j, k, 3) = f2 - omega * (f2 - weights(3) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = -vel_x
                    lattice_out(j, k, 4) = f3 - omega * (f3 - weights(4) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = -vel_y
                    lattice_out(j, k, 5) = f4 - omega * (f4 - weights(5) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = vel_x + vel_y
                    lattice_out(j, k, 6) = f5 - omega * (f5 - weights(6) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = -vel_x + vel_y
                    lattice_out(j, k, 7) = f6 - omega * (f6 - weights(7) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = -vel_x - vel_y
                    lattice_out(j, k, 8) = f7 - omega * (f7 - weights(8) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = vel_x - vel_y
                    lattice_out(j, k, 9) = f8 - omega * (f8 - weights(9) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                end do
            end do
        end subroutine stream_and_collide

        ! Boundary conditions for sliding lid simulations !
        subroutine sliding_lid_boundary(lattice, coords, dims)

            ! Input variable declarations !
            real(8), intent(inout) :: lattice(instance_width, instance_height, directions)
            integer, intent(in) :: coords(2), dims(2)

            ! Iterator !
            integer :: j

            ! Left side bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            if (coords(1) == 0) then
                !$acc kernels present(lattice)
                lattice(1,:,2) = lattice(2,:,4)
                lattice(1,2:instance_height,6) = lattice(2,1:instance_height-1,8)
                lattice(1,1:instance_height-1,9) = lattice(2,2:instance_height,7)
                !$acc end kernels
            end if

            ! Right side bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            if (coords(1) == dims(1) - 1) then
                !$acc kernels present(lattice)
                lattice(instance_width,:,4) = lattice(instance_width-1,:,2)
                lattice(instance_width,2:instance_height,7) = lattice(instance_width-1,1:instance_height-1,9)
                lattice(instance_width,1:instance_height-1,8) = lattice(instance_width-1,2:instance_height,6)
                !$acc end kernels
            end if

            ! Bottom bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            if (coords(2) == dims(2) - 1) then
                !$acc kernels present(lattice)
                lattice(:,instance_height,3) = lattice(:,instance_height-1,5)
                lattice(1:instance_width-1,instance_height,6) = lattice(2:instance_width,instance_height-1,8)
                lattice(2:instance_width,instance_height,7) = lattice(1:instance_width-1,instance_height-1,9)
                !$acc end kernels
            end if

            ! Top moving boundary (y velocity is 0, so nothing is added to direction 4 (indexed as 5)) !
            if (coords(2) == 0) then
                !$acc parallel loop present(lattice)
                do j = 1, instance_width
                    lattice(j, 1, 5) = lattice(j, 2, 3)

                    if (j < instance_width) then
                        lattice(j+1, 1, 8) = lattice(j, 2, 6) - 2.0_8 * weights(6) * (lattice(j, 2, 1) + lattice(j, 2, 2) + lattice(j, 2, 3) + lattice(j, 2, 4) + lattice(j, 2, 5) + lattice(j, 2, 6) + lattice(j, 2, 7) + lattice(j, 2, 8) + lattice(j, 2, 9)) * ((shift_directions_x(6) * wall_speed%x) / (1.0_8 / 3.0_8))
                        lattice(j, 1, 9) = lattice(j+1, 2, 7) - 2.0_8 * weights(7) * (lattice(j+1, 2, 1) + lattice(j+1, 2, 2) + lattice(j+1, 2, 3) + lattice(j+1, 2, 4) + lattice(j+1, 2, 5) + lattice(j+1, 2, 6) + lattice(j+1, 2, 7) + lattice(j+1, 2, 8) + lattice(j+1, 2, 9)) * ((shift_directions_x(7) * wall_speed%x) / (1.0_8 / 3.0_8))
                    end if
                end do
            end if

        end subroutine sliding_lid_boundary

        ! Initialize lattice for sliding lid simulation !
        subroutine populate_lattice_sliding_lid(lattice)

            ! Input variable declarations !
            real(8), intent(inout) :: lattice(instance_width, instance_height, directions)

            ! Initial density and velocity array values !
            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)

            ! Populate arrays to 1.0 density and 0.0 velocity !
            density_arr = 1.0_8
            velocity_arr = velocity(0.0_8, 0.0_8)

            ! Set lid velocity given desired reynold's number !
            call set_lid_velocity_given_reynolds(1000.0_8)

            ! Set input lattice to new lattice with 1.0 density and 0.0 velocity !
            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_sliding_lid

        ! Helper function to set lid speed based on width and desired reynold's number !
        ! IMPORTANT: only works with top or bottom walls, or side walls if lattice is a square (this code only takes width into account) !
        subroutine set_lid_velocity_given_reynolds(reynolds_number)

            ! Input variable declaration !
            real(8), intent(in) ::  reynolds_number

            ! Calculated viscosity variable !
            real(8) :: analytical_viscosity

            ! Calculate viscosity !
            analytical_viscosity = (1.0_8 / 3.0_8) * ((1.0_8 / omega) - 0.5_8)
            
            ! Update wall speed to fit desired reynold's number !
            wall_speed = velocity((reynolds_number * analytical_viscosity) / (coarray_dimensions * (instance_width - 2.0_8)), 0.0_8)
        
        end subroutine set_lid_velocity_given_reynolds

        ! Helper function to calculate equilibrium !
        function calculate_equilibrium(density_arr, velocity_arr)

            ! Input variable declarations !
            real(8), intent(in) :: density_arr(instance_width, instance_height)
            type(velocity), intent(in) :: velocity_arr(instance_width, instance_height)

            ! Output variable declaration !
            real(8) :: calculate_equilibrium(instance_width, instance_height, directions)

            ! Interator !
            integer :: i

            ! Calculation variables !
            real(8) :: u_squared(instance_width, instance_height), velocity_sum(instance_width, instance_height)

            ! Populate the whole array with 0.0 values !
            calculate_equilibrium = 0.0_8

            ! Total velocity value used for equilibrium calculation !
            u_squared(2:instance_width-1, 2:instance_height-1) = velocity_arr(2:instance_width-1, 2:instance_height-1)%x**2 + velocity_arr(2:instance_width-1, 2:instance_height-1)%y**2

            do i=1, directions
                ! Calculation values used in equilibrium equation !
                velocity_sum(2:instance_width-1, 2:instance_height-1) = shift_directions_x(i) * velocity_arr(2:instance_width-1, 2:instance_height-1)%x + shift_directions_y(i) * velocity_arr(2:instance_width-1, 2:instance_height-1)%y

                ! Calculate equilibrium lattice !
                calculate_equilibrium(2:instance_width-1, 2:instance_height-1, i) = weights(i) * density_arr(2:instance_width-1, 2:instance_height-1) * (1.0_8 + 3.0_8 * velocity_sum(2:instance_width-1, 2:instance_height-1) + (4.5_8 * velocity_sum(2:instance_width-1, 2:instance_height-1)**2) - 1.5_8 * u_squared(2:instance_width-1, 2:instance_height-1))
            end do

        end function calculate_equilibrium

	    subroutine init_mpi_buffers()

            ! Allocate buffers to CPU memory !
            allocate(send_w_buf(instance_height, directions))
            allocate(recv_w_buf(instance_height, directions))
            allocate(send_e_buf(instance_height, directions))
            allocate(recv_e_buf(instance_height, directions))
            
            allocate(send_n_buf(instance_width, directions))
            allocate(recv_n_buf(instance_width, directions))
            allocate(send_s_buf(instance_width, directions))
            allocate(recv_s_buf(instance_width, directions))

            ! Map buffers to GPU !
            !$acc enter data create(send_w_buf, recv_w_buf, send_e_buf, recv_e_buf)
            !$acc enter data create(send_n_buf, recv_n_buf, send_s_buf, recv_s_buf)

        end subroutine init_mpi_buffers

        subroutine cleanup_mpi_buffers()

            ! Remove buffers from GPU memory !
            !$acc exit data delete(send_w_buf, recv_w_buf, send_e_buf, recv_e_buf)
            !$acc exit data delete(send_n_buf, recv_n_buf, send_s_buf, recv_s_buf)

            ! Deallocate buffers on CPU memory !
            deallocate(send_w_buf, recv_w_buf, send_e_buf, recv_e_buf)
            deallocate(send_n_buf, recv_n_buf, send_s_buf, recv_s_buf)

        end subroutine cleanup_mpi_buffers

end module lattice_boltzmann_OpenACC

program LatticeBoltzmannOpenACCMain
    use lattice_boltzmann_OpenACC
    use mpi
    use openacc
    implicit none

    integer :: ierr, rank, comm_size
    integer :: dims(2), coords(2)
    integer :: comm_cart
    integer :: neighbor_n, neighbor_s, neighbor_e, neighbor_w
    integer :: args_count, io_status, interval_len, interval_num
    integer :: norm_width, norm_height, mod_width, mod_height
    character(len=20) :: width_arg, height_arg, coarray_dim_arg, interval_length_arg, num_intervals_arg
    logical :: periods(2)

    integer :: local_comm, local_rank
    integer :: num_devices, my_device

    ! Lattice arrays !
    real(8), allocatable :: lattice_even(:,:,:), lattice_odd(:,:,:)

    ! Initialize MPI !
    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, comm_size, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! Split global communicator to local communicator for processes on the same node !
    call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, local_comm, ierr)

    ! Determine the rank of this process in its node !
    call MPI_Comm_rank(local_comm, local_rank, ierr)
    call MPI_Comm_free(local_comm, ierr)

    ! Get total available NVIDIA GPUs on this node !
    num_devices = acc_get_num_devices(acc_device_nvidia)

    if (num_devices > 0) then

        ! Map the local rank to a unique GPU index !
        my_device = mod(local_rank, num_devices)
        call acc_set_device_num(my_device, acc_device_nvidia)
    else
        if (rank == 0) print *, "CRITICAL ERROR: No NVIDIA GPUs detected."
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    ! Get command line arguments !
    args_count = command_argument_count()
    if (args_count < 5) then
        if (rank == 0) print *, "Missing arguments: [width height grid_x interval_length num_intervals]"
        call MPI_Finalize(ierr)
        stop
    end if

    ! Parse arguments !
    call get_command_argument(1, width_arg); read(width_arg, *) global_width
    call get_command_argument(2, height_arg); read(height_arg, *) global_height
    call get_command_argument(3, coarray_dim_arg); read(coarray_dim_arg, *) coarray_dimensions
    call get_command_argument(4, interval_length_arg); read(interval_length_arg, *) interval_len
    call get_command_argument(5, num_intervals_arg); read(num_intervals_arg, *) interval_num

    ! Define the MPI grid !
    dims(1) = coarray_dimensions
    dims(2) = comm_size / coarray_dimensions

    ! Don't use built-in periodic bounds !
    periods = [.false., .false.]

    ! Create MPI grid !
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .false., comm_cart, ierr)
    call MPI_Comm_rank(comm_cart, rank, ierr)
    call MPI_Cart_coords(comm_cart, rank, 2, coords, ierr)

    ! Get neighboring image ranks !
    call MPI_Cart_shift(comm_cart, 0, 1, neighbor_w, neighbor_e, ierr)
    call MPI_Cart_shift(comm_cart, 1, 1, neighbor_n, neighbor_s, ierr)

    ! Determine width of image (and get remainder if not evenly divisible) !
    norm_width = global_width / dims(1)
    mod_width = mod(global_width, dims(1))

    ! If current image is on the border, then add remainder to its width !
    if (coords(1) == dims(1) - 1) then
        instance_width = norm_width + mod_width + 2
    else
        instance_width = norm_width + 2
    end if

    ! Do the same operation to height !
    norm_height = global_height / dims(2)
    mod_height = mod(global_height, dims(2))
    if (coords(2) == dims(2) - 1) then
        instance_height = norm_height + mod_height + 2
    else
        instance_height = norm_height + 2
    end if

    ! Allocate local image arrays !
    allocate(lattice_even(instance_width, instance_height, directions))
    allocate(lattice_odd(instance_width, instance_height, directions))

    ! Intialize buffers on the GPU and CPU !
    call init_mpi_buffers()

    ! MPI + GPU version only supports sliding lid !
    call populate_lattice_sliding_lid(lattice_even)
    lattice_odd = lattice_even
    
    ! Sync processes before performance test !
    call MPI_Barrier(comm_cart, ierr)

    ! Maps arrays onto the GPU !
    !$acc enter data copyin(lattice_even, lattice_odd)

    ! Perform performance test !
    call do_parallel_performance_test(lattice_even, lattice_odd, comm_cart, neighbor_n, neighbor_s, neighbor_e, neighbor_w, coords, dims)

    ! Clean up GPU memory !
    !$acc exit data copyout(lattice_even)
    !$acc exit data delete(lattice_odd)

    ! Clean up buffers !
    call cleanup_mpi_buffers()

    ! Shut down MPI !
    call MPI_Finalize(ierr)

    contains

        ! Runs the simulation for interval_len * interval_num time steps, outputting the lattice state every lattice_len steps interval_num times !
        subroutine output_results(lattice_even, lattice_odd, interval_length, num_intervals, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
            use mpi

            ! Input variable declarations !
            real(8), intent(inout) :: lattice_even(instance_width, instance_height, directions)
            real(8), intent(inout) :: lattice_odd(instance_width, instance_height, directions)
            integer, intent(in) :: interval_length, num_intervals
            integer, intent(in) :: comm_cart, n_n, n_s, n_e, n_w, coords(2), dims(2)

            ! Iterators !
            integer :: interval, sub_interval, ierr
            integer(8) :: total_steps

            ! Initial lattice output (step 0) !
            call gather_and_write(lattice_even, 0, comm_cart, dims)
            total_steps = 0

            ! Interval handling loop !
            do interval = 1, num_intervals
                do sub_interval = 1, interval_length
                    total_steps = total_steps + 1
                    
                    ! If total_steps is odd, use lattice_even as input data and lattice_odd as output lattice, if total_steps is even then vice versa !
                    if (mod(total_steps, 2) /= 0) then
                        call perform_one_time_step(lattice_even, lattice_odd, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
                    else 
                        call perform_one_time_step(lattice_odd, lattice_even, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
                    end if
                end do

                ! Ensure all images finish calculating before image 1 starts reading !
                !$acc wait
                call MPI_Barrier(comm_cart, ierr)

                ! If total_steps is odd, output lattice_odd, otherwise output lattice_even !
                if (mod(total_steps, 2) /= 0) then
                    call gather_and_write(lattice_odd, interval * interval_length, comm_cart, dims)
                else 
                    call gather_and_write(lattice_even, interval * interval_length, comm_cart, dims)
                end if

            end do

        end subroutine output_results

        ! Helper subroutine to collect and write lattice data to output file !
        subroutine gather_and_write(lattice, step_num, comm, dims)
            use mpi

            ! Input variable declarations !
            real(8), intent(inout) :: lattice(instance_width, instance_height, directions)
            integer, intent(in) :: step_num, comm, dims(2)
            
            ! Global lattice aggregates !
            real(8), allocatable :: global_density(:,:)
            type(velocity), allocatable :: global_velocity(:,:)
            
            ! Image variables !
            integer :: rank, comm_size, remote_rank, ierr, max_w, max_h
            integer :: remote_coords(2), rx, ry, rw, rh, rx_start, ry_start
            integer :: stat(MPI_STATUS_SIZE)

            ! File I/O variables !
            character(len=50) :: file_name
            character(len=20) :: interval_str

            ! Iterators !
            integer :: d, j, k
            
            ! Get information about image grid !
            call MPI_Comm_rank(comm, rank, ierr)
            call MPI_Comm_size(comm, comm_size, ierr)

            ! Update the CPU copy with the data from the GPU !
            !$acc update host(lattice)
            
            ! Only run on first image (prevents race conditions / duplicate outputs) !
            if (rank == 0) then
                allocate(global_density(global_width, global_height))
                allocate(global_velocity(global_width, global_height))
                
                ! Loop over all images !
                do remote_rank = 0, comm_size - 1

                    ! Determine image position in image grid !
                    call MPI_Cart_coords(comm, remote_rank, 2, remote_coords, ierr)
                    rx = remote_coords(1)
                    ry = remote_coords(2) 
                    
                    ! Figure out size of the lattice on this image !
                    if (rx == dims(1) - 1) then
                        rw = (global_width / dims(1)) + mod(global_width, dims(1))
                    else
                        rw = global_width / dims(1)
                    end if
                    
                    if (ry == dims(2) - 1) then
                        rh = (global_height / dims(2)) + mod(global_height, dims(2))
                    else
                        rh = global_height / dims(2)
                    end if
                    
                    ! Map (image) local lattice coords onto global lattice !
                    rx_start = (rx) * (global_width / dims(1)) + 1
                    ry_start = (ry) * (global_height / dims(2)) + 1
                    
                    block
                        ! Temp variables !
                        real(8) :: exact_buffer(rw + 2, rh + 2, directions)
                        real(8) :: local_density(rw + 2, rh + 2)
                        type(velocity) :: local_velocity(rw + 2, rh + 2)
                    
                        ! Copy if local rank, otherwise receive from network !
                        if (remote_rank == 0) then
                            exact_buffer = lattice(1:rw+2, 1:rh+2, :)
                        else
                            call MPI_Recv(exact_buffer, directions * (rw+2) * (rh+2), MPI_DOUBLE_PRECISION, &
                                          remote_rank, 0, comm, stat, ierr)
                        end if
                        
                        ! Calculate density if image's lattice !
                        local_density(1:rw+2, 1:rh+2) = 1.0_8
                        local_density(2:rw+1, 2:rh+1) = sum(exact_buffer(2:rw+1, 2:rh+1, :), dim=3)
                        
                        ! Intialize velocity arrays to 0 (prevents any NaN errors for ghost nodes) !
                        local_velocity(1:rw+2, 1:rh+2)%x = 0.0_8
                        local_velocity(1:rw+2, 1:rh+2)%y = 0.0_8
                        
                        ! Calculate lattice momentum !
                        do d = 1, directions
                            local_velocity(2:rw+1, 2:rh+1)%x = local_velocity(2:rw+1, 2:rh+1)%x + exact_buffer(2:rw+1, 2:rh+1, d) * shift_directions_x(d)
                            local_velocity(2:rw+1, 2:rh+1)%y = local_velocity(2:rw+1, 2:rh+1)%y + exact_buffer(2:rw+1, 2:rh+1, d) * shift_directions_y(d)
                        end do
                        
                        ! Calculate lattice velocity !
                        local_velocity(2:rw+1, 2:rh+1)%x = local_velocity(2:rw+1, 2:rh+1)%x / local_density(2:rw+1, 2:rh+1)
                        local_velocity(2:rw+1, 2:rh+1)%y = local_velocity(2:rw+1, 2:rh+1)%y / local_density(2:rw+1, 2:rh+1)
                        
                        ! Add image date to the global array !
                        do j = 1, rw
                            do k = 1, rh
                                global_density(rx_start + j - 1, ry_start + k - 1) = local_density(j + 1, k + 1)
                                global_velocity(rx_start + j - 1, ry_start + k - 1) = local_velocity(j + 1, k + 1)
                            end do
                        end do
                    end block
                end do
                
                ! Write the global arrays to file (NOTE: files are currently appended with "-sliding-lid", change this for different simulations so the python code knows which files to use) !
                write(interval_str, '(I0)') step_num
                file_name = "./Visualization/output-" // trim(adjustl(interval_str)) // "-sliding-lid.txt"
                open(1, file=file_name, status="replace", action="write")
                do j = 1, global_width
                    do k = 1, global_height
                        write(1, '(I0, A, I0, A, G15.7, A, G15.7, A, G15.7)') &
                            j, ", ", k, ", ", global_density(j,k), ", ", global_velocity(j,k)%x, ", ", global_velocity(j,k)%y
                    end do
                end do
                close(1)
                
                deallocate(global_density, global_velocity)
                
            else
                ! Other ranks just send their lattices !
                call MPI_Send(lattice, directions * instance_width * instance_height, MPI_DOUBLE_PRECISION, 0, 0, comm, ierr)
            end if
            
        end subroutine gather_and_write

        ! Runs the simulation for 1000 time steps without any file I/O and calculates the MLUPS, printing the results !
        subroutine do_parallel_performance_test(lattice_even, lattice_odd, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
            use mpi

            ! Input variable declarations !
            real(8), intent(inout) :: lattice_even(instance_width, instance_height, directions)
            real(8), intent(inout) :: lattice_odd(instance_width, instance_height, directions)
            integer, intent(in) :: comm_cart, n_n, n_s, n_e, n_w
            integer, intent(in) :: coords(2), dims(2)

            ! MPI variables !
            integer :: rank, ierr
            
            ! Iterator !
            integer :: step_idx

            ! Timer variables !
            real(8) :: start_time, end_time, t_elapsed, max_t_elapsed
            
            ! MLUPS calculation variables !
            real (8) :: mlups
            integer(8) :: total_nodes, total_steps

            ! Get current rank !
            call MPI_Comm_rank(comm_cart, rank, ierr)

            total_steps = 1000

            ! Total fluid nodes !
            total_nodes = int(global_width, 8) * int(global_height, 8)

            ! Make sure all images are synced before starting timer !
            !$acc wait
            call MPI_Barrier(comm_cart, ierr)
            
            ! Start timer !
            start_time = MPI_Wtime()

            ! Run for 1000 steps !
            do step_idx = 1, total_steps

                ! If total_steps is odd, use lattice_even as input data and lattice_odd as output lattice, if total_steps is even then vice versa !
                if (mod(step_idx, 2) /= 0) then
                    call perform_one_time_step(lattice_even, lattice_odd, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
                else
                    call perform_one_time_step(lattice_odd, lattice_even, comm_cart, n_n, n_s, n_e, n_w, coords, dims)
                end if
            end do

            !$acc wait
            
            ! Stop timer !
            end_time = MPI_Wtime()

            ! Calculate time elapsed on current image !
            t_elapsed = end_time - start_time

            ! Gets maximum time elapsed on all images !
            call MPI_Reduce(t_elapsed, max_t_elapsed, 1, MPI_REAL8, MPI_MAX, 0, comm_cart, ierr)

            ! Calculate elapsed time and MLUPS (only on image 1 to avoid duplicate printing) !
            if (rank == 0) then
                mlups = (real(total_nodes, 8) * real(total_steps, 8)) / (max_t_elapsed * 1.0e6_8)
                
                print *, "Total time (s): ", max_t_elapsed
                print *, "MLUPS: ", mlups
            end if

        end subroutine do_parallel_performance_test

end program LatticeBoltzmannOpenACCMain
