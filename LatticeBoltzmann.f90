module lattice_boltzmann
    implicit none

    ! Type used for velocity vectors !
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
    real(8), parameter :: pi = 4.0_8 * atan(1.0_8)
    real(8), parameter :: epsilon = 0.01
    
    ! Wall velocity, used when simulation has a moving wall !
    type(velocity) :: wall_speed = velocity(0.1_8, 0.0_8)

    ! Input and output pressure used for Poiseuille flow !
    real(8), parameter :: poiseuille_in_pressure = 0.3, poiseuille_out_pressure = 0.29

    ! Width of coarray grid !
    integer :: coarray_dimensions = 1

    ! Encoding for boundary configuration. 0 is normal periodic boundaries, 1 is Couette flow, 2 is Poiseuille flow, 3 is sliding lid !
    integer :: boundary_configuration = 0

    contains

        ! Perform one time step of the LBM simulation !
        subroutine perform_one_time_step(lattice_in, lattice_out)

            ! Input variable declarations ! 
            real(8), intent(inout) :: lattice_in(directions, instance_width, instance_height)[coarray_dimensions,*]
            real(8), intent(inout) :: lattice_out(directions, instance_width, instance_height)[coarray_dimensions,*]

            ! Update all ghost nodes in the input array !
            call sync_ghost_nodes(lattice_in)

            ! Apply boundary conditions to ghost nodes !
            if (boundary_configuration == 0) then
                call periodic_boundary(lattice_in)
            else if (boundary_configuration == 1) then
                call couette_boundary(lattice_in)
            else if (boundary_configuration == 2) then
                call poiseuille_boundary(lattice_in)
            else if (boundary_configuration == 3) then
                call sliding_lid_boundary(lattice_in)
            end if

            ! Combined streaming and collision (pulls from lattice_in, computes, writes to lattice_out) !
            call stream_and_collide(lattice_in, lattice_out)

        end subroutine perform_one_time_step

        ! Synchronize ghost nodes across all images !
        subroutine sync_ghost_nodes(lattice)
            
            ! Input variable declaration !
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]
            
            ! Image iteration variables !
            integer :: img(2), total_images, img_grid_y

            ! Get information about the image grid !
            img = this_image(lattice)
            total_images = num_images()
            img_grid_y = (total_images + coarray_dimensions - 1) / coarray_dimensions
            
            ! Send left inner edge to Western neighbor's right ghost edge !
            if (img(1) /= 1) then
                lattice(:, instance_width, :)[img(1) - 1, img(2)] = lattice(:, 2, :)
            end if

            ! Send right inner edge to Eastern neighbor's left ghost edge !
            if (img(1) /= coarray_dimensions) then
                lattice(:, 1, :)[img(1) + 1, img(2)] = lattice(:, instance_width - 1, :)
            end if

            ! Sync so that corner nodes are not lost !
            sync all
            
            ! Send bottom inner edge to Southern neighbor's top ghost edge !
            if (img(2) /= img_grid_y) then
                lattice(:, :, 1)[img(1), img(2) + 1] = lattice(:, :, instance_height - 1)
            end if

            ! Send top inner edge to Northern neighbor's bottom ghost edge !
            if (img(2) /= 1) then
                lattice(:, :, instance_height)[img(1), img(2) - 1] = lattice(:, :, 2)
            end if

            ! Final sync between images !
            sync all

        end subroutine sync_ghost_nodes

        ! Combined streaming and collision step !
        subroutine stream_and_collide(lattice_in, lattice_out)

            ! Input variable declarations !
            real(8), intent(in) :: lattice_in(directions, instance_width, instance_height)
            real(8), intent(inout) :: lattice_out(directions, instance_width, instance_height)

            ! Iterators and data variables !
            integer :: j, k
            real(8) :: f0, f1, f2, f3, f4, f5, f6, f7, f8
            real(8) :: density, vel_x, vel_y, u_sq, c_dot_u

            ! Loop through the lattice !
            do k = 2, instance_height - 1
                do j = 2, instance_width - 1
                    
                    ! Take data that will "stream in" from neighbors in input lattice !
                    f0 = lattice_in(1, j, k)
                    f1 = lattice_in(2, j - 1, k)     ! x: 1, y: 0   -> pulled from left node !
                    f2 = lattice_in(3, j, k + 1)     ! x: 0, y: 1   -> pulled from bottom node !
                    f3 = lattice_in(4, j + 1, k)     ! x:-1, y: 0   -> pulled from right node !
                    f4 = lattice_in(5, j, k - 1)     ! x: 0, y:-1   -> pulled from top node !
                    f5 = lattice_in(6, j - 1, k + 1) ! x: 1, y: 1   -> pulled from bottom left node !
                    f6 = lattice_in(7, j + 1, k + 1) ! x:-1, y: 1   -> pulled from bottom rght node !
                    f7 = lattice_in(8, j + 1, k - 1) ! x:-1, y:-1   -> pulled from top right node !
                    f8 = lattice_in(9, j - 1, k - 1) ! x: 1, y:-1   -> pulled from top left node !

                    ! Values used for collision calculations !
                    density = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8
                    vel_x = (f1 - f3 + f5 - f6 - f7 + f8) / density
                    vel_y = (f2 - f4 + f5 + f6 - f7 - f8) / density
                    u_sq = vel_x**2 + vel_y**2

                    ! Do collision calculations and save to output lattice !
                    lattice_out(1, j, k) = f0 - omega * (f0 - weights(1) * density * (1.0_8 - 1.5_8 * u_sq))
                    
                    c_dot_u = vel_x
                    lattice_out(2, j, k) = f1 - omega * (f1 - weights(2) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = vel_y
                    lattice_out(3, j, k) = f2 - omega * (f2 - weights(3) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = -vel_x
                    lattice_out(4, j, k) = f3 - omega * (f3 - weights(4) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = -vel_y
                    lattice_out(5, j, k) = f4 - omega * (f4 - weights(5) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = vel_x + vel_y
                    lattice_out(6, j, k) = f5 - omega * (f5 - weights(6) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = -vel_x + vel_y
                    lattice_out(7, j, k) = f6 - omega * (f6 - weights(7) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = -vel_x - vel_y
                    lattice_out(8, j, k) = f7 - omega * (f7 - weights(8) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))
                    
                    c_dot_u = vel_x - vel_y
                    lattice_out(9, j, k) = f8 - omega * (f8 - weights(9) * density * (1.0_8 + 3.0_8 * c_dot_u + 4.5_8 * c_dot_u**2 - 1.5_8 * u_sq))

                end do
            end do

        end subroutine stream_and_collide

        ! Standard periodic boundary conditions !
        subroutine periodic_boundary(lattice)

            ! Input variable declariation !
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]
            
            ! Image iteration variables !
            integer :: img(2), total_images, img_grid_y

            ! Get information about image grid !
            img = this_image(lattice)
            total_images = num_images()
            img_grid_y = (total_images + coarray_dimensions - 1) / coarray_dimensions
            
            ! If on the left edge of the global grid, pull from the right edge !
            if (img(1) == 1) then
                lattice(:, 1, :) = lattice(:, instance_width - 1, :)[coarray_dimensions, img(2)]
            end if

            ! If on the right edge of the global grid, pull from the left edge !
            if (img(1) == coarray_dimensions) then
                lattice(:, instance_width, :) = lattice(:, 2, :)[1, img(2)]
            end if

            ! Sync so that corner nodes are not lost !
            sync all
            
            ! If on the top edge of the global grid, pull from the bottom edge !
            if (img(2) == 1) then
                lattice(:, :, 1) = lattice(:, :, instance_height - 1)[img(1), img_grid_y]
            end if

            ! If on the bottom edge of the global grid, pull from the top edge !
            if (img(2) == img_grid_y) then
                lattice(:, :, instance_height) = lattice(:, :, 2)[img(1), 1]
            end if

            ! Final sync between images !
            sync all

        end subroutine periodic_boundary

        ! Boundary conditions for couette flow simulations !
        subroutine couette_boundary(lattice)

            ! Input variable declaration !
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]
            
            ! Image iteration variables !
            integer :: img(2), total_images, img_grid_y
            real(8) :: average_density

            ! Get information about image grid !
            img = this_image(lattice)
            total_images = num_images()
            img_grid_y = (total_images + coarray_dimensions - 1) / coarray_dimensions
            
            ! Calculate average density !
            average_density = sum(lattice) / (instance_width * instance_height)

            ! Handle periodic boundaries on left and right !
            if (img(1) == 1) then
                lattice(:,1,:) = lattice(:, instance_width-1, :)[coarray_dimensions, img(2)]
            end if
            if (img(1) == coarray_dimensions) then
                lattice(:,instance_width,:) = lattice(:, 2, :)[1, img(2)]
            end if
            
            ! Ensure periodic bounds are finished !
            sync all 

            ! Bottom bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            if (img(2) == img_grid_y) then
                lattice(3,:,instance_height) = lattice(5,:,instance_height-1)
                lattice(6,1:instance_width-1,instance_height) = lattice(8,2:instance_width,instance_height-1)
                lattice(7,2:instance_width,instance_height) = lattice(9,1:instance_width-1,instance_height-1)
            end if

            ! Top moving boundary (indexes are +1 since fortran is 1-indexed) !
            if (img(2) == 1) then
                lattice(5,:,1) = lattice(3,:,2)
                lattice(8,2:instance_width,1) = lattice(6,1:instance_width-1,2) - 2.0_8 * weights(6) * average_density * ((shift_directions_x(6) * wall_speed%x) / (1.0_8 / 3.0_8))
                lattice(9,1:instance_width-1,1) = lattice(7,2:instance_width,2) - 2.0_8 * weights(7) * average_density * ((shift_directions_x(7) * wall_speed%x) / (1.0_8 / 3.0_8))
            end if
        end subroutine couette_boundary

        ! Boundary conditions for poiseuille flow simulations !
        subroutine poiseuille_boundary(lattice)

            ! Input variable declaration !
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]
            
            ! Image iteration variables !
            integer :: img(2), total_images, img_grid_y

            ! Calculation variables !
            real(8) :: density_in, density_out
            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)
            real(8) :: f_star_east_minus_feq(directions, 1, instance_height), f_star_west_minus_feq(directions, 1, instance_height)
            real(8) :: equilibriums(directions, instance_width, instance_height)

            ! Get information about image grid !
            img = this_image(lattice)
            total_images = num_images()
            img_grid_y = (total_images + coarray_dimensions - 1) / coarray_dimensions

            ! Values used in pressure differential calculation !
            density_in = poiseuille_in_pressure / (1.0_8 / 3.0_8)
            density_out = poiseuille_out_pressure / (1.0_8 / 3.0_8)

            ! Initial measurements !
            density_arr = calculate_density_array(lattice)
            velocity_arr = calculate_average_velocity_array(lattice, density_arr)
            equilibriums = calculate_equilibrium(density_arr, velocity_arr)

            ! Calculate f_star (current array) - initial equilibrium !
            if (img(1) == coarray_dimensions) then
                f_star_east_minus_feq(:,1,:) = lattice(:, instance_width-1, :) - equilibriums(:, instance_width-1, :)
            end if
            if (img(1) == 1) then
                f_star_west_minus_feq(:,1,:) = lattice(:, 2, :) - equilibriums(:, 2, :)
            end if

            ! Calculate target velocities and densities on the boundaries !
            if (img(1) == coarray_dimensions) then
                velocity_arr(instance_width-1,:)%x = velocity_arr(instance_width-1,:)%x * (density_arr(instance_width-1,:) / density_in)
                velocity_arr(instance_width-1,:)%y = velocity_arr(instance_width-1,:)%y * (density_arr(instance_width-1,:) / density_in)
                density_arr(instance_width-1,:) = density_in
            end if
            if (img(1) == 1) then
                velocity_arr(2,:)%x = velocity_arr(2,:)%x * (density_arr(2,:) / density_out)
                velocity_arr(2,:)%y = velocity_arr(2,:)%y * (density_arr(2,:) / density_out)
                density_arr(2,:) = density_out
            end if

            ! Re-evaluate equilibrium with updated target densities !
            equilibriums = calculate_equilibrium(density_arr, velocity_arr)

            ! Apply periodic boundary with pressure differential !
            if (img(1) == coarray_dimensions) then
                ! Push right boundary to the left-most image's left ghost node
                lattice(:, 1, :)[1, img(2)] = equilibriums(:, instance_width-1, :) + f_star_east_minus_feq(:,1,:)
            end if

            if (img(1) == 1) then
                ! Push left boundary to the right-most image's right ghost node
                lattice(:, instance_width, :)[coarray_dimensions, img(2)] = equilibriums(:, 2, :) + f_star_west_minus_feq(:,1,:)
            end if
            
            sync all

            ! Top bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            if (img(2) == 1) then
                lattice(5,:,1) = lattice(3,:,2)
                lattice(8,2:instance_width,1) = lattice(6,1:instance_width-1,2)
                lattice(9,1:instance_width-1,1) = lattice(7,2:instance_width,2)
            end if

            ! Bottom bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            if (img(2) == img_grid_y) then
                lattice(3,:,instance_height) = lattice(5,:,instance_height-1)
                lattice(6,1:instance_width-1,instance_height) = lattice(8,2:instance_width,instance_height-1)
                lattice(7,2:instance_width,instance_height) = lattice(9,1:instance_width-1,instance_height-1)
            end if
        end subroutine poiseuille_boundary

        ! Boundary conditions for sliding lid simulations !
        subroutine sliding_lid_boundary(lattice)

            ! Input variable declaration !
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]

            ! Image iteration variables !
            integer :: img(2), total_images, img_grid_y

            ! Get information about image grid !
            img = this_image(lattice)
            total_images = num_images()
            img_grid_y = (total_images + coarray_dimensions - 1) / coarray_dimensions

            ! Left side bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            if (img(1) == 1) then
                lattice(2,1,:) = lattice(4,2,:)
                lattice(6,1,2:instance_height) = lattice(8,2,1:instance_height-1)
                lattice(9,1,1:instance_height-1) = lattice(7,2,2:instance_height)
            end if

            ! Right side bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            if (img(1) == coarray_dimensions) then
                lattice(4,instance_width,:) = lattice(2, instance_width-1,:)
                lattice(7,instance_width,2:instance_height) = lattice(9,instance_width-1,1:instance_height-1)
                lattice(8,instance_width,1:instance_height-1) = lattice(6,instance_width-1,2:instance_height)
            
            end if

            ! Bottom bounce-back boundary (indexes are +1 since fortran is 1-indexed) !
            if (img(2) == img_grid_y) then
                lattice(3,:,instance_height) = lattice(5,:,instance_height-1)
                lattice(6,1:instance_width-1,instance_height) = lattice(8,2:instance_width,instance_height-1)
                lattice(7,2:instance_width,instance_height) = lattice(9,1:instance_width-1,instance_height-1)

            end if

            ! Top moving boundary (y velocity is 0, so nothing is added to direction 4 (indexed as 5)) !
            if (img(2) == 1) then
                lattice(5,:,1) = lattice(3,:,2)
                lattice(8,2:instance_width,1) = lattice(6,1:instance_width-1,2) - 2.0_8 * weights(6) * sum(lattice(:, 1:instance_width-1, 2), dim=1) * ((shift_directions_x(6) * wall_speed%x) / (1.0_8 / 3.0_8))
                lattice(9,1:instance_width-1,1) = lattice(7,2:instance_width,2) - 2.0_8 * weights(7) * sum(lattice(:, 2:instance_width, 2), dim=1) * ((shift_directions_x(7) * wall_speed%x) / (1.0_8 / 3.0_8))

            end if

        end subroutine sliding_lid_boundary

        ! Helper function to calculate density !
        function calculate_density_array(lattice)

            ! Input variable declaration !
            real(8), intent(in) :: lattice(directions, instance_width, instance_height)

            ! Output variable declaration !
            real(8) :: calculate_density_array(instance_width, instance_height)

            ! Populate the whole array with 1.0 values !
            calculate_density_array = 1.0_8

            ! Calculate the density of each node and save them to a 2d array !
            calculate_density_array(2:instance_width-1, 2:instance_height-1) = sum(lattice(:, 2:instance_width-1, 2:instance_height-1), dim=1)

        end function calculate_density_array

        ! Helper function for calculating average velocity !
        function calculate_average_velocity_array(lattice, density_arr)

            ! Input variable declarations !
            real(8), intent(in) :: lattice(directions, instance_width, instance_height)
            real(8), intent(in) :: density_arr(instance_width, instance_height)

            ! Output variable declaration !
            type(velocity) :: calculate_average_velocity_array(instance_width, instance_height)

            ! Iterator !
            integer :: i

            ! Populate the whole array with 0 velocities !
            calculate_average_velocity_array%x = 0.0_8
            calculate_average_velocity_array%y = 0.0_8

            ! Loop through directions and calculate momentum per node !
            do i=1, directions
                calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%x = calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%x + lattice(i,2:instance_width-1, 2:instance_height-1) * shift_directions_x(i)
                calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%y = calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%y + lattice(i,2:instance_width-1, 2:instance_height-1) * shift_directions_y(i)
            end do

            ! Calculate the average velocity of each node !
            calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%x = calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%x / density_arr(2:instance_width-1, 2:instance_height-1)
            calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%y = calculate_average_velocity_array(2:instance_width-1, 2:instance_height-1)%y / density_arr(2:instance_width-1, 2:instance_height-1)

        end function calculate_average_velocity_array

        ! Helper function to calculate equilibrium !
        function calculate_equilibrium(density_arr, velocity_arr)

            ! Input variable declarations !
            real(8), intent(in) :: density_arr(instance_width, instance_height)
            type(velocity), intent(in) :: velocity_arr(instance_width, instance_height)

            ! Output variable declaration !
            real(8) :: calculate_equilibrium(directions, instance_width, instance_height)

            ! Iterator !
            integer :: i

            ! Calculation variables !
            real(8) :: u_squared(instance_width, instance_height), velocity_sum(instance_width, instance_height)

            ! Populate the whole array with 0.0 values !
            calculate_equilibrium = 0.0_8

            ! Total velocity value used for equilibrium calculation !
            u_squared(2:instance_width-1, 2:instance_height-1) = velocity_arr(2:instance_width-1, 2:instance_height-1)%x**2 + velocity_arr(2:instance_width-1, 2:instance_height-1)%y**2

            do i=1, directions
                ! Calculate values used in equilibrium equation !
                velocity_sum(2:instance_width-1, 2:instance_height-1) = shift_directions_x(i) * velocity_arr(2:instance_width-1, 2:instance_height-1)%x + shift_directions_y(i) * velocity_arr(2:instance_width-1, 2:instance_height-1)%y

                ! Calculate equilibrium lattice !
                calculate_equilibrium(i,2:instance_width-1, 2:instance_height-1) = weights(i) * density_arr(2:instance_width-1, 2:instance_height-1) * (1.0_8 + 3.0_8 * velocity_sum(2:instance_width-1, 2:instance_height-1) + (4.5_8 * velocity_sum(2:instance_width-1, 2:instance_height-1)**2) - 1.5_8 * u_squared(2:instance_width-1, 2:instance_height-1))
            end do

        end function calculate_equilibrium

        ! Initialize lattice with (semi)-random values !
        subroutine populate_lattice_random(lattice)

            ! Input variable declarations !
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            ! Iterators !
            integer :: i, j, k

            ! Random number !
            real(8) :: u

            call random_seed()

            ! populate with random numbers !
            do j=1, instance_width
                do k=1, instance_height
                    do i=1, directions
                        call random_number(u)
                        lattice(i,j,k) = 0.07 + (u * 0.02)
                    end do
                end do
            end do

        end subroutine populate_lattice_random

        ! Initialize lattice with dense center !
        subroutine populate_lattice_dense_center(lattice)

            ! Input variable declaration !
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]

            ! Iterators !
            integer :: j, k, global_j, global_k

            ! Image coordinates !
            integer :: img(2)

            ! Small value used for calulations !
            real(8) :: epsilon1 = 0.01

            ! Number of total nodes !
            real(8) :: global_nodes

            ! Get coordinates of current image !
            img = this_image(lattice)

            ! At one point set f(0) to 4/9 + 4/9 * epsilon, f(1,2,3,4) to 1/9 + 1/9 * epsilon, f(5,6,7,8) to 1/36 + 1/36 + epsilon. Now roh will be 1+epsilon !
            ! At all other points set f(0) to 4/9 - 4/9 * 1/(m*n-1) * epsilon, f(1,2,3,4) to 1/9 - 1/9 * 1/(m*n-1) * epsilon, f(5,6,7,8) to 1/36 - 1/36 * 1/(m*n-1) * epsilon !

            ! Get total number of nodes !
            global_nodes = real(global_width * global_height, 8)

            ! Loop through lattice !
            do j=1, instance_width
                do k=1, instance_height
                    
                    ! Map local ghost/inner nodes to global coordinates !
                    global_j = (img(1) - 1) * (instance_width - 2) + (j - 1)
                    global_k = (img(2) - 1) * (instance_height - 2) + (k - 1)

                    ! Populate global center with stated values, otherwise use other values !
                    if (global_j == global_width/2 .and. global_k == global_height/2) then
                        lattice(1,j,k) = weights(1) + weights(1) * epsilon1
                        lattice(2:5,j,k) = weights(2) + weights(2) * epsilon1
                        lattice(6:,j,k) = weights(6) + weights(6) * epsilon1
                    else
                        lattice(1,j,k) = weights(1) - weights(1) * epsilon1 * (1.0_8 / (global_nodes - 1.0_8))
                        lattice(2:5,j,k) = weights(2) - weights(2) * epsilon1 * (1.0_8 / (global_nodes - 1.0_8))
                        lattice(6:,j,k) = weights(6) - weights(6) * epsilon1  * (1.0_8 / (global_nodes - 1.0_8))
                    end if
                    
                end do
            end do

        end subroutine populate_lattice_dense_center

        ! Initialize lattice with shear wave !
        subroutine populate_lattice_shear_wave(lattice)

            ! Input variable declaration !
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)[coarray_dimensions,*]

            ! Initial density and velocity array values !
            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)

            ! Y-coordinate variables !
            integer :: y_pos, global_k

            ! Image coordinates !
            integer :: img(2)

            ! Get coordinates of current image !
            img = this_image(lattice)

            ! Populate arrays with 1.0 density and 0.0 velocity !
            density_arr = 1.0_8
            velocity_arr%x = 0.0_8
            velocity_arr%y = 0.0_8

            ! Loop through instance y-positions !
            do y_pos=2, instance_height-1
                ! Get global y-position !
                global_k = (img(2) - 1) * (instance_height - 2) + y_pos

                ! Populate x velocity with sinusoidal wave based on y-position !
                velocity_arr(:,y_pos)%x = epsilon * sin((2.0_8 * pi * (global_k-1)) / (global_height))
            end do

            ! Set input lattice to new lattice using this velocity !
            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_shear_wave

        ! Initialize lattice for Couette flow !
        subroutine populate_lattice_couette(lattice)

            ! Input variable declaration !
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            ! Initial density and velocity array values !
            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)

            ! Set boundary to 1 (couette flow) !
            boundary_configuration = 1

            ! Populate arrays with 1.0 density and 0.0 velocity !
            density_arr = 1.0_8
            velocity_arr = velocity(0.0_8, 0.0_8)

            ! Set input lattice to new lattice with 1.0 density and 0.0 velocity !
            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_couette

        ! Intialize lattice for poiseuille flow !
        subroutine populate_lattice_poiseuille(lattice)

            ! Input variable declariation !
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            ! Initial density and velocity array values !
            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)

            ! Set boundary to 2 (poiseuille flow) !
            boundary_configuration = 2

            ! Populate arrays to 1.0 density and 0.0 velocity !
            density_arr = 1.0_8
            velocity_arr = velocity(0.0_8, 0.0_8)

            ! Set input lattice to new lattice with 1.0 density and 0.0 velocity !
            lattice = calculate_equilibrium(density_arr, velocity_arr)

        end subroutine populate_lattice_poiseuille

        ! Intialize lattice for sliding lid simulation !
        subroutine populate_lattice_sliding_lid(lattice)

            ! Input variable declaration !
            real(8), intent(inout) :: lattice(directions, instance_width, instance_height)

            ! Initial density and velocity array values !
            real(8) :: density_arr(instance_width, instance_height)
            type(velocity) :: velocity_arr(instance_width, instance_height)

            ! Set boundary to 3 (sliding lid) !
            boundary_configuration = 3

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

            ! Calculated viscocity variable!
            real(8) :: analytical_viscosity

            ! Calculate viscocity !
            analytical_viscosity = (1.0_8 / 3.0_8) * ((1.0_8 / omega) - 0.5_8)
            
            ! Update wall speed to fit desired reynold's number !
            wall_speed = velocity((reynolds_number * analytical_viscosity) / (coarray_dimensions * (instance_width - 2.0_8)), 0.0_8)

        end subroutine set_lid_velocity_given_reynolds
        
end module lattice_boltzmann