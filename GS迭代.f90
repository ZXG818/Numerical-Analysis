!  GS_iteration.f90 
!  Author: Zhou Xingguang 3120103311
!  Organization: Xi'an Jiaotong University - NuTHeL

    module GS
    implicit none
        real(8) :: stop_step = 10
        real(8) :: epsilon = 1.0D-6
        
    contains
        subroutine IterSolver(matrix, b, x)
        implicit none
            real(8), intent(inout) :: matrix(:, :)
            real(8), intent(inout) :: b(:)
            real(8), intent(inout) :: x(:)
            real(8), allocatable   :: x_front(:)
            real(8), allocatable   :: x_temp(:)
            integer                :: N
            integer                :: i, j
            N = size(matrix, 1)
            allocate(x_front(N))
            allocate(x_temp(N))
            
            x = 0
            x_front = x
            x_temp = x
            !begin the iteration
            do i=1, stop_step
                write(*, *) "step", i
                write(*, *) x
                do j=1, N
                    x(j) = (b(j)-sum(matrix(j, :) * x_temp) + matrix(j, j) * x_temp(j)) / matrix(j, j)
                    x_temp(j) = x(j)
                end do
                ! Stop condition
                ! Try to judge the 2-Norm of the subtract vector 
                ! between the solution vector and its front vector.
                if(sqrt(sum((x - x_front) * (x - x_front))) .LT. epsilon) then
                    write(*, *) "iteration stop:", x
                    exit
                end if
                x_front = x
            end do
            deallocate(x_front)
            deallocate(x_temp)
            return 
        end subroutine
    end MODULE
    
    ! Main Program
    program main
    use GS
    USE Jacobi
    implicit none
        real(8) :: matrix(3, 3);
        real(8) :: b(3) = (/ 1, 3, 5 /)
        real(8) :: x(3)
        matrix(1, :) = (/ 1, 2, -2 /)
        matrix(2, :) = (/ 1, 1,  1 /)
        matrix(3, :) = (/ 2, 2,  1 /)
        call IterSolver(matrix, b, x)
        write(*, *) "The final result:", x
    end program
    