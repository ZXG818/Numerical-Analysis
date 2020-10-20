! Jacobi Iteration
    module JACOBI
    implicit none
        real(8), parameter :: epsilon = 1.0D-6
        integer, parameter :: stop_step = 1e3
    contains
        ! Iteration solver
        subroutine IterSolver(matrix, b, x)
        implicit none
            real(8), intent(inout) :: matrix(:, :)
            real(8), intent(inout) :: b(:)
            real(8), intent(inout) :: x(:)
            real(8), allocatable   :: x_front(:)    ! try to record the value of the front calculation.
            integer                :: N
            integer                :: i, j
            N = size(matrix, 1)
            allocate(x_front(N))
            ! We should to try to calculate the initial value of x
            ! and, we try to initial the x = 0 here.
            x = 0
            x_front = x
            ! begin to iteration
            do i=1, stop_step
                write(*, *) "step", i
                write(*, *) x
                do j=1, N
                    x(j) = (b(j)-sum(matrix(j, :)*x_front)+matrix(j, j)*x_front(j)) / matrix(j, j)
                end do
                ! Stop condition
                ! Try to judge the 2-Norm of the subtract vector 
                ! between the solution vector and its front vector.
                if(sqrt(sum((x - x_front) * (x - x_front))) .LT. epsilon) then
                    write(*, *) "iteration stop:", x
                    exit
                end if
                ! give the value to x_front
                x_front = x
            end do
            deallocate(x_front)
            return
        end subroutine
    end module
    !
    ! MAIN
    !
    program main
    use JACOBI
    implicit none
        real(8) :: matrix(3, 3)
        real(8) :: b(3) = (/ 72, 83, 42 /)
        real(8) :: x(3)
        matrix(1, :) = (/ 10, -1, -2 /)
        matrix(2, :) = (/ -1, 10, -2 /)
        matrix(3, :) = (/ -1, -1, 5 /)
        call IterSolver(matrix, b, x)
        write(*, *) "The final result:", x
    end program
 