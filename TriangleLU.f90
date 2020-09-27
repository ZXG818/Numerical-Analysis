    ! 追赶法分解三对角矩阵
    module TridiagonalMatrix
    implicit none
        ! 本结构体在该程序中用处不大，值得一提的是，
        ! 在稀疏矩阵的运算过程中，可以使用该结构体来记录那些不为零的元素。
        type :: point
            integer :: i
            integer :: j
            real    :: value
        end type
    contains
        ! 为了使用矩阵的列相关操作，举例：定义matrix(3, 4)为4行3列矩阵
        subroutine PrintMatrix(matrix)
        implicit none
            real, intent(in) :: matrix(:, :)
            integer          :: i
            integer          :: row
            row = size(matrix, 2)
            ! 按行打印矩阵
            do i=1, row
                write(*, *) matrix(:, i)
            end do
            return
        end subroutine
        ! 开始进行追赶法进行LU分解
        subroutine Run_Chase(matrix, L, U, b)
        implicit none
            real, intent(in)    :: matrix(:, :)
            real, intent(out)   :: L(:, :)
            real, intent(out)   :: U(:, :)
            real, intent(inout) :: b(:)   ! 方程等号右边的列矩阵
            real, allocatable   :: y(:)
            integer             :: i, j
            integer             :: row
            integer             :: col
            row = size(matrix, 2)
            ! 为中间方程结果y申请内存空间
            allocate(y(row))
            ! 先计算矩阵的第一个数
            U(1, 1) = matrix(1, 1)
            L(1, 1) = 1.0
            y(1) = b(1)
            ! 开始进行L和U的求解
            do i=2, row
                L(i-1, i) = matrix(i-1, i) / U(i-1, i-1)
                L(i, i) = 1
                U(i, i) = matrix(i, i) - L(i-1, i) * matrix(i, i-1)
                U(i, i-1) = matrix(i, i-1)
                y(i) = b(i) - L(i-1, i) * y(i-1)
            end do
            ! 将y计算的值返回给b
            b = y
            ! 释放y申请的空间
            deallocate(y)
            return
        end subroutine
        ! 利用高斯消元法中的回代过程完成方程组的最终求解
        subroutine Solution(U, b, x)
        implicit none
            real, intent(inout) :: U(:, :)
            real, intent(inout) :: b(:)
            real, intent(out)   :: x(:)  ! 用于存放最终结果
            integer             :: row, col, i, j
            row = size(U, 2)
            col = size(U, 1)
            x(row) = b(row) / U(row, row)
            do i=row-1, 1, -1
                do j=col, i+1, -1
                    b(i) = b(i) - U(j, i) * x(j)
                end do
                x(i) = b(i) / U(i, i)
            end do
            return
        end subroutine
    end module
    ! 开始主程序
    program main
    use TridiagonalMatrix
    implicit none
        integer :: i, j
        real :: matrix(5, 5) = (/ 1, 2, 0, 0, 0, 2, 3, 1, 0, 0, 0, -3, 4, 2, 0, 0, 0, 4, 7, 1, 0, 0, 0, -5, 6 /)
        real :: L(5, 5) = 0
        real :: U(5, 5) = 0
        real :: b(5) = (/ 5, 9, 2, 19, -4 /)
        real :: x(5) = 0
        integer :: row = 5
        call Run_Chase(matrix, L, U, b)
        
        write(*, *) "三对角矩阵的LU分解如下所示："
        write(*, *) "==============================L矩阵=============================="
        call PrintMatrix(L)
        write(*, *) "==============================U矩阵=============================="
        call PrintMatrix(U)
        write(*, *) "============================中间b向量============================="
        write(*, *) b
        ! 解方程
        call Solution(U, b, x)
        write(*, *) "方程的最终结果："
        write(*, *) x
    end program
    