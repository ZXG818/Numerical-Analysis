    !作者：ZXG
    !组织：Xi`an Jiaotong University - NuTHeL
    !日期：2020-9-29
    !版本：version 1.0
    !高斯消元法求解矩阵
    !要求：该矩阵为方阵且可逆
    !特殊声明：A(M, N)代表M行N列矩阵
    module GaussMethod
    implicit none
        integer, parameter :: N = 3
    contains
        !打印矩阵
        subroutine PrintMatrix(matrix)
        implicit none
            real, intent(in) :: matrix(N, N)
            integer          :: i
            !integer :: row
            !row = size(matrix, 1)
            do i=1, N
                write(*, "(3F10.4)") matrix(i,:)
            end do
            return 
        end subroutine
        ! 进行高斯消元法的运算，并解方程组
        subroutine Gauss(matrix, x, b)
        implicit none
            real, intent(inout) :: matrix(:, :)
            real, intent(out)   :: x(:)
            real, intent(inout) :: b(:)
            real                :: coeff !用于计算每一次进行行变换时的系数
            integer             :: i, j
            !integer :: row
            !row = size(matrix, 1)
            !开始进行高斯消元法
            !由于Fortran的数组可以直接进行行列操作，所以可以少写一个循环
            do i=1, N-1
                do j=i+1, N
                    !计算每次行变换的系数
                    coeff = matrix(j, i) / matrix(i, i)
                    !开始进行行变换
                    matrix(j, :) = matrix(j, :) - matrix(i, :) * coeff
                    b(j) = b(j) - b(i) * coeff
                end do
            end do
            !打印一下变换后的矩阵
            write(*, *) "变换后的系数矩阵："
            call PrintMatrix(matrix)
            write(*, *) "变换后，等号右边的列向量："
            write(*, "(3F10.4)") b
            !TODO:开始进行回代过程来解方程组
            !先求解最后一个元素的x
            x(N) = b(N) / matrix(N, N)
            do i=N-1, 1, -1
                do j=N, i+1, -1
                    b(i) = b(i) - matrix(i, j) * x(j)
                end do
                x(i) = b(i) / matrix(i, i)
            end do
            return
        end subroutine
    end module
    !主函数
    program main
    use GaussMethod
    implicit none
        real :: matrix(N, N)
        real :: x(N)
        real :: b(N) = (/ 0, 3, 2 /)
        !开始给矩阵赋值
        matrix(1, :) = (/  1,  2,  1 /)
        matrix(2, :) = (/  2,  2,  3 /)
        matrix(3, :) = (/ -1, -3,  0 /)
        call Gauss(matrix, x, b)
        write(*, *) "求解的x向量："
        write(*, "(3F10.4)") x
    end program
    
            