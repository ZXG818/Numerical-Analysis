    !作者：ZXG
    !组织：Xi`an Jiaotong University - NuTHeL
    !日期：2020-9-29
    !版本：version 1.0
    !列主元高斯消元法求解矩阵
    !要求：该矩阵为方阵且可逆
    !特殊声明：A(M, N)代表M行N列矩阵
    module COL_MAIN
    implicit none
        integer, parameter :: N = 3
    contains
        !打印矩阵
        subroutine PrintMatrix(matrix)
        implicit none
            real, intent(in) :: matrix(:, :)
            integer          :: i
            do i=1, N
                write(*, "(3F12.4)") matrix(i, :)
            end do
            return
        end subroutine
        !列主元高斯消元法
        subroutine Col_Main_Gauss(matrix, b, x)
        implicit none
            real, intent(inout) :: matrix(:, :)
            real, intent(out)   :: x(:)
            real, intent(inout) :: b(:)
            real                :: temp(N)            !用于交换主元所在行时的临时变量
            real                :: big
            real                :: coeff              !记录高斯消元时的系数
            integer             :: i, j
            integer             :: index              !用于记录主元所在的行号
            logical             :: flag = .false.     !看是否需要交换主元行
            real, intrinsic     :: abs
            do i=1, N-1
                !选取列主元
                big = abs(matrix(i,i))
                do j=i+1, N
                    if(big < abs(matrix(j, i))) then
                        index = j
                        big = matrix(j, i)
                        flag = .true.
                    end if
                end do
                !交换主元所在的行
                if(flag == .true.) then
                    temp(:) = matrix(i, :)
                    matrix(i, :) = matrix(index, :)
                    matrix(index, :) = temp(:)
                    !交换等号右边列向量的行
                    temp(1) = b(i)
                    b(i) = b(index)
                    b(index) = temp(1)
                    !恢复flag的状态
                    flag = .false.                    
                end if
                !接下来进行高斯消元
                do j=i+1, N
                    coeff = matrix(j, i) / matrix(i, i)
                    !进行行变换
                    matrix(j, :) = matrix(j, :) - matrix(i, :) * coeff
                    !同时进行等号右边列向量的行变换
                    b(j) = b(j) - b(i) * coeff
                end do
            end do
            !打印一下中间矩阵
            write(*, *) "变换后的系数矩阵："
            call PrintMatrix(matrix)
            write(*, *) "变换后，等号右边的列向量："
            write(*, "(3F12.4)") b
            !接下来进行高斯消元的回代过程
            !先求解最后一个元素的值
            x(N) = b(N) / matrix(N, N)
            do i=N-1, 1, -1
                do j=N, i+1, -1
                    b(i) = b(i) - matrix(i, j)*x(j)
                end do
                x(i) = b(i) / matrix(i, i)
            end do
            return
        end subroutine
    end module
    !主函数
    program main
    use COL_MAIN
    implicit none
        real :: matrix(N, N)
        real :: b(N) = (/ 0, 3, 2 /)
        real :: x(N)
        matrix(1, :) = (/ 1, 2, 1 /)
        matrix(2, :) = (/ 2, 2, 3 /)
        matrix(3, :) = (/ -1, -3, 0 /)
        call Col_Main_Gauss(matrix, b, x)
        write(*, *) "方程的解为："
        write(*, "(3F12.4)") x
    end program
    