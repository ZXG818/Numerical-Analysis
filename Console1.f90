! 使用梯形法计算积分
    module INTEGRAL
    implicit none 
        real, parameter :: PI = 3.1415926
    contains
        ! 生成离散的点数据
        subroutine GenerateData(sample, lower, upper, func)
        implicit none
            real, intent(inout) :: sample(:)
            real, intent(in) :: lower, upper ! 积分上下界
            real, external :: func           ! 被积函数
            real :: r, width
            integer :: i, n
            n = size(sample, 1)
            width = (upper-lower) / (n-1)
            r = lower
            do i=1, n
                sample(i) = func(r)
                r = r + width
            end do
            return
        end subroutine
        ! 使用梯形法计算积分
        real function Trape_Integral(sample, lower, upper)
        implicit none
            real, intent(inout) :: sample(:)
            real, intent(in) :: lower, upper
            real :: width
            real :: sum = 0
            integer :: i, n
            n = size(sample, 1)
            width = (upper - lower) / (n - 1)
            do i=1, n-1
                sum = sum + (sample(i) + sample(i+1))*width/2 ! 计算梯形面积
            end do
            Trape_Integral = sum
            return 
        end function
    end module
    
    program main
    use INTEGRAL
    implicit none
        real :: sample(10000001) ! 计划采取10000001个点
        real :: upper = PI
        real :: lower = 0
        real :: ans
        real, intrinsic :: sin
        call GenerateData(sample, lower, upper, sin) !计算0到PI的sin的积分
        ans = Trape_Integral(sample, lower, upper)
        write(*, "(A6, F20.12)") 'ans=', ans
    end program
    
    ! 计算结果：  ans=      2.000000000000