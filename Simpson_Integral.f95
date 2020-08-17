    ! Simpson辛普森积分法
    ! 采样点的个数必须为奇数
    ! 计算公式：Integeral(fx) = (f0 + 4*f1 + 2*f2 + 4*f3 + 2*f4 + ... + 4*fn-1 + fn) * h / 3 
    ! 其中h为程序中的width变量
    module MY_INTEGRAL
    implicit none
        real, parameter :: PI = 3.1415926
    contains
        ! 获取采样点
        subroutine Sample(s, lower, upper, func)
        implicit none
            real, intent(inout) :: s(:)
            real, intent(in) :: lower, upper
            real, external :: func
            integer :: i, N
            real :: width, r
            r = lower
            N = size(s, 1)
            width = (upper - lower) / (N-1)
            do i=1, N
                s(i) = func(r)
                r = r + width
            end do
            return
        end subroutine
        ! 利用公式进行计算积分
        real function GetIntegeral(s, lower, upper)
        implicit none
            real, intent(in) :: s(:)
            real, intent(in) :: lower, upper
            integer :: N, i
            real :: ans = 0.0
            real :: width
            N = size(s, 1)
            width = (upper-lower)/(N-1)
            ans = ans + s(1) + s(N)
            do i=2, N-2, 2
                ans = ans + 4*s(i) + 2*s(i+1)
            end do
            GetIntegeral = ans * width / 3
            return 
        end function
    end module
    
    program main
    use MY_INTEGRAL
    implicit none
        real :: s(10001)
        real :: lower=0
        real :: upper=PI
        real, intrinsic :: sin
        real :: ans = 0.0
        call Sample(s, lower, upper, sin)
        ans = GetIntegeral(s, lower, upper)
        write(*, *) ans
    end program
    ! 运行结果：1.999978 此为sin(x)在0到PI上的积分
