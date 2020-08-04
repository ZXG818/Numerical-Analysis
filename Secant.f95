    ! 此为Fortran95语法的，利用割线法求解函数根的FORTRAN程序；
    ! 该程序由trace-shadow / watermelon / ZXG进行编写，实验函数为sin(x)，实验结果良好。
    module NUMERICAL
    implicit none
        private zero        ! 此为私有变量
        private GetPoint    ! 此函数为私有函数，进行封装
        real, parameter :: zero = 0.00001
    contains
    
        real function GetPoint(a, b, func)
        implicit none
            real, intent(in) :: a
            real, intent(in) :: b
            real, external :: func
            
            real :: FA
            real :: FB
            real :: k
            
            FA = func(a)
            FB = func(b)
            k = (FB-FA) / (b-a)
            
            GetPoint = a - FA/k
            return 
        end function
    

        real function Solver(a, b, func)
        implicit none
            ! 声明参数
            real, intent(inout) :: a
            real, intent(inout) :: b
            real, external :: func  ! 声明外部函数，进行调用
            ! 声明函数内临时变量
            real :: c
            real :: FA
            real :: FB
            real :: FC
            
            c = GetPoint(a, b, func)
            FA = func(a)
            FB = func(b)
            FC = func(c)
            
            do while(abs(FC - 0) >= zero)
                a = b
                FA = FB
                b = c
                FB = FC
                c = GetPoint(a, b, func)
                FC = func(c)
            end do
            
            Solver = c
            return
        end function
        
        real function func(X)
        implicit none
            real, intent(in) :: X
            func = sin(X)
            return
        end function
    end module
    
    program main
    use NUMERICAL
    implicit none
        real :: a, b
        real :: ans
        write(*, *) "输入两个猜测的数值:"
        read(*, *) a, b
        ans = Solver(a, b, func)
        write(*, *) ans
    end program
