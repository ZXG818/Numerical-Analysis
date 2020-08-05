    module NUMERICAL
    implicit none
        private zero
        real, parameter :: zero = 0.0001
        
    contains
        real function Newton(a, f, df)
        implicit none   
            real, intent(inout) :: a
            real, external :: f         ! 需要求解的函数
            real, external :: df        ! 需要求解的函数的导数
            
            real :: b
            real :: FB
            
            b = a - f(a)/df(a)
            FB = f(b)
            
            do while(abs(FB) > zero)
                a = b
                b = a - f(a)/df(a)
                FB = f(b)
            end do
            
            Newton = b
            return
        end function
        
        ! 求值函数
        real function func(x)
        implicit none
            real, intent(in) :: x
            func = sin(x)
            return
        end function
        
        ! 求值函数的导函数
        real function dfunc(x)
        implicit none
            real, intent(in) :: x
            dfunc = cos(x)
            return
        end function
    end module
    
    program main
    use NUMERICAL
    implicit none
        real :: a
        write(*, *) "请输入起始值:"
        read(*, *) a
        write(*, *) Newton(a, func, dfunc)
    end program