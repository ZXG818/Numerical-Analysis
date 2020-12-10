    module QR
    implicit none
        private ERR
        real(8) :: ERR = 0.000001
    contains
        !矩阵单位化
        subroutine IdentityMatrix(matrix)
        implicit none
            real(8), intent(inout) :: matrix(:, :)
            integer                :: N
            integer                :: i, j
            N = size(matrix, 1)
            do i=1, N
                do j=1, N
                    if(i .eq. j) then
                        matrix(i, j) = 1.D0
                    else
                        matrix(i, j) = 0.D0
                    end if
                end do
            end do
        end subroutine
        !将向量单位化
        subroutine IdentityVector(vector)
        implicit none
            real(8), intent(inout) :: vector(:)
            integer                :: i, j
            integer                :: N
            N = size(vector, 1)
            vector(1) = 1
            do i=2, N
                vector(i) = 0
            end do
        end subroutine
        !向量转置乘法
        subroutine VectorInverseMul(matrix, vector)
        implicit none
            real(8), intent(inout) :: matrix(:, :)
            real(8), intent(in)    :: vector(:)
            integer                :: N
            integer                :: i, j
            N = size(vector, 1)
            do i=1, N
                do j=1, N
                    matrix(i, j) = vector(i) * vector(j)
                end do
            end do
        end subroutine
        !开始豪斯霍尔德变换的QR分解
        subroutine HouseHolder(matrix, Q, R)
        implicit none
            real(8), target, intent(inout) :: matrix(:, :)
            real(8), intent(inout)         :: Q(:, :)
            real(8), target, intent(inout) :: R(:, :)
            real(8), target, allocatable   :: H(:, :)
            real(8), target, allocatable   :: I(:, :)
            real(8), allocatable           :: temp(:, :)   !临时矩阵
            real(8), target, allocatable   :: u(:)
            real(8), target, allocatable   :: e(:)
            real(8), pointer               :: a_ptr(:) => NULL()
            real(8), pointer               :: u_ptr(:) => NULL()
            real(8), pointer               :: e_ptr(:) => NULL()
            real(8), pointer               :: I_ptr(:, :) => NULL()
            real(8), pointer               :: H_ptr(:, :) => NULL()
            real                           :: delta
            integer                        :: k
            integer                        :: N
            N = size(matrix, 1)
            allocate(u(N))  
            allocate(e(N))
            allocate(I(N, N))
            allocate(H(N, N))
            allocate(temp(N, N))
            R = matrix
            call IdentityMatrix(Q)
            !开始进行变换，进行N-1次变换
            do k=1, N-1
                call IdentityMatrix(H)
                a_ptr => R(k:N, k)
                e_ptr => e(k:N)
                u_ptr => u(k:N)
                I_ptr => I(k:N, k:N)
                H_ptr => H(k:N, k:N)
                !向量单位化
                call IdentityVector(e_ptr)
                !矩阵单位化
                call IdentityMatrix(I_ptr)
                !计算
                delta = sqrt(sum(a_ptr*a_ptr))
                e_ptr = delta * e_ptr
                u_ptr = (a_ptr - e_ptr) / sqrt(sum((a_ptr-e_ptr)*(a_ptr-e_ptr)))
                call VectorInverseMul(temp, u_ptr)
                H_ptr = I_ptr - 2*temp
                R = matmul(H, R)
                Q = matmul(H, Q)
            end do
            deallocate(u)
            deallocate(e)
            deallocate(I)
            deallocate(H)
            deallocate(temp)
            !将Q转置
            Q = transpose(Q)
        end subroutine
    end module
    !
    !MAIN PROGRAM
    !
    program main
    use QR
    implicit none
        real(8) :: matrix(3, 3)
        real(8) :: Q(3, 3) = 0
        real(8) :: R(3, 3) = 0
        integer :: i
        matrix(1, :) = (/ 3, 14, 9 /)
        matrix(2, :) = (/ 6, 43, 3 /)
        matrix(3, :) = (/ 6, 22, 15 /)
        call HouseHolder(matrix, Q, R)
        do i=1, 3
            write(*, *) Q(i, :)
        end do
        write(*, *) "***************"
        do i=1, 3
            write(*, *) R(i, :)
        end do
    end program
    