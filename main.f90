module config
    implicit none
    public

    !Segment boundaries
    real(8), parameter :: a = 0.d0, b = 1.d0
    !Integration step
    real(8), parameter :: step = 1d-3
    !Threads number
    integer, parameter :: thr = 4
end module config

module com
    implicit none
    public
    real(8) :: xx
end module com

program output
    use config
    use omp_lib
    use com
    implicit none

    abstract interface
        real(8) function func(x)
            implicit none
            real(8), intent(in) :: x
        end function func
    end interface

    procedure(func), pointer :: fptr
    real(8), parameter :: pi = 4*atan(1.d0)

    fptr => intgrCircle
    !print *, "Circle:    ", 2*QG10(fptr, 0.d0, 1.d0, step), "cf.", pi/24

    fptr => ftst
    !print *, "Circle(1x):", 4*QG10(fptr, 0.d0, 1.d0, step)/3, "cf.", pi/24
    !print *, 2*ftst(25.d-2)/3, intgrCircle(25.d-2)

    fptr => intgrVert
    !print *, "Vertical:  ", 2*QG10(fptr, 0.d0, 1.d0, step), "cf.", pi/24-pi/192

    fptr => intgrHor
    !print *, "Horizontal:", 2*QG10(fptr, 0.d0, 1.d0, step), "cf.", pi/24-pi/192

    fptr => polIntgrTil
    print *, "Tilted:    ", QG10(fptr, pi/4, 3*pi/4, step), "cf.", pi/24-43*pi/3072

contains

    pure real(8) function ftst(x)
        implicit none
        real(8), intent(in) :: x
        ftst = x**2*sqrt(1-x**2)**3
        !ftst = x**2
    end function

    real(8) function QG2(fptr, min, max, step)
        implicit none

        procedure(func), pointer, intent(in) :: fptr
        real(8), intent(in) :: min, max, step

        real(8) :: x, half, mid
        real(8) :: tmp, er, next
        real(8), dimension(2) :: a
        integer(kind=8) :: i, j, n

        data a /0.577350269189626d0,-0.577350269189626d0/

        n = (max-min)/step
        half = step/2
        QG2=0
        !$omp parallel default(none) private(x, mid, next, er, tmp) shared(half, n, min, step, a, QG2, fptr) num_threads(thr)
            er = 0
            !$omp do reduction(+: QG2)
            do i = 1, n
                x = min + (i-1)*step
                mid = x + half
                next = half*fptr(a(1)*half+mid) - er
                tmp = QG2 + next
                er = (tmp - QG2) - next
                QG2 = tmp
                next = half*fptr(a(2)*half+mid) - er
                tmp = QG2 + next
                er = (tmp - QG2) - next
                QG2 = tmp
            end do
            !$omp end do
        !$omp end parallel

    end function QG2

    real(8) function QG6(fptr, min, max, step)
        implicit none

        procedure(func), pointer, intent(in) :: fptr
        real(8), intent(in) :: min, max, step

        real(8) :: x, half, mid
        real(8) :: tmp, er, next
        real(8), dimension(6) :: a, w
        integer(kind=8) :: i, j, n

        data w /0.467913934572691d0, 0.467913934572691d0, 0.360761573048139d0, &
                0.360761573048139d0, 0.171324492379170d0, 0.171324492379170d0/
        data a /0.238619186083197d0,-0.238619186083197d0, 0.661209386466265d0, &
                -0.661209386466265d0, 0.932469514203152d0, -0.932469514203152d0/

        n = (max-min)/step
        half = step/2
        QG6=0
        !$omp parallel default(none) private(x, mid, next, er, tmp) shared(half, n, min, step, a, w, QG6, fptr) num_threads(thr)
            er = 0
            !$omp do reduction(+: QG6)
            do i = 1, n
                x = min + (i-1)*step
                mid = x + half
                do j = 1, 6
                    next = half*w(j)*fptr(a(j)*half+mid) - er
                    tmp = QG6 + next
                    er = (tmp - QG6) - next
                    QG6 = tmp
                end do
            end do
            !$omp end do
        !$omp end parallel

    end function QG6

    real(8) function QG8(fptr, min, max, step)
        implicit none

        procedure(func), pointer, intent(in) :: fptr
        real(8), intent(in) :: min, max, step

        real(8) :: x, half, mid
        real(8) :: tmp, er, next
        real(8), dimension(8) :: a, w
        integer(kind=8) :: i, j, n

        data w /0.362683783378362d0, 0.362683783378362d0, 0.313706645877887d0, 0.313706645877887d0, &
                0.222381034453374d0, 0.222381034453374d0, 0.101228536290376d0, 0.101228536290376d0/
        data a /0.183434642495650d0, -0.183434642495650d0, 0.525532409916329d0, -0.525532409916329d0, &
                0.796666477413627d0, -0.796666477413627d0, 0.960289856497536d0, -0.960289856497536d0/

        n = (max-min)/step
        half = step/2
        QG8=0
        !$omp parallel default(none) private(x, mid, next, er, tmp) shared(half, n, min, step, a, w, QG8, fptr) num_threads(thr)
            er = 0
            !$omp do reduction(+: QG8)
            do i = 1, n
                x = min + (i-1)*step
                mid = x + half
                do j = 1, 8
                    next = half*w(j)*fptr(a(j)*half+mid) - er
                    tmp = QG8 + next
                    er = (tmp - QG8) - next
                    QG8 = tmp
                end do
            end do
            !$omp end do
        !$omp end parallel

    end function QG8

    real(8) function QG10(fptr, min, max, step)
        implicit none

        procedure(func), pointer, intent(in) :: fptr
        real(8), intent(in) :: min, max, step

        real(8) :: x, half, mid
        real(8) :: tmp, er, next
        real(8), dimension(10) :: a, w
        integer(kind=8) :: i, j, n

        data w /0.066671344308688d0, 0.149451349150581d0, 0.219086362515982d0, &
                0.269266719309996d0, 0.295524224714753d0, 0.295524224714753d0, &
                0.269266719309996d0, 0.219086362515982d0, 0.149451349150581d0, &
                0.066671344308688d0/
        data a /0.97390652851717174d0, 0.86506336668898454d0, 0.67940956829902444d0, &
                0.43339539412924721d0, 0.14887433898163122d0, -0.14887433898163119d0, &
                -0.43339539412924721d0, -0.67940956829902444d0, -0.86506336668898454d0, &
                -0.97390652851717174d0/

        n = (max-min)/step
        half = step/2
        QG10=0
        !$omp parallel default(none) private(x, mid, next, er, tmp) shared(half, n, min, step, a, w, QG10, fptr) num_threads(thr)
            er = 0
            !$omp do reduction(+: QG10)
            do i = 1, n
                x = min + (i-1)*step
                mid = x + half
                do j = 1, 10
                    next = half*w(j)*fptr(a(j)*half+mid) - er
                    tmp = QG10 + next
                    er = (tmp - QG10) - next
                    QG10 = tmp
                end do
            end do
            !$omp end do
        !$omp end parallel

    end function QG10

    real(8) function innerFunc(y)
        use com
        implicit none
        real(8), intent(in) :: y
        innerFunc = xx**2*y**2
    end function innerFunc

    pure real(8) function arc(x)
        implicit none
        real(8), intent(in) :: x
        arc = sqrt(1-x**2)
    end function arc

    real(8) function intgrCircle(x)
        use com
        implicit none
        procedure(func), pointer :: fptr
        real(8), intent(in) :: x

        fptr => innerFunc
        !$omp critical
        xx = x
        intgrCircle = QG10(fptr, 0.d0, arc(x), step)
        !$omp end critical
        intgrCircle = 2*intgrCircle
    end function intgrCircle

    pure real(8) function vertArc(x)
        implicit none
        real(8), intent(in) :: x
        vertArc = sqrt(1-4*x**2)
    end function vertArc

    real(8) function intgrVert(x)
        use com
        implicit none
        procedure(func), pointer :: fptr
        real(8), intent(in) :: x

        fptr => innerFunc
        !$omp critical
        xx = x
        if (x - 0.5d0 .lt. 1d-8) then
            intgrVert = QG10(fptr, vertArc(x), arc(x), step)
        else
            intgrVert = QG10(fptr, 0.d0, arc(x), step)
        end if
        !$omp end critical
        intgrVert = 2*intgrVert
    end function intgrVert

    real(8) function intgrHor(x)
        use com
        implicit none
        procedure(func), pointer :: fptr
        real(8), intent(in) :: x

        fptr => innerFunc
        !$omp critical
        xx = x
        intgrHor = QG10(fptr, arc(x)/2, arc(x), step)
        !$omp end critical
        intgrHor = 2*intgrHor
    end function intgrHor

    real(8) function polFunc(rho)
        use com
        implicit none
        real(8), intent(in) :: rho
        polFunc = rho**5*(sin(2*xx)**2)
    end function polFunc

    pure real(8) function polArc(x)
        implicit none
        real(8), intent(in) :: x
        polArc = sqrt(2.d0/(5-3*sin(2*x)))
    end function polArc

    real(8) function polIntgrTil(phi)
        use com
        implicit none
        procedure(func), pointer :: fptr
        real(8), intent(in) :: phi

        fptr => polFunc
        !$omp critical
        xx = phi
        polIntgrTil = QG10(fptr, polArc(phi), 1.d0, step)
        !$omp end critical
    end function polIntgrTil

end program output
