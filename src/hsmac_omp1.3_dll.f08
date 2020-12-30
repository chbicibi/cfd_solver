include "hsmac_omp1.2.f08"


module io
    use iso_c_binding
    use hsmac2d
    implicit none

    contains

    subroutine io_read_inputfile() bind(c, name="f_read_inputfile")
        call read_inputfile("in2d.txt")
    end subroutine io_read_inputfile


    subroutine io_initialize() bind(c, name="f_initialize")
        call initialize("grid.csv")
    end subroutine io_initialize


    subroutine io_advance(uin, vin, tin, uout, vout, tout, nx, ny) bind(c, name="f_advance")
        integer(c_int), intent(in) :: nx, ny
        real(c_double), intent(in) :: uin(0:nx, 0:ny+1), vin(0:nx+1, 0:ny), tin(0:nx+1, 0:ny+1)
        real(c_double), intent(out) :: uout(0:nx, 0:ny+1), vout(0:nx+1, 0:ny), tout(0:nx+1, 0:ny+1)

        call advance(uin, vin, tin, uout, vout, tout)
    end subroutine io_advance


    subroutine io_calc_velociry(uin, vin, pin, tin, uout, vout, nx, ny) bind(c, name="f_calc_velociry")
        integer(c_int), intent(in) :: nx, ny
        real(c_double), intent(in) :: uin(0:nx, 0:ny+1), vin(0:nx+1, 0:ny), pin(0:nx+1, 0:ny+1), tin(0:nx+1, 0:ny+1)
        real(c_double), intent(out) :: uout(0:nx, 0:ny+1), vout(0:nx+1, 0:ny)

        call calc_velociry(uin, vin, pin, tin, uout, vout)
    end subroutine io_calc_velociry


    subroutine io_calc_pressure(u, v, p, itr, flag, nx, ny, msg) bind(c, name="f_calc_pressure")
        integer(c_int), intent(in) :: nx, ny
        integer(c_int), intent(in) :: itr
        real(c_double), intent(inout) :: u(0:nx, 0:ny+1), v(0:nx+1, 0:ny), p(0:nx+1, 0:ny+1)
        integer(c_int), intent(out) :: flag
        character(c_char), intent(out) :: msg(120)

        if (msg(1) == "x") then
            call calc_pressure(itr, u, v, p, flag)
        else
            call calc_pressure(itr, u, v, p, flag, msg)
        end if
    end subroutine io_calc_pressure


    subroutine io_calc_temperature(uin, vin, tin, tout, nx, ny) bind(c, name="f_calc_temperature")
        integer(c_int), intent(in) :: nx, ny
        real(c_double), intent(in) :: uin(0:nx, 0:ny+1), vin(0:nx+1, 0:ny), tin(0:nx+1, 0:ny+1)
        real(c_double), intent(out) :: tout(0:nx+1, 0:ny+1)

        call calc_temperature(uin, vin, tin, tout)
    end subroutine io_calc_temperature


    subroutine io_bind_velocity(u, v, field, nx, ny) bind(c, name="f_bind_velocity")
        integer(c_int), intent(in) :: nx, ny
        real(c_double), intent(inout) :: u(0:nx, 0:ny+1), v(0:nx+1, 0:ny), field(0:nx+1, 0:ny+1)

        call bind_velocity(u, v, field)
    end subroutine io_bind_velocity


    ! subroutine io_output_raw(filename) bind(c, name="output_raw")
    !     character(kind=c_char), intent(in) :: filename(*)

    !     call output_raw(filename)
    ! end subroutine io_output_raw


    subroutine testsub(a, n) bind(c, name="testsub")
        integer(c_int), intent(in) :: n
        character(c_char), intent(in) :: a(n)
        ! integer(c_int), intent(inout) :: a(3, 3)
        ! integer(c_int), intent(in), value :: n
        ! real(c_double), intent(in), value :: d
        ! real(c_double), intent(in) :: elevator_in(n)
        print *, a
        print *, nx, ny
        ! a = 2
        ! print *, a

        ! print *, filename(1:3)
    end subroutine testsub
end module io
