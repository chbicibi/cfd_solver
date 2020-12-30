! 2020.12.09 version 1.2.0
! 2020.12.10 version 1.2.1
! 2020.12.12 version 1.2.2


module string
    implicit none

    contains

    function to_str(n)
        character(:), allocatable :: to_str
        integer,      intent(in)  :: n
        integer                   :: len

        if (n == 0) then
            len = 1
        else if (n > 0) then
            len = int(log10(dble(n))) + 1
        else
            len = int(log10((dble(-n)))) + 2
        end if

        allocate(character(len) :: to_str)
        write(to_str, "(i0)") n
    end function to_str


    integer function to_i(str)
        character(*), intent(in) :: str

        read(str, *) to_i
    end function to_i
end module string


module hsmac2d
    use string

    implicit none

    ! ここで定義した変数は program~ から end program~ までの範囲で共有される
    ! integer, parameter   :: nx0 = 120, ny0 = 80 ! parameter属性がついていると変更不可になる
    integer              :: nx, ny
    integer              :: itype, icycle, nitr, ncycle, itr, iflg, irelp, method, iter
    integer              :: iout
    integer              :: stridex, stridey ! 2020.12.9 added
    real(8)              :: dx, dy, dt
    real(8)              :: vis, alp, buo
    real(8)              :: re, pr, gr, time, omg, epsp
    real(8)              :: dmax
    real(8), allocatable :: uo(:, :), un(:, :), vo(:, :), vn(:, :)
    real(8), allocatable :: to(:, :), tn(:, :), po(:, :), fi(:, :), div(:, :)
    real(8), allocatable :: f2x(:, :), f2y(:, :)
    real(8), allocatable :: delpa(:, :)
    integer, allocatable :: n_pitr(:)
    character(100) :: dirname


    contains

    subroutine main
        character(:), allocatable :: dir_out, file_out, file_grid, index
        allocate(character(10) :: index)

        ! パラメータファイルの読み込み
        call read_inputfile("in2d.txt")
        ! 出力ディレクトリ設定
        dir_out  = trim(dirname)
        ! 出力ファイル名設定
        file_out = "out_"
        ! グリッドファイル名設定
        file_grid = "grid.csv"
        ! 出力ディレクトリ作成
        ! call system("mkdir "// dir_out // " > nul 2>&1")
        if (access(dir_out, " ") > 0) call system("mkdir "// dir_out)
        ! 初期値の設定
        call initialize(file_grid)

        if (iter > 0) then
            write(index, "(i0.4)") iter - 1
            call load_value(dir_out // "/" // file_out // trim(index) // ".dat")
        end if

        ! 時間進行の反復
        ! 時間進行カウンタ(icycle)がncycleより小さい時
        do while (icycle < ncycle)
            ! 時間進行
            call advance(un, vn, tn, uo, vo, to) ! xn -> xo

            ! 速度場の計算
            call calc_velociry(uo, vo, po, un, vn)

            ! 圧力反復
            do itr = 1, nitr
                ! 圧力場の計算
                call calc_pressure(itr, un, vn, po)

                ! iflg -> 0:収束 1:収束せず
                ! newton法による圧力場の計算が収束したとき圧力反復終了
                if (iflg == 0) exit

                ! 圧力場の計算が収束していないとき
                ! 圧力計算の反復回数がnitrとなったら発散とみなして計算終了(→メッセージを表示して続行)
                if (itr == nitr) then
                    write (6, *) 'not converge !'
                    ! call exit(1)
                end if
            end do
            n_pitr(icycle) = itr ! 2020.12.9 test

            ! 非等温場計算の場合は温度場を計算
            if (itype == 2) call calc_temperature(uo, vo, to, tn)

            if (mod(icycle, iout) == 0) then
                write(index, "(i0.4)") iter
                call output_raw(dir_out // "/" // file_out // trim(index) // ".dat")
                iter = iter + 1

                block ! 2020.12.09 test
                    integer :: i, unit
                    open(newunit=unit, file=dir_out // "/nitr.csv")
                    write(unit, "(a)") "icycle,nitr"
                    do i = 1, ncycle
                        if (n_pitr(i) > 0) then
                            write(unit, "(i0,a,i0)") i, ",", n_pitr(i)
                        end if
                    end do
                    close(unit)
                end block

            end if
        end do
    end subroutine main


    ! *********************************************************************
    ! *                        パラメータ読み込み
    ! *********************************************************************
    subroutine read_inputfile(filename)
        character(*), intent(in) :: filename
        real(8)                  :: dlx, dly
        integer                  :: i

        ! パラメータファイルのオープン
        open (10, file = filename)
            ! in2d.txt中のコメント行(1-15行目)のスキップ
            do i = 1, 15
                read (10, *)
            end do
            ! itype  -> 1:等温場 2:非等温場
            ! icycle -> 時間進行カウンタの初期値(普通は0)
            ! nitr   -> 圧力反復回数の上限
            ! ncycle -> 時間進行回数の上限
            read (10, *) itype, icycle, nitr, ncycle
            ! in2d.txt中のコメント行(17-18行目)のスキップ
            do i = 1, 2
                read (10, *)
            end do
            ! epsp -> 許容誤差
            ! omg  -> 重み付け
            read (10, *) epsp, omg
            ! in2d.txt中のコメント行(20-21行目)のスキップ
            do i = 1, 2
                read (10, *)
            end do
            ! dt -> 時間刻み幅
            ! re -> レイノルズ数
            ! pr -> プラントル数
            ! gr -> グラスホフ数
            read (10, *) dt, re, pr, gr
            ! in2d.txt中のコメント行(23-24行目)のスキップ
            do i = 1, 2
                read (10, *)
            end do
            read (10, *) nx, ny, stridex, stridey
            ! in2d.txt中のコメント行(26-27行目)のスキップ
            do i = 1, 2
                read (10, *)
            end do
            ! dlx    -> 解析領域の横幅
            ! dly    -> 解析領域の高さ
            ! irelp  -> 1:圧力を正規化する
            ! method -> 使用しない
            read (10, *) dlx, dly, irelp, method, iter
            ! in2d.txt中のコメント行(29-30行目)のスキップ
            do i = 1, 2
                read (10, *)
            end do
            read (10, *) iout, dirname
        close (10)

        if (mod(nx, stridex) > 0 .or. mod(ny, stridey) > 0) then
            print *, "Error: Both nx and ny must be a number that is divisible by stride."
            call exit(1)
        end if

        ! x方向の格子幅
        dx  = dlx / dble(nx)
        ! y方向の格子幅
        dy  = dly / dble(ny)
        ! 運動方程式中の拡散項の係数(ここではre)
        vis = 1.0d0 / re
        ! エネルギー方程式中の拡散項の係数(ここでは1)
        alp = 1.0d0
        ! 等温場なら浮力項の係数はゼロに設定
        ! 非等温場なら浮力項の係数を計算(ここでは gr * pr ** 2)
        if (itype == 1) then
            buo = 0.0d0
        else
            buo = gr * pr ** 2
        end if
    end subroutine read_inputfile


    ! *********************************************************************
    ! *                        初期設定
    ! *********************************************************************
    subroutine initialize(filename)
        character(*), intent(in) :: filename
        integer :: ix, iy, unit

        ! 配列割り当て
        allocate(uo(0:nx  , 0:ny+1), un(0:nx ,  0:ny+1  ), &
                 vo(0:nx+1, 0:ny  ), vn(0:nx+1, 0:ny    ), &
                 to(0:nx+1, 0:ny+1), tn(0:nx+1, 0:ny+1), &
                 po(0:nx+1, 0:ny+1), fi(0:nx+1, 0:ny+1), div(1:nx, 1:ny))
        allocate(f2x(0:nx, 0:ny), source=0d0)
        allocate(f2y(0:nx, 0:ny), source=0d0)
        allocate(delpa(1:nx+1, 1:ny+1), source=0d0)
        allocate(n_pitr(ncycle), source=-1)

        ! 計算条件の設定
        ! uの初期値設定
        un(0:nx, 0:ny+1) = 1.0d0
        ! vの初期値設定
        vn(0:nx+1, 0:ny) = 0.0d0
        ! pの初期値設定
        po(0:nx+1, 0:ny+1) = 0.0d0

        ! 格子情報の初期値設定
        fi(0:nx+1, 0:ny+1) = 1.0d0

        open(newunit=unit, file=filename, status="old")
            fi(0:nx+1, 0:ny+1) = 1.0d0
            do iy = 1, ny
                read (unit, *) fi(1:nx, iy)
            end do
        close(unit)

        ! ----------------------------------------------------------------------
        ! （注意）浮力項の計算で温度の配列を使用しているので等温場でもt=0として
        ! 初期条件だけは設定する必要がある．ゼロ以外の値を入れると浮力項が計算
        ! される可能性があるので注意．
        ! ----------------------------------------------------------------------
        ! tの初期値設定(領域内は高温(+0.5)と低温(-0.5)の中間温度)
        tn(0:nx+1, 0:ny+1) = 0.0d0
        ! tの境界：右面（冷却）t=-0.5
        tn(nx+1  , 0:ny+1) = -1.0d0 - tn(nx, 0:ny+1) ! 2.0d0 * -0.5d0
        ! tの境界：左面（加熱）t=+0.5
        tn(0     , 0:ny+1) = 1.0d0 - tn(1, 0:ny+1) ! 2.0d0 * +0.5d0
        ! tの境界：上面（断熱）
        tn(1:nx  , ny+1) = tn(1:nx, ny)
        ! tの境界：下面（断熱）
        tn(1:nx  , 0) = tn(1:nx, 1)
    end subroutine initialize


    subroutine load_value(filename)
        character(*), intent(in) :: filename
        integer :: ix, iy, unit, info(3)
        real(8) :: dummy
        real :: temp(0:nx+1, 0:ny+1)

        if (access(filename, " ") > 0) then
            print *, "Error: no file"
            print *, filename
            call exit(1)
        end if

        ! open(newunit=unit, file=filename, status="old")
        !     read(unit, *)
        !     read(unit, *)
        !     do iy = 1, ny
        !         do ix = 1, nx
        !             read(unit, *) dummy, dummy, un(ix, iy), vn(ix, iy), po(ix, iy), fi(ix, iy)
        !         end do
        !     end do
        ! close(unit)

        ! 2020.12.10 changed
        open(newunit=unit, file=filename, access="stream", status="old")
            read(unit) info
            call assert("un", info, [nx+1, ny+2, 4])
            read(unit) temp(0:nx, 0:ny+1)
            un = dble(temp(0:nx, 0:ny+1))

            read(unit) info
            call assert("vn", info, [nx+2, ny+1, 4])
            read(unit) temp(0:nx+1, 0:ny)
            vn = dble(temp(0:nx+1, 0:ny))

            read(unit) info
            call assert("po", info, [nx+2, ny+2, 4])
            read(unit) temp(0:nx+1, 0:ny+1)
            po = dble(temp(0:nx+1, 0:ny+1))
        close(unit)

        ! un(1:nx - 1, :) = (un(1:nx - 1, :) + un(2:nx, :)) * 0.5d0
        ! vn(:, 1:ny - 1) = (vn(:, 1:ny - 1) + vn(:, 2:ny)) * 0.5d0
        ! un(0   , :   ) = un(1 , : )
        ! un(:   , 0   ) = un(: , 0 )
        ! un(:   , ny+1) = un(: , ny)
        ! vn(0   , :   ) = vn(1 , : )
        ! vn(:   , 0   ) = vn(: , 0 )
        ! vn(nx+1, :   ) = vn(nx, : )

        contains

        subroutine assert(name, info, collect)
            character(2), intent(in) :: name
            integer, intent(in) :: info(3), collect(3)
            if (.not. all(info == collect)) then
                print *,  "Invalid data shape of '" // name // "'"
                call exit(1)
            end if
        end subroutine assert
    end subroutine load_value


    ! *********************************************************************
    ! *                         時間進行
    ! *********************************************************************
    subroutine advance(uin, vin, tin, uout, vout, tout)
        real(8), intent(in) :: uin(0:, 0:), vin(0:, 0:), tin(0:, 0:)
        real(8), intent(out) :: uout(0:, 0:), vout(0:, 0:), tout(0:, 0:)
        integer :: ix, iy

        time   = dt * dble(icycle)
        icycle = icycle + 1
        ! 時間進行カウンタ(icycle)を100回毎に表示
        if (mod(icycle, 100) == 0) write (6, "('A:cycle=', i8, ' itr=', i8)") icycle, itr
        ! --------------------------------------------------------------------
        ! un -> uo:必要なら入れ替える前にunとuoから変動量を求める
        ! un:前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される
        ! uo:新しい時間ステップでの初期値．unを保存．
        ! --------------------------------------------------------------------
        uout(0:nx, 0:ny+1) = uin(0:nx, 0:ny+1)
        ! --------------------------------------------------------------------
        ! vn -> vo:必要なら入れ替える前にvnとvoから変動量を求める
        ! vn:前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される
        ! vo:新しい時間ステップでの初期値．vnを保存．
        ! --------------------------------------------------------------------
        vout(0:nx+1, 0:ny) = vin(0:nx+1, 0:ny)
        ! --------------------------------------------------------------------
        ! tn -> to:必要なら入れ替える前にtnとtoから変動量を求める
        ! tn:前の時間ステップでの値
        ! to:新しい時間ステップでの初期値．tnを保存．
        ! --------------------------------------------------------------------
        tout(0:nx+1, 0:ny+1) = tin(0:nx+1, 0:ny+1)
    end subroutine advance


    ! *********************************************************************
    ! *                    速度場の計算
    ! *********************************************************************
    subroutine calc_velociry(uin, vin, pin, uout, vout)
        real(8), intent(in) :: uin(0:, 0:), vin(0:, 0:), pin(0:, 0:)
        real(8), intent(out) :: uout(0:, 0:), vout(0:, 0:)
        integer :: ix, iy

        ! --------------------------------------------------------------------
        !                      u(ix, iy)の計算
        ! --------------------------------------------------------------------
        !$omp parallel do private(iy)
        do iy = 1, ny
            do ix = 1, nx - 1
                call calcu(ix, iy)
            end do
            ! call calcu_vec(iy)
        end do
        !$omp end parallel do

        ! --------------------------------------------------------------------
        ! v(ix, iy)の計算
        ! --------------------------------------------------------------------
        !$omp parallel do private(iy)
        do iy = 1, ny - 1
            do ix = 1, nx
                call calcv(ix, iy)
            end do
            ! call calcv_vec(iy)
        end do
        !$omp end parallel do

        ! 速度の境界条件の処理
        call bind_velocity(uout, vout)

        contains
        subroutine calcu(ix, iy)
            integer, intent(in) :: ix, iy
            real(8) :: vv, cnvux, cnvuy, tu, buou, difu

            ! vvはu(ix, iy)におけるvの補間値
            vv = (vin(ix, iy) + vin(ix+1 , iy) + vin(ix, iy-1) + vin(ix+1 , iy-1)) * 0.25d0
            ! 対流項(cnvux, cnvuy)を１次精度風上差分にて計算
            if (uin(ix, iy) >= 0.0d0) then
                cnvux = uin(ix, iy) * (uin(ix, iy) - uin(ix-1, iy)) / dx
            else
                cnvux = uin(ix, iy) * (uin(ix+1, iy) - uin(ix, iy)) / dx
            end if

            if (vv >= 0.0d0) then
                cnvuy = vv * (uin(ix, iy) - uin(ix, iy-1)) / dy
            else
                cnvuy = vv * (uin(ix, iy+1) - uin(ix, iy)) / dy
            end if

            ! x方向の浮力項(buou)はゼロ
            tu   = 0.0d0
            buou = buo * tu
            ! 拡散項(difu)の計算
            difu = vis * ((uin(ix-1, iy) - 2.0d0 * uin(ix, iy) + uin(ix+1, iy)) / dx ** 2 &
                        + (uin(ix, iy-1) - 2.0d0 * uin(ix, iy) + uin(ix, iy+1)) / dy ** 2)
            ! 仮の速度(u)の計算
            uout(ix, iy) = uin(ix, iy) + dt * (-cnvux - cnvuy + difu + buou + (pin(ix, iy) - pin(ix+1, iy)) / dx)
        end subroutine calcu

        subroutine calcv(ix, iy)
            integer, intent(in) :: ix, iy
            real(8) :: uu, cnvvx, cnvvy, tv, buov, difv

            ! uuはv(ix, iy)におけるuの補間値
            uu = (uin(ix-1, iy) + uin(ix, iy) + uin(ix-1, iy+1) + uin(ix, iy+1)) * 0.25d0
            ! 対流項(cnvvx, cnvvy)を１次精度風上差分にて計算
            if (uu >= 0.0d0) then
                cnvvx = uu * (vin(ix, iy) - vin(ix-1, iy)) / dx
            else
                cnvvx = uu * (vin(ix+1, iy) - vin(ix, iy)) / dx
            end if

            if (vin(ix, iy) >= 0.0d0) then
                cnvvy = vin(ix, iy) * (vin(ix, iy) - vin(ix, iy-1)) / dy
            else
                cnvvy = vin(ix, iy) * (vin(ix, iy+1) - vin(ix, iy)) / dy
            end if
            ! 浮力項(buov)の計算
            tv   = (to(ix, iy) + to(ix, iy+1)) * 0.5d0
            buov = buo * tv
            ! 拡散項(difv)の計算
            difv = vis * ((vin(ix-1, iy) - 2.0d0 * vin(ix, iy) + vin(ix+1, iy)) / dx ** 2 &
                        + (vin(ix, iy-1) - 2.0d0 * vin(ix, iy) + vin(ix, iy+1)) / dy ** 2)
            ! 仮の速度(v)の計算
            vout(ix, iy) = vin(ix, iy) + dt * (-cnvvx - cnvvy + difv + buov + (pin(ix, iy) - pin(ix, iy+1)) / dy)
        end subroutine calcv

        subroutine calcu_vec(iy)
            integer, intent(in) :: iy
            real(8) :: vv(nx-1), cnvux(nx-1), cnvuy(nx-1), tu(nx-1), buou(nx-1), difu(nx-1)

            ! vvはu(ix, iy)におけるvの補間値
            vv = (vin(1:nx-1, iy) + vin(2:nx , iy) + vin(1:nx-1, iy-1) + vin(2:nx , iy-1)) * 0.25d0
            ! 対流項(cnvux, cnvuy)を１次精度風上差分にて計算
            where (uin(1:nx-1, iy) >= 0.0d0)
                cnvux = uin(1:nx-1, iy) * (uin(1:nx-1, iy) - uin(0:nx-2, iy)) / dx
            else where
                cnvux = uin(1:nx-1, iy) * (uin(2:nx , iy) - uin(1:nx-1, iy)) / dx
            end where

            where (vv >= 0.0d0)
                cnvuy = vv * (uin(1:nx-1, iy) - uin(1:nx-1, iy-1)) / dy
            else where
                cnvuy = vv * (uin(1:nx-1, iy+1) - uin(1:nx-1, iy)) / dy
            end where

            ! x方向の浮力項(buou)はゼロ
            tu   = 0.0d0
            buou = buo * tu
            ! 拡散項(difu)の計算
            difu = vis * ((uin(0:nx-2, iy) - 2.0d0 * uin(1:nx-1, iy) + uin(2:nx , iy)) / dx ** 2 &
                        + (uin(1:nx-1, iy-1) - 2.0d0 * uin(1:nx-1, iy) + uin(1:nx-1, iy+1)) / dy ** 2)
            ! 仮の速度(u)の計算
            uout(1:nx-1, iy) = uin(1:nx-1, iy) + dt * (-cnvux - cnvuy + difu + buou &
                                                           + (pin(1:nx-1, iy) - pin(2:nx , iy)) / dx)
        end subroutine calcu_vec

        subroutine calcv_vec(iy)
            integer, intent(in) :: iy
            real(8) :: uu(nx), cnvvx(nx), cnvvy(nx), tv(nx), buov(nx), difv(nx)

            ! uuはv(ix, iy)におけるuの補間値
            uu = (uin(0:nx-1, iy) + uin(1:nx, iy) + uin(0:nx-1, iy+1) + uin(1:nx, iy+1)) * 0.25d0
            ! 対流項(cnvvx, cnvvy)を１次精度風上差分にて計算
            where (uu >= 0.0d0)
                cnvvx = uu * (vin(1:nx, iy) - vin(0:nx-1, iy)) / dx
            else where
                cnvvx = uu * (vin(2:nx+1, iy) - vin(1:nx, iy)) / dx
            end where

            where (vin(1:nx, iy) >= 0.0d0)
                cnvvy = vin(1:nx, iy) * (vin(1:nx, iy) - vin(1:nx, iy-1)) / dy
            else where
                cnvvy = vin(1:nx, iy) * (vin(1:nx, iy+1) - vin(1:nx, iy)) / dy
            end where
            ! 浮力項(buov)の計算
            tv   = (to(1:nx, iy) + to(1:nx, iy+1)) * 0.5d0
            buov = buo * tv
            ! 拡散項(difv)の計算
            difv = vis * ((vin(0:nx-1, iy) - 2.0d0 * vin(1:nx, iy) + vin(2:nx+1, iy)) / dx ** 2 &
                        + (vin(1:nx, iy-1) - 2.0d0 * vin(1:nx, iy) + vin(1:nx, iy+1)) / dy ** 2)
            ! 仮の速度(v)の計算
            vout(1:nx, iy) = vin(1:nx, iy) + dt * (-cnvvx - cnvvy + difv + buov &
                                                           + (pin(1:nx, iy) - pin(1:nx, iy+1)) / dy)
        end subroutine calcv_vec

        subroutine calcu_mat
            real(8) :: vv(nx-1, ny), cnvux(nx-1, ny), cnvuy(nx-1, ny)
            real(8) :: tu(nx-1, ny), buou(nx-1, ny), difu(nx-1, ny)

            ! vvはu(ix, iy)におけるvの補間値
            vv = (vin(1:nx-1, 1:ny) + vin(2:nx , 1:ny) + vin(1:nx-1, 0:ny-1) + vin(2:nx , 0:ny-1)) * 0.25d0
            ! 対流項(cnvux, cnvuy)を１次精度風上差分にて計算
            where (uin(1:nx-1, 1:ny) >= 0.0d0)
                cnvux = uin(1:nx-1, 1:ny) * (uin(1:nx-1, 1:ny) - uin(0:nx-2, 1:ny)) / dx
            else where
                cnvux = uin(1:nx-1, 1:ny) * (uin(2:nx , 1:ny) - uin(1:nx-1, 1:ny)) / dx
            end where

            where (vv >= 0.0d0)
                cnvuy = vv * (uin(1:nx-1, 1:ny) - uin(1:nx-1, 0:ny-1)) / dy
            else where
                cnvuy = vv * (uin(1:nx-1, 2:ny+1) - uin(1:nx-1, 1:ny)) / dy
            end where

            ! x方向の浮力項(buou)はゼロ
            tu   = 0.0d0
            buou = buo * tu
            ! 拡散項(difu)の計算
            difu = vis * ((uin(0:nx-2, 1:ny) - 2.0d0 * uin(1:nx-1, 1:ny) + uin(2:nx , 1:ny)) / dx ** 2 &
                        + (uin(1:nx-1, 0:ny-1) - 2.0d0 * uin(1:nx-1, 1:ny) + uin(1:nx-1, 2:ny+1)) / dy ** 2)
            ! 仮の速度(u)の計算
            uout(1:nx-1, 1:ny) = uin(1:nx-1, 1:ny) + dt * (-cnvux - cnvuy + difu + buou &
                                                           + (pin(1:nx-1, 1:ny) - pin(2:nx , 1:ny)) / dx)
        end subroutine calcu_mat

        subroutine calcv_mat
            real(8) :: uu(nx, ny-1), cnvvx(nx, ny-1), cnvvy(nx, ny-1)
            real(8) :: tv(nx, ny-1), buov(nx, ny-1), difv(nx, ny-1)

            ! uuはv(ix, iy)におけるuの補間値
            uu = (uin(0:nx-1, 1:ny-1) + uin(1:nx, 1:ny-1) + uin(0:nx-1, 2:ny) + uin(1:nx, 2:ny)) * 0.25d0
            ! 対流項(cnvvx, cnvvy)を１次精度風上差分にて計算
            where (uu >= 0.0d0)
                cnvvx = uu * (vin(1:nx, 1:ny-1) - vin(0:nx-1, 1:ny-1)) / dx
            else where
                cnvvx = uu * (vin(2:nx+1, 1:ny-1) - vin(1:nx, 1:ny-1)) / dx
            end where

            where (vin(1:nx, 1:ny-1) >= 0.0d0)
                cnvvy = vin(1:nx, 1:ny-1) * (vin(1:nx, 1:ny-1) - vin(1:nx, 0:ny-2)) / dy
            else where
                cnvvy = vin(1:nx, 1:ny-1) * (vin(1:nx, 2:ny) - vin(1:nx, 1:ny-1)) / dy
            end where
            ! 浮力項(buov)の計算
            tv   = (to(1:nx, 1:ny-1) + to(1:nx, 2:ny)) * 0.5d0
            buov = buo * tv
            ! 拡散項(difv)の計算
            difv = vis * ((vin(0:nx-1, 1:ny-1) - 2.0d0 * vin(1:nx, 1:ny-1) + vin(2:nx+1, 1:ny-1)) / dx ** 2 &
                        + (vin(1:nx, 0:ny-2) - 2.0d0 * vin(1:nx, 1:ny-1) + vin(1:nx, 2:ny)) / dy ** 2)
            ! 仮の速度(v)の計算
            vout(1:nx, 1:ny-1) = vin(1:nx, 1:ny-1) + dt * (-cnvvx - cnvvy + difv + buov &
                                                           + (pin(1:nx, 1:ny-1) - pin(1:nx, 2:ny)) / dy)
        end subroutine calcv_mat
    end subroutine calc_velociry


    ! *********************************************************************
    ! *                        圧力場の計算
    ! *********************************************************************
    subroutine calc_pressure(itr, u, v, p)
        integer, intent(in) :: itr
        real(8), intent(inout) :: u(0:, 0:), v(0:, 0:), p(0:, 0:)
        real(8) :: del, div_max!, delp
        real(8), save :: div_max_last, div_max_diff, omg1 = 1.0d0, epsp1
        integer :: pos_max(2)
        integer :: ix, iy
        ! real(8) :: uadd(nx, ny), usub(nx, ny), vadd(nx, ny), vsub(nx, ny)

        ! p(ix, iy)の計算
        del = dt * (2.0d0 / dx ** 2 + 2.0d0 / dy ** 2)
        ! omg1 = max(min(1.0d0 + (itr - 1) * 0.1d0, omg) - itr * 1d-3, 1.0d0)
        ! omg1 = 1.0d0 + (omg - 1) * 3.591d-2 * (itr - 1) / (1 + exp(0.01d0 * (itr - 1)))

        ! if (itr <= 10) then
            ! omg1 = min(1 + (omg - 1) * (itr - 1) / 200d0, omg)
            ! omg1 = 1.0d0
        ! else if (itr <= 100) then
            ! omg1 = min(1 + (omg - 1) * (itr - 1) / 200d0, omg)
            ! omg1 = min(max(omg * (div_max_last / epsp - 1), 0.5d0), 1.0d0)
            ! omg1 = min(max(0.99d0 * omg1 + &
            !                0.01d0 * omg * 1.4d0 ** log10(div_max_last / epsp - 1), 0.5d0), 1.0d0)
        if (itr == 1) then
            omg1 = max(0.9999d0 * omg1, 1.0d0)
            epsp1 = epsp
        ! else if (div_max_diff / epsp < -0.1d0) then
        !     omg1 = min(omg1 * 1.01d0, omg)
        ! else if (div_max_diff / epsp > 0.1d0) then
        !     omg1 = max(omg1 / 1.01d0, 1.0d0)
        else if (itr == 100) then
            omg1 = min(1.1d0 * omg1, omg)
        else if (itr >= 1000) then
            epsp1 = epsp * 10d0 * (itr - 1000) / (nitr - 1000)
        end if

        !$omp parallel do private(ix, iy)
        do iy = 1, ny, stridey
            do ix = 1, nx, stridex
                call subloop(ix, iy)
            end do
        end do
        !$omp end parallel do

        ! div(1:nx, 1:ny)   = (u(1:nx, 1:ny) - u(0:nx-1, 1:ny)) / dx &
        !                   + (v(1:nx, 1:ny) - v(1:nx, 0:ny-1)) / dy
        ! delpa(1:nx, 1:ny) = -omg * div(1:nx, 1:ny) / del / 2
        ! p(1:nx, 1:ny)     = p(1:nx, 1:ny) + delpa(1:nx, 1:ny)
        ! u(1:nx, 1:ny)     = u(1:nx, 1:ny) + dt / dx * (delpa(1:nx, 1:ny) - delpa(2:nx+1, 1:ny))
        ! v(1:nx, 1:ny)     = v(1:nx, 1:ny) + dt / dy * (delpa(1:nx, 1:ny) - delpa(1:nx, 2:ny+1))

        ! 圧力の相対性に関する処理(irelp=1なら以下の処理を行う)
        if (irelp == 1) p(1:nx, 1:ny) = p(1:nx, 1:ny) - p(1, 1)
        ! 発散の最大値を調べる
        ! div(1, :) = 0 ! 2020.12.07 temp
        pos_max = maxloc(abs(div))
        div_max = abs(div(pos_max(1), pos_max(2)))
        div_max_diff = div_max - div_max_last
        div_max_last = div_max

        ! iflg=1なら，連続の式を満たしていないと判定し再び圧力計算を行う．
        if (div_max >= epsp1) then
        ! if (count(abs(div) >= epsp) >= 10) then ! 2020.12.07 temp
            iflg = 1
        else
            iflg = 0
        end if

        ! 圧力計算の回数を100回ごとに表示
        if (mod(itr, 100) == 0) then
            write (6, "('P:cycle=', i8, ' itr=', i8, '   div(max)(', 2i6, ')=', 1pe13.5, '   omg=', g11.4)") &
                icycle, itr, pos_max, div_max, omg1
        end if

        ! 新たに得られた速度を用いて境界条件を処理する
        call bind_velocity(u, v)

        contains

        subroutine calcp(ix, iy)
            integer, intent(in) :: ix, iy
            real(8) :: delp

            div(ix, iy) = (u(ix, iy) - u(ix-1, iy)) / dx &
                        + (v(ix, iy) - v(ix, iy-1)) / dy
            delp          = -omg1 * div(ix, iy) / del
            p(ix,   iy  ) = p(ix,   iy  ) + delp
            u(ix,   iy  ) = u(ix,   iy  ) + dt / dx * delp
            u(ix-1, iy  ) = u(ix-1, iy  ) - dt / dx * delp
            v(ix,   iy  ) = v(ix,   iy  ) + dt / dy * delp
            v(ix,   iy-1) = v(ix,   iy-1) - dt / dy * delp
        end subroutine calcp

        subroutine subloop(ix, iy)
            integer, intent(in) :: ix, iy
            integer :: i, j

            do j = 0, stridey - 1
                do i = 0, stridex - 1
                    call calcp(ix+i, iy+j)
                end do
            end do
        end subroutine subloop
    end subroutine calc_pressure

    ! *********************************************************************
    ! *                      温度場の計算
    ! *********************************************************************
    subroutine calc_temperature(uin, vin, tin, tout)
        real(8), intent(in) :: uin(0:, 0:), vin(0:, 0:), tin(0:, 0:)
        real(8), intent(out) :: tout(0:, 0:)
        integer :: ix, iy

        ! t(ix, iy)の計算
        !$omp parallel do private(ix, iy)
        do iy = 1, ny
            do ix = 1, nx
                call calct(ix, iy)
            end do
        end do
        !$omp end parallel do

        ! 境界条件の処理
        call bind_temperature(tout)

        contains

        subroutine calct(ix, iy)
            integer, intent(in) :: ix, iy
            real(8) :: uut, vvt, cnvtx, cnvty, dift

            ! uut, vvtはそれぞれt(ix, iy)におけるu, vの補間値
            uut = (uin(ix, iy) + uin(ix-1, iy    )) * 0.5d0
            vvt = (vin(ix, iy) + vin(ix,     iy-1)) * 0.5d0
            ! 対流項(cnvtx, cnvty)を１次精度風上差分にて計算
            if (uut >= 0.0d0) then
                cnvtx = uut * (tin(ix, iy) - tin(ix-1, iy)) / dx
            else
                cnvtx = uut * (tin(ix+1, iy) - tin(ix, iy)) / dx
            end if

            if (vvt >= 0.0d0) then
                cnvty = vvt * (tin(ix, iy) - tin(ix, iy-1)) / dy
            else
                cnvty = vvt * (tin(ix, iy+1) - tin(ix, iy)) / dy
            end if

            ! 拡散項(dift)の計算
            dift = alp * ((tin(ix-1, iy  ) - 2.0d0 * tin(ix, iy) + tin(ix+1, iy  )) / dx ** 2 &
                        + (tin(ix,   iy-1) - 2.0d0 * tin(ix, iy) + tin(ix,   iy+1)) / dy ** 2)
            ! 次の時間のtの計算
            tout(ix, iy) = tin(ix, iy) + dt * (-cnvtx - cnvty + dift)
        end subroutine calct
    end subroutine calc_temperature


    ! *********************************************************************
    ! *                    速度の境界条件の処理
    ! *********************************************************************
    subroutine bind_velocity(u, v)
        real(8), intent(inout) :: u(0:, 0:), v(0:, 0:)
        integer :: ix, iy

        ! u（右面）
        u(nx, 1:ny) = u(nx-1, 1:ny)
        ! u（左面）
        u(0, 1:ny) = 1.0d0
        ! u（上面）
        ! u(0:nx, ny + 1) = -u(0:nx, ny)
        u(0:nx, ny+1) = u(0:nx, ny)
        ! u（下面）
        ! u(0:nx, 0) = -u(0:nx, 1)
        u(0:nx, 0) = u(0:nx, 1)
        ! v（右面）
        v(nx+1, 1:ny-1) = v(nx, 1:ny-1)
        ! v（左面）
        v(0, 1:ny-1) = 0.0d0
        ! v（上面）
        v(0:nx+1, ny) = 0.0d0
        ! v(0:nx+1, ny) =  v(0:nx+1, ny-1)
        ! v（下面）
        v(0:nx+1, 0) = 0.0d0
        ! v(0:nx+1, 0) = v(0:nx+1, 1)
        ! 物体内部の処理
        u(0:nx, 0:ny) = u(0:nx, 0:ny) * ((fi(0:nx, 0:ny) + fi(1:nx+1, 0:ny  )) * 0.5d0)
        v(0:nx, 0:ny) = v(0:nx, 0:ny) * ((fi(0:nx, 0:ny) + fi(0:nx,   1:ny+1)) * 0.5d0)
    end subroutine bind_velocity


    ! *********************************************************************
    ! *                  温度の境界条件の処理
    ! *********************************************************************
    subroutine bind_temperature(t)
        real(8), intent(inout) :: t(0:, 0:)
        integer :: ix, iy

        ! 右面
        t(nx+1, 0:ny+1) = -1.0d0 - t(nx, 0:ny+1)
        ! 左面
        t(0, 0:ny+1) = 1.0d0 - t(1, 0:ny+1)
        ! 上面
        t(1:nx, ny+1) = t(1:nx, ny)
        ! 下面
        t(1:nx, 0) = t(1:nx, 1)
    end subroutine bind_temperature


    ! ********************************************************************
    ! tecplot用データ出力
    ! ********************************************************************
    subroutine output_tecplt(filename)
        character(*), intent(in) :: filename
        real(8) :: x, y, u, v, p, f, w
        integer :: ix, iy, unit

        ! write(*, *) 'output_tecplt'

        open(newunit=unit, file=filename)
            write (unit, *) 'variables = "x", "y", "u", "v", "p", "f", "w"'
            write (unit, "(1h , 'zone i=', i5, ', j=', i5, ', f=point')") nx, ny
            do iy = 1, ny
                do ix = 1, nx
                    ! write(*, *) 'iy=', iy
                    x = dx * dble(ix)
                    y = dy * dble(iy)
                    ! u = un(ix, iy)
                    ! v = vn(ix, iy)
                    ! uとvの位置をずらす
                    u = (un(ix, iy) + un(ix-1, iy  )) * 0.5d0
                    v = (vn(ix, iy) + vn(ix  , iy-1)) * 0.5d0
                    p = po(ix, iy)
                    f = fi(ix, iy)
                    w = (un(ix-1, iy-1) + un(ix    , iy-1) - un(ix-1, iy+1) - un(ix  , iy+1) &
                       + vn(ix+1, iy-1) + vn(ix+1, iy    ) - vn(ix-1, iy  ) - vn(ix-1, iy-1)) * 0.25d0
                    write(unit, "(7es11.3)") x, y, u, v, p, f, w
                end do
            end do
        close(unit)
    end subroutine output_tecplt


    ! ********************************************************************
    ! vtk binary用データ出力
    ! ********************************************************************
    subroutine output_vtk
        integer      :: nx1, ny1, i, k
        integer(1)   :: kai(2)
        character*30 :: str0g30
        character*20 :: fname

        kai(1) = 13 ! 改行記号の設定
        kai(2) = 10

        nx1 = nx + 1
        ny1 = ny + 1

        open(99, file="res.vtk", form="unformatted", convert="big_endian", status="replace", access='stream')
            ! ヘッダを書き込む
            ! 文字列の書き込みができない。
            write(99) "# vtk datafile version 3.0", kai
            write(99) "vtk output", kai
            write(99) "binary", kai
            write(99) "dataset structured_grid", kai
            write(str0g30, "(i6, i6, i6)") nx1, ny1 , 1
            write(99) "dimensions ", str0g30, kai
            write(str0g30, "(i6)") nx1 * ny1
            write(99) "points ", str0g30, " float", kai

            ! 格子点を記録します。
            do k = 0, ny
                do i = 0, nx ! vtk出力ではループをiから回す必要がある。
                    write(99) real(dx * dble(i)), real(dy * dble(k)), real(0)
                end do
            end do

            ! じゃんじゃん書き込みます
            write(str0g30, "(a, i8)") "point_data ", nx1 * ny1
            write(99) str0g30, kai
            write(99) "vectors velocity float", kai

            do k = 0, ny
                do i = 0, nx ! vtk出力ではループをiから回す必要がある。
                    write(99) real((un(i, k) + un(i  , k+1)) / 2)
                    write(99) real((vn(i, k) + vn(i+1, k  )) / 2)
                    write(99) real(0)
                end do
            end do

            write(99) "scalars t float", kai
            write(99) "lookup_table default", kai

            do k = 0, ny
                do i = 0, nx ! vtk出力ではループをiから回す必要がある。
                    write(99) real((fi(i, i)  + fi(i+1, k) + fi(i, i+1) + fi(i+1, k+1)) / 4)
                end do
            end do
        close(99)
    end subroutine output_vtk


    ! ********************************************************************
    ! 生バイナリデータ出力
    ! ********************************************************************
    subroutine output_raw(filename)
        character(*), intent(in) :: filename
        integer :: unit

        open(newunit=unit, file=filename, access="stream", status="replace")
        write(unit) shape(un), 4, real(un)
        write(unit) shape(vn), 4, real(vn)
        write(unit) shape(po), 4, real(po)
        close(unit)
    end subroutine output_raw
end module hsmac2d


program main_program
    use hsmac2d
    implicit none

    call main
end program main_program
