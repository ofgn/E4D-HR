module mod_con

    use vars
    use input
    use reorder_mesh

    implicit none
    integer, dimension(:), allocatable :: wrows, wcols, par_map
    real, dimension(:), allocatable :: Wm, sigma_par, rsigma_par, sigi, v, Wmw
    real, dimension(:, :, :), allocatable :: WA
    real, dimension(3) :: CV, CK, CJ
    integer :: ccount, npar, ji_count, cur_el, ncgrad, nnz_cgrad
    integer, dimension(:, :), allocatable :: neighbors
    logical :: itst = .true.
    logical :: iitst = .false.
    logical, dimension(:), allocatable :: use4cgrad
contains

    !_____________________________________________________________________
    subroutine build_WmII
        implicit none
        integer :: n_homo, i, j, k, nzmax, n_homo_elem, count, row, col, i1, i2, lzone
        integer :: znn, nbr, nbr2, x1, x2, x3, z1, z2, z3, ii, nzn, nelem_tmp, cnt, jnk
        integer, dimension(:), allocatable :: homo_par
        integer, dimension(4) :: neighbors_temp
        real :: midxi, midyi, midzi, midxn, midyn, midzn, rx, ry, rz, r, x, rho, eps
        logical :: ffound, extfil
        character*80 :: extfile

        if (invi .and. itst) then
            itst = .false.
            iitst = .true.
            if (allocated(Wm)) deallocate (Wm)

        elseif (iitst .and. .not. invi) then
            itst = .true.
            iitst = .false.
            if (allocated(Wm)) deallocate (Wm)
        end if

        if (allocated(Wm)) then
            do i = 1, ccount
                call comp_Wm(i, Wm(i))
            end do
            return
        else

            if (allocated(neighbors)) deallocate (neighbors)
            if (allocated(rblock)) deallocate (rblock)

       !!Check to see if there are NN regularized zones. If so, read neighbors
            !if(sum(smetric(:,2))>nrz) then
            allocate (neighbors(n_elements, 4))
            do i = 1, 80
                if (cfg_filename(i:i) == '.') then
                    open (21, file=cfg_filename(1:i + 1)//".neigh", status='old', action='read')

                    if (allocated(element_map)) then
                        !if element_map is allocated then there are inactive elements
                        !and we need to re_map the element neighbors
                        read (21, *) nelem_tmp
                        cnt = 0
                        do j = 1, nelem_tmp
                            read (21, *) jnk, neighbors_temp(1:4)
                            if (element_map(j) .ne. 0) then
                                cnt = cnt + 1
                                do k = 1, 4
                                    if (neighbors_temp(k) > 0) then
                                        neighbors(cnt, k) = element_map(neighbors_temp(k))
                                    else
                                        neighbors(cnt, k) = neighbors_temp(k)
                                    end if
                                end do
                            end if
                        end do
                        close (21)

                    else
                        read (21, *)
                        do j = 1, n_elements
                            read (21, *) k, neighbors(j, 1:4)
                        end do
                        close (21)
                    end if
                    exit
                end if
            end do
            nzn = maxval(zones)

            !count the number of constraint equations
            row = 0

            do i = 1, n_elements
                do j = 1, nrz
             !!12 and 13 are joint inversion constraints, skip them for now
                    if (smetric(j, 2) .eq. 12 .or. smetric(j, 2) .eq. 13) goto 10
                    if (smetric(j, 1) == nzn + 1 .and. i == 1) then
                        !this constraint block is an external constraint block
                        !read the number of constraint equations but only do this once
                        open (11, file='external_files.txt', status='old', action='read', IOSTAT=io_stat)
                        if (io_stat .ne. 0) then
                            write (*, *) 'Did not find the external constraint list file external_confiles.txt'
                            write (*, *) 'Not implementing external constraint'
                            close (11)
                        else
                            !loop over the list file
                            ffound = .false.
                            do while (.not. ffound)
                                read (11, *, IOSTAT=io_stat) ii, extfile
                                if (ii == j) then
                                    ffound = .true.
                                    close (11)
                                end if
                            end do
                            if (ffound) then
                                open (11, file=extfile, status='old', action='read')
                                read (11, *) ii
                                ccount = ccount + ii
                                !write(*,*) 'adding ',ii,' constraints for block ',j
                            end if
                        end if
                    else
                        if (smetric(j, 2) > 0) then
                            if (smetric(j, 1) == zones(i)) then
                                if (smetric(j, 2) == 3 .or. smetric(j, 2) == 4) then
                                    ccount = ccount + 1

                                else
                                    do k = 1, 4
                                        nbr = neighbors(i, k)

                                        if (zones(nbr) .ne. zones(i)) then
                                            do ii = 2, zone_links(j, 1) + 1
                                                if (zone_links(j, ii) == zones(nbr)) then
                                                    ccount = ccount + 1
                                                end if
                                            end do
                                        else
                                            ccount = ccount + 1
                                        end if

                                    end do
                                end if
                            end if
                        end if
                    end if
10                  continue
                end do
            end do

            allocate (Wm(ccount), rblock(ccount, 3))
            Wm = 0.0
            ccount = 0
            ji_count = 0
            row = 0
            col = 0

            rblock = 0
            do i = 1, n_elements
                do j = 1, nrz
                    if (smetric(j, 2) .eq. 12 .or. smetric(j, 2) .eq. 13) goto 20

                    if (smetric(j, 1) == nzn + 1 .and. i == 1) then
                        !this constraint block is an external constraint block
                        !read the constraints only once
                        open (11, file='external_files.txt', status='old', action='read', IOSTAT=io_stat)
                        if (io_stat .ne. 0) then
                            close (11)
                        else
                            ffound = .false.
                            do while (.not. ffound)
                                read (11, *, IOSTAT=io_stat) ii, extfile
                                if (ii == j) then
                                    ffound = .true.
                                    close (11)
                                    !write(*,*) 'Adding external constraints for block: ',j
                                end if
                            end do
                            if (ffound) then
                                open (11, file=extfile, status='old', action='read')
                                read (11, *) ii
                                do k = 1, ii
                                    ccount = ccount + 1
                                    read (11, *) rblock(ccount, 1:2)
                                    rblock(ccount, 3) = j
                                end do
                                close (11)
                                write (*, *) ' ADDING EXTERNAL CONSTRAINTS FROM FILE: ', trim(extfile), ' FOR REG BLOCK: ', j
                            end if
                        end if

                    else

                        if (smetric(j, 2) > 0) then
                            if (smetric(j, 1) == zones(i)) then
                                if (smetric(j, 2) == 3 .or. smetric(j, 2) == 4) then
                                    ccount = ccount + 1
                                    rblock(ccount, 1) = i
                                    rblock(ccount, 3) = j

                                else
                                    do k = 1, 4
                                        nbr = neighbors(i, k)

                                        if (zones(nbr) .ne. zones(i)) then
                                            do ii = 2, zone_links(j, 1) + 1
                                                if (zone_links(j, ii) == zones(nbr)) then
                                                    ccount = ccount + 1
                                                    rblock(ccount, 1) = i
                                                    rblock(ccount, 2) = nbr
                                                    rblock(ccount, 3) = j
                                                end if
                                            end do
                                        else
                                            if (smetric(j, 2) < 9 .or. smetric(j, 2) > 10) then
                                                ccount = ccount + 1
                                                rblock(ccount, 1) = i
                                                rblock(ccount, 2) = nbr
                                                rblock(ccount, 3) = j
                                            end if
                                        end if

                                    end do
                                end if

                            end if
                        end if
                    end if
20                  continue
                end do
            end do
        end if

        !count the joint inversion constraints if necessary
!!$    if(cgmin_flag(1)) then
!!$       do i=1,nelem
!!$          do j=1,nrz
!!$             if(smetric(j,2).eq.12 .or. smetric(j,2).eq.13) then
!!$                if(smetric(j,1) .eq. zones(i)) then
!!$                   !one constraint for each component of the
!!$                   !cross-gradient vector
!!$                   do k=1,3
!!$                      ccount=ccount+1
!!$                      ji_count = ji_count+1
!!$                      rblock(ccount,1) = i
!!$                      rblock(ccount,2) = ji_count
!!$                      rblock(ccount,3) = j
!!$                   end do
!!$                end if
!!$             end if
!!$          end do
!!$       end do
!!$    end if

        !set the on/off flag for element estimation
        if (.not. allocated(J_on_off)) then
            allocate (J_on_off(n_elements))
            J_on_off = .true.
            do i = 1, nrz

                if (smetric(i, 2) == 0 .and. smetric(i, 1) .ne. nzn + 1) then
                    do j = 1, n_elements
                        if (smetric(i, 1) == zones(j)) then
                            J_on_off(j) = .false.
                        end if
                    end do
                end if

                if (smetric(i, 2) == 0 .and. smetric(i, 1) == nzn + 1) then
                    !these are specified in and external constraint file
                    open (11, file='external_files.tmp', status='old', action='read')
                    ffound = .false.
                    do while (.not. ffound)
                        read (11, *, IOSTAT=io_stat) ii, extfile
                        if (ii == i) then
                            ffound = .true.
                            close (11)
                        end if
                    end do
                    if (ffound) then
                        open (11, file=extfile, status='old', action='read')
                        read (11, *) ii
                        do k = 1, ii
                            read (11, *) j
                            J_on_off(j) = .false.
                        end do
                        close (11)
                    end if
                end if

            end do
        end if

        do i = 1, ccount
            call comp_Wm(i, Wm(i))
        end do

    end subroutine build_WmII
    !_____________________________________________________________________

    !_____________________________________________________________________
    subroutine comp_Wm(indx, wm)
        implicit none
        integer, intent(in) :: indx
        real, intent(out) :: wm
        integer :: rbi, eli, eln, i
        real :: X, mn, sd
        real :: midxi, midyi, midzi, midxn, midyn, midzn, rx, ry, rz, r, rho, eps

        rbi = abs(rblock(indx, 3))
        select case (smetric(rbi, 2))

        case (1)
            if (invi) then
                X = (log(phase(rblock(indx, 1))) - log(phase(rblock(indx, 2))))
            else
                X = (log(sigma_re(rblock(indx, 1))) - log(sigma_re(rblock(indx, 2))))
            end if

        case (2)
            if (invi) then
                !X=abs(log(sigmai(rblock(indx,1)))-log(sigmai(rblock(indx,2))))
                X = abs(log(phase(rblock(indx, 1))) - log(phase(rblock(indx, 2))))
            else
                X = abs(log(sigma_re(rblock(indx, 1))) - log(sigma_re(rblock(indx, 2))))
            end if

        case (3)
            select case (smetric(rbi, 3))
            case (0)
                if (invi) then
                    X = log(phase(rblock(indx, 1))) - (C_targ(rbi))
                else
                    X = log(sigma_re(rblock(indx, 1))) - (C_targ(rbi))
                end if
            case (1)
                if (invi) then
                    X = log(phase(rblock(indx, 1))) - log(refsig(rblock(indx, 1)))
                else
                    X = log(sigma_re(rblock(indx, 1))) - log(refsig(rblock(indx, 1)))
                end if
            case (2)
                if (invi) then
                    X = log(phase(rblock(indx, 1))) - log(prefsig(rblock(indx, 1)))
                else
                    X = log(sigma_re(rblock(indx, 1))) - log(prefsig(rblock(indx, 1)))
                end if
            end select

        case (4)
            select case (smetric(rbi, 3))
            case (0)
                if (invi) then
                    X = abs(log(phase(rblock(indx, 1))) - (C_targ(rbi)))
                else
                    X = abs(log(sigma_re(rblock(indx, 1))) - (C_targ(rbi)))
                end if
            case (1)
                if (invi) then
                    X = abs(log(phase(rblock(indx, 1))) - log(refsig(rblock(indx, 1))))
                else
                    X = abs(log(sigma_re(rblock(indx, 1))) - log(refsig(rblock(indx, 1))))
                end if
            case (2)
                if (invi) then
                    X = abs(log(phase(rblock(indx, 1))) - log(prefsig(rblock(indx, 1))))
                else
                    X = abs(log(sigma_re(rblock(indx, 1))) - log(prefsig(rblock(indx, 1))))
                end if
            end select

        case (5)
            if (invi) then
                X = (log(phase(rblock(indx, 1))) - log(phase(rblock(indx, 2))))
            else
                X = (log(sigma_re(rblock(indx, 1))) - log(sigma_re(rblock(indx, 2))))
            end if

            eli = rblock(indx, 1)
            eln = rblock(indx, 2)
            midxi = 0.25*sum(nodes(elements(eli, 1:4), 1))
            midyi = 0.25*sum(nodes(elements(eli, 1:4), 2))
            midzi = 0.25*sum(nodes(elements(eli, 1:4), 3))
            midxn = 0.25*sum(nodes(elements(eln, 1:4), 1))
            midyn = 0.25*sum(nodes(elements(eln, 1:4), 2))
            midzn = 0.25*sum(nodes(elements(eln, 1:4), 3))
            !rx=((midxi-midxn)**2)
            !ry=((midyi-midyn)**2)
            !rz=((midzi-midzn)**2)
            !r=(rx+ry+rz)
            rx = ((midxi - midxn))
            ry = ((midyi - midyn))
            rz = ((midzi - midzn))
            r = sqrt((rx**2 + ry**2 + rz**2))
            rx = rx/r
            ry = ry/r
            rz = rz/r

        case (6)
            if (invi) then
                X = abs(log(phase(rblock(indx, 1))) - log(phase(rblock(indx, 2))))
            else
                X = abs(log(sigma_re(rblock(indx, 1))) - log(sigma_re(rblock(indx, 2))))
            end if

            eli = rblock(indx, 1)
            eln = rblock(indx, 2)
            midxi = 0.25*sum(nodes(elements(eli, 1:4), 1))
            midyi = 0.25*sum(nodes(elements(eli, 1:4), 2))
            midzi = 0.25*sum(nodes(elements(eli, 1:4), 3))
            midxn = 0.25*sum(nodes(elements(eln, 1:4), 1))
            midyn = 0.25*sum(nodes(elements(eln, 1:4), 2))
            midzn = 0.25*sum(nodes(elements(eln, 1:4), 3))

            rx = ((midxi - midxn))
            ry = ((midyi - midyn))
            rz = ((midzi - midzn))
            r = sqrt((rx**2 + ry**2 + rz**2))
            rx = rx/r
            ry = ry/r
            rz = rz/r

       !!extra for fracture testing
            !rx=abs((midxi-midxn))
            !ry=abs((midyi-midyn))
            !rz=abs((midzi-midzn))
            !rx = rx/r
            !ry = ry/r
            !rz = rz/r

            !r=sqrt(midxi**2 + midyi**2 + midzi**2)
            !midxi=midxi/r
            !midyi=midyi/r
            !midzi=midzi/r
            !r=rx*midxi + ry*midyi + rz*midzi

        case (7)
            select case (smetric(rbi, 3))
            case (0)
                if (invi) then
                    X = (log(phase(rblock(indx, 1))) - C_targ(rbi)) - &
                        (log(phase(rblock(indx, 2))) - C_targ(rbi))
                else
                    X = (log(sigma_re(rblock(indx, 1))) - C_targ(rbi)) - &
                        (log(sigma_re(rblock(indx, 2))) - C_targ(rbi))
                end if

            case (1)
                if (invi) then
                    X = (log(phase(rblock(indx, 1))) - log(refsig(rblock(indx, 1)))) - &
                        (log(phase(rblock(indx, 2))) - log(refsig(rblock(indx, 2))))

                else
                    X = (log(sigma_re(rblock(indx, 1))) - log(refsig(rblock(indx, 1)))) - &
                        (log(sigma_re(rblock(indx, 2))) - log(refsig(rblock(indx, 2))))
                end if
            case (2)
                if (invi) then
                    X = (log(phase(rblock(indx, 1))) - log(prefsig(rblock(indx, 1)))) - &
                        (log(phase(rblock(indx, 2))) - log(prefsig(rblock(indx, 2))))
                else
                    X = (log(sigma_re(rblock(indx, 1))) - log(prefsig(rblock(indx, 1)))) - &
                        (log(sigma_re(rblock(indx, 2))) - log(prefsig(rblock(indx, 2))))
                end if
            end select

        case (8)
            select case (smetric(rbi, 3))
            case (0)
                if (invi) then
                    X = abs((log(phase(rblock(indx, 1))) - C_targ(rbi)) - &
                            (log(phase(rblock(indx, 2))) - C_targ(rbi)))
                else
                    X = abs((log(sigma_re(rblock(indx, 1))) - C_targ(rbi)) - &
                            (log(sigma_re(rblock(indx, 2))) - C_targ(rbi)))
                end if
            case (1)
                if (invi) then
                    X = abs((log(phase(rblock(indx, 1))) - log(refsig(rblock(indx, 1)))) - &
                            (log(phase(rblock(indx, 2))) - log(refsig(rblock(indx, 2)))))
                else
                    X = abs((log(sigma_re(rblock(indx, 1))) - log(refsig(rblock(indx, 1)))) - &
                            (log(sigma_re(rblock(indx, 2))) - log(refsig(rblock(indx, 2)))))
                end if
            case (2)
                if (invi) then
                    X = abs((log(phase(rblock(indx, 1))) - log(prefsig(rblock(indx, 1)))) - &
                            (log(phase(rblock(indx, 2))) - log(prefsig(rblock(indx, 2)))))
                else
                    X = abs((log(sigma_re(rblock(indx, 1))) - log(prefsig(rblock(indx, 1)))) - &
                            (log(sigma_re(rblock(indx, 2))) - log(prefsig(rblock(indx, 2)))))
                end if
            end select

        case (9)
            if (invi) then
                X = (log(phase(rblock(indx, 1))) - log(phase(rblock(indx, 2))))
            else
                X = (log(sigma_re(rblock(indx, 1))) - log(sigma_re(rblock(indx, 2))))
            end if

        case (10)
            if (invi) then
                X = abs(log(phase(rblock(indx, 1))) - log(phase(rblock(indx, 2))))
            else
                X = abs(log(sigma_re(rblock(indx, 1))) - log(sigma_re(rblock(indx, 2))))
            end if

        case (11) !test case for horizontal radial regularization TCJ 11/04/15
            if (invi) then
                X = abs(log(phase(rblock(indx, 1))) - log(phase(rblock(indx, 2))))
            else
                X = abs(log(sigma_re(rblock(indx, 1))) - log(sigma_re(rblock(indx, 2))))
            end if
            eli = rblock(indx, 1)
            eln = rblock(indx, 2)

            midxi = 0.25*sum(nodes(elements(eli, 1:4), 1))!-zwts(rbi,2)
            midyi = 0.25*sum(nodes(elements(eli, 1:4), 2))!-zwts(rbi,3)
            !midzi=0.25*sum(nodes(elements(eli,1:4),3))-zwts(rbi,4)
            midxn = 0.25*sum(nodes(elements(eln, 1:4), 1))!-zwts(rbi,2)
            midyn = 0.25*sum(nodes(elements(eln, 1:4), 2))!-zwts(rbi,3)
            !midzn=0.25*sum(nodes(elements(eln,1:4),3))-zwts(rbi,4)
            !get the normal from neighbor to me
            midxn = midxn - midxi
            midyn = midyn - midyi
            r = sqrt(midxn**2 + midyn**2)
            if (r == 0) then
                r = 1
            else
                midxn = midxn/r
                midyn = midyn/r

                r = sqrt(midxi**2 + midyi**2); 
                if (r == 0) then
                    r = 1
                else
                    midxi = midxi/r
                    midyi = midyi/r
                    r = (midxn*midxi + midyn*midyi)**2
                end if
            end if

        case DEFAULT

        end select

        mn = Fw_parm(rbi, 1)
        sd = Fw_parm(rbi, 2)

        select case (Fw_type(rbi))

        case (1)
            wm = .5*(1 - erf((X - mn)/sqrt(2*sd**2)))
        case (2)
            wm = .5*(1 + erf((X - mn)/sqrt(2*sd**2)))

        case (3)
            wm = 1 - exp(-((X - mn)**2)/(2*sd**2))

        case (4)
            wm = exp(-((X - mn)**2)/(2*sd**2))

        case (5)
            if ((X - mn) < 0) then
                wm = sd**(-2)
            else
                wm = sd**(2)*((X - mn)**2 + sd**2)**(-2)
            end if

        case (6)
            if ((X - mn) > 0) then
                wm = sd**(-2)
            else
                wm = sd**(2)*((X - mn)**2 + sd**2)**(-2)
            end if

        case DEFAULT
        end select

        if (smetric(rbi, 2) == 5 .or. smetric(rbi, 2) == 6) then
            wm = (1 - abs(rx*zwts(rbi, 2) + ry*zwts(rbi, 3) + rz*zwts(rbi, 4)))**2

        end if
        if (smetric(rbi, 2) == 11) then
            wm = r*wm
        end if
        wm = wm*zwts(rbi, 1)

    !! test adjustment with face area
        !wm=wm*Wmw(indx)

    end subroutine comp_Wm
    !_____________________________________________________________________

    !_____________________________________________________________________
    subroutine build_BM
        implicit none
        write (*, *) my_rank, 'Building BM'
    end subroutine build_BM
    !_____________________________________________________________________

    !_____________________________________________________________________
    subroutine build_V
        implicit none
        write (*, *) my_rank, 'Building V'
    end subroutine build_V
    !_____________________________________________________________________

    !_____________________________________________________________________
    !_____________________________________________________________________

    !_____________________________________________________________________
    subroutine build_evol
        !build the vector of element volumes
        implicit none
        integer :: i, j, k
        real*8 :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
        real*8, dimension(4, 4) :: a

        if (allocated(evol)) deallocate (evol)
        allocate (evol(n_elements))

        do k = 1, n_elements
            x1 = dble(nodes(elements(k, 1), 1))
            y1 = dble(nodes(elements(k, 1), 2))
            z1 = dble(nodes(elements(k, 1), 3))

            x2 = dble(nodes(elements(k, 2), 1))
            y2 = dble(nodes(elements(k, 2), 2))
            z2 = dble(nodes(elements(k, 2), 3))

            x3 = dble(nodes(elements(k, 3), 1))
            y3 = dble(nodes(elements(k, 3), 2))
            z3 = dble(nodes(elements(k, 3), 3))

            x4 = dble(nodes(elements(k, 4), 1))
            y4 = dble(nodes(elements(k, 4), 2))
            z4 = dble(nodes(elements(k, 4), 3))

            a(1:4, 1) = (/x1, x2, x3, x4/)
            a(1:4, 2) = (/y1, y2, y3, y4/)
            a(1:4, 3) = (/z1, z2, z3, z4/)
            a(1:4, 4) = (/1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00/)

            evol(k) = abs(det4(a))/6.0E+00
        end do

    end subroutine build_evol
    !_____________________________________________________________________
    !___________________________________________________________________________
    function det4(a1)

        !computes the determinant of a 4x4 matrix
        implicit none
        real*8 :: a1(4, 4)
        real*8 :: a(4, 4)
        real*8 :: det4
        integer ::i
        real*8 :: c1, c2, c3, c4

        a = dble(a1)
        c1 = a(1, 1)*(a(2, 2)*(a(3, 3)*a(4, 4) - a(3, 4)*a(4, 3)) - a(2, 3)*&
             &(a(3, 2)*a(4, 4) - a(3, 4)*a(4, 2)) + a(2, 4)*(a(3, 2)*a(4, 3) - a(3, 3)*a(4, 2)))

        c2 = a(1, 2)*(a(2, 1)*(a(3, 3)*a(4, 4) - a(3, 4)*a(4, 3)) - a(2, 3)*&
             &(a(3, 1)*a(4, 4) - a(3, 4)*a(4, 1)) + a(2, 4)*(a(3, 1)*a(4, 3) - a(3, 3)*a(4, 1)))

        c3 = a(1, 3)*(a(2, 1)*(a(3, 2)*a(4, 4) - a(3, 4)*a(4, 2)) - a(2, 2)*&
             &(a(3, 1)*a(4, 4) - a(3, 4)*a(4, 1)) + a(2, 4)*(a(3, 1)*a(4, 2) - a(3, 2)*a(4, 1)))

        c4 = a(1, 4)*(a(2, 1)*(a(3, 2)*a(4, 3) - a(3, 3)*a(4, 2)) - a(2, 2)*&
             &(a(3, 1)*a(4, 3) - a(3, 3)*a(4, 1)) + a(2, 3)*(a(3, 1)*a(4, 2) - a(3, 2)*a(4, 1)))

        det4 = real(c1 - c2 + c3 - c4)

        return
    end function det4
    !__________________________________________________________________________
    !_____________________________________________________________________

    !_____________________________________________________________________
    subroutine get_areas_normal_vol(el, A, N, vol)
        implicit none
        integer, intent(in) :: el
        real*8, dimension(4), intent(out) :: A
        real*8, dimension(4, 3), intent(out) :: N
        real*8, intent(out) :: vol
        real*8, dimension(3) :: u, v, x, elmid, pmid
        real*8 :: vlen
        integer :: i

        do i = 1, 3
            elmid(i) = .25*sum(nodes(elements(el, :), i))
        end do

        !face 1,2,3
        u = nodes(elements(el, 2), :) - nodes(elements(el, 1), :)
        v = nodes(elements(el, 3), :) - nodes(elements(el, 1), :)
        x = cross_prod(u, v)
        vlen = sqrt(dot_product(x, x))
        A(1) = .5*vlen
        N(1, :) = x/vlen

        !make sure the normal is pointing outward
        do i = 1, 3
            pmid(i) = (nodes(elements(el, 1), i) + nodes(elements(el, 2), i) + &
                       nodes(elements(el, 3), i))/3
        end do
        u = pmid - elmid
        if (dot_product(u, N(1, :)) < 0) N(1, :) = -N(1, :)

        !get the volume while we're here
        x = nodes(elements(el, 4), :) - nodes(elements(el, 1), :)
        vol = sqrt(dot_product(u, cross_prod(v, x)))/6

        !face 1,2,4
        u = nodes(elements(el, 2), :) - nodes(elements(el, 1), :)
        v = nodes(elements(el, 4), :) - nodes(elements(el, 1), :)
        x = cross_prod(u, v)
        vlen = sqrt(dot_product(x, x))
        A(2) = .5*vlen
        N(2, :) = x/vlen

        !make sure the normal is pointing outward
        do i = 1, 3
            pmid(i) = (nodes(elements(el, 1), i) + nodes(elements(el, 2), i) + &
                       nodes(elements(el, 4), i))/3
        end do
        u = pmid - elmid
        if (dot_product(u, N(2, :)) < 0) N(2, :) = -N(2, :)

        !face 1,3,4
        u = nodes(elements(el, 3), :) - nodes(elements(el, 1), :)
        v = nodes(elements(el, 4), :) - nodes(elements(el, 1), :)
        x = cross_prod(u, v)
        vlen = sqrt(dot_product(x, x))
        A(3) = .5*vlen
        N(3, :) = x/vlen

        !make sure the normal is pointing outward
        do i = 1, 3
            pmid(i) = (nodes(elements(el, 1), i) + nodes(elements(el, 3), i) + &
                       nodes(elements(el, 4), i))/3
        end do
        u = pmid - elmid
        if (dot_product(u, N(3, :)) < 0) N(3, :) = -N(3, :)

        !face 2,3,4
        u = nodes(elements(el, 3), :) - nodes(elements(el, 2), :)
        v = nodes(elements(el, 4), :) - nodes(elements(el, 2), :)
        x = cross_prod(u, v)
        vlen = sqrt(dot_product(x, x))
        A(4) = .5*vlen
        N(4, :) = x/vlen

        !make sure the normal is pointing outward
        do i = 1, 3
            pmid(i) = (nodes(elements(el, 2), i) + nodes(elements(el, 3), i) + &
                       nodes(elements(el, 4), i))/3
        end do
        u = pmid - elmid
        if (dot_product(u, N(4, :)) < 0) N(4, :) = -N(4, :)

    end subroutine get_areas_normal_vol
    !_____________________________________________________________________

    !_____________________________________________________________________
    subroutine comp_Wmw
        implicit none
        !computes the area of the face between neighbor tets and assigns to Wmw
        integer :: i, j, k, icnt
        integer, dimension(3) :: n_indexes
        real, dimension(3) :: u, v, cp

        do i = 1, ccount
            !compute only if this is a spatial constraing (i.e. has a neighbor specified)
            if (rblock(i, 2) .ne. 0) then
                icnt = 0
                do j = 1, 4
                    do k = 1, 4
                        if (elements(rblock(i, 1), j) == elements(rblock(i, 2), k)) then
                            icnt = icnt + 1
                            n_indexes(icnt) = elements(rblock(i, 1), j)
                            exit
                        end if
                    end do
                    if (icnt == 3) then
                        exit
                    end if
                end do
                if (icnt .ge. 3) then
                    u = nodes(n_indexes(2), :) - nodes(n_indexes(1), :)
                    v = nodes(n_indexes(3), :) - nodes(n_indexes(1), :)
                    cp = 0.5*cross_prod(dble(u), dble(v))
                    Wmw(i) = 0.5*sqrt(cp(1)**2 + cp(2)**2 + cp(3)**2)
                end if
            end if
        end do
        !adjust Wmw so that the average value is 1.0
        Wmw = Wmw*ccount/sum(Wmw)

    end subroutine comp_Wmw
    !_____________________________________________________________________
    !_____________________________________________________________________
    function cross_prod(u, v)
        real*8, dimension(3) :: cross_prod, u, v
        cross_prod(1) = u(2)*v(3) - u(3)*v(2)
        cross_prod(2) = u(3)*v(1) - u(1)*v(3)
        cross_prod(3) = u(1)*v(2) - u(2)*v(1)
    end function cross_prod
    !______________________________________________________________________

    !_____________________________________________________________________
    subroutine build_Wm

        implicit none
        integer :: n_homo, i, j, k, nzmax, n_homo_elem, count, row, col, i1, i2, lzone
        integer :: znn, nbr, nbr2, x1, x2, x3, z1, z2, z3
        integer, dimension(:), allocatable :: homo_par
        integer, dimension(:, :), allocatable :: neighbors
        real :: midxi, midyi, midzi, midxn, midyn, midzn, rx, ry, rz, r, x, rho, eps

        if (gs_flag) then
            if (allocated(Wm)) then
                return
            end if
            allocate (Wm(n_elements), wrows(n_elements), wcols(n_elements))
            ccount = n_elements
            do i = 1, n_elements
                wrows(i) = i
                wcols(i) = i
                Wm(i) = 1
            end do
            return
        end if

    !!If Wm is already constructed then there is nothing to do
        if (allocated(Wm)) then
            !return
            deallocate (par_map)
            deallocate (sigma_par, rsigma_par)
            deallocate (wrows, wcols, Wm)
        end if

    !!check to see if there are heterogeneous zones. if so, read in the neighbors
        if (sum(reg_opt) > 0) then
            allocate (neighbors(n_elements, 4))
            do i = 1, 80
                if (cfg_filename(i:i) == '.') then
                    open (21, file=cfg_filename(1:i + 1)//".neigh", status='old', action='read'); 
                    read (21, *)
                    do j = 1, n_elements
                        read (21, *) k, neighbors(j, 1:4)
                    end do
                    goto 25
                end if
            end do
        end if
25      continue
        close (21)

    !!send a warning if there is a discrepancy be the number of zones
    !!and the number of zones which are regularized
        nzmax = int(maxval(zones))

        if (nzmax .ne. nrz) then
            if (.not. gs_flag) then
                open (51, file='e4d.log', status='old', action='write', position='append')
                write (51, *) '!!!WARNING!!! THERE ARE ', nzmax, ' ZONES IN THE ELEMENTS FILE AND ', nrz, ' ZONES REGULARIZED'
                close (51)
            else
                open (51, file='e4d.log', status='old', action='write', position='append')
                write (51, *) '!!!WARNING!!! THERE AR', nzmax, ' ZONES AND ', nrz, ' ARE REGULARIZED'
                write (51, *) 'MAKES SURE THAT ZONES ', nrz + 1, ' THROUGH ', nzmax, ' HAVE SEMIVARIOGRAM CONSTRAINTS'
                write (51, *) 'THIS IS NOT CHECKED INTERNALLY'
                close (51)
            end if

        end if

    !!determine how many parameters we will estimate
        n_homo = 0
        allocate (homo_par(nrz))
        do i = 1, nrz
            if (reg_opt(i) == 0) then
                n_homo = n_homo + 1
                homo_par(n_homo) = i
            end if
        end do

    !!count the number of elements in homogeneous zones
        n_homo_elem = 0
        do i = 1, n_elements
            if (reg_opt(zones(i)) == 0) then
                n_homo_elem = n_homo_elem + 1
            end if
        end do

    !!allocate the new parameter vector
        npar = n_elements - n_homo_elem + n_homo
        allocate (sigma_par(npar), rsigma_par(npar))

    !!build the vector mapping the conductivities to the inversion parameters
        allocate (par_map(n_elements))
        count = 0
        do i = 1, n_elements
            if (reg_opt(zones(i)) == 0) then
          !!this is an element in a homogeneous zone ... find which zone
                do j = 1, n_homo
                    if (homo_par(j) == zones(i)) then
                        par_map(i) = j
                        if (invi) then
                            sigma_par(j) = sigma_im(i)
                        else
                            sigma_par(j) = sigma_re(i)
                        end if
                    end if
                end do
            else
                count = count + 1
                par_map(i) = n_homo + count
                if (invi) then
                    sigma_par(par_map(i)) = sigma_im(i)
                    !rsigma_par(par_map(i))=refsigi(i)
                else
                    sigma_par(par_map(i)) = sigma_re(i)
                    rsigma_par(par_map(i)) = refsig(i)
                end if
            end if

        end do

    !!now map sigma_par back to sigma to make sure they're consistent
        if (invi) then
            do i = 1, n_elements
                sigma_im(i) = sigma_par(par_map(i))
            end do
        else
            do i = 1, n_elements
                sigma_re(i) = sigma_par(par_map(i))
                refsig(i) = rsigma_par(par_map(i))
            end do
        end if

    !!Now we're ready to build the sparse regularization matrix between inversion
    !!parameters. We start with the homogeneous zones and then move to the
    !!heterogeneous zones. Since we don't know how many elements are in the matrix
    !!apriori, we'll do a false run to count, allocate the matrix, and then do the
    !!computations to fill in the matrix.

        ccount = 0
    !!Homogeneous zone regularization count
        do i = 1, nrz
       !!check to see if this is a homogeneous zone and is linked to other zones
            if (reg_opt(i) == 0 .and. zone_links(i, 1) > 0) then
                do j = 1, zone_links(i, 1)
                    if (reg_opt(zone_links(i, j + 1)) == 0) ccount = ccount + 2
                end do
            end if
        end do

    !!smallest model regularization count
        do i = 1, n_elements
            if (reg_opt(zones(i)) .ge. 1 .and. reg_opt(zones(i)) .le. 4) then
          !!this element belongs to a zone with smallest model smoothing
                do j = 1, 4
                    nbr = neighbors(i, j)                           !!nbr = element index of the neighbor
                    if (nbr <= 0) goto 102                                !!there is no neighbor on face j so skip this neighbor
                                                          !!or we have already covered this pair in nbr<i.
                    znn = zones(nbr)                             !!zone of the neighbor
                    if (zones(i) == znn) goto 101                   !!elements i and nbr are in the same zone

                    do k = 1, zone_links(zones(i), 1)
                        if (zone_links(zones(i), k + 1) == znn) goto 101          !!elements i and nbr are linked neighbors
                    end do
                    goto 102

101                 continue
                    ccount = ccount + 2

102                 continue
                end do
            end if
        end do

    !!if there are homogeneous zone, external constraints are not allowed
    !!check to make sure we don't have both and allocate
        if (allocated(ex_vals) .and. n_homo > 0) then
            open (51, file='e4d.log', status='old', action='write', position='append')
            write (51, *) "Homogeneous zones are specified. External constraints will be ignored"
            close (51)
            nec = 0
        end if
        allocate (Wm(ccount), wrows(ccount), wcols(ccount))

    !!build the regularization matrix
        ccount = 0
        row = 0
        col = 0
    !!build the homogeneous zone regularization
        do i = 1, nrz
       !!check to see if this is a homogeneous zone and is linked to another zones
            if (reg_opt(i) == 0 .and. zone_links(i, 1) > 0) then
                do j = 1, zone_links(i, 1)

             !!find which parameters hold the linked zones
                    lzone = zone_links(i, j + 1)
                    if (reg_opt(lzone) == 0) then

                        do k = 1, n_homo
                            if (homo_par(k) == lzone) i2 = k
                            if (homo_par(k) == i) i1 = k
                        end do

                        row = row + 1; col = i1; ccount = ccount + 1
                        wrows(ccount) = row; wcols(ccount) = col; Wm(ccount) = zwts(i, 1)

                        col = i2; ccount = ccount + 1
                        wrows(ccount) = row; wcols(ccount) = col; Wm(ccount) = -zwts(i, 1)

                    end if
                end do
            end if
        end do

        do i = 1, n_elements
       !!smallest model regularization
            if (reg_opt(zones(i)) .ge. 1 .and. reg_opt(zones(i)) .le. 4) then

                do j = 1, 4
                    nbr = neighbors(i, j)                        !!nbr = element index of the neighbor
                    if (nbr <= 0) goto 111            !!there is no neighbor on face j so skip this neighbor
                                                       !!or we have already covered this pair
                    znn = zones(nbr)                          !!zone of the neighbor
                    if (zones(i) == znn) goto 110                !!elements i and nbr are in the same zone

                    do k = 1, zone_links(zones(i), 1)
                        if (zone_links(zones(i), k + 1) == znn) goto 110   !!elements i and nbr are linked neighbors
                    end do
                    goto 111

110                 continue

                    if (reg_opt(zones(i)) == 1 .or. reg_opt(zones(i)) == 3) then              !!smallest model smoothing
                        row = row + 1; col = par_map(i); ccount = ccount + 1
                        wrows(ccount) = row; wcols(ccount) = col; Wm(ccount) = zwts(zones(i), 1)

                        col = par_map(nbr); ccount = ccount + 1
                        wrows(ccount) = row; wcols(ccount) = col; Wm(ccount) = -zwts(zones(i), 1)
                    end if

                    if (reg_opt(zones(i)) == 4) then
                        eps = zwts(zones(i), 2)
                        x = log(sigma_par(par_map(i)))*zwts(zones(i), 1) - log(sigma_par(par_map(nbr)))*zwts(zones(i), 1)
                        rho = (eps**2)/(x**2 + eps**2)
                        !write(*,*) eps,x,rho
                        row = row + 1; col = par_map(i); ccount = ccount + 1
                        wrows(ccount) = row; wcols(ccount) = col; Wm(ccount) = rho*zwts(zones(i), 1)
                        !write(*,*) wrows(ccount),wcols(ccount), Wm(ccount)
                        col = par_map(nbr); ccount = ccount + 1
                        wrows(ccount) = row; wcols(ccount) = col; Wm(ccount) = -rho*zwts(zones(i), 1)
                        !write(*,*) wrows(ccount),wcols(ccount), Wm(ccount)

                    end if

                    if (reg_opt(zones(i)) == 2) then              !!first order weighted smoothing

                        midxi = 0.25*sum(nodes(elements(i, 1:4), 1))
                        midyi = 0.25*sum(nodes(elements(i, 1:4), 2))
                        midzi = 0.25*sum(nodes(elements(i, 1:4), 3))
                        midxn = 0.25*sum(nodes(elements(nbr, 1:4), 1))
                        midyn = 0.25*sum(nodes(elements(nbr, 1:4), 2))
                        midzn = 0.25*sum(nodes(elements(nbr, 1:4), 3))
                        rx = (midxi - midxn)**2
                        ry = (midyi - midyn)**2
                        rz = (midzi - midzn)**2
                        r = rx + ry + rz

                        row = row + 1; col = par_map(i); ccount = ccount + 1
                        wrows(ccount) = row; wcols(ccount) = col; 
                        Wm(ccount) = zwts(zones(i), 2)*rx/r + zwts(zones(i), 3)*ry/r + zwts(zones(i), 4)*rz/r

                        col = par_map(nbr); ccount = ccount + 1
                        wrows(ccount) = row; wcols(ccount) = col; Wm(ccount) = -Wm(ccount - 1)!zwts(zones(i),1)

                        !row = row+1; col=par_map(i); ccount=ccount+1
                        !wrows(ccount)=row; wcols(ccount)=col;
                        !Wm(ccount) =1! zwts(zones(i),2)*rx/r + zwts(zones(i),3)*ry/r + zwts(zones(i),4)*rz/r

                        !col=par_map(nbr); ccount=count+1
                        !wrows(ccount)=row; wcols(ccount)=col;
                        !Wm(ccount)=-1!Wm(ccount-1)

                    end if

111                 continue
                end do

            end if
        end do

        return

1001    continue

        allocate (wrows(n_elements), wcols(n_elements), par_map(n_elements))
        allocate (Wm(n_elements), sigma_par(n_elements), rsigma_par(n_elements))
        ccount = n_elements
        npar = n_elements
        do i = 1, n_elements
            wrows(i) = i
            wcols(i) = i
            Wm(i) = 1
            par_map(i) = i
            if (invi) then
                sigma_par(i) = sigma_im(i)
            else
                sigma_par(i) = sigma_re(i)
            end if
        end do

   !!now map sigma_par back to sigma to make sure they're consistent
        if (invi) then
            do i = 1, n_elements
                sigma_im(i) = sigma_par(par_map(i))
            end do
        else
            do i = 1, n_elements
                sigma_re(i) = sigma_par(par_map(i))
            end do
        end if

    end subroutine build_Wm
    !_____________________________________________________________________

end module mod_con
