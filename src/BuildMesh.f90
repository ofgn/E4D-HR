! Changelog
! 27/4/23 - OFGN
! Add format specifier to force single line output with the Intel Fortran Compiler. - OFGN 27/4/23

module mesh

    use vars
    use input

    use iso_fortran_env, only: real128

contains

!__________________________________________________________________________________________________
    subroutine build_tetgen4(i_flag)

        !!This subroutine builds the computational mesh using triangle (for the surface) and
        !!tetgen (for the body)

        implicit none

        integer :: cfg_unit, io_status
        integer :: translate_flag
        real(real128) :: x_origin, y_origin, z_origin
        character(:), allocatable :: vis_loc

        logical :: i_flag

        integer :: m, d1, d2, d3, el1, el2, el3, el4, el5, i, j, k, tcount, nedge, ecount, bv, nsp2, ii, fstat
        integer :: n_points, count, npre, nplc, ns_points, n_surfac, dum1, dum2, dum3, ncplc, splc_count, nedge0, n2pts
        integer :: pt1, pt2, tetflag, vtk_flag, nh, nz, n, nedge_seq, nedge_flag, writit
        integer, dimension(:), allocatable ::  nc_plc, edge_flag, edge_seq, tflag, tf2, np_bounds
        integer, dimension(:), allocatable :: pmap, tmp_flag, tedge_flag
        integer, dimension(:, :), allocatable :: zon_labs, surfac
        integer, dimension(10000) :: temp_int
        real*8 :: min_elev, r, miny_ap, maxy_ap
        real*8, dimension(:, :), allocatable :: s_points, all_points, sp2, tall_points
        real*8, dimension(:), allocatable :: dist
        real*8, dimension(:, :), allocatable :: edges, all_points2, hole_xyz, zone_xyz
        real*8, dimension(:), allocatable :: edge_sep, min_vols, z_sig, z_isig
        character*20 :: qual, max_vol
        character*80 :: tetcom, exocom, tricom, exocom2d
        logical :: sdone = .false.
        integer, dimension(1) :: mind
        real*8 :: w1, w2, w3
        logical :: file_exists, mtran
        character(len=100) :: fmt

        open (52, file="mesh_build.log", action='write', status='replace')
        close(52)

        !!OPEN AND READ MESH FILE
        inquire (file=trim(cfg_filename), exist=file_exists)
        if (.not. file_exists) call check_minp(1, -1)

        open (newunit=cfg_unit, file=cfg_filename, status='old', action='read')
        read(cfg_unit, *, IOSTAT=io_status) qual, max_vol; call check_minp(2, io_status)
        read(cfg_unit, *, IOSTAT=io_status) min_elev; call check_minp(3, io_status)
        read(cfg_unit, *, IOSTAT=io_status) tetflag; call check_minp(4, io_status)
        read(cfg_unit, *, IOSTAT=io_status) tetcom; call check_minp(5, io_status)
        read(cfg_unit, *, IOSTAT=io_status) tricom; call check_minp(6, io_status)
        read(cfg_unit, *, IOSTAT=io_status) n_points; call check_minp(7, io_status)

        call check_minp(8, n_points)
        allocate (all_points(n_points, 3), tflag(n_points), edge_flag(n_points))
        allocate (tall_points(n_points, 3), tmp_flag(n_points), pmap(n_points))

        !!Read in the points. The surface points (including the edges) must
        !!be specified first in the file.
        all_points = 0
        tflag = 0
        nedge = 0
        ns_points = 0
        edge_flag = 0
        tall_points = 0
        tmp_flag = 0
        pmap = 0
        n2pts = 0
        open (52, file="mesh_build.log", action='write', status='old', position='append')
        do i = 1, n_points

            read(cfg_unit, *, IOSTAT=io_status) j, tall_points(i, 1:3), tmp_flag(i)
            if (io_status .ne. 0) then
                close(52)
                call check_minp(9, i)
            end if
            write (52, *) "Control point ", i, " reads: ", j, tall_points(i, 1:3), tmp_flag(i)

            !count the number of boundary points
            if (tmp_flag(i) == 2) n2pts = n2pts + 1;
        end do
        close(52)

        if (n2pts < 3) call check_minp(27, n2pts)

        !!reorder the points ... 1 first, then 2, then others
        do i = 1, n_points
            if (tmp_flag(i) == 1) then
                ns_points = ns_points + 1
                all_points(ns_points, 1:3) = tall_points(i, 1:3)
                tflag(ns_points) = 1
                pmap(i) = ns_points
            end if
        end do
        do i = 1, n_points
            if (tmp_flag(i) == 2) then
                nedge = nedge + 1
                ns_points = ns_points + 1
                all_points(ns_points, 1:3) = tall_points(i, 1:3)
                edge_flag(ns_points) = 1
                tflag(ns_points) = 2
                pmap(i) = ns_points
            end if
        end do
        m = ns_points
        do i = 1, n_points
            if (tmp_flag(i) .ne. 1 .and. tmp_flag(i) .ne. 2) then
                m = m + 1
                all_points(m, 1:3) = tall_points(i, 1:3)
                tflag(m) = tmp_flag(i)
                pmap(i) = m
            end if
        end do
        deallocate (tall_points, tmp_flag)

        !!build a sequential list of edge pnts defining the surface boundary
        call check_minp(11, io_status)
        allocate (edge_sep(nedge), edges(nedge, 2), edge_seq(nedge))
        ecount = 0
        do i = 1, ns_points
            if (edge_flag(i) == 1) then
                ecount = ecount + 1
                edges(ecount, 1:2) = all_points(i, 1:2)
            end if
        end do

        edge_seq(1) = 1

        do i = 1, nedge - 1
            edge_sep = sqrt((edges(:, 1) - edges(edge_seq(i), 1))**2 + &
                (edges(:, 2) - edges(edge_seq(i), 2))**2)
            edge_sep(edge_seq(1:i)) = maxval(edge_sep)
            edge_seq(i + 1) = minloc(edge_sep(:), 1)
        end do

        do i = 1, nedge
            do j = 1, ns_points
                if (edges(edge_seq(i), 1) == (all_points(j, 1)) .and. &
                    edges(edge_seq(i), 2) == (all_points(j, 2))) then
                    edge_seq(i) = j
                    exit
                end if
            end do
        end do

        !!If there are other PLCs defined besides the boundary,
        !!store them in a scratch file for now
        read(cfg_unit, *, IOSTAT=io_status) nplc; call check_minp(12, io_status)
        call check_minp(13, nplc)

        open (52, file="mesh_build.log", action='write', status='old', position='append')
        if (nplc > 0) then

            open (13, status='scratch', form='unformatted')
            allocate (nc_plc(nplc), np_bounds(nplc))

            do i = 1, nplc
                read(cfg_unit, *, IOSTAT=io_status) nc_plc(i), np_bounds(i)

                if (io_status .ne. 0) then
                    close(52)
                    call check_minp(14, i)
                else
                    write (52, *) "Number of control points and boundary number for: ", i, " are ", nc_plc(i), np_bounds(i)
                end if

                read(cfg_unit, *, IOSTAT=io_status) temp_int(1:nc_plc(i))
                if (io_status .ne. 0) then
                    close(52)
                    call check_minp(15, i)
                end if
                do j = 1, nc_plc(i)
                    if (temp_int(j) < 0 .or. temp_int(j) > n_points) then
                        close(52)
                        call check_minp(26, temp_int(j))
                    end if
                end do
                write (13) pmap(temp_int(1:nc_plc(i)))
            end do

        end if
        close(52)

        !!Read the hole and region options
        read(cfg_unit, *, IOSTAT=io_status) nh; call check_minp(16, io_status)
        call check_minp(17, nh)

        if (nh > 0) then
            allocate (hole_xyz(nh, 3))
            do i = 1, nh
                read(cfg_unit, *, IOSTAT=io_status) j, hole_xyz(i, 1:3)
                if (io_status .ne. 0) call check_minp(18, i)
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) "Coordinates for hole: ", i, " are ", hole_xyz(i, 1:3)
                close(52)
            end do

        end if

        read(cfg_unit, *, IOSTAT=io_status) nz; call check_minp(19, io_status)
        call check_minp(20, nz)

        if (nz > 0) then
            allocate (zone_xyz(nz, 3), zon_labs(nz, 2), min_vols(nz), z_sig(nz))
            if (i_flag) allocate (z_isig(nz))
            do i = 1, nz
                if (i_flag) then
                    read(cfg_unit, *, IOSTAT=io_status) zon_labs(i, 1), zone_xyz(i, 1:3), min_vols(i), z_sig(i), z_isig(i)
                else
                    read(cfg_unit, *, IOSTAT=io_status) zon_labs(i, 1), zone_xyz(i, 1:3), min_vols(i), z_sig(i)
                end if
                if (io_status .ne. 0) call check_minp(21, i)
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *)
                write (52, *) "Config info for zone: ", i, " ....... "
                write (52, *) 'Zone Number: ', zon_labs(i, 1)
                write (52, *) 'Point in zone: ', zone_xyz(i, 1:3)
                write (52, *) 'Maximum volume: ', min_vols(i)
                write (52, *) 'Conductivity: ', z_sig(i)
                if (i_flag) write (52, *) "Complex Conductivity: ", z_isig(i)
                close(52)
                zon_labs(i, 2) = zon_labs(i, 1)
            end do
        end if

        !! Read the VTK output flag
        read(cfg_unit, *, IOSTAT=io_status) vtk_flag
        call check_minp(22, io_status)

        read(cfg_unit, *, IOSTAT=io_status) vis_loc
        read(cfg_unit, *, IOSTAT=io_status) translate_flag
        if (translate_flag == 0) then
            mtran = .false.
        else
            mtran = .true.
        end if

        close(cfg_unit)

        if (mode == 0) then
            return
        end if

        !!BUILD THE MATRIX OF All POINTS USED TO BUILD THE MESH
        !!Open a scratch file to put the point positions into
        open (11, status='scratch', form='unformatted')

        !!Record the surface points
        do i = 1, n_points
            write (11) (all_points(i, 1:3))
        end do

        !!Record the points defining the bottom boundary
        !!there will be nedge of these points
        do i = 1, nedge
            write (11) (all_points(edge_seq(i), 1:2)), min_elev
        end do
        count = n_points + nedge
        nedge0 = nedge

        !!ALLOCATE AND INITIALIZE THE FULL POINTS ARRAY
        if (nplc > 0) then
            tcount = count
            rewind (13)
            rewind (11)
        else
            tcount = count
            rewind (11)
        end if

        ii = 0
        x_origin = 0
        y_origin = 0
        z_origin = 0

        if (mtran) then
            do i = 1, n_points
                if (tflag(i) .ne. 2) then
                    x_origin = x_origin + all_points(i, 1)
                    y_origin = y_origin + all_points(i, 2)
                    z_origin = z_origin + all_points(i, 3)
                    ii = ii + 1
                end if
            end do
            if (ii > 0) then
                x_origin = x_origin/real(ii)
                y_origin = y_origin/real(ii)
                z_origin = z_origin/real(ii)
            else
                x_origin = sum(all_points(:, 1))/real(n_points)
                y_origin = sum(all_points(:, 2))/real(n_points)
                z_origin = sum(all_points(:, 3))/real(n_points)
            end if
        end if

        deallocate (all_points)
        allocate (all_points2(tcount, 3))

        do i = 1, tcount
            read (11) all_points2(i, 1:3)
        end do
        close(11)

        !translate the coordinates to optimize precision
        all_points2(:, 1) = all_points2(:, 1) - x_origin
        all_points2(:, 2) = all_points2(:, 2) - y_origin
        all_points2(:, 3) = all_points2(:, 3) - z_origin

        min_elev = min_elev - z_origin

        if (nh > 0) then
            hole_xyz(:, 1) = hole_xyz(:, 1) - x_origin
            hole_xyz(:, 2) = hole_xyz(:, 2) - y_origin
            hole_xyz(:, 3) = hole_xyz(:, 3) - z_origin
        end if

        if (nz > 0) then
            zone_xyz(:, 1) = zone_xyz(:, 1) - x_origin
            zone_xyz(:, 2) = zone_xyz(:, 2) - y_origin
            zone_xyz(:, 3) = zone_xyz(:, 3) - z_origin
        end if

        !!record the translation coordinates and write
        !!to file for runs in mode > 1

        open (22, file=trim(mesh_prefix)//".trn", status='replace')
        write (22, "(3E20.12)") x_origin, y_origin, z_origin
        close(22)

        !!CONSTRUCT THE TRIANGLE .poly INPUT FILE
        !if (msh_opt.eq.1) then
        open (33, file="surface.poly", status='replace', action='write')
        write (33, *) "#GENERATED BY FERM3D BUILD_TETGEN4 DURING CALL TO TRIANGLE"
        WRITE (33, *) "#THERE ARE ", ns_points, " MESH CONFIGURATION POINTS"
        maxy_ap = maxval(all_points2(1:ns_points, 2))
        miny_ap = minval(all_points2(1:ns_points, 2))
        write (33, 111) ns_points, "        2        0        1  # Total nodes, dimensions, attributes, boundary marker flag"
        do i = 1, ns_points
            write (33, 119) i, all_points2(i, 1:2), tflag(i)
        end do

        !!THIS SECTIONS COUNTS THE NUMBER OF PLCS ON THE SURFACE
        splc_count = 0
        if (nplc > 0) then
            do i = 1, nplc
                read (13) temp_int(1:nc_plc(i))
                do m = 1, nc_plc(i) - 1
                    if (tflag(temp_int(m)) > 0 .and. tflag(temp_int(m + 1)) > 0 .and. &
                        tflag(temp_int(m)) < 3 .and. tflag(temp_int(m + 1)) < 3) then
                        splc_count = splc_count + 1
                    end if
                end do

                if (tflag(temp_int(1)) > 0 .and. tflag(temp_int(nc_plc(i))) > 0 .and. &
                    tflag(temp_int(1)) < 3 .and. tflag(temp_int(nc_plc(i))) < 3) then
                    splc_count = splc_count + 1
                end if

            end do
            rewind (13)
        end if

        !!now do the boundary segments...
        write (33, *) nedge + splc_count, "        1   #THERE ARE", nedge + splc_count, "BOUNDARY SEGMENTS"
        do i = 1, nedge - 1
            write (33, *) i, edge_seq(i), edge_seq(i + 1), 2
        end do
        write (33, *) nedge, edge_seq(nedge), edge_seq(1), 2

        !!NOW FIND THE SURFACE PLCS AND WRITE OUT TO THE .POLY FILE
        if (nplc > 0) then
            write (33, *) "#THERE ARE", splc_count, "USER DEFINED SURFACE PLCS"
            splc_count = 0
            do i = 1, nplc

                read (13) temp_int(1:nc_plc(i))
                do m = 1, nc_plc(i) - 1
                    if (tflag(temp_int(m)) > 0 .and. tflag(temp_int(m + 1)) > 0 .and. &
                        tflag(temp_int(m)) < 3 .and. tflag(temp_int(m + 1)) < 3) then
                        splc_count = splc_count + 1
                        write (33, *) nedge + splc_count, temp_int(m), temp_int(m + 1), 3 + i
                    end if
                end do
                if (tflag(temp_int(1)) > 0 .and. tflag(temp_int(nc_plc(i))) > 0 .and. &
                    tflag(temp_int(1)) < 3 .and. tflag(temp_int(nc_plc(i))) < 3) then
                    splc_count = splc_count + 1
                    write (33, *) nedge + splc_count, temp_int(nc_plc(i)), temp_int(1), 6 + i
                end if
            end do
            rewind (13)
        end if
        write (33, *) "#THESE ARE THE HOLES ZONES"
        write (33, *) "        0"
        write (33, *) "#THESE ARE THE DEFINED ZONES FOR TRIANGLE"
        write (33, *) "        0"

        close(33)

        !!Delete previous surface mesh files if they exist
        inquire (file='surface.1.node', exist=file_exists)
        if (file_exists) call system('rm surface.1.node');
        inquire (file='surface.1.ele', exist=file_exists)
        if (file_exists) call system('rm surface.1.ele');
        inquire (file='surface.1.neigh', exist=file_exists)
        if (file_exists) call system('rm surface.1.neigh');
        inquire (file='surface.1.edge', exist=file_exists)
        if (file_exists) call system('rm surface.1.edge');
        write (*, *) "Calling Triangle"
        call system(trim(tricom)//" -pnq30e surface.poly")
        call check_minp(24, 0)
        call build_sig_triangle

        !!THIS SAVES A COPY OF THE .ELE FILE, CHANGES THE .ELE FILE SO THAT IT HAS ATTRIBUTES
        !!WHICH ARE NECESSARY TO BUILD THE TRIANGLE 2D.EXO FILE
        call system("cp surface.1.ele surface.1.ele.old")
        open (77, file="surface.1.ele.old", status='old', action='read')
        open (88, file="surface.1.ele", status='replace', action='write')
        read (77, *) d1, d2, d3
        write (88, *) d1, 3, 1
        do i = 1, d1
            read (77, *) el1, el2, el3, el4
            write (88, *) el1, el2, el3, el4, 1
        end do
        close(88)
        close(77)

119     format(I10, 2ES22.12, I10)
120     format(4I10)

        write (*, *) "DONE WITH SURFACE TRIANGULATION"

!____________________________________________________________________________
        open (23, file='surface.1.node', status='old', action='read')
        read (23, *) nsp2, d1, d2, d3
        allocate (sp2(nsp2, 3), tf2(nsp2))
        allocate (dist(ns_points))
        sp2 = 0
        tf2 = 0
        do i = 1, nsp2
            read (23, *) d1, sp2(i, 1:2), tf2(i)
        end do
        close(23)

        do i = 1, nsp2

            dist = sqrt((dble(all_points2(1:ns_points, 1)) - sp2(i, 1))**2 + (dble(all_points2(1:ns_points, 2)) - sp2(i, 2))**2)
            mind = minloc(dist)
            d1 = mind(1)
            w1 = dist(d1)
            dist(d1) = 1e30
            mind = minloc(dist)
            d2 = mind(1)
            w2 = dist(d2)
            dist(d2) = 1e30
            mind = minloc(dist)
            w3 = dist(d3)
            d3 = mind(1)

            if (w1 == 0.0) then
                sp2(i, 3) = dble(all_points2(d1, 3))
                goto 10
            elseif (w2 == 0.0) then
                sp2(i, 3) = dble(all_points2(d2, 3))
                goto 10
            elseif (w3 == 0.0) then
                sp2(i, 3) = dble(all_points2(d3, 3))
                goto 10
            end if

            sp2(i, 3) = 1/(1/w1 + 1/w2 + 1/w3)*(dble(all_points2(d1, 3))/w1 + dble(all_points2(d2, 3))/w2 + dble(all_points2(d3, 3))/w3)
10          continue
        end do

        !__________________________________________________________________________________
        nedge = 0
        do i = 1, nsp2
            if (tf2(i) == 2) then
                nedge = nedge + 1
            elseif (tf2(i) == 0) then
                tf2(i) = 1
            end if
        end do
        if (allocated(edge_sep)) deallocate (edge_sep)
        if (allocated(edges)) deallocate (edges)
        if (allocated(edge_seq)) deallocate (edge_seq)
        allocate (edge_sep(nedge), edges(nedge, 2), edge_seq(nedge))
        ecount = 0
        do i = 1, nsp2
            if (tf2(i) == 2) then
                ecount = ecount + 1
                edges(ecount, 1:2) = sp2(i, 1:2)
            end if
        end do

        edge_seq(1) = 1

        do i = 1, nedge - 1
            edge_sep = sqrt((edges(:, 1) - edges(edge_seq(i), 1))**2 + &
                (edges(:, 2) - edges(edge_seq(i), 2))**2)
            edge_sep(edge_seq(1:i)) = maxval(edge_sep)
            edge_seq(i + 1) = minloc(edge_sep(:), 1)
        end do

        do i = 1, nedge
            do j = 1, nsp2
                if (edges(edge_seq(i), 1) == (sp2(j, 1)) .and. &
                    edges(edge_seq(i), 2) == (sp2(j, 2))) then
                    edge_seq(i) = j
                    exit
                end if
            end do
        end do

!!!!!!!!!!!!!!!!!!NOW WE ARE READY TO CONSTRUCT THE TETGEN INPUT FILE!!!!!!!!!!!!!!!!!!!!!!

        write (*, *) "BEGIN CONSTRUCTING TETGEN INPUT .poly FILE"
        !!Change the extension of mesh_file to .poly
        k = len_trim(cfg_filename)
        do i = 1, k
            if (cfg_filename(i:i) == '.') then
                cfg_filename(i + 1:i + 4) = 'poly'
                npre = i
                exit
            end if
        end do

        tcount = 0
        do i = 1, n_points
            if (tflag(i) .ne. 1 .and. tflag(i) .ne. 2) then
                tcount = tcount + 1
            end if
        end do

        open (11, file=cfg_filename, status='replace', action='write')
        write (11, 111) tcount + nedge + nsp2, ", 3, 1, 1,  # Total nodes, dimensions,         attributes, boundary marker flag"

        !!Surface points
        write (11, *) "# The next ", nsp2, " points contains the surface points created by triangle"
        write (11, *) "# including electrodes"
        write (11, *) "# Boundary marker 1 indicates the surface (Positive z         normal) boundary"

        do i = 1, nsp2
            write (11, 113) i, sp2(i, :), 1, tf2(i)
        end do

        !!Lower Boundary points
        write (11, 117) nedge
        write (11, *) "# Boundary 2 indicates the lower (negative z normal)         boundary"
        do i = 1, nedge
            write (11, 113) i + nsp2, sp2(edge_seq(i), 1:2), min_elev, 1, 2
        end do

        !!Internal refine points
        write (11, 114) "# The next ", tcount, " points are internal refine points"
        write (11, *) "# A boundary flag of zero indicates the point is not on         the boundary"
        do i = 1, tcount
            write (11, 113) i + nsp2 + nedge, all_points2(i + ns_points, 1:3), 1, tflag(i + ns_points)
        end do

        !!Write the 'waterproof' facets that define the mesh boundaries
        !!Define the surface boundary,this read in the .ele file created by triangle
        open (44, file="surface.1.ele", status='old', action='read')
        read (44, *) n_surfac, dum1, dum2
        write (11, *)
        write (11, *)
        write (11, *) n_surfac + nplc + nedge + 1, 1, "  #Total number of facets and         boundary marker flag"
        !write(11,*) n_surfac+nedge+1,1,"  #Total number of facets and         boundary marker flag"
        write (11, *) "# 1 = surface boundary"
        write (11, *) "# 2 = side and bottom boundary"

121     format(I10)
122     format(4I10)
123     format(3I10)
124     format(5I10)
        !!WRITE OUT THE FACETS READ FROM THE TRIANGLE.ELE FILE INTO THE TETGEN.POLY FILE
        !!THESE ARE THE FACETS THAT DEFINE THE SURFACE
        allocate (surfac(n_surfac, 3))
        write (11, *) "# The next", n_surfac, " facets define the surface,         boundary 1"
        do i = 1, n_surfac
            write (11, *) "1 0 1"
            read (44, *) dum3, surfac(i, 1:3)
            write (11, 122) 3, surfac(i, 1:3)
        end do

        !!Define the side (vertical) boundaries
        write (11, *)
        write (11, 118) nedge
        do i = 1, nedge - 1
            write (11, *) "1 0 2"
            write (11, "(5I10)") 4, edge_seq(i), nsp2 + i, nsp2 + i + 1, edge_seq(i + 1)
        end do

        write (11, *) "1 0 2"
        write (11, "(5I10)") 4, edge_seq(nedge), nsp2 + nedge, nsp2 + 1, edge_seq(1)

        !!Define the lower boundary
        write (11, *)
        write (11, *) " # This facet defines the lower boundary"
        write (11, *) "1 0 2"
        write (11, "(I6)", advance='no') nedge
        do i = 1, nedge
            write (11, "(I7)", advance='no') i + nsp2
        end do
        write (11, *)
        write (11, *)

        !!Define the non boundary PLC's if necessary
        if (nplc > 0) then
            rewind (13)

            write (11, *) '# ', nplc, "USER DEFINED PLC'S"
            do i = 1, nplc
                d1 = 3 + i
                read (13) temp_int(1:nc_plc(i))
                do j = 1, nc_plc(i)
                    if (tflag(temp_int(j)) .ne. 1 .and. tflag(temp_int(j)) .ne. 2) then
                        temp_int(j) = nsp2 + nedge + (temp_int(j) - ns_points)
                    end if
                end do

                do j = 1, nsp2
                    if (tf2(j) == d1) then
                        nc_plc(i) = nc_plc(i) + 1
                        temp_int(nc_plc(i)) = j
                    end if
                end do

                deallocate (edge_sep, edges, edge_seq)
                allocate (edge_sep(nc_plc(i)), edges(nc_plc(i), 3), edge_seq(nc_plc(i)))

                do j = 1, nc_plc(i)
                    if (temp_int(j) > nsp2) then
                        edges(j, :) = all_points2(ns_points + (temp_int(j) - nsp2 - nedge), :)
                    else
                        edges(j, :) = (sp2(temp_int(j), :))
                    end if
                end do

                edge_seq(1) = 1
                do j = 1, nc_plc(i) - 1
                    edge_sep = sqrt((edges(:, 1) - edges(edge_seq(j), 1))**2 + &
                        (edges(:, 2) - edges(edge_seq(j), 2))**2 + &
                        (edges(:, 3) - edges(edge_seq(j), 3))**2)
                    edge_sep(edge_seq(1:j)) = 1e30 !maxval(edge_sep)
                    edge_seq(j + 1) = minloc(edge_sep(:), 1)
                end do

                write (11, "(3I10)") 1, 0, np_bounds(i)
                write (fmt, '(a, I0, a)') '(I8,', nc_plc(i), 'I8)'
                write (11, trim(fmt)) nc_plc(i), temp_int(edge_seq)

            end do
            close(13)
        end if

11      continue

        write (11, *)
        write (11, *) nh, "   #Number of holes"
        if (nh > 0) then
            fmt = "(I0, X, G0.15, X, G0.15, X, G0.15)" ! Add format specifier to force single line output with the Intel Fortran Compiler. - OFGN 27/4/23
            do i = 1, nh
                write (11, trim(fmt)) i, hole_xyz(i, :)
            end do
        end if

        write (11, *) nz, "   #Number of zones"
        if (nz > 0) then
            fmt = "(I0, X, G0.15, X, G0.15, X, G0.15, X, I0, X, G0.15)" ! Add format specifier to force single line output with the Intel Fortran Compiler. - OFGN 27/4/23
            do i = 1, nz
                write (11, trim(fmt)) zon_labs(i, 1), zone_xyz(i, :), zon_labs(i, 2), min_vols(i)
            end do
        end if

        close(11)
        deallocate (all_points2)
        write (*, *) "DONE BUILDING ", trim(cfg_filename)

1001    continue

        !!Delete previous mesh files if they exist
        inquire (file=cfg_filename(1:npre)//'1.node', exist=file_exists)
        if (file_exists) call system('rm '//cfg_filename(1:npre)//'1.node')
        inquire (file=cfg_filename(1:npre)//'1.ele', exist=file_exists)
        if (file_exists) call system('rm '//cfg_filename(1:npre)//'1.ele')
        inquire (file=cfg_filename(1:npre)//'1.face', exist=file_exists)
        if (file_exists) call system('rm '//cfg_filename(1:npre)//'1.face')
        inquire (file=cfg_filename(1:npre)//'1.neigh', exist=file_exists)
        if (file_exists) call system('rm '//cfg_filename(1:npre)//'1.neigh')

        ! Call TetGen
        ! (ofgn - 11/12/24): Added additional TetGen command line switches
        write (*, *) "    CALLING TETGEN         "
        if (nz > 0) then
            write (*, *) trim(tetcom) // " -pnq" // trim(qual) // "a" // &
                trim(max_vol) // "aAAnnefV " // trim(cfg_filename)
            call system(trim(tetcom) // " -pnq" // trim(qual) // "a" // &
                trim(max_vol) // "aAAnnefV " // trim(cfg_filename))
        else
            write (*, *) trim(tetcom) // " -pnq" // trim(qual) // "a" // &
                trim(max_vol) // "AAnnefV " // trim(cfg_filename)
            call system(trim(tetcom) // " -pnq" // trim(qual) // "a" // &
                trim(max_vol) // "AAnnefV " // trim(cfg_filename))
        end if
        !call system('rm surface.*')

        do i = 1, 40
            if (cfg_filename(i:i) == '.') then
                npre = i
                exit
            end if
        end do

        !!Make sure the mesh files exist to see if tetgen ran succesfully
        inquire (file=cfg_filename(1:npre)//'1.node', exist=file_exists)
        if (.not. file_exists) call check_minp(25, 0)
        inquire (file=cfg_filename(1:npre)//'1.ele', exist=file_exists)
        if (.not. file_exists) call check_minp(25, 0)
        inquire (file=cfg_filename(1:npre)//'1.face', exist=file_exists)
        if (.not. file_exists) call check_minp(25, 0)
        inquire (file=cfg_filename(1:npre)//'1.neigh', exist=file_exists)
        if (.not. file_exists) call check_minp(25, 0)

        !!make sure the boundary flags are correct and rebuild if not
        !if(.not. tank_flag) then
        write (*, *) "Checking node boundary flags"
        call check_nodes(cfg_filename)
        !end if

        !! Build the conductivity file if necessary
        if (nz > 0) then
            call build_sig(cfg_filename, npre, nz, z_sig, z_isig, zon_labs, i_flag)
        end if

        !! Export the mesh to VTK format if specified
        if (vtk_flag .eq. 1) then
            inquire (file=trim(mesh_prefix)//".1.node", exist=file_exists)
            if (file_exists) then
                inquire (file=trim(mesh_prefix)//".1.ele", exist=file_exists)
                if (file_exists) then
                    inquire (file=trim(mesh_prefix)//".sigma", exist=file_exists)
                    if (file_exists) then
                        ! call export_mesh_model(trim(mesh_prefix), nz)
                    else
                        call mesh_error_log(106, "Cannot find "//trim(mesh_prefix)//".sigma")
                    end if
                else
                    call mesh_error_log(102, "Cannot find "//trim(mesh_prefix)//".1.ele")
                end if
            else
                call mesh_error_log(101, "Cannot find "//trim(mesh_prefix)//".1.node")
            end if
        end if

111     format(BN, (I10), A71)
112     format(A12, BN, (I10), A63)
113     format(I10, 3ES22.12, 2I10)
114     format(A11, I10, A61)
115     format("# The next ", I10, " points contain nodes for PCL ", I10)
116     format("# The next ", I10, "facets define the boundary of PLC ", I10)
117     format(" # The next ", I10, " points define the lower boundary")
118     format(" # The next ", I10, " facets define the vertical boundaries,         boundary 3")

    end subroutine build_tetgen4
!___________________________________________________________________________________________________

!___________________________________________________________________________________________________
    subroutine check_nodes(mpre)
        !!checks to make sure nodes have correct boundary markers
        character*40 :: mpre
        logical :: exst
        integer :: i, j, k, l, ecount, npre                          !counters
        integer :: ncpts                                        !number of control points
        integer :: nplc
        integer :: nnods                                        !number of nodes
        integer :: nedge                                        !number of edge points

        integer, dimension(:), allocatable :: bpts, cbpts, fbpts  !boundary flags
        integer, dimension(:), allocatable :: edge_seq          !boundary flags

        integer, dimension(:, :), allocatable :: plc             !plc definitions
        real*8 :: x, y, z, pd, a1, a2, a3, dp, dist, lb                  !points
        real*8 :: xp, yp, zp, t, mx, my, mz, tl
        real*8, dimension(3) :: p1, p2, p3, N, pt, v1, v2, v3          !points and normal
        real*8, dimension(:, :), allocatable :: cpts             !control points
        real*8, dimension(:, :), allocatable :: nods             !node points
        real*8, dimension(:, :), allocatable :: edges            !node points
        real*8, dimension(:), allocatable :: edge_sep            !node points
        real*8, dimension(3, 3) :: M                             !matrix
        real*8, parameter :: dtl = 1e-3

        k = len_trim(mpre)
        do i = 1, k
            if (mpre(i:i) == '.') then
                mpre(i + 1:i + 4) = 'poly'
                npre = i
                exit
            end if
        end do

        !!read the mesh translation
        open (20, file=mpre(1:npre)//'trn', status='old', action='read')
        read (20, *) pt
        close(20)

        !!read the config file
        open (20, file=mpre(1:npre)//'cfg', status='old', action='read')
        read (20, *)
        read (20, *) lb
        read (20, *)
        read (20, *)
        read (20, *)

        !read and translate the control points
        nedge = 0
        read (20, *) ncpts
        allocate (cpts(ncpts, 3), cbpts(ncpts))
        do i = 1, ncpts
            read (20, *) j, cpts(i, 1:3), cbpts(i)
            if (cbpts(i) == 2) nedge = nedge + 1
            do j = 1, 3
                cpts(i, j) = cpts(i, j) - pt(j);
            end do
        end do
        lb = lb - pt(3)

        !read the plcs
        read (20, *) nplc
        allocate (plc(nplc, 100))
        do i = 1, nplc
            read (20, *) plc(i, 1:2)
            read (20, *) plc(i, 3:2 + plc(i, 1))
        end do
        close(20)

        !read the nodes
        open (20, file=mpre(1:npre)//'1.node', status='old', action='read')
        read (20, *) nnods
        allocate (nods(nnods, 3), bpts(nnods), fbpts(nnods))
        do i = 1, nnods
            read (20, *) j, nods(i, 1:3), k, bpts(i)
        end do
        close(20)
        fbpts = bpts

        !compute the normal for each plc
        do i = 1, nplc

            do k = 1, nnods
                x = nods(k, 1)
                y = nods(k, 2)
                z = nods(k, 3)

                !first see if this point is near to the segement connecting two control points
                !defining this plc
                do l = 1, plc(i, 1) - 1
                    xp = x - cpts(plc(i, 2 + l), 1)
                    yp = y - cpts(plc(i, 2 + l), 2)
                    zp = z - cpts(plc(i, 2 + l), 3)

                    mx = cpts(plc(i, 3 + l), 1) - cpts(plc(i, 2 + l), 1)
                    my = cpts(plc(i, 3 + l), 2) - cpts(plc(i, 2 + l), 2)
                    mz = cpts(plc(i, 3 + l), 3) - cpts(plc(i, 2 + l), 3)

                    dist = sqrt(mx**2 + my**2 + mz**2)
                    tl = dtl*dist

                    t = (xp*mx + yp*my + zp*mz)/(mx**2 + my**2 + mz**2)
                    dist = sqrt((xp - t*mx)**2 + (yp - t*my)**2 + (zp - t*mz)**2)
                    if (dist < tl .and. (t) .ge. 0 .and. (t) .le. 1) then
                        fbpts(k) = plc(i, 2)

                        goto 100
                    end if
                end do

                !do the last segment here
                xp = x - cpts(plc(i, 2 + plc(i, 1)), 1)
                yp = y - cpts(plc(i, 2 + plc(i, 1)), 2)
                zp = z - cpts(plc(i, 2 + plc(i, 1)), 3)

                mx = cpts(plc(i, 3), 1) - cpts(plc(i, 2 + plc(i, 1)), 1)
                my = cpts(plc(i, 3), 2) - cpts(plc(i, 2 + plc(i, 1)), 2)
                mz = cpts(plc(i, 3), 3) - cpts(plc(i, 2 + plc(i, 1)), 3)

                dist = sqrt(mx**2 + my**2 + mz**2)
                tl = dtl*dist

                t = (xp*mx + yp*my + zp*mz)/(mx**2 + my**2 + mz**2)
                dist = sqrt((xp - t*mx)**2 + (yp - t*my)**2 + (zp - t*mz)**2)
                if (dist < tl .and. (t) .ge. 0 .and. (t) .le. 1) then
                    fbpts(k) = plc(i, 2)

                    goto 100
                end if

                if (plc(i, 1) > 2) then

                    !build the normal to this plc
                    do j = 1, 4
                        M(1, :) = cpts(plc(i, 3), :);
                        M(2, :) = cpts(plc(i, 4), :);
                        M(3, :) = cpts(plc(i, 5), :);
                        if (j < 4) then
                            M(:, j) = 1
                            N(j) = getdet(M, 3)
                        else
                            D = -getdet(M, 3)
                        end if
                    end do

                    !set the tolerance
                    dist = sqrt((M(1, 1) - M(2, 1))**2 + (M(1, 2) - M(2, 2))**2 + (M(1, 3) - M(2, 3))**2)
                    tl = dtl*dist
                    dist = sqrt((M(3, 1) - M(2, 1))**2 + (M(3, 2) - M(2, 2))**2 + (M(3, 3) - M(2, 3))**2)
                    if (dist < (tl/dtl)) tl = dtl*dist

                    !if we're here then this point is in the plane and
                    !doesn't have boundary marker specified in the config
                    !file. We need to determine if the point is bounded by
                    !by the points defining the plc and if so, change the
                    !boundary marker

                    !The center piont is M(2,:) defined above
                    !a1 is the angle between the vectors (M(1,:)-M(2,:)) and (M(3,:)-M(2,:))
                    !a2 is the angle between the vectors (nods(k,:)-M(2,:)) and (M(1,:)-M(2,:))
                    !a3 is the angle between the veoctors(nods(k,:)-M(2,:)) and (M(3,:)-M(2,:))
                    !(a2+a3)>=a1 for every set of points on the boundary if point nods(k,:) is
                    !within the boundary

                    !pd is the distance from this point to the plane
                    pd = abs(N(1)*x + N(2)*y + N(3)*z + D)/sqrt(dot_product(N, N))

                    if (pd < tl .and. plc(i, 2) .ne. bpts(k)) then
                        v1 = cpts(plc(i, plc(i, 1) + 2), :) - cpts(plc(i, 3), :);
                        v2 = cpts(plc(i, 4), :) - cpts(plc(i, 3), :);
                        v3 = nods(k, :) - cpts(plc(i, 3), :);
                        a1 = acos(dot_product(v2, v1)/(sqrt(dot_product(v2, v2))*sqrt(dot_product(v1, v1))))
                        a2 = acos(dot_product(v3, v1)/(sqrt(dot_product(v3, v3))*sqrt(dot_product(v1, v1))))
                        a3 = acos(dot_product(v3, v2)/(sqrt(dot_product(v3, v3))*sqrt(dot_product(v2, v2))))

                        if ((a2 + a3) .gt. (a1)) then
                            goto 100
                        end if

                        do l = 2, plc(i, 1) - 1

                            v1 = cpts(plc(i, l + 1), :) - cpts(plc(i, l + 2), :)
                            v2 = cpts(plc(i, l + 3), :) - cpts(plc(i, l + 2), :)
                            v3 = nods(k, :) - cpts(plc(i, l + 2), :)

                            a1 = acos(dot_product(v2, v1)/(sqrt(dot_product(v2, v2))*sqrt(dot_product(v1, v1))))
                            a2 = acos(dot_product(v3, v1)/(sqrt(dot_product(v3, v3))*sqrt(dot_product(v1, v1))))
                            a3 = acos(dot_product(v3, v2)/(sqrt(dot_product(v3, v3))*sqrt(dot_product(v2, v2))))

                            if ((a2 + a3) .gt. (a1)) then
                                goto 100
                            end if
                        end do

                        !if we're here then the point is on the plane and within the plane boundary
                        !so fix the flag
                        fbpts(k) = plc(i, 2)

                    end if
                end if
100             continue
            end do

        end do

        !now we need to check the outer boundaries to make the outer boundary points
        !have the correct boundary markers.
        !!build a sequential list of edge pnts defining the surface boundary
        allocate (edge_sep(nedge), edges(nedge, 3), edge_seq(nedge))
        ecount = 0
        do i = 1, ncpts
            if (cbpts(i) == 2) then
                ecount = ecount + 1
                edges(ecount, 1:3) = cpts(i, 1:3)
            end if
        end do
        edge_seq(1) = 1
        do i = 1, nedge - 1
            edge_sep = sqrt((edges(:, 1) - edges(edge_seq(i), 1))**2 + &
                (edges(:, 2) - edges(edge_seq(i), 2))**2)
            edge_sep(edge_seq(1:i)) = maxval(edge_sep)
            edge_seq(i + 1) = minloc(edge_sep(:), 1)
        end do

        !compute the normal for each plc
        do i = 1, nedge

            if (i < nedge) then
                do j = 1, 4
                    M = 0
                    M(1, 1:3) = edges(edge_seq(i), 1:3)
                    M(2, 1:3) = edges(edge_seq(i + 1), 1:3)
                    M(3, 1:3) = M(2, 1:3)
                    M(3, 3) = lb
                    if (j < 4) then
                        M(:, j) = 1
                        N(j) = getdet(M, 3)
                    else
                        D = -getdet(M, 3)
                    end if
                end do
            else
                do j = 1, 4
                    M = 0
                    M(1, 1:3) = edges(edge_seq(nedge), 1:3)
                    M(2, 1:3) = edges(edge_seq(1), 1:3)
                    M(3, 1:3) = M(2, 1:3)
                    M(3, 3) = lb
                    if (j < 4) then
                        M(:, j) = 1
                        N(j) = getdet(M, 3)
                    else
                        D = -getdet(M, 3)
                    end if
                end do

            end if

            !set the tolerance
            dist = sqrt((M(1, 1) - M(2, 1))**2 + (M(1, 2) - M(2, 2))**2 + (M(1, 3) - M(2, 3))**2)
            tl = dtl*dist
            dist = sqrt((M(3, 1) - M(2, 1))**2 + (M(3, 2) - M(2, 2))**2 + (M(3, 3) - M(2, 3))**2)
            if (dist < tl/dtl) tl = dtl*dist

            do k = 1, nnods
                x = nods(k, 1)
                y = nods(k, 2)
                z = nods(k, 3)

                if ((z) == (lb)) then
                    fbpts(k) = 2
                else
                    !pd is the distance from this point to the plane
                    pd = abs(N(1)*x + N(2)*y + N(3)*z + D)/sqrt(dot_product(N, N))

                    if (pd < tl) then
                        fbpts(k) = 2
                    end if
                end if
            end do

        end do

        !rewrite the node file with the corrected boundary markers
        call system("cp "//mpre(1:npre)//"1.node "//trim(mpre)//".1_orig.node")
        open (20, file=mpre(1:npre)//'1.node', status='replace', action='write')
        write (20, *) nnods, 3, 1, 1
        do i = 1, nnods
            write (20, "(I10,G20.10,G20.10,G20.10,I3,I10)") i, nods(i, 1), nods(i, 2), nods(i, 3), 1, fbpts(i)
        end do
        write (20, *) "THIS NODE FILE WAS MODIFIED BY E4D. SEE ORIGINAL TETGEN NODE FILE IN: ", trim(mpre)//".1_orig.node"
        close(20)
    end subroutine check_nodes
!___________________________________________________________________________________________________

!___________________________________________________________________________________________________
    subroutine check_minp(indx, ios)
        implicit none
        integer :: indx, ios
        logical :: exst

        select case (indx)

          case (1)
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *) "Cannot find the mesh configuration file: ", trim(cfg_filename)
            close(52)
            write (*, *) "Cannot find the mesh configuration file: ", trim(cfg_filename)
            call crash_exit

          case (2)
            if (ios .ne. 0) then
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) "There was a problem reading the quality factor max element volume in ", trim(cfg_filename)
                close(52)
                write (*, *) "There was a problem reading the quality factor max element volume in ", trim(cfg_filename)
                call crash_exit
            end if

          case (3)
            if (ios .ne. 0) then
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) "There was a problem reading the minimum elevation in ", trim(cfg_filename)
                close(52)
                write (*, *) "There was a problem reading the minimum elevation in ", trim(cfg_filename)
                call crash_exit
            end if

          case (4)
            if (ios .ne. 0) then
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) "There was a problem reading the build tetgen flag ", trim(cfg_filename)
                close(52)
                write (*, *) "There was a problem reading the build tetgen flag ", trim(cfg_filename)
                call crash_exit
            end if

          case (5)
            if (ios .ne. 0) then
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) "There was a problem reading the tetgen executable name in ", trim(cfg_filename)
                close(52)
                write (*, *) "There was a problem reading the tetgen executable name in ", trim(cfg_filename)
                call crash_exit
            end if

          case (6)
            if (ios .ne. 0) then
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) "There was a problem reading the triangle executable name in ", trim(cfg_filename)
                close(52)
                write (*, *) "There was a problem reading the triangle executable name in ", trim(cfg_filename)
                call crash_exit
            end if

          case (7)
            if (ios .ne. 0) then
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) "There was a problem reading the number of control points ", trim(cfg_filename)
                close(52)
                write (*, *) "There was a problem reading the number of control points ", trim(cfg_filename)
                call crash_exit
            end if

          case (8)
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *) "Allocating arrays for ", ios, ' control points'
            close(52)

          case (9)
            if (ios .ne. 0) then
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) "There was a problem reading control point number: ", ios
                close(52)
                write (*, *) "There was a problem reading control point number: ", ios
                call crash_exit
            end if

          case (10)
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *) 'SURFACE POINTS (flag 1 or 2) MUST BE DEFINED BEFORE INTERNAL POINTS (flag 0)'
            write (52, *) 'SEE POINT ', ios, ' IN MESH CONFIGURATION FILE ', trim(cfg_filename)
            close(52)
            write (*, *) 'SURFACE POINTS (flag 1 or 2) MUST BE DEFINED BEFORE INTERNAL POINTS (flag 0)'
            write (*, *) 'SEE POINT ', ios, ' IN MESH CONFIGURATION FILE ', trim(cfg_filename)
            call crash_exit

          case (11)
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *) "Building sequetial list of control points for outer boundary"
            close(52)
          case (12)
            if (ios .ne. 0) then
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) 'There was a problem reading the number of piecewise linear complexes in', trim(cfg_filename)
                close(52)
                write (*, *) 'There was a problem reading the number of piecewise linear complexes in', trim(cfg_filename)
                call crash_exit
            end if

          case (13)
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *)
            write (52, *) 'There are: ', ios, 'piecewise linear complexes'
            close(52)

          case (14)
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *) 'There was a problem reading the number of points and boundary number'
            write (52, *) 'for piecewise linear complex: ', ios, ' in ', trim(cfg_filename)
            close(52)
            write (*, *) 'There was a problem reading the number of points and boundary number'
            write (*, *) 'for piecewise linear complex: ', ios, ' in ', trim(cfg_filename)
            call crash_exit

          case (15)
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *) 'There was a problem reading the control points for plc number: ', ios
            close(52)
            write (*, *) 'There was a problem reading the control points for plc number: ', ios
            call crash_exit

          case (16)
            if (ios .ne. 0) then
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) 'There was a problem reading the number of holes', trim(cfg_filename)
                close(52)
                write (*, *) 'There was a problem reading the number of holes', trim(cfg_filename)
                call crash_exit
            end if

          case (17)
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *)
            write (52, *) 'There are: ', ios, 'holes'
            close(52)

          case (18)
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *) 'There was a problem reading the hole coordinates for hole', ios
            close(52)
            write (*, *) 'There was a problem reading the hole coordinates for hole', ios
            call crash_exit

          case (19)
            if (ios .ne. 0) then
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) 'There was a problem reading the number of zones in ', trim(cfg_filename)
                close(52)
                write (*, *) 'There was a problem reading the number of zones in ', trim(cfg_filename)
                call crash_exit
            end if

          case (20)
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *)
            write (52, *) 'There are: ', ios, 'zones'
            close(52)

          case (21)
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *) 'There was a problem reading config. info for zone ', ios
            close(52)
            write (*, *) 'There was a problem reading config. info for zone ', ios
            call crash_exit

          case (22)
            if (ios .ne. 0) then
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) 'There was a problem reading the exodus build flag'
                close(52)
                write (*, *) 'There was a problem reading the exodus build flag'
                call crash_exit
            end if

          case (23)
            if (ios .ne. 0) then
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) 'There was a problem reading the location of the exodus executable'
                close(52)
                write (*, *) 'There was a problem reading the location of the exodus executable'
                call crash_exit
            end if

          case (24)
            inquire (file='surface.1.node', exist=exst)
            if (.not. exst) goto 10
            inquire (file='surface.1.ele', exist=exst)
            if (.not. exst) goto 10
            inquire (file='surface.1.neigh', exist=exst)
            if (.not. exst) goto 10
            inquire (file='surface.1.edge', exist=exst)
            if (.not. exst) goto 10
10          continue
            if (.not. exst) then
                write (*, *) 'It appears there was a problem building the surface mesh'
                write (*, *) 'Cannot find one of the surface mesh files surface.1.*'
                write (*, *) 'Aborting'
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) 'It appears there was a problem building the surface mesh'
                write (52, *) 'Cannot find one of the surface mesh files surface.1.*'
                write (52, *) 'Aborting'
                close(52)
                call crash_exit
            end if

          case (25)
            write (*, *) 'It appears tetgen failed, please check for errors in'
            write (*, *) 'the mesh configuration file'
            write (*, *) 'Aborting'
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *) 'It appears tetgen failed, please check for errors in'
            write (52, *) 'the mesh configuration file'
            write (52, *) 'Aborting'
            close(52)
            call crash_exit

          case (26)
            write (*, *) 'A PLC uses control point ', ios, 'which is out of range'
            write (*, *) 'Aborting'
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *) 'A PLC uses control point ', ios, 'which is out of range'
            write (52, *) 'Aborting'
            close(52)
            call crash_exit

          case (27)

            write (*, *)
            write (*, *) 'ERROR'
            write (*, *) 'There are only ', ios, ' boundary points specified in: ', trim(cfg_filename)
            write (*, *) 'There must be at least 3 boundary points specified (boundary flag = 2)'
            write (*, *) 'Aborting ...'
            write (*, *)

            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *)
            write (52, *) 'ERROR'
            write (52, *) 'There are only ', ios, ' boundary points specified in: ', trim(cfg_filename)
            write (52, *) 'There must be at least 3 boundary points specified (boundary flag = 2)'
            write (52, *) 'Aborting ...'
            write (52, *)
            close(52)

            open (52, file="e4d.log", action='write', status='old', position='append')
            write (*, *)
            write (52, *) 'ERROR'
            write (52, *) 'There are only ', ios, ' boundary points specified in: ', trim(cfg_filename)
            write (52, *) 'There must be at least 3 boundary points specified (boundary flag = 2)'
            write (52, *) 'Aborting ...'
            write (52, *)
            close(52)

            call crash_exit

          case (28)
            if (ios .ne. 0) then
                open (52, file="mesh_build.log", action='write', status='old', position='append')
                write (52, *) 'There was a problem reading the mesh translation option'
                close(52)
                write (*, *) 'There was a problem reading the mesh translation option'
                call crash_exit
            end if

          case (29)
            open (52, file="mesh_build.log", action='write', status='old', position='append')
            write (52, *) 'The mesh translation option must be 0 for no translation'
            write (52, *) 'or 1 for translation.'
            write (52, *) 'You entered: ', ios
            write (51, *) 'Aborting ...'
            close(52)

            write (*, *) 'The mesh translation option must be 0 for no translation'
            write (*, *) 'or 1 for translation.'
            write (*, *) 'You entered: ', ios
            write (*, *) 'Aborting ...'
            call crash_exit
        end select
    end subroutine check_minp

    subroutine mesh_error_log(mode, error_message)
        implicit none

        integer, intent(in) :: mode
        character*255, intent(in) :: error_message
        integer :: io_status
        logical :: file_exists

        select case (mode)

          case (101)
            inquire (file='mesh_build.log', exist=file_exists)
            if (file_exists) then
                open (io_status, file='mesh_build.log', status='old', action='write', position='append')
            else
                open (io_status, file='mesh_build.log', status='new', action='write')
            end if
            write (io_status, "(a)") error_message
            close(io_status)
            write (*, "(a)") error_message

          case (102)
            inquire (file='mesh_build.log', exist=file_exists)
            if (file_exists) then
                open (io_status, file='mesh_build.log', status='old', action='write', position='append')
            else
                open (io_status, file='mesh_build.log', status='new', action='write')
            end if
            write (io_status, "(a)") error_message
            close(io_status)
            write (*, "(a)") error_message
        end select

    end subroutine mesh_error_log
!__________________________________________________________________________________________________

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine build_sig(cfg_filename, npre, nzon, zsigs, zsigsi, zlabs, ifg)
        implicit none
        character*40 :: cfg_filename
        integer :: npre, nzon
        real*8, dimension(nzon) :: zsigs, zsigsi
        integer, dimension(nzon, 2) :: zlabs
        integer :: i, j, nsig, a, b, m, n
        integer, dimension(:), allocatable :: zones
        logical :: ifg

        write (*, *) "BUILDING "//cfg_filename(1:npre)//"sigma"
        open (11, file=cfg_filename(1:npre)//"1.ele", status='old', action='read')
        open (12, file=cfg_filename(1:npre)//"sigma", status='replace', action='write')

        read (11, *) nsig
        allocate (zones(nsig))
        do i = 1, nsig
            read (11, *) j, a, b, m, n, zones(i)
        end do
        close(11)

        if (ifg) then
            write (12, *) nsig, 2
        else
            write (12, *) nsig, 1
        end if
        do i = 1, nsig
            do j = 1, nzon
                if (zones(i) == zlabs(j, 2)) then
                    if (ifg) then
                        write (12, *) zsigs(j), zsigsi(j)
                    else
                        write (12, *) zsigs(j)
                    end if
                    goto 100
                end if
            end do
            write (*, *) "While trying to build sigma, "
            write (*, *) "could not find a zone match for element ", i
            write (*, *) "setting too 1"
            write (12, *) "1.0 1.0 1.0"
100         continue
        end do
        close(12)

    end subroutine build_sig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine build_sig_triangle
        implicit none
        integer :: nsigs, i
        real*8 :: sigval
        open (66, file='surface.sigma', status='replace', action='write')
        open (55, file='surface.1.ele', status='old', action='read')
        read (55, *) nsigs
        close(55)
        sigval = 0.1
        write (66, *) nsigs
        do i = 1, nsigs
            write (66, *) sigval
        end do
        close(66)
    end subroutine build_sig_triangle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !_________________________________________________________________
    real*8 function getdet(matrix, n)
        implicit none
        real*8, dimension(n, n) :: matrix
        integer, intent(IN) :: n
        real*8 :: m, temp
        integer :: i, j, k, l
        logical :: DetExists = .true.
        l = 1
        !Convert to upper triangular form
        do k = 1, n - 1
            if (matrix(k, k) == 0) then
                DetExists = .false.
                do i = k + 1, n
                    if (matrix(i, k) .ne. 0) then
                        do j = 1, n
                            temp = matrix(i, j)
                            matrix(i, j) = matrix(k, j)
                            matrix(k, j) = temp
                        end do
                        DetExists = .true.
                        l = -l
                        exit
                    end if
                end do
                if (DetExists .eqv. .false.) then
                    getdet = 0
                    return
                end if
            end if
            do j = k + 1, n
                m = matrix(j, k)/matrix(k, k)
                do i = k + 1, n
                    matrix(j, i) = matrix(j, i) - m*matrix(k, i)
                end do
            end do
        end do

        !Calculate determinant by finding product of diagonal elements
        getdet = l
        do i = 1, n
            getdet = getdet*matrix(i, i)
        end do

    end function getdet
    !_____________________________________________________________________
    ! !> Export the mesh as a Multiblock VTK file
    ! !! @param[in] prefix The prefix for input/output files.
    ! !! @param[in] max_blocks The number of blocks to write.
    ! subroutine export_mesh_model(prefix, max_blocks)
    !     implicit none
    !     character(len=*), intent(in) :: prefix
    !     integer, intent(in) :: max_blocks
    !     character(len=255) :: node_file, ele_file, mod_filename, vtm_file, vtu_file
    !     real(8), allocatable :: points(:, :)
    !     integer, allocatable :: cells(:, :)
    !     integer, allocatable :: block_id(:, :)
    !     real(8), allocatable :: cell_data(:, :), output_cell_data(:, :)
    !     character(len=255), allocatable :: cell_data_labels(:)
    !     integer :: i, block, n_cells_used
    !     integer, allocatable :: selected_cells(:)
    !     character(len=10) :: block_str

    !     call system('mkdir -p '//trim("mesh"))

    !     node_file = trim(prefix)//".1.node"
    !     ele_file = trim(prefix)//".1.ele"
    !     mod_filename = trim(prefix)//".sigma"
    !     vtm_file = trim("mesh/")//trim(prefix)//".vtm"

    !     call read_node_file(node_file, points)
    !     call read_ele_file(ele_file, cells, block_id)
    !     call read_cdata_file(mod_filename, cell_data)

    !     if ((size(cell_data, dim=1) .eq. 1) .and. (mode .eq. 1) .or. (mode .eq. 31)) then
    !         allocate (output_cell_data(2, size(cell_data, dim=2)))
    !         output_cell_data(1, :) = cell_data(1, :)
    !         output_cell_data(2, :) = 1.0d0/cell_data(1, :)
    !         allocate (cell_data_labels(2))
    !         cell_data_labels = ["Conductivity [S/m]   ", "Resistivity [ohm-m]  "]
    !     else if ((size(cell_data, dim=1) .eq. 2)  .and. (mode .eq. 21) .or. (mode .eq. 41)) then
    !         allocate (output_cell_data(3, size(cell_data, dim=2)))
    !         output_cell_data(1, :) = cell_data(1, :)
    !         output_cell_data(2, :) = 1.0d0/cell_data(1, :)
    !         output_cell_data(3, :) = atan2(cell_data(2, :), cell_data(1, :))
    !         allocate (cell_data_labels(3))
    !         cell_data_labels = ["Conductivity [S/m]   ", "Resistivity [ohm-m]  ", &
    !                             "Phase [mrad]         "]
    !     end if

    !     write (*, '(a)') new_line('a')
    !     write (*, '(a)') "Writing VTK multiblock files..."
    !     do block = 1, max_blocks
    !         n_cells_used = count(block_id(1, :) == block)
    !         allocate (selected_cells(n_cells_used))
    !         selected_cells = pack([(i, i=1, size(cells, dim=2))], block_id(1, :) == block)
    !         write (block_str, '(I0)') block
    !         vtu_file = trim("mesh/")//trim(prefix)//"_zone_"//trim(block_str)//".vtu"
    !         call write_vtu_block(vtu_file, points, cells(:, selected_cells), &
    !                              output_cell_data(:, selected_cells), &
    !                              cell_data_labels, block)
    !         deallocate (selected_cells)
    !     end do

    !     call write_vtm_file(vtm_file, prefix, max_blocks)
    !     write (*, '(a)') vtm_file//" written."
    ! end subroutine export_mesh_model

end module mesh
