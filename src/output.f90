! -----------------------------------------------------------------------------
! @brief Module to handle output of results.
! @author ofgn, Tim C. Johnson
! @date 2024-07-01
! @note Refactored several routines - ofgn: 2024-07-01
! -----------------------------------------------------------------------------
module output
    use iso_fortran_env, only: int32
    use vars
    use e4d_report
    use reorder_mesh
    use vtk

    implicit none

    integer, dimension(:, :), allocatable :: ipot
    real, dimension(:, :), allocatable :: pot

contains

    ! --------------------------------------------------------------------------
    ! @brief Write the simulated data to a .dat file.
    ! --------------------------------------------------------------------------
    subroutine export_simulated_data()

        implicit none

        integer(int32) :: simulated_data_flag                                   ! Flag for simulated data (0 = no, 1 = yes)
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: i                                                     ! Loop counter
        logical :: exists                                                       ! File existence check
        character(len=4096) :: file_path                                        ! File path for output

        inquire(file=trim(out_file), exist=exists)

        if (.not. exists) then
            call log_error("Cannot find the output options file: " &
                // trim(out_file))
            return
        end if

        open(newunit=unit, file=trim(out_file), status="old", action="read", &
            iostat=io_status)

        if (io_status .ne. 0) then
            call log_error("Problem reading the output options file: " &
                // trim(out_file))
            return
        end if

        read(unit, *, iostat=io_status) simulated_data_flag

        if (io_status .ne. 0) then
            call log_error("Problem reading the output flag in: " &
                // trim(out_file))
            close(unit)
            return
        end if

        call execute_command_line("mkdir -p data", exitstat=io_status)

        if (io_status .ne. 0) then
            call log_error("Failed to create data directory")
            return
        end if

        if (simulated_data_flag .eq. 1) then

            if (invi .and. (mode .eq. 3)) then
                write(file_path, "(A,'/IP_',I3.3,'.dat')") "data", iter
            else if (mode .eq. 3) then
                write(file_path, "(A,'/DC_',I3.3,'.dat')") "data", iter
            end if

            open(newunit=unit, file=file_path, status="replace", &
                action="write", iostat=io_status)

            if (io_status .ne. 0) then
                call log_error("Failed to open file: " // trim(file_path))
                return
            end if

            if (invi .and. (mode .eq. 3)) then
                write(unit, "(A8, 4A12, 4A16)") &
                    "Index", "C1", "C2", "P1", "P2", &
                    "DC_Observed", "DC_Simulated", &
                    "IP_Observed", "IP_Simulated"
            else if (mode .eq. 3) then
                write(unit, "(A8, 4A12, 2A16)") &
                    "Index", "C1", "C2", "P1", "P2", &
                    "DC_Observed", "DC_Simulated"
            end if

            do i = 1, nm
                if (invi .and. (mode .eq. 3)) then
                    write(unit, "(I8, 4I12, 4F16.6)") &
                        i, s_conf(i, 1:4), &
                        dobs(i), dpred(i), &
                        dobsi(i), dpredi(i)
                else if (mode .eq. 3) then
                    write(unit, "(I8, 4I12, 2F16.6)") &
                        i, s_conf(i, 1:4), &
                        dobs(i), dpred(i)
                end if
            end do

            close(unit)
        end if
        
        return
    end subroutine export_simulated_data

    ! --------------------------------------------------------------------------
    ! @brief Exports the simulated data to a .srv file.
    ! --------------------------------------------------------------------------
    subroutine export_simulated_survey()
        implicit none

        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: i                                                     ! Loop counter
        logical :: exists                                                       ! File existence check
        character(len=4096) :: file_path                                        ! File path for output

        write(file_path, "(A,'.srv')") trim(sig_filename)

        call nreport(69)

        open(newunit=unit, file=file_path, status="replace", action="write", &
            iostat=io_status)

        if (io_status .ne. 0) then
            call log_error("Failed to open file: " // trim(file_path))
            return
        end if

        write(unit, *) ne

        do i = 1, ne
            write(unit, "(I10,3F15.5,I10)") i, &
                e_pos(i, 1) + xorig, &
                e_pos(i, 2) + yorig, &
                e_pos(i, 3) + zorig, &
                int(e_pos(i, 4))
        end do

        write(unit, *)
        write(unit, *) nm

        if (i_flag) then
            do i = 1, nm
                write(unit, "(I8,4I10,4G15.5)") i, s_conf(i, 1:4), &
                    dpred(i), 0.05 * abs(dpred(i)), &
                    dpredi(i), 0.05 * abs(dpredi(i))
            end do
        else
            do i = 1, nm
                write(unit, "(I8,4I10,2G15.5)") i, s_conf(i, 1:4), &
                    dpred(i), 0.05 * abs(dpred(i))
            end do
        end if

        close(unit)
    end subroutine export_simulated_survey

    ! --------------------------------------------------------------------------
    ! @brief Writes the real conductivity model to a .sig file.
    ! @param[in] model_name The name of the model file.
    ! --------------------------------------------------------------------------
    subroutine export_real_conductivity_model(model_name)
        implicit none

        character(len=*), intent(in), optional :: model_name                    ! The name of the model file.

        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: i                                                     ! Loop counter
        character(4096) :: file_path                                            ! File path for output

        call execute_command_line("mkdir -p models", exitstat=io_status)

        if (io_status .ne. 0) then
            call log_error("Failed to create models directory")
            return
        end if

        if (present(model_name)) then
            file_path = "models/" // trim(model_name) // ".sig"
        else
            write(file_path, "(A,'/DC_',I3.3,'.sig')") "models", iter
        end if

        open(newunit=unit, file=trim(file_path), status="replace", action="write", &
            iostat=io_status)

        if (io_status .ne. 0) then
            call log_error("Failed to open file: " // trim(file_path))
            return
        end if

        if (allocated(element_map)) then
            write(unit, *) size(element_map), merge(2, 1, i_flag), chi2
            do i = 1, size(element_map)
                if (element_map(i) .ne. 0) then
                    if (i_flag) then
                        write(unit, *) sigma_re(element_map(i)), &
                            real(sigma_im(element_map(i)))
                    else
                        write(unit, *) sigma_re(element_map(i))
                    end if
                else
                    if (i_flag) then
                        write(unit, *) -999, -999
                    else
                        write(unit, *) -999
                    end if
                end if
            end do
        else
            write(unit, *) n_elements, merge(2, 1, i_flag), chi2
            if (i_flag) then
                do i = 1, n_elements
                    write(unit, *) sigma_re(i), real(sigma_im(i))
                end do
            else
                do i = 1, n_elements
                    write(unit, *) sigma_re(i)
                end do
            end if
        end if

        close(unit)
    end subroutine export_real_conductivity_model

    ! --------------------------------------------------------------------------
    ! @brief Writes the complex conductivity model to a .sig file.
    ! @param[in] model_name The name of the model file.
    ! --------------------------------------------------------------------------
    subroutine export_complex_conductivity_model(model_name)
        implicit none

        character(len=*), intent(in), optional :: model_name                    ! The name of the model file.

        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: i                                                     ! Loop counter
        character(4096) :: file_path                                            ! File path for output

        call execute_command_line("mkdir -p models", exitstat=io_status)

        if (io_status .ne. 0) then
            call log_error("Failed to create models directory")
            return
        end if

        if (present(model_name)) then
            file_path = "models/" // trim(model_name) // ".sig"
        else
            write(file_path, "(A,'/IP_',I3.3,'.sig')") "models", iter
        end if

        open(newunit=unit, file=trim(file_path), status="replace", action="write", &
            iostat=io_status)

        if (io_status .ne. 0) then
            call log_error("Failed to open file: " // trim(file_path))
            return
        end if

        if (allocated(element_map)) then
            write(unit, *) size(element_map), 2
            do i = 1, size(element_map)
                if (element_map(i) .ne. 0) then
                    write(unit, *) sigma_re(element_map(i)), &
                        real(sigma_im(element_map(i)))
                else
                    write(unit, *) -999, -999
                end if
            end do
        else
            write(unit, *) n_elements, 2
            do i = 1, n_elements
                write(unit, *) sigma_re(i), sigma_im(i)
            end do
        end if

        close(unit)
    end subroutine export_complex_conductivity_model

    ! --------------------------------------------------------------------------
    ! @brief Exports the model as a VTK file.
    ! --------------------------------------------------------------------------
    subroutine export_vtk(node_file_path, ele_file_path, sigma_file_path, &
            vtk_file_path, n_zones)
            
        use vtk

        implicit none

        character(len=*) :: node_file_path                                      ! File path for node file
        character(len=*) :: ele_file_path                                       ! File path for element file
        character(len=*) :: sigma_file_path                                     ! File path for sigma file
        character(len=*) :: vtk_file_path                                       ! File path for VTK file
        integer(int32) :: n_zones                                               ! Number of zones

        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: i                                                     ! Loop counter
        logical :: exists                                                       ! File existence check
        logical, allocatable :: cell_mask(:)                                    ! Cell mask
        type(VtkUnstructuredGrid) :: u_grid                                     ! VTK unstructured grid
        type(VtkData) :: point_data                                             ! VTK point data                         
        type(VtkData) :: cell_data                                              ! VTK cell data   
                           
        call execute_command_line("mkdir -p vtk", exitstat=io_status)

        if (io_status .ne. 0) then
            call log_error("Failed to create vtk directory")
            return
        end if

        call u_grid%read_tetgen_node(node_file_path)
        call u_grid%read_tetgen_ele(ele_file_path)
        call cell_data%read_scalar_real64(sigma_file_path)

        cell_data%scalar_real64(3, :) = atan2(cell_data%scalar_real64(3, :), &
            cell_data%scalar_real64(2, :)) * 1000.0
        cell_data%scalar_real64(2, :) = 1.0 / cell_data%scalar_real64(2, :)

        cell_data%scalar_real64_labels(1) = 'Zone'
        cell_data%scalar_real64_labels(2) = 'Resistivity[Ohm-m]'
        cell_data%scalar_real64_labels(3) = 'Phase[mrad]'

        allocate(cell_mask(u_grid%n_cells))
        cell_mask(:) = (cell_data%scalar_real64(1, :) .lt. real(n_zones, real64))

        call u_grid%mask_cells(cell_mask, cell_data=cell_data)
        call u_grid%write_legacy_vtk(vtk_file_path, cell_data=cell_data)
        
        return
    end subroutine export_vtk


    ! !> @brief Exports the model as a Multiblock VTK file.
    ! !!
    ! !! This subroutine writes the mesh data as VTK files with multiblock structure. 
    ! !! The files are named using the given prefix and are saved in the "mesh" directory, 
    ! !! which is created if it doesn't already exist.
    ! !!
    ! !! @param[in] node_fname The filename for the node file.
    ! !! @param[in] ele_fname The filename for the element file.
    ! !! @param[in] mod_fname The filename for the model file.
    ! !! @param[in] max_blocks The number of blocks to write.
    ! subroutine export_3d(node_fname, ele_fname, mod_fname, max_blocks)
    !     implicit none

    !     character(len=*), intent(in) :: node_fname, ele_fname, mod_fname
    !     integer, intent(in) :: max_blocks
    !     character(255) :: mod_prefix, vtm_fname, vtu_fname
    !     real(8), allocatable :: points(:, :)
    !     integer, allocatable :: cells(:, :)
    !     integer, allocatable :: block_id(:, :)
    !     real(8), allocatable :: cell_data(:, :), output_cell_data(:, :)
    !     character(255), allocatable :: cell_data_labels(:)
    !     integer :: i, block, n_cells_used
    !     integer, allocatable :: selected_cells(:)
    !     character(10) :: block_str
    !     integer :: io_status

    !     ! Create the mesh directory if it doesn't exist
    !     call execute_command_line('mkdir -p 3d', exitstat=io_status)
    !     if (io_status .ne. 0) then
    !         call log_error("Failed to create mesh directory")
    !         return
    !     end if
    !     if (i_flag) then
    !         write(mod_prefix, "(A,'/IP_',I0,'.sig')") "3d", iter
    !     else
    !         write(mod_prefix, "(A,'/DC_',I0,'.sig')") "3d", iter
    !     end if

    !     vtm_fname = trim(mod_prefix) // ".vtm"
        

    !     print *, node_fname
    !     print *, ele_fname
    !     print *, mod_fname
    !     print *, vtm_fname

    !     ! Read the node, element, and cell data
    !     call read_node_file(trim(node_fname), points)
    !     call read_ele_file(trim(ele_fname), cells, block_id)
    !     call read_cdata_file(trim(mod_fname), cell_data)

    !     ! Process cell data based on mode and data dimensions
    !     if ((size(cell_data, dim=1) == 1) .and. (mode == 1 .or. mode == 31)) then
    !         allocate(output_cell_data(2, size(cell_data, dim=2)))
    !         output_cell_data(1, :) = cell_data(1, :)
    !         output_cell_data(2, :) = 1.0d0 / cell_data(1, :)
    !         allocate(cell_data_labels(2))
    !         cell_data_labels = ["Conductivity [S/m]", "Resistivity [ohm-m]"]
    !     else if ((size(cell_data, dim=1) == 2) .and. (mode == 21 .or. mode == 41)) then
    !         allocate(output_cell_data(3, size(cell_data, dim=2)))
    !         output_cell_data(1, :) = cell_data(1, :)
    !         output_cell_data(2, :) = 1.0d0 / cell_data(1, :)
    !         output_cell_data(3, :) = atan2(cell_data(2, :), cell_data(1, :))
    !         allocate(cell_data_labels(3))
    !         cell_data_labels = ["Conductivity [S/m]", "Resistivity [ohm-m]", "Phase [mrad]"]
    !     end if

    !     ! Inform the user about the file writing process
    !     write(*, '(a)') new_line('a')
    !     write(*, '(a)') "Writing VTK multiblock files..."

    !     ! Write each block as a VTK file
    !     do block = 1, max_blocks
    !         n_cells_used = count(block_id(1, :) == block)
    !         allocate(selected_cells(n_cells_used))
    !         selected_cells = pack([(i, i=1, size(cells, dim=2))], block_id(1, :) == block)
    !         write(block_str, '(I0)') block
    !         vtu_fname = trim(mod_prefix) // "_zone_" // trim(block_str) // ".vtu"
    !         call write_vtu_block(vtu_fname, points, cells(:, selected_cells), &
    !                             output_cell_data(:, selected_cells), &
    !                             cell_data_labels, block)
    !         deallocate(selected_cells)
    !     end do

    !     ! Write the VTM file that links all VTU blocks
    !     call write_vtm_file(vtm_fname, prefix, max_blocks)
    !     write(*, '(a)') trim(vtm_fname) // " written."

    ! end subroutine export_3d

    !_________________________________________________________________________________

    !_______________________________________________________________________________________
    subroutine write_sigma_rttl(filename)
        implicit none
        character*40 :: filename, oname
        integer(int32) :: i, np

      !!find the location of the . in the file name
        np = 40
        do i = 1, np
            if (filename(i:i) == ".") then
                np = i
                exit
            end if
        end do

        oname = filename
        oname(np:np + 3) = ".sigma"
        open (12, file=trim(oname)//".part", status="replace", action="write")
        if (i_flag) then
            write (12, *) n_elements, 2, chi2
        else
            write (12, *) n_elements, 1, chi2
        end if

        if (i_flag) then
            do i = 1, n_elements
                write (12, *) sigma_re(i), real(sigma_im(i))
            end do
        else
            do i = 1, n_elements
                write (12, *) sigma_re(i)
            end do
        end if
        close (12)
        call system("mv "//trim(oname)//".part "//trim(oname))

    end subroutine write_sigma_rttl
    !_______________________________________________________________________________________
    subroutine write_sigma_tl(tm)
        implicit none
        integer :: i
        real :: tm
        character(8) ::ts
        character(20) :: filename = ""
        write (ts, "(f8.3)") tm
        write (filename, "(A6,A8)") "tl_sigma_", adjustl(ts)
        open (12, file=filename, status="replace", action="write")

        if (allocated(element_map)) then
            if (i_flag) then
                write (12, *) size(element_map), 2, chi2
                do i = 1, size(element_map)
                    if (element_map(i) .ne. 0) then
                        write (12, *) sigma_re(element_map(i)), real(sigma_im(element_map(i)))
                    else
                        write (12, *) - 999, -999
                    end if
                end do
            else
                write (12, *) size(element_map), 1, chi2
                do i = 1, size(element_map)
                    if (element_map(i) .ne. 0) then
                        write (12, *) sigma_re(element_map(i))
                    else
                        write (12, *) - 999
                    end if
                end do
            end if
        else
            if (i_flag) then
                write (12, *) n_elements, 2, chi2
            else
                write (12, *) n_elements, 1, chi2
            end if

            if (i_flag) then
                do i = 1, n_elements
                    write (12, *) sigma_re(i), real(sigma_im(i))
                end do
            else
                do i = 1, n_elements
                    write (12, *) sigma_re(i)
                end do
            end if
        end if
        close (12)
    end subroutine write_sigma_tl
    !_______________________________________________________________________________________

    !_______________________________________________________________________________________
    subroutine write_pots
        implicit none
        integer :: resid_flag, pot_flag, npot, o_opt, io_status, jflag
        logical :: fcheck
        character*80 :: resid_prefix
        character*20 :: filename, jformat
        integer :: i, a, b, j, emin, emax, ra, rb
        integer, dimension(2) :: spack
        real, dimension(nnodes) :: pa, pb, rp, cp
        integer ::  status(MPI_STATUS_SIZE)

        inquire (file=trim(out_file), exist=fcheck); if (.not. fcheck) goto 10

        call nreport(21)
        open (15, file=out_file, status="old", action="read")
        read (15, *, IOSTAT=io_status) resid_flag; if (io_status .ne. 0) goto 11
        read (15, *, IOSTAT=io_status) resid_prefix; if (io_status .ne. 0) goto 12
        read (15, *, IOSTAT=io_status) npot; if (io_status .ne. 0) goto 13

        if (npot > 0) then
            allocate (ipot(npot, 2))
            do i = 1, npot
                read (15, *, IOSTAT=io_status) ipot(i, 1); if (io_status .ne. 0) goto 14
            end do
        end if

      !!read the jacobian output flag
        read (15, *, IOSTAT=io_status) jflag
        if (io_status .ne. 0) then
            open (51, file="e4d.log", status="old", action="write", position="append")
            write (51, *) " There was a problem Jacobian matrix output option: ", i, " in: ", trim(out_file)
            write (51, *) " Not printing the Jacobian matrix."
            close (51)
            write (*, *)
            write (*, *) " There was a problem Jacobian matrix output option: ", i, " in: ", trim(out_file)
            write (*, *) " Not printing the Jacobian matrix."
            goto 9
        end if
        read (15, *, IOSTAT=io_status) jformat
        if (io_status .ne. 0) then
            open (51, file="e4d.log", status="old", action="write", position="append")
            write (51, *) " There was a problem Jacobian output format option: ", i, " in: ", trim(out_file)
            write (51, *) " Printing in binary format"
            close (51)
            write (*, *)
            write (*, *) " There was a problem Jacobian matrix output option: ", i, " in: ", trim(out_file)
            write (*, *) " Printing in binary format"

        end if

        if (jflag == 1) then
            jaco_out_opt = .true.
            jaco_ascii_opt = .true.
            open (51, file="e4d.log", status="old", action="write", position="append")
            write (51, *) " Printing Jacobian matrix in "
            write (*, *) " Printing Jacobian matrix in "

            if (trim(jformat) == "ASCII" .or. trim(jformat) == "ascii") then
                jaco_ascii_opt = .true.
                write (51, *) " ascii format"
                write (*, *) " ascii format"
            else
                jaco_ascii_opt = .false.
                write (51, *) " binary format"
                write (*, *) " binary format"
            end if
            close (51)
        else
            jaco_out_opt = .false.
            open (51, file="e4d.log", status="old", action="write", position="append")
            write (51, *) " Not printing Jacobian matrix"
            write (*, *) " Not printing Jacobian matrix"
            close (51)
        end if

9       continue
        close (15)
        do i = 1, npot
            pa = 0
            pb = 0
            if (ipot(i, 1) > nm) goto 100
            a = s_conf(ipot(i, 1), 1)
            b = s_conf(ipot(i, 1), 2)
            do j = 1, n_rank - 1
                emin = eind(j, 1); emax = eind(j, 2)
                if ((emin .le. a) .and. (emax .ge. a)) ra = j
                if ((emin .le. b) .and. (emax .ge. b)) rb = j
            end do

            if (a .ne. 0) then
                spack(1) = ra
                spack(2) = a
                call send_commando(23)
                call MPI_BCAST(spack, 2, MPI_INTEGER, 0, E4D_COMM, ierr)
                call MPI_RECV(pa, nnodes, MPI_DOUBLE, ra, 0, E4D_COMM, status, ierr)
            end if

            if (b .ne. 0) then
                spack(1) = rb
                spack(2) = b
                call send_commando(23)
                call MPI_BCAST(spack, 2, MPI_INTEGER, 0, E4D_COMM, ierr)
                call MPI_RECV(pb, nnodes, MPI_DOUBLE, rb, 0, E4D_COMM, status, ierr)
            end if

            do j = 1, nnodes
                rp(j) = pa(j) - pb(j)
            end do

            if (i_flag) then
                pa = 0
                pb = 0
                if (a .ne. 0) then
                    spack(1) = ra
                    spack(2) = a
                    call send_commando(123)
                    call MPI_BCAST(spack, 2, MPI_INTEGER, 0, E4D_COMM, ierr)
                    call MPI_RECV(pa, nnodes, MPI_DOUBLE, ra, 0, E4D_COMM, status, ierr)
                end if

                if (b .ne. 0) then
                    spack(1) = rb
                    spack(2) = b
                    call send_commando(123)
                    call MPI_BCAST(spack, 2, MPI_INTEGER, 0, E4D_COMM, ierr)
                    call MPI_RECV(pb, nnodes, MPI_DOUBLE, rb, 0, E4D_COMM, status, ierr)
                end if

                do j = 1, nnodes
                    cp(j) = pa(j) - pb(j)
                end do

            end if

            write (filename, "(A,I0)") "potential.", ipot(i, 1)
            open (27, file=filename, status="replace", action="write")

            if (allocated(node_map)) then

                if (i_flag) then
                    write (27, *) size(node_map), 2, ipot(i, 1)
                    do j = 1, size(node_map)
                        if (node_map(i) .ne. 0) then
                            write (27, *) rp(node_map(j)), cp(node_map(j))
                        else
                            write (27, *) - 999, -999
                        end if
                    end do
                else
                    write (27, *) size(node_map), 1, ipot(i, 1)
                    do j = 1, size(node_map)
                        if (node_map(i) .ne. 0) then
                            write (27, *) rp(node_map(j))
                        else
                            write (27, *) - 999
                        end if
                    end do
                end if

            else
                if (i_flag) then
                    write (27, *) nnodes, 2, ipot(i, 1)
                    do j = 1, nnodes
                        write (27, *) rp(j), cp(j)
                    end do
                else
                    write (27, *) nnodes, 1, ipot(i, 1)
                    do j = 1, nnodes
                        write (27, *) rp(j)
                    end do
                end if

            end if
            close (27)

100         continue
        end do
        return

10      continue
        open (51, file="e4d.log", status="old", action="write", position="append")
        write (51, *)
        write (51, *) " Cannot find the output options file: ", trim(out_file)
        close (51)
        write (*, *)
        write (*, *) " Cannot find the output options file: ", trim(out_file)
        return

11      continue
        open (51, file="e4d.log", status="old", action="write", position="append")
        write (51, *)
        write (51, *) " The was a problem reading the first line in the output file: ", trim(out_file)
        close (51)
        write (*, *)
        write (*, *) " There was a problem reading the first line in the output file: ", trim(out_file)
        return

12      continue
        open (51, file="e4d.log", status="old", action="write", position="append")
        write (51, *)
        write (51, *) "There was a problem reading the simulated data file name in: ", trim(out_file)
        close (51)
        write (*, *)
        write (*, *) "The was a problem reading the simulated data file in: ", trim(out_file)
        return

13      continue
        open (51, file="e4d.log", status="old", action="write", position="append")
        write (51, *)
        write (51, *) " There was a problem reading the number of potential fields to write in: ", trim(out_file)
        close (51)
        write (*, *)
        write (*, *) " There was a problem reading the number of potential fields to write in: ", trim(out_file)
        return

14      continue
        open (51, file="e4d.log", status="old", action="write", position="append")
        write (51, *) " There was a problem reading potential field index: ", i, " in: ", trim(out_file)
        close (51)
        write (*, *)
        write (*, *) " There was a problem reading potential field index: ", i, " in: ", trim(out_file)
        return

    end subroutine write_pots
    !_______________________________________________________________________________________

    !_______________________________________________________________________________________
    subroutine write_sigiter
        implicit none
        character*40 :: fiter
        integer :: i

        write (fiter, "(A,I0)") "si.", iter
        open (23, file=trim(fiter), status="replace", action="write")
        if (i_flag) then
            write (23, *) n_elements, 2, chi2
            do i = 1, n_elements
                write (23, *) sigma_re(i), sigma_im(i)
            end do

        else
            write (23, *) n_elements, 1, chi2
            do i = 1, n_elements
                write (23, *) sigma_re(i)
            end do

        end if
        close (23)

    end subroutine write_sigiter
    !_______________________________________________________________________________________

    !_____________________________________________________________________
    subroutine send_commando(com)
      !!Send a general command to the slaves
        integer :: com

        call MPI_BCAST(com, 1, MPI_INTEGER, 0, E4D_COMM, ierr)

    end subroutine send_commando
    !____________________________________________________________________

end module output

