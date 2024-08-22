module vtk
    implicit none
contains
    !> Read .node file and extract node coordinates
    !! @param[in] filename Path to the .node file
    !! @param[out] points Array of coordinates of points
    subroutine read_node_file(filename, points)
        implicit none

        ! Input parameters
        character(len=*) :: filename                                        !< Path to the .node file

        ! Output parameters
        real(8), allocatable, intent(out) :: points(:, :)                    !< Array of coordinates of points

        ! Local variables
        integer :: n_points                                                 !< Number of points
        integer :: i, unit, io_stat                                         !< Loop index, file unit, I/O status
        integer :: id_flag                                                  !< Flag indicating presence of boundary ids
        integer :: index                                                    !< Index of the point
        character(len=255) :: line                                          !< Line buffer

        ! Open the file
        open (newunit=unit, file=filename, status='old', action='read')

        ! Read the header
        read (unit, '(A)', IOSTAT=io_stat) line
        if (io_stat .ne. 0) then
            print *, "Error reading header from node file: IOSTAT=", io_stat
            stop
        end if
        read (line, *, IOSTAT=io_stat) n_points
        if (io_stat .ne. 0) then
            print *, "Error parsing header from node file: IOSTAT=", io_stat
            stop
        end if

        ! Allocate arrays
        allocate (points(3, n_points))

        ! Read point data
        do i = 1, n_points
            read (unit, '(A)', IOSTAT=io_stat) line
            if (io_stat .ne. 0) then
                print *, "Error reading data for point ", i, ": IOSTAT=", io_stat
                stop
            end if
            read (line, *, IOSTAT=io_stat) index, points(1:3, i)
            if (io_stat .ne. 0) then
                print *, "Error parsing data for point ", i, ": IOSTAT=", io_stat
                stop
            end if
        end do

        ! Close the file
        close (unit)
    end subroutine read_node_file

    !> Read .ele file and extract tetrahedra connectivity
    !! @param[in] filename Path to the .ele file
    !! @param[out] cells Connectivity array of tetrahedra
    !! @param[out] block_id Array of block IDs associated with cells
    subroutine read_ele_file(filename, cells, block_id)
        implicit none

        ! Input parameters
        character(len=*) :: filename                                    !< Path to the .ele file

        ! Output parameters
        integer, allocatable, intent(out) :: cells(:, :)       !< Connectivity array of tetrahedra
        integer, allocatable, intent(out) :: block_id(:, :)    !< Array of block IDs associated with cells

        ! Local variables
        integer :: n_cells                      !< Number of tetrahedra
        integer :: n_points_per_cell, id_flag   !< Number of nodes per cell, flag indicating presence of block ids
        integer :: i, unit, io_stat             !< Loop index, file unit, I/O status
        character(len=255) :: line              !< Line buffer
        integer :: index                        !< Index of the cell

        ! Open the file
        open (newunit=unit, file=filename, status='old', action='read')

        ! Read the header
        read (unit, '(A)', IOSTAT=io_stat) line
        if (io_stat .ne. 0) then
            print *, "Error reading header from ele file: IOSTAT=", io_stat
            stop
        end if
        read (line, *, IOSTAT=io_stat) n_cells, n_points_per_cell, id_flag
        if (io_stat .ne. 0) then
            print *, "Error parsing header from ele file: IOSTAT=", io_stat
            stop
        end if

        ! Check cell type
        if (n_points_per_cell .ne. 4) then
            print *, "Error: Cells must be tetrahedral."
            stop
        end if

        ! Allocate arrays
        allocate (cells(n_points_per_cell, n_cells))
        allocate (block_id(id_flag, n_cells))

        ! Read cell data
        do i = 1, n_cells
            read (unit, '(A)', IOSTAT=io_stat) line
            if (io_stat .ne. 0) then
                print *, "Failed to read data for cell ", i, ": IOSTAT=", io_stat
                stop
            end if
            read (line, *, IOSTAT=io_stat) index, cells(1:n_points_per_cell, i), block_id(1:id_flag, i)
            if (io_stat .ne. 0) then
                print *, "Error parsing data for cell ", i, ": IOSTAT=", io_stat
                stop
            end if
        end do

        ! Close the file
        close (unit)
    end subroutine read_ele_file

    subroutine read_cell_data(filename, cell_data)
        implicit none

        ! Input parameters
        character(len=*) :: filename                                  !< Path to the cell data file

        ! Output parameters
        real(8), allocatable, intent(out) :: cell_data(:, :)          !< Cell data arrays

        ! Local variables
        integer :: n_cells                  !< Number of tetrahedra
        integer :: n_attr                   !< Number of attributes
        integer :: i, unit, io_stat         !< Loop index, file unit, I/O status
        character(len=255) :: line          !< Line buffer

        ! Open the file
        open (newunit=unit, file=filename, status='old', action='read')

        ! Read the header
        read (unit, '(A)', IOSTAT=io_stat) line
        if (io_stat .ne. 0) then
            print *, "Error reading header from file: IOSTAT=", io_stat
            stop
        end if
        read (line, *, IOSTAT=io_stat) n_cells, n_attr
        if (io_stat .ne. 0) then
            print *, "Error parsing header from file: IOSTAT=", io_stat
            stop
        end if


        allocate (cell_data(n_attr, n_cells))

        ! Read cell data
        do i = 1, n_cells
            read (unit, '(A)', IOSTAT=io_stat) line
            if (io_stat .ne. 0) then
                print *, "Failed to read data for cell ", i, ": IOSTAT=", io_stat
                stop
            end if
            read (line, *, IOSTAT=io_stat) cell_data(1:n_attr, i)
            if (io_stat .ne. 0) then
                print *, "Error parsing data for cell ", i, ": IOSTAT=", io_stat
                stop
            end if
        end do

        ! Close the file
        close (unit)
    end subroutine read_cell_data

    !> Writes mesh data to a VTK file.
    !! This subroutine writes point and cell data to a VTK formatted file.
    !! It supports filtering cells based on block IDs.
    !!
    !! @param[in] filename The name of the VTK file to write.
    !! @param[in] points Array containing point coordinates.
    !! @param[in] cells Array containing cell connectivity information.
    !! @param[in] cell_data Array containing additional cell data (e.g., conductivity).
    !! @param[in] cell_data_labels Array containing labels for the cell data.
    !! @param[in] block_id Array containing block IDs associated with cells.
    !! @param[in] block_list Optional filter to specify which cells to include based on their block ID.
    subroutine write_vtk(filename, points, cells, cell_data, cell_data_labels, block_id, block_list)
        implicit none

        ! Input parameters
        character(len=*) :: filename                                                !< Path to the output VTK file
        real(8), dimension(:, :), intent(in) :: points                              !< Array of coordinates of points
        integer, dimension(:, :), intent(in) :: cells                               !< Connectivity array of tetrahedra
        real(8), dimension(:, :), intent(in) :: cell_data                           !< Array of data associated with cells
        character(len=255), dimension(:), intent(in) :: cell_data_labels            !< Array of labels for the cell data
        integer, dimension(:, :), intent(in) :: block_id                            !< Optional array of block IDs associated with cells
        integer, dimension(:), intent(in), optional :: block_list                   !< Optional array of block IDs to include.

        ! Local variables
        integer :: n_points, n_cells, n_cells_used                                  !< Number of points, cells, filtered cell count
        integer :: i, j, unit                                                       !< Loop index, file unit

        open (newunit=unit, file=filename, status='replace', action='write')

        write (unit, '(a)') "# vtk DataFile Version 3.0"
        write (unit, '(a)') "Converted from TetGen mesh to VTK format"
        write (unit, '(a)') "ASCII"
        write (unit, '(a)') "DATASET UNSTRUCTURED_GRID"

        n_points = size(points, dim=2)
        n_cells = size(cells, dim=2)

        ! Check cell data array
        if (size(cell_data, dim=2) .ne. n_cells) then
            print *, "Error: Cell data array does not match the number of cells."
            stop
        end if

        ! Check cell data labels
        if (size(cell_data_labels) .ne. size(cell_data, dim=1)) then
            print *, "Error: Cell data labels array does not match the number of cell data arrays."
            print *, "Expected ", size(cell_data, dim=1), " labels, got ", size(cell_data_labels)
            stop
        end if

        ! Check block ID array
        if (size(block_id, dim=2) .ne. n_cells) then
            print *, "Error: Block ID array does not match the number of cells."
            stop
        end if

        ! Write point coordinates
        write (unit, '("POINTS ", i0, " double")') n_points
        do i = 1, n_points
            write (unit, '(3F17.6)') points(:, i)
        end do

        ! Check if block filter is provided
        if (present(block_list)) then
            ! Initialise count for filtered cells
            n_cells_used = 0
            do i = 1, n_cells
                if (all(block_id(1, i) .ne. block_list)) cycle
                n_cells_used = n_cells_used + 1
            end do
        else
            n_cells_used = n_cells
        end if

        write (unit, '("CELLS ", i0, " ", i0)') n_cells_used, 5*n_cells_used
        do i = 1, n_cells
            if (present(block_list)) then
                if (all(block_id(1, i) .ne. block_list)) cycle
                write (unit, '(i0, 1x, i0, 1x, i0, 1x, i0, 1x, i0, 1x)') 4, cells(1:4, i) - 1
            else
                write (unit, '(i0, 1x, i0, 1x, i0, 1x, i0, 1x, i0, 1x)') 4, cells(1:4, i) - 1
            end if
        end do

        write (unit, '("CELL_TYPES ", i0)') n_cells_used
        do i = 1, n_cells
            if (present(block_list)) then
                if (all(block_id(1, i) .ne. block_list)) cycle
            end if
            write (unit, '(i0)') 10
        end do

        write (unit, '("CELL_DATA ", i0)') n_cells_used
        write (unit, '(a)') "SCALARS Block_Id int 1"
        write (unit, '(a)') "LOOKUP_TABLE default"
        do i = 1, n_cells
            if (present(block_list)) then
                if (all(block_id(1, i) .ne. block_list)) cycle
            end if
            write (unit, '(i0)') block_id(1, i)
        end do

        do i = 1, size(cell_data, dim=1)
            write (unit, '(a)') "SCALARS "//trim(cell_data_labels(i))//" double 1"
            write (unit, '(a)') "LOOKUP_TABLE default"
            do j = 1, n_cells
                if (present(block_list)) then
                    if (all(block_id(1, j) .ne. block_list)) cycle
                end if
                write (unit, '(f17.6)') cell_data(i, j)
            end do
        end do

        close (unit)
    end subroutine write_vtk

    subroutine convert_to_vtk(prefix, max_blocks)
        implicit none
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: max_blocks
        character(len=255) :: node_file, ele_file, sig_file, vtk_file
        real(8), allocatable :: points(:, :)
        integer, allocatable :: cells(:, :)
        integer, allocatable :: block_id(:, :)
        real(8), allocatable :: cell_data(:, :), output_cell_data(:, :)
        character(len=255), allocatable :: cell_data_labels(:)
        integer, allocatable :: block_list(:)
        integer :: i

        ! Construct file names
        node_file = trim(prefix)//".1.node"
        ele_file = trim(prefix)//".1.ele"
        sig_file = trim(prefix)//".sig"
        vtk_file = trim(prefix)//".vtk"

        ! Allocate and initialize block list
        allocate (block_list(max_blocks))
        do i = 1, max_blocks
            block_list(i) = i
        end do

        ! Read the .node file
        call read_node_file(node_file, points)

        ! Read the .ele file
        call read_ele_file(ele_file, cells, block_id)

        ! Read the .sig file
        call read_cell_data(sig_file, cell_data)

        if (size(cell_data, dim=1) .eq. 1) then
            allocate(output_cell_data(2, size(cell_data, dim=2)))
            output_cell_data(1, :) = cell_data(1, :)
            output_cell_data(2, :) = 1.0d0 / cell_data(1, :)
            allocate(cell_data_labels(2))
            cell_data_labels = ["Conductivity_[S/m]    ", "Resistivity_[ohm-m]   "]
        else if (size(cell_data, dim=1) .eq. 2) then
            allocate(output_cell_data(3, size(cell_data, dim=2)))
            output_cell_data(1, :) = cell_data(1, :)
            output_cell_data(2, :) = 1.0d0 / cell_data(1, :)
            output_cell_data(3, :) = atan2(cell_data(2, :), cell_data(1, :))
            allocate(cell_data_labels(3))
            cell_data_labels = ["Conductivity_[S/m]    ", "Resistivity_[ohm-m]   ", &
                               "Phase_[mrad]          "]
        end if

        ! Write to VTK file
        write (*, '(a)') new_line('a')
        write (*, '(a)') "Writing VTK file..."
        call write_vtk(vtk_file, points, cells, output_cell_data, &
        cell_data_labels, block_id, block_list)
        write (*, '(a)') vtk_file//" written."

    end subroutine convert_to_vtk

end module vtk
