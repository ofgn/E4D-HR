module SparseMatrixCSB
    implicit none
  
    ! Define the CSBMatrix type
    type :: CSBMatrix
      integer :: n_rows          ! Total number of rows
      integer :: n_cols          ! Total number of columns
      integer :: block_size_row  ! Block size in rows
      integer :: block_size_col  ! Block size in columns
      integer :: n_blocks_row    ! Number of blocks in rows
      integer :: n_blocks_col    ! Number of blocks in columns
      integer :: nnz_local       ! Number of local non-zero elements
      integer, allocatable :: row_idx(:)  ! Local row indices within block
      integer, allocatable :: col_idx(:)  ! Local column indices within block
      real(kind=8), allocatable :: data(:) ! Local non-zero values
    end type CSBMatrix
  
  contains
  
    ! Initialize the CSB matrix and distribute among MPI processes
    subroutine initialize_csb_matrix(matrix, n_rows, n_cols, block_size_row, block_size_col, &
                                     row_indices, col_indices, values, comm)
      type(CSBMatrix), intent(out) :: matrix
      integer, intent(in) :: n_rows, n_cols, block_size_row, block_size_col
      integer, intent(in) :: row_indices(:), col_indices(:)
      real(kind=8), intent(in) :: values(:)
      integer, intent(in) :: comm
  
      integer :: num_procs, rank, ierr
      integer :: i, nnz
      integer :: n_blocks_row, n_blocks_col, block_row, block_col
      integer, allocatable :: counts(:), displs(:)
      integer :: nnz_local
      integer :: procs_per_row, procs_per_col
      integer :: owner_rank
  
      call MPI_Comm_size(comm, num_procs, ierr)
      call MPI_Comm_rank(comm, rank, ierr)
  
      ! Initialize matrix metadata
      matrix%n_rows = n_rows
      matrix%n_cols = n_cols
      matrix%block_size_row = block_size_row
      matrix%block_size_col = block_size_col
      n_blocks_row = ceiling(real(n_rows) / block_size_row)
      n_blocks_col = ceiling(real(n_cols) / block_size_col)
      matrix%n_blocks_row = n_blocks_row
      matrix%n_blocks_col = n_blocks_col
  
      ! Total number of non-zero elements
      nnz = size(values)
  
      ! Determine process grid dimensions
      procs_per_row = floor(sqrt(real(num_procs * n_blocks_row / n_blocks_col)))
      if (procs_per_row == 0) procs_per_row = 1
      procs_per_col = num_procs / procs_per_row
  
      if (procs_per_row * procs_per_col /= num_procs) then
        if (rank == 0) then
          print *, 'Number of MPI processes must match the process grid dimensions.'
        end if
        call MPI_Abort(comm, 1, ierr)
      end if
  
      ! Allocate counts and displacements for MPI_Scatterv
      allocate(counts(num_procs))
      allocate(displs(num_procs))
      counts = 0
      displs = 0
  
      ! Determine the owner of each non-zero element
      do i = 1, nnz
        block_row = (row_indices(i) - 1) / block_size_row
        block_col = (col_indices(i) - 1) / block_size_col
        owner_rank = mod(block_row, procs_per_row) * procs_per_col + mod(block_col, procs_per_col)
        counts(owner_rank+1) = counts(owner_rank+1) + 1
      end do
  
      displs(1) = 0
      do i = 2, num_procs
        displs(i) = displs(i-1) + counts(i-1)
      end do
  
      ! Prepare data for scattering
      integer, allocatable :: send_row(:), send_col(:)
      real(kind=8), allocatable :: send_data(:)
      allocate(send_row(nnz))
      allocate(send_col(nnz))
      allocate(send_data(nnz))
      counts = 0  ! Reuse counts as indices
  
      do i = 1, nnz
        block_row = (row_indices(i) - 1) / block_size_row
        block_col = (col_indices(i) - 1) / block_size_col
        owner_rank = mod(block_row, procs_per_row) * procs_per_col + mod(block_col, procs_per_col)
        integer :: idx
        idx = displs(owner_rank+1) + counts(owner_rank+1)
        send_row(idx+1) = row_indices(i)
        send_col(idx+1) = col_indices(i)
        send_data(idx+1) = values(i)
        counts(owner_rank+1) = counts(owner_rank+1) + 1
      end do
  
      ! Allocate local arrays
      nnz_local = counts(rank+1)
      allocate(matrix%row_idx(nnz_local))
      allocate(matrix%col_idx(nnz_local))
      allocate(matrix%data(nnz_local))
  
      ! Scatter the data
      call MPI_Scatterv(send_row, counts, displs, MPI_INTEGER, &
                        matrix%row_idx, nnz_local, MPI_INTEGER, 0, comm, ierr)
      call MPI_Scatterv(send_col, counts, displs, MPI_INTEGER, &
                        matrix%col_idx, nnz_local, MPI_INTEGER, 0, comm, ierr)
      call MPI_Scatterv(send_data, counts, displs, MPI_DOUBLE_PRECISION, &
                        matrix%data, nnz_local, MPI_DOUBLE_PRECISION, 0, comm, ierr)
  
      ! Adjust local indices to local block coordinates
      do i = 1, nnz_local
        matrix%row_idx(i) = mod(matrix%row_idx(i)-1, block_size_row) + 1
        matrix%col_idx(i) = mod(matrix%col_idx(i)-1, block_size_col) + 1
      end do
  
      ! Clean up
      if (rank == 0) then
        deallocate(send_row)
        deallocate(send_col)
        deallocate(send_data)
      end if
      deallocate(counts)
      deallocate(displs)
    end subroutine initialize_csb_matrix
  
    ! Parallel matrix-vector multiplication y = A * x
    subroutine spmv_csb(matrix, x_global, y_global, comm)
      type(CSBMatrix), intent(in) :: matrix
      real(kind=8), intent(in) :: x_global(:)
      real(kind=8), intent(out) :: y_global(:)
      integer, intent(in) :: comm
  
      integer :: ierr, rank, num_procs
      integer :: i
      integer :: procs_per_row, procs_per_col
      integer :: local_block_size_row, local_block_size_col
      real(kind=8), allocatable :: x_local(:), y_local(:)
      integer :: nnz_local
      integer :: start_col, end_col
      integer :: start_row, end_row
  
      call MPI_Comm_size(comm, num_procs, ierr)
      call MPI_Comm_rank(comm, rank, ierr)
  
      ! Determine process grid dimensions
      procs_per_row = floor(sqrt(real(num_procs * matrix%n_blocks_row / matrix%n_blocks_col)))
      if (procs_per_row == 0) procs_per_row = 1
      procs_per_col = num_procs / procs_per_row
  
      ! Determine the process's position in the grid
      integer :: proc_row, proc_col
      proc_row = rank / procs_per_col
      proc_col = mod(rank, procs_per_col)
  
      ! Determine local block ranges
      start_row = proc_row * matrix%block_size_row + 1
      end_row = min(start_row + matrix%block_size_row - 1, matrix%n_rows)
      start_col = proc_col * matrix%block_size_col + 1
      end_col = min(start_col + matrix%block_size_col - 1, matrix%n_cols)
      local_block_size_row = end_row - start_row + 1
      local_block_size_col = end_col - start_col + 1
  
      ! Allocate local vectors
      allocate(x_local(local_block_size_col))
      allocate(y_local(local_block_size_row))
      y_local = 0.0_8
  
      ! Scatter relevant portions of x_global to x_local
      call scatter_vector(x_global, x_local, start_col, end_col, comm)
  
      ! Perform local SpMV
      nnz_local = size(matrix%data)
      do i = 1, nnz_local
        y_local(matrix%row_idx(i)) = y_local(matrix%row_idx(i)) + &
                                     matrix%data(i) * x_local(matrix%col_idx(i))
      end do
  
      ! Gather y_local contributions to y_global
      call gather_vector(y_local, y_global, start_row, end_row, comm)
  
      deallocate(x_local)
      deallocate(y_local)
    end subroutine spmv_csb
  
    ! Parallel transpose matrix-vector multiplication y = Aᵗ * x
    subroutine spmv_csb_transpose(matrix, x_global, y_global, comm)
      type(CSBMatrix), intent(in) :: matrix
      real(kind=8), intent(in) :: x_global(:)
      real(kind=8), intent(out) :: y_global(:)
      integer, intent(in) :: comm
  
      integer :: ierr, rank, num_procs
      integer :: i
      integer :: procs_per_row, procs_per_col
      integer :: local_block_size_row, local_block_size_col
      real(kind=8), allocatable :: x_local(:), y_local(:)
      integer :: nnz_local
      integer :: start_col, end_col
      integer :: start_row, end_row
  
      call MPI_Comm_size(comm, num_procs, ierr)
      call MPI_Comm_rank(comm, rank, ierr)
  
      ! Determine process grid dimensions
      procs_per_row = floor(sqrt(real(num_procs * matrix%n_blocks_row / matrix%n_blocks_col)))
      if (procs_per_row == 0) procs_per_row = 1
      procs_per_col = num_procs / procs_per_row
  
      ! Determine the process's position in the grid
      integer :: proc_row, proc_col
      proc_row = rank / procs_per_col
      proc_col = mod(rank, procs_per_col)
  
      ! Determine local block ranges
      start_row = proc_row * matrix%block_size_row + 1
      end_row = min(start_row + matrix%block_size_row - 1, matrix%n_rows)
      start_col = proc_col * matrix%block_size_col + 1
      end_col = min(start_col + matrix%block_size_col - 1, matrix%n_cols)
      local_block_size_row = end_row - start_row + 1
      local_block_size_col = end_col - start_col + 1
  
      ! Allocate local vectors
      allocate(x_local(local_block_size_row))
      allocate(y_local(local_block_size_col))
      y_local = 0.0_8
  
      ! Scatter relevant portions of x_global to x_local
      call scatter_vector(x_global, x_local, start_row, end_row, comm)
  
      ! Perform local SpMV transpose
      nnz_local = size(matrix%data)
      do i = 1, nnz_local
        y_local(matrix%col_idx(i)) = y_local(matrix%col_idx(i)) + &
                                     matrix%data(i) * x_local(matrix%row_idx(i))
      end do
  
      ! Gather y_local contributions to y_global
      call gather_vector(y_local, y_global, start_col, end_col, comm)
  
      deallocate(x_local)
      deallocate(y_local)
    end subroutine spmv_csb_transpose
  
    ! Subroutine to scatter a global vector to local portions
    subroutine scatter_vector(global_vec, local_vec, start_idx, end_idx, comm)
      real(kind=8), intent(in) :: global_vec(:)
      real(kind=8), intent(out) :: local_vec(:)
      integer, intent(in) :: start_idx, end_idx
      integer, intent(in) :: comm
  
      integer :: ierr, rank, num_procs
      integer :: counts, displs
      integer :: local_size
  
      call MPI_Comm_size(comm, num_procs, ierr)
      call MPI_Comm_rank(comm, rank, ierr)
  
      local_size = end_idx - start_idx + 1
      allocate(counts(num_procs))
      allocate(displs(num_procs))
  
      ! Gather counts and displacements from all processes
      counts = 0
      displs = 0
      counts(rank+1) = local_size
      displs(rank+1) = start_idx - 1
  
      ! Allgather counts and displacements
      call MPI_Allgather(counts(rank+1), 1, MPI_INTEGER, counts, 1, MPI_INTEGER, comm, ierr)
      call MPI_Allgather(displs(rank+1), 1, MPI_INTEGER, displs, 1, MPI_INTEGER, comm, ierr)
  
      ! Scatter the vector
      call MPI_Scatterv(global_vec, counts, displs, MPI_DOUBLE_PRECISION, &
                        local_vec, local_size, MPI_DOUBLE_PRECISION, 0, comm, ierr)
  
      deallocate(counts)
      deallocate(displs)
    end subroutine scatter_vector
  
    ! Subroutine to gather local vector contributions to a global vector
    subroutine gather_vector(local_vec, global_vec, start_idx, end_idx, comm)
      real(kind=8), intent(in) :: local_vec(:)
      real(kind=8), intent(inout) :: global_vec(:)
      integer, intent(in) :: start_idx, end_idx
      integer, intent(in) :: comm
  
      integer :: ierr, rank, num_procs
      integer :: counts, displs
      integer :: local_size
  
      call MPI_Comm_size(comm, num_procs, ierr)
      call MPI_Comm_rank(comm, rank, ierr)
  
      local_size = end_idx - start_idx + 1
      allocate(counts(num_procs))
      allocate(displs(num_procs))
  
      ! Gather counts and displacements from all processes
      counts = 0
      displs = 0
      counts(rank+1) = local_size
      displs(rank+1) = start_idx - 1
  
      ! Allgather counts and displacements
      call MPI_Allgather(counts(rank+1), 1, MPI_INTEGER, counts, 1, MPI_INTEGER, comm, ierr)
      call MPI_Allgather(displs(rank+1), 1, MPI_INTEGER, displs, 1, MPI_INTEGER, comm, ierr)
  
      ! Initialize global vector to zero
      if (rank == 0) then
        global_vec = 0.0_8
      end if
  
      ! Reduce the local vectors into the global vector
      call MPI_Reduce_scatter_block(local_vec, global_vec, local_size, MPI_DOUBLE_PRECISION, &
                                    MPI_SUM, comm, ierr)
  
      deallocate(counts)
      deallocate(displs)
    end subroutine gather_vector
  
  end module SparseMatrixCSB
  