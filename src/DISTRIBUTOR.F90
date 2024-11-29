module distributor
  !!This module divides the processes and sets corresponding variables
  !!depending on command line input
  !!Yilin Fang and Tim Johnson, 12/18/2015
    use vars
    implicit none
    integer :: nargs, iarg, stat
    integer :: my_color, my_key, my_comm, my_group, my_grank, my_commsize
    integer :: my_wrank, n_wrank
    character*4 :: arg, tstr
    logical :: petsc_initialized

contains

    subroutine distribute
        implicit none

        call MPI_COMM_DUP(mpi_comm_world, E4D_COMM, ierr)
        call MPI_COMM_GROUP(E4D_COMM, my_group, ierr)
        call MPI_COMM_RANK(E4D_COMM, my_grank, ierr)
        call MPI_COMM_SIZE(E4D_COMM, my_commsize, ierr)
    end subroutine distribute

    subroutine dist_abort()
        implicit none

        call PetscInitialized(petsc_initialized, ierr)
        if (petsc_initialized) then
            call PetscFinalize(ierr)
        end if
        stop
    end subroutine dist_abort
end module distributor
