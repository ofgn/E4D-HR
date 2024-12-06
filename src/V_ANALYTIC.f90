! Changelog 
! 7/4/23 - OFGN
! Use phase rather than phase/resistance
! Max phase set to the current measurement rather than the running total
! Min phase set to the current measurement rather than the running total
! 1/10/23 - OFGN
! Increase witdth in format descriptor.

module v_analytic

    use vars
    use input
    use master
    use output

    implicit none
    real :: se,s_ave
    integer :: resid_flag,n_pot
    character*80 :: resid_prefix
    integer, dimension(:), allocatable :: poti
    real, dimension(:), allocatable :: an_pot
    contains

    !______________________________________________________________________________
    subroutine compute_analytic
        implicit none
        logical :: stat
        integer :: i
        
        n_pot=0
        call get_out_opts(stat)
        if(.not.stat) return
        
        !call write_an_header

        if(resid_flag .eq.1 .or. n_pot>1) then

            call read_nodes(stat)
            if(.not.stat) call crash_exit
        
            call read_survey
            call translate_electrodes
            call read_conductivity

            call get_level
            
            if(resid_flag==1) then
            call build_an_dpred
            if(resid_flag==1) then
                call write_an_dpred
                call write_an_srv
            end if
            end if

            if(n_pot>0) then
            do i=1,n_pot
                if(poti(i) .le. nm) then
                    call build_pot_an(poti(i))
                    call write_pot_an(poti(i))
                end if
            end do
            end if
        end if

        if(allocated(poti)) deallocate(poti)
        if(allocated(an_pot)) deallocate(an_pot)

        open(51,file='e4d.log',status='old',action='write',position="append")
        write(51,*)
        write(51,*) "Computation complete."
        close(51)
        return
    end subroutine compute_analytic
    !______________________________________________________________________________


    ! @brief Subroutine to calculate the weighted mean model for apparent conductivity and phase
    ! @details Based on E4D V_ANALYTIC.f90 subroutine get_ave_sig()
    !   ! Uses phase rather than phase/resistance and fixed the calculation of the running total
    ! @author ofgn
    subroutine calculate_mean_model()
        implicit none

        integer :: i, j                                     ! Loop index
        real(8) :: pi = 4 * atan(1.0d0)                  ! Pi constant
        real(8) :: gfam, gfan, gfbm, gfbn                ! Terms for geometric factor calculation
        real(8) :: geometric_factor                      ! Geometric factor
        real(8) :: value_sum, weight_sum                 ! Variables for calculating the mean
        real(8) :: min_app_sigma                         ! Minimum apparent conductivity
        real(8) :: max_app_sigma                         ! Maximum apparent conductivity
        real(8) :: min_app_phase                         ! Minimum apparent phase
        real(8) :: max_app_phase                         ! Maximum apparent phase
        real(8), allocatable :: app_sigma(:)             ! Array of apparent conductivities
        real(8), allocatable :: app_phase(:)             ! Array of apparent phases
        logical :: status                                ! Status flag for reading data

        call read_nodes(status)

        if (.not. status) then
            call crash_exit()
        end if

        call read_elements(status)
        if (.not. status) then
            call crash_exit()
        end if

        if (.not. allocated(sigma_re)) then
            allocate(sigma_re(n_elements))
            model_size = n_elements
        end if

        if (i_flag .and. .not. allocated(sigma_im)) then
            allocate(sigma_im(n_elements))
        end if

        call get_level()

        allocate(app_sigma(nm))

        ! Calculate the geometric factor and apparent conductivity for each measurement
        ! Update the minimum and maximum apparent conductivities
        do i = 1, nm
            call get_gf(gfam, s_conf(i, 1), s_conf(i, 3))
            call get_gf(gfan, s_conf(i, 1), s_conf(i, 4))
            call get_gf(gfbm, s_conf(i, 2), s_conf(i, 3))
            call get_gf(gfbn, s_conf(i, 2), s_conf(i, 4))
            geometric_factor = (gfam - gfan - (gfbm - gfbn)) / (4.0d0 * pi)
            app_sigma(i) = geometric_factor / dobs(i)
        end do

        ! Calculate the weighted mean apparent conductivity
        min_app_sigma = minval(pack(app_sigma, app_sigma > 0.0d0))
        max_app_sigma = maxval(pack(app_sigma, app_sigma > 0.0d0))
        sigma_0 = sum(pack(app_sigma / Wd, app_sigma > 0.0d0)) / sum(pack(1.0d0 / Wd, app_sigma > 0.0d0))
        sigma_re = sigma_0
        deallocate(app_sigma)

        ! For complex conductivity, update the minimum and maximum apparent phases
        ! Calculate the weighted mean apparent phase
        ! Calculate the corresponding imaginary conductivity values
        if (i_flag) then
            allocate(app_phase(nm))
            app_phase = -dobsi
            min_app_phase = minval(pack(app_phase, app_phase > 0.0d0))
            max_app_phase = maxval(pack(app_phase, app_phase > 0.0d0))
            phase_0 = sum(pack(app_phase / Wdi, app_phase > 0.0d0)) / sum(pack(1.0d0 / Wdi, app_phase > 0.0d0))
            sigma_im = tan(phase_0) * sigma_re
            deallocate(app_phase)
        end if


        call export_real_conductivity_model(model_name="mean")

        open(51, file='e4d.log', status='old', action='write', position='append')
        write(51, "(A36, ES12.5)") "  Minimum apparent resistivity:     ", 1.0d0 / max_app_sigma
        write(51, "(A36, ES12.5)") "  Maximum apparent resistivity:     ", 1.0d0 / min_app_sigma
        write(51, "(A36, ES12.5)") "  Mean apparent resistivity:        ", 1.0d0 / sigma_0
        write(51, "(A36, ES12.5)") "  Minimum apparent conductivity:    ", min_app_sigma
        write(51, "(A36, ES12.5)") "  Maximum apparent conductivity:    ", max_app_sigma
        write(51, "(A36, ES12.5)") "  Mean apparent conductivity:       ", sigma_0
        close(51)

        if (i_flag) then
            open(51, file='e4d.log', status='old', action='write', position='append')
            write(51, "(A36, ES12.5)") "  Minimum apparent phase:           ", min_app_phase
            write(51, "(A36, ES12.5)") "  Maximum apparent phase:           ", max_app_phase
            write(51, "(A36, ES12.5)") "  Mean apparent phase:              ", phase_0
            close(51)
        end if
    end subroutine calculate_mean_model

    ! @brief Subroutine to calculate the median model for apparent conductivity and phase.
    ! @author ofgn
    subroutine calculate_median_model()
        implicit none

        integer :: i                                     ! Loop index
        real(8) :: pi = 4 * atan(1.0d0)                  ! Pi constant
        real(8) :: gfam, gfan, gfbm, gfbn                ! Terms for geometric factor calculation
        real(8) :: geometric_factor                      ! Geometric factor
        real(8) :: min_app_sigma                         ! Minimum apparent conductivity
        real(8) :: max_app_sigma                         ! Maximum apparent conductivity
        real(8) :: min_app_phase                         ! Minimum apparent phase
        real(8) :: max_app_phase                         ! Maximum apparent phase
        real(8), allocatable :: app_sigma(:)             ! Array of apparent conductivities
        real(8), allocatable :: app_phase(:)             ! Array of apparent phases
        logical :: status                                ! Status flag for reading data

        call read_nodes(status)
        if (.not. status) then
            call crash_exit()
        end if

        call read_elements(status)
        if (.not. status) then
            call crash_exit()
        end if

        if (.not. allocated(sigma_re)) then
            allocate(sigma_re(n_elements))
            model_size = n_elements
        end if

        if (i_flag .and. .not. allocated(sigma_im)) then
            allocate(sigma_im(n_elements))
        end if

        call get_level()

        min_app_sigma = huge(0.0d0)
        max_app_sigma = 0.0d0

        allocate(app_sigma(nm))

        ! Calculate the geometric factor and apparent conductivity for each measurement
        ! Update the minimum and maximum apparent conductivities
        do i = 1, nm
            call get_gf(gfam, s_conf(i, 1), s_conf(i, 3))
            call get_gf(gfan, s_conf(i, 1), s_conf(i, 4))
            call get_gf(gfbm, s_conf(i, 2), s_conf(i, 3))
            call get_gf(gfbn, s_conf(i, 2), s_conf(i, 4))
            geometric_factor = (gfam - gfan - (gfbm - gfbn)) / (4.0d0 * pi)
            app_sigma(i) = geometric_factor / dobs(i)
            if ((app_sigma(i) > 0.0d0) .and. (app_sigma(i) > max_app_sigma)) then
                max_app_sigma = app_sigma(i)
            else if ((app_sigma(i) > 0.0d0) .and. (app_sigma(i) < min_app_sigma)) then
                min_app_sigma = app_sigma(i)
            end if
        end do

        sigma_0 = median(app_sigma)
        sigma_re = sigma_0

        ! For complex conductivity, update the minimum and maximum apparent phases
        ! Calculate the median apparent phase
        ! Calculate the corresponding imaginary conductivity values
        if (i_flag) then
            app_phase = pack((-dobsi), (-dobsi) > 0.0d0)
            min_app_phase = minval(app_phase)
            max_app_phase = maxval(app_phase)
            phase_0 = median(app_phase)
            sigma_im = tan(phase_0) * sigma_re
            deallocate(app_phase)
        end if

        call export_real_conductivity_model(model_name="median")

        open(51, file='e4d.log', status='old', action='write', position='append')
        write(51, "(A36, ES12.5)") "  Minimum apparent resistivity:     ", 1.0d0 / max_app_sigma
        write(51, "(A36, ES12.5)") "  Maximum apparent resistivity:     ", 1.0d0 / min_app_sigma
        write(51, "(A36, ES12.5)") "  Median apparent resistivity:      ", 1.0d0 / sigma_0
        write(51, "(A36, ES12.5)") "  Minimum apparent conductivity:    ", min_app_sigma
        write(51, "(A36, ES12.5)") "  Maximum apparent conductivity:    ", max_app_sigma
        write(51, "(A36, ES12.5)") "  Median apparent conductivity:     ", sigma_0
        close(51)

        if (i_flag) then
            open(51, file='e4d.log', status='old', action='write', position='append')
            write(51, "(A36, ES12.5)") "  Minimum apparent phase:           ", min_app_phase
            write(51, "(A36, ES12.5)") "  Maximum apparent phase:           ", max_app_phase
            write(51, "(A36, ES12.5)") "  Median apparent phase:            ", phase_0
            close(51)
        end if
    end subroutine calculate_median_model


    ! ---------------------------------------------------------------------------------------------------
    ! @brief Function to calculate the median of an array using quickselect algorithm
    ! @param[in] arr The input array for which the median is to be calculated
    ! @return The median value of the array
    ! @author ofgn
    ! ---------------------------------------------------------------------------------------------------
    function median(arr)
        implicit none
        real(8), intent(in) :: arr(:)                                           ! Input array
        real(8) :: median                                                       ! Median value
        integer :: n                                                            ! Size of the input array
        real(8), allocatable :: arr_copy(:)                                     ! Local copy of the array for manipulation

        n = size(arr)
        allocate(arr_copy(n))
        arr_copy = arr

        if (mod(n, 2) == 1) then
            median = quickselect(arr_copy, 1, n, (n+1)/2)
        else
            median = 0.5 * (quickselect(arr_copy, 1, n, n/2) + quickselect(arr_copy, 1, n, n/2+1))
        end if
    end function median

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Quickselect algorithm to find the k-th smallest element
    ! @param[inout] arr The array to be processed
    ! @param[in] left Left index of the subarray
    ! @param[in] right Right index of the subarray
    ! @param[in] k The position of the k-th smallest element to find
    ! @return The k-th smallest element in the array
    ! @author ofgn
    ! ---------------------------------------------------------------------------------------------------
    recursive function quickselect(arr, left, right, k) result(kth_smallest)
        implicit none
        real(8), intent(inout) :: arr(:)                                        ! Input array
        integer, intent(in) :: left, right, k                                   ! Left, right, and k indices
        real(8) :: kth_smallest                                                 ! The k-th smallest element
        integer :: pivot_index                                                  ! Index of the pivot element

        if (left == right) then
            kth_smallest = arr(left)
        else
            ! Partition the array and find the pivot index
            pivot_index = partition(arr, left, right)
            
            ! Recur based on pivot position relative to k
            if (k == pivot_index) then
                kth_smallest = arr(k)
            else if (k < pivot_index) then
                kth_smallest = quickselect(arr, left, pivot_index - 1, k)
            else
                kth_smallest = quickselect(arr, pivot_index + 1, right, k)
            end if
        end if
    end function quickselect

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Partition function for quickselect
    ! @param[inout] arr The array to be partitioned
    ! @param[in] left Left index of the subarray
    ! @param[in] right Right index of the subarray
    ! @return The index of the pivot element after partitioning
    ! @author ofgn
    ! ---------------------------------------------------------------------------------------------------
    function partition(arr, left, right) result(pivot_index)
        implicit none
        real(8), intent(inout) :: arr(:)                                        ! Input array
        integer, intent(in) :: left, right                                      ! Left and right indices
        integer :: pivot_index, i, j                                            ! Pivot and loop variables
        real(8) :: pivot_value, temp                                            ! Pivot value

        pivot_value = arr(right)
        pivot_index = left
        do i = left, right - 1
            if (arr(i) < pivot_value) then
                ! Swap arr(i) and arr(pivot_index)
                temp = arr(i)
                arr(i) = arr(pivot_index)
                arr(pivot_index) = temp
                pivot_index = pivot_index + 1
            end if
        end do

        ! Swap arr(pivot_index) and arr(right)
        temp = arr(pivot_index)
        arr(pivot_index) = arr(right)
        arr(right) = temp
    end function partition
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine write_pot_an(id)
        implicit none
        integer :: id
        integer :: i
        character*80 :: fstr
        
        write(fstr,"(A,I0)") 'an_potential.',id
        
        open(51,file='e4d.log',status='old',action='write',position="append")
        write(51,*) "  Writing analytic potential field ",trim(fstr)
        close(51)

        open(10,file=trim(fstr),status='replace',action='write')
        write(10,*) nnodes,1
        do i=1,nnodes
            write(10,*) an_pot(i)
        end do
        close(10)

    end subroutine write_pot_an
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine build_pot_an(id)
        implicit none
        integer :: id
        integer :: i,ei
        real :: ixp,iyp,izp,ixn,iyn,izn
        real, dimension(nnodes) :: r
        real, parameter :: pi = 3.14159265359
        
        if(.not. allocated(an_pot)) allocate(an_pot(nnodes))

        if(s_conf(id,1)>0) then
            ixp=e_pos(s_conf(id,1),1)
            iyp=e_pos(s_conf(id,1),2)
            izp=e_pos(s_conf(id,1),3)
        end if
        if(s_conf(id,2)>0) then
            ixn=e_pos(s_conf(id,2),1)
            iyn=e_pos(s_conf(id,2),2)
            izn=e_pos(s_conf(id,2),3)
        end if
        an_pot=0
        if(s_conf(id,1)>0) then
            r=sqrt( (ixp-nodes(:,1))**2 + (iyp-nodes(:,2))**2 + (izp-nodes(:,3))**2)
            r=r+1e-15
            an_pot = an_pot + 1/(4*pi*r*s_ave);

            izp=2*se-izp
            r=sqrt( (ixp-nodes(:,1))**2 + (iyp-nodes(:,2))**2 + (izp-nodes(:,3))**2)
            r=r+1e-15
            an_pot = an_pot + 1/(4*pi*r*s_ave);
        end if

        if(s_conf(id,2)>0) then
            r=sqrt( (ixn-nodes(:,1))**2 + (iyn-nodes(:,2))**2 + (izn-nodes(:,3))**2)
            r=r+1e-15


            an_pot = an_pot - 1/(4*pi*r*s_ave);

            izn=2*se-izn
            r=sqrt( (ixn-nodes(:,1))**2 + (iyn-nodes(:,2))**2 + (izn-nodes(:,3))**2)
            r=r+1e-15
            an_pot = an_pot - 1/(4*pi*r*s_ave);
        end if

    end subroutine build_pot_an
    !______________________________________________________________________________
    !______________________________________________________________________________
    subroutine write_an_dpred
        implicit none
        integer :: i
        open(51,file='e4d.log',status='old',action='write',position='append')
        write(51,*) 
        write(51,*) "  Writing analytic measurements to ",trim(resid_prefix)
        write(*,*) "  Writing analytic measurements to ",trim(resid_prefix)
        close(51)

        open(10,file=trim(resid_prefix),status='replace',action='write')
        write(10,*) nm
        do i=1,nm
            write(10,*) i,s_conf(i,1),s_conf(i,2),s_conf(i,3),s_conf(i,4),dobs(i),dpred(i)
        end do
        close(10)
    end subroutine write_an_dpred
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine write_an_srv
        implicit none
        integer :: i
        open(51,file='e4d.log',status='old',action='write',position='append')
        write(51,*) 
        write(51,*) "  Writing analytic survey to ",trim(sig_filename),".srv"
        write(*,*) "  Writing analytic survey to ",trim(sig_filename),".srv"
        
        close(51)

        open(10,file=trim(sig_filename)//".srv",status='replace',action='write')
        write(10,*) ne
        do i=1,ne
            write(10,"(I7,4F10.3,I5)") i,e_pos(i,1)+xorig,e_pos(i,2)+yorig,e_pos(i,3)+zorig,e_pos(i,4)
        end do
        write(10,*)
        write(10,*) nm
        do i=1,nm
            write(10,"(5I7,2g15.6)") i,s_conf(i,1),s_conf(i,2),s_conf(i,3),s_conf(i,4),dpred(i),0.05*abs(dpred(i))
        end do
        close(10)
    end subroutine write_an_srv
    !______________________________________________________________________________


    !______________________________________________________________________________
    subroutine build_an_dpred
        implicit none
        real :: gfam,gfan,gfbm,gfbn,gft
        integer :: i
        real, parameter :: pi = 3.14159265359
        
        
        if(allocated(dpred)) deallocate(dpred)
        allocate(dpred(nm))
        
        do i=1,nm
            call get_gf(gfam,s_conf(i,1),s_conf(i,3))
            call get_gf(gfan,s_conf(i,1),s_conf(i,4))
            call get_gf(gfbm,s_conf(i,2),s_conf(i,3))
            call get_gf(gfbn,s_conf(i,2),s_conf(i,4))

            gft = (gfam- gfan - (gfbm - gfbn))/(4*pi)
            dpred(i) = gft/s_ave
        end do
        
        
    end subroutine build_an_dpred
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine get_gf(gf,ei,ep)
        implicit none
        real :: gf
        integer :: ei,ep
        real :: xi,yi,zi,zii,xp,yp,zp,r,ri
        
        gf=0;
        if(ei==0 .or. ep==0) return

        xi=e_pos(ei,1); yi=e_pos(ei,2); zi=e_pos(ei,3); zii=2*se-e_pos(ei,3)
        
        xp=e_pos(ep,1); yp=e_pos(ep,2); zp=e_pos(ep,3)
        
        r=sqrt(  (xi-xp)**2 + (yi-yp)**2 + (zi-zp)**2);
        ri=sqrt( (xi-xp)**2 + (yi-yp)**2 + (zii-zp)**2);
        
        gf = 1/r + 1/ri;
        
    end subroutine get_gf
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine get_level
        implicit none
        integer :: i,cc
        logical :: flat = .true.
        real :: mx,mn
        cc = 0
        se = 0

        mx=nodes(1,3);
        mn=nodes(1,3);
        
        do i=2,nnodes
            if(nbounds(i)==1) then
            if(nodes(i,3)>mx) mx=nodes(i,3)
            if(nodes(i,3)<mn) mn=nodes(i,3)
            end if
        end do
        !se=.5*(mx+mn)
        se = mx
        
        open(51,file='e4d.log',status='old',action='write',position='append');
        write(51,"(A,I10.10)") "  Number of nodes:                  ",nnodes
        write(51,*) 
        write(51,"(A36, G12.5)") " Maximum surface node elevation:   ",mx+zorig ! Increase witdth in format descriptor. - OFGN 1/10/23
        write(51,"(A36, G12.5)") " Minimum surface node elevation:   ",mn+zorig ! Increase witdth in format descriptor. - OFGN 1/10/23
        write(51,"(A36, G12.5)") " Using surface elevation of    :   ",se+zorig ! Increase witdth in format descriptor. - OFGN 1/10/23
        if(mx .ne. se .or. mn .ne. se) then
            write(51,*) "  !!! WARNING: It appears there may be some surface variability."
            write(51,*) "  !!! The analytic solution requires a flat surface."
        end if
        close(51)
        
        if(use_mean .or. use_median) return

        mx=sigma_re(1);
        mn=sigma_re(1);
        
        do i=2,n_elements
            if(sigma_re(i)>mx) mx=sigma_re(i)
            if(sigma_re(i)<mn) mn=sigma_re(i)
        end do
        s_ave=.5*(mx+mn)

        open(51,file='e4d.log',status='old',action='write',position='append');
        write(51,*) 
        write(51,*) " Maximum conductivity:             ",mx
        write(51,*) " Minimum conductivity:             ",mn
        write(51,*) " Using conductivity  :             ",s_ave
        if(mx .ne. s_ave .or. mn .ne. s_ave) then
            write(51,*) "  !!!WARNING: it appears there may be some conductive variability"
            write(51,*) "  !!!The analytic solution requires homogeneous conductivity"
        end if
        close(51)

    end subroutine get_level      
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine write_an_header
        implicit none
            
        open(51,file='e4d.log',status='old',action='write')
        write(51,*) "RUNNING IN ANALYTIC FORWARD SOLUTION MODE"
        write(51,*) "MESH FILE = ",cfg_filename
        write(51,*) "SURVEY FILE = ",efile
        write(51,*) "CONDUCTIVITY FILE = ",sig_filename
        write(51,*) "OUTPUT OPTIONS FILE = ",out_file
        close(51);
    end subroutine write_an_header
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine get_out_opts(chck)
        implicit none
        logical :: chck
        logical :: exst
        integer :: ist,i

        chck=.true.
        inquire(file=trim(out_file),exist=exst)
        if(.not. exst) goto 10

        open(10,file=trim(out_file),status='old',action='read')
        read(10,*,IOSTAT=ist) resid_flag; if(ist.ne.0) goto 11
        read(10,*,IOSTAT=ist) resid_prefix; if(ist.ne.0) goto 12
        read(10,*,IOSTAT=ist) n_pot;   if(ist.ne.0) goto 13
        if(n_pot>0) then
            allocate(poti(n_pot))
            do i=1,n_pot
            read(10,*,IOSTAT=ist) poti(i); if(ist.ne.0) goto 14
            end do
        end if
        close(10)
        
        return

    10    continue
        chck = .false.
        open(51,file='e4d.log',status='old',action='write',position='append')
        write(51,*) ' Cannot find the output options file: ',trim(out_file)
        write(51,*) ' Aborting ...'
        close(51)
        write(*,*) ' Cannot find the output options file: ',trim(out_file)
        write(*,*) ' Aborting ...'
        return 
        
    11    continue
        chck = .false.
        open(51,file='e4d.log',status='old',action='write',position='append')
        write(51,*) ' The was a problem reading the first line in the output file: ',trim(out_file)
        write(51,*) ' Aborting ...'
        close(51)
        write(*,*) ' The was a problem reading the first line in the output file: ',trim(out_file)
        write(*,*) ' Aborting ...'
        return

    12    continue
        chck = .false.
        open(51,file='e4d.log',status='old',action='write',position='append')
        write(51,*) ' The was a problem reading the predicted data file in: ',trim(out_file)
        write(51,*) ' Aborting ...'
        close(51)
        write(*,*) ' The was a problem reading the predicted data file in: ',trim(out_file)
        write(*,*) ' Aborting'
        return

    13    continue
        chck = .false.
        open(51,file='e4d.log',status='old',action='write',position='append')
        write(51,*) ' There was a problem reading the number of potential fields to write in: ',trim(out_file)
        write(51,*) ' Aborting ...'
        close(51)
        write(*,*) ' The was a problem reading the number of potential fields to write in: ',trim(out_file)
        write(*,*) ' Aborting ...'
        return

    14    continue
        chck = .false.
        open(51,file='e4d.log',status='old',action='write',position='append')
        write(51,*) ' The was a problem reading potential index: ',i,' in: ',trim(out_file)
        write(51,*) ' Aborting ...'
        close(51)
        write(*,*) ' The was a problem reading potential index: ',i,' in: ',trim(out_file)
        write(*,*) ' Aborting ...'

        return


    end subroutine get_out_opts
    !______________________________________________________________________________

end module v_analytic
