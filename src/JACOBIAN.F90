module jacob
  
  use vars
  
  implicit none
  
  
contains 

  !_______________________________________________________________________________
  subroutine build_jaco
    implicit none
    
    integer :: di,a,b,m,n,i,ii,jj,emin,emax,ra,rb,rm,rn,mown,my_jrow,irow,count
    integer :: a_old,b_old,m_old,n_old,a_olds,b_olds,m_olds,n_olds
    integer :: tag, ierr , status(MPI_STATUS_SIZE)
    integer, dimension(nm) :: sreqa, sreqb,sreqm,sreqn
    integer :: sca,scb,scm,scn
    integer :: rreqa, rreqb,rreqm, rreqn,aflg,bflg,mflg,nflg
    real, dimension(nnodes) :: pa,pb,pm,pn
    logical :: ta,tb,tm,tn
    integer :: Cstart, Rstart, Cend, Rend,scount,rcount,req,cs
    logical :: fd = .false.
    logical :: rflag
    real :: tt1,tt2,ett
    real :: rtarg = 0.05 
    tag=0
   
    call cpu_time(tt1)

    n_data_assigned = data_assignments(my_rank,2) - data_assignments(my_rank,1) + 1
    if(.not. allocated(Jaco)) then
       allocate(Jaco(n_data_assigned, n_elements))
    elseif(res_flag .and. allocated(Jaco)) then
       deallocate(Jaco) 
       allocate(Jaco(n_data_assigned,n_elements))
    end if



    pa=0; pb=0; pm=0; pn=0
  
    a_olds=-1
    b_olds=-1
    m_olds=-1
    n_olds=-1

    !!check for pole pole surveys and send poles to save comp. time if appropriate
    do di=1,4

       if(pcheck(di)) then
      
          !determine who owns the pole
          if(s_conf(1,di)==0) goto 3
          do i=1,n_rank-1 
             emin=eind(i,1); emax=eind(i,2)
             if((emin .le. s_conf(1,di)) .and. (emax .ge. s_conf(1,di))) then
                mown = i
                exit
             end if
          end do
          
          if(mown==my_rank) then

             select case(di)
             case(1)
                !if(invi) then
                !   pa=polesi(:,s_conf(1,di)-eind(my_rank,1)+1)
                !else
                   pa=poles(:,s_conf(1,di)-eind(my_rank,1)+1)
                !end if
                call MPI_BCAST(pa,nnodes,MPI_DOUBLE,mown-1,SCOMM,ierr)
             case(2)
                !if(invi) then
                !   pb=polesi(:,s_conf(1,di)-eind(my_rank,1)+1)
                !else
                   pb=poles(:,s_conf(1,di)-eind(my_rank,1)+1)
                !end if
                call MPI_BCAST(pb,nnodes,MPI_DOUBLE,mown-1,SCOMM,ierr)
             case(3)
                pm=poles(:,s_conf(1,di)-eind(my_rank,1)+1)
                call MPI_BCAST(pm,nnodes,MPI_DOUBLE,mown-1,SCOMM,ierr)
             case(4)
                pn=poles(:,s_conf(1,di)-eind(my_rank,1)+1)
                call MPI_BCAST(pn,nnodes,MPI_DOUBLE,mown-1,SCOMM,ierr)
             end select
       
          else

             select case(di)
             case(1)
                call MPI_BCAST(pa,nnodes,MPI_DOUBLE,mown-1,SCOMM,ierr)
             case(2)
                call MPI_BCAST(pb,nnodes,MPI_DOUBLE,mown-1,SCOMM,ierr)
             case(3)
                call MPI_BCAST(pm,nnodes,MPI_DOUBLE,mown-1,SCOMM,ierr)
             case(4)
                call MPI_BCAST(pn,nnodes,MPI_DOUBLE,mown-1,SCOMM,ierr)
             end select
                
          end if
        
          select case(di)
          case(1)
             a_olds=s_conf(1,di)
          case(2)
             b_olds=s_conf(1,di)
          case(3)
             m_olds=s_conf(1,di)
          case(4)
             n_olds=s_conf(1,di)
          end select
       end if
3      continue
    end do

 
    !!do the non-blocking sends first
    scount=0
    fd=.false.
    a_old=a_olds
    b_old=b_olds
    m_old=m_olds
    n_old=n_olds

   
    call MPI_BARRIER(SCOMM,ierr)
    call cpu_time(tt2)

    
    do ii=1,nm
   
       di=rr_seq(ii)
       a = s_conf(di,1)
       b = s_conf(di,2)
       m = s_conf(di,3)
       n = s_conf(di,4)
       
       !!determine which rank has each pole in this measurement
       ra=0; rb=0; rm=0; rn=0
       do i=1,n_rank-1
          emin=eind(i,1); emax=eind(i,2)
          if((emin .le. a) .and. (emax .ge. a)) ra = i	
          if((emin .le. b) .and. (emax .ge. b)) rb = i 
          if((emin .le. m) .and. (emax .ge. m)) rm = i
          if((emin .le. n) .and. (emax .ge. n)) rn = i
          if((data_assignments(i,1) .le. di) .and. (data_assignments(i,2) .ge. di)) mown = i
       end do
       
       if(mown==my_rank) then          
    
          !!I own this measurement, I'll either recieve poles or use my own
          !if(a .ne. a_old) then 
          if(ra==my_rank) then

             !if(invi) then
             !   pa=polesi(:,a-eind(my_rank,1)+1)
             !else
                pa=poles(:,a-eind(my_rank,1)+1)
             !end if

          else

             if(ra .ne. 0) then
                call MPI_RECV(pa, nnodes, MPI_DOUBLE, ra,a, E4D_COMM, status, ierr)
             end if

          end if
          !end if
         
          !if(b .ne. b_old ) then
          if(rb==my_rank) then

             !if(invi) then
             !   pb=polesi(:,b-eind(my_rank,1)+1)
             !else
                pb=poles(:,b-eind(my_rank,1)+1)
             !end if

          else

             if(rb .ne. 0) then
                call MPI_RECV(pb, nnodes, MPI_DOUBLE, rb,b, E4D_COMM, status, ierr) 
             end if
          end if
          !end if
 

          !if(m .ne. m_old) then
          if(rm==my_rank) then
             pm=poles(:,m-eind(my_rank,1)+1)
          else
             if(rm .ne. 0) then
                call MPI_RECV(pm, nnodes, MPI_DOUBLE, rm,m, E4D_COMM, status, ierr)
             end if
          end if
          !end if

          !if(n .ne. n_old) then
          if(rn==my_rank) then
             pn = poles(:,n-eind(my_rank,1)+1)
          else 
             if(rn .ne. 0) then
                call MPI_RECV(pn, nnodes, MPI_DOUBLE, rn,n, E4D_COMM, status, ierr)
             end if
          end if
          !end if
         
          my_jrow = di - data_assignments(my_rank,1) + 1
     
          call build_my_jaco(my_jrow,pa-pb,pm-pn)
          call MPI_SEND(my_rank,1,MPI_INTEGER,0,1,E4D_COMM,ierr)
          
       else
          
          !!I don't own this measurement, so I'll send poles if they belong to me
          !if(a .ne. a_old) then
          if(ra==my_rank) then

             !if(invi) then
             !   pa = polesi(:,a-eind(my_rank,1)+1)
             !else
                pa = poles(:,a-eind(my_rank,1)+1)
             !end if
             call MPI_SEND(pa,nnodes,MPI_DOUBLE,mown,a,E4D_COMM, ierr)
             
          end if
          !end if
          
          !if(b .ne. b_old) then
          if(rb==my_rank) then
             !if(invi) then
             !   pb = polesi(:,b-eind(my_rank,1)+1)
             !else
                pb = poles(:,b-eind(my_rank,1)+1)
             !end if
             call MPI_SEND(pb,nnodes,MPI_DOUBLE,mown,b,E4D_COMM, ierr)
          end if
          !end if

          !if(m .ne. m_old) then
          if(rm==my_rank) then
             pm = poles(:,m-eind(my_rank,1)+1)
             call MPI_SEND(pm,nnodes,MPI_DOUBLE,mown,m,E4D_COMM, ierr)
          end if
          !end if

          !if(n .ne. n_old) then
          if(rn==my_rank) then
             pn = poles(:,n-eind(my_rank,1)+1)
             call MPI_SEND(pn,nnodes,MPI_DOUBLE,mown,n,E4D_COMM, ierr)
          end if
          !end if
       end if
 
       a_old=a
       b_old=b
       m_old=m
       n_old=n

    end do

    

    !!turn elements off if necessary
    if(allocated(J_on_off)) then
       do i=1,n_elements
          if(.not. J_on_off(i)) then
             Jaco(:,i)=0
          end if;
       end do
    end if

  end subroutine build_jaco
  !_______________________________________________________________________________

 
  !_______________________________________________________________________________
  subroutine build_my_jaco(irow,phi_s,phi_r)
    implicit none
    integer :: irow
    real, dimension(nnodes) ::phi_s, phi_r
    integer :: mcount,i,j,k,row,col,ind
    real :: C2
    integer :: js
    
    mcount=0  
    do j=1,n_elements

       C2 = 0   
       do k=1,10
          mcount=mcount+1
          ind=A_map(mcount)
          row=rows(ind)
          col=cols(ind)
          if(row .ne. col) then
             !!acount for upper and lower triangle
             C2=C2+real(delA(mcount))*(phi_s(col)*phi_r(row) + phi_s(row)*phi_r(col))
          else
             C2=C2+real(delA(mcount))*phi_s(col)*phi_r(row)
          end if
       end do
       
       !!Use this for estimating log(delta_sigma)
       if(invi) then
          Jaco(irow,j) = -sigma_im(j)*C2
       else
          Jaco(irow,j) = -sigma_re(j)*C2
       end if



    end do

    
   
  end subroutine build_my_jaco
  !_______________________________________________________________________________



end module jacob
