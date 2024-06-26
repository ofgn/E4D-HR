module assemble

!!____________________________________________________________________________________________
!!  E4D
!!  Copyright N) 2014, Battelle Memorial Institute
!!  All rights reserved.
!!
!!  1 . Battelle Memorial Institute (hereinafter Battelle) hereby grants permission to any
!!  person or entity lawfully obtaining a copy of this software and associated documentation
!!  files (hereinafter $B!H(Bthe Software$B!I(B) to redistribute and use the Software in source and
!!  binary forms, with or without modification.  Such person or entity may use, copy, modify,
!!  merge, publish, distribute, sublicense, and/or sell copies of the Software, and may permit
!!  others to do so, subject to the following conditions:
!!
!!   -Redistributions of source code must retain the above copyright notice, this list of
!!    conditions and the following disclaimers.
!!   -Redistributions in binary form must reproduce the above copyright notice, this list
!!    of conditions and the following disclaimer in the documentation and/or other materials
!!    provided with the distribution.
!!   -Other than as used herein, neither the name Battelle Memorial Institute or Battelle may
!!    be used in any form whatsoever without the express written consent of Battelle.
!!   -Redistributions of the software in any form, and publications based on work performed
!!    using the software should include the following citation as a reference:
!!    [Cite published manuscript.  If this does not apply, delete this bulleted item].
!!
!!  2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
!!  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
!!  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
!!  BATTELLE OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
!!  OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!!  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!!  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
!!  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!  POSSIBILITY OF SUCH DAMAGE.
!!
!!  DISCLAIMER
!!  The Software was produced by Battelle under Contract No. DE-AC05-76RL01830 with the
!!  Department of Energy.  The U.S. Government is granted for itself and others acting on its
!!  behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce,
!!  prepare derivative works, distribute copies to the public, perform publicly and display
!!  publicly, and to permit others to do so.  The specific term of the license can be identified
!!  by inquiry made to Battelle or DOE.  Neither the United States nor the United States
!!  Department of Energy, nor any of their employees, makes any warranty, express or implied,
!!  or assumes any legal liability or responsibility for the accuracy, completeness or
!!  usefulness of any data, apparatus, product or process disclosed, or represents that its use
!!  would not infringe privately owned rights.
!!______________________________________________________________________________________________
!!______________________________________________________________________________________________
!! Author: Tim Johnson
!!
!! email: e4d@pnnl.gov
!!
!! Description: The assemble module contains subroutines for collecting
!! the potentials from each slave processors that are necessary to
!! to reconstruct a simulated ERT or SIP survey
!!
!! Input variables: see individual subroutines
!!
!! Output variables: see individual subroutines
!!____________________________________________________________________

  use vars
  implicit none

contains

  !___________________________________________________________________
  subroutine assemble_data(flg)

    implicit none
    integer :: flg                    !Input variable. 1=real; 2=complex
    integer :: i,j                      !indexing variable
    integer :: a,b                    !current source electrode indexes
    integer :: m,n                    !potential electrode indexes
    integer :: e1,e2                  !first and last poles stored by this process
    integer :: indx                   !current index in my_drows and my_dvals
    integer :: p                      !pole solution index


    !! set the pole solution ownership bounds for this process
    e1=eind(my_rank,1)
    e2=eind(my_rank,2)

    !! find the number of measurements using a pole owned by this process
    !! and allocate my_drows, my_dvals
    if(allocated(my_drows).and.res_flag) then
       deallocate(my_drows)
       deallocate(my_dvals)
    end if
    if(.not.allocated(my_drows)) then
       nmy_drows=0
       do i=1,nm
          a=s_conf(i,1)
          b=s_conf(i,2)
          if((a>=e1 .and. a<=e2) .or. (b>=e1.and. b<=e2)) then
             nmy_drows=nmy_drows+1
          end if
       end do
       allocate(my_drows(nmy_drows),my_dvals(nmy_drows))
    end if

    indx=0
    my_dvals=0


    if(flg==1) then
       !! Assemble the part of the simulated data owned by this process
       do i=1,nm
          !! Set a,b,m and n for this measurement
          a = s_conf(i,1)
          b = s_conf(i,2)
          m = s_conf(i,3)
          n = s_conf(i,4)

          if((a>=e1 .and. a<=e2) .or. (b>=e1.and. b<=e2)) then
             !! If this process owns a source or sink electrode pole for this
             !! measurement, then compute the part of the measurement depending on
             !! on this pole
             indx=indx+1
             my_drows(indx)=i

             do p=e1,e2
                !! Loop of the poles owned by this process, find the ones used by this
                !! measurement, and add (or subtract) the simulated potential at the
                !! potential electrodes (m and n)
                if(p==a) then 
                   if(m.ne.0) then
                      do j=1,4
                         my_dvals(indx)=my_dvals(indx) + real(poles(source_nodes(m,j),a-e1+1)*source_currents(m,j))
                      end do
                      !my_dvals(indx)=my_dvals(indx) + real(poles(e_nods(m),a-e1+1)
                   end if
                   if(n.ne.0) then
                      do j=1,4
                         my_dvals(indx)=my_dvals(indx) - real(poles(source_nodes(n,j),a-e1+1)*source_currents(n,j))
                      end do
                      !my_dvals(indx)=my_dvals(indx) - real(poles(e_nods(n),a-e1+1))
                   end if
                end if
                
                if(p==b) then
                   if(m.ne.0) then
                      do j=1,4
                         my_dvals(indx)=my_dvals(indx) - real(poles(source_nodes(m,j),b-e1+1)*source_currents(m,j))
                      end do
                      !my_dvals(indx)=my_dvals(indx) - real(poles(e_nods(m),b-e1+1))
                   end if
                   if(n.ne.0) then
                      do j=1,4
                         my_dvals(indx)=my_dvals(indx) + real(poles(source_nodes(n,j),b-e1+1)*source_currents(n,j))
                      end do
                      !my_dvals(indx)=my_dvals(indx) + real(poles(e_nods(n),b-e1+1))
                   end if
                end if

             end do
          end if
       end do
    end if

    if(flg==2) then
       !! Assemble the data for the complex potential (see comments above for details)
       do i=1,nm
          a = s_conf(i,1)
          b = s_conf(i,2)
          m = s_conf(i,3)
          n = s_conf(i,4)

          if((a>=e1 .and. a<=e2) .or. (b>=e1.and. b<=e2)) then
             indx=indx+1
             my_drows(indx)=i

             do p=e1,e2
                
                if(p==a) then
                   if(m.ne.0) then
                      do j=1,4
                         my_dvals(indx)=my_dvals(indx) + real(polesi(source_nodes(m,j),a-e1+1)*source_currents(m,j))
                      end do
                   end if
                   if(n.ne.0) then
                      do j=1,4
                         my_dvals(indx)=my_dvals(indx) - real(polesi(source_nodes(n,j),a-e1+1)*source_currents(n,j))
                      end do
                   end if
                   !if(m.ne.0) my_dvals(indx)=my_dvals(indx) + real(polesi(e_nods(m),a-e1+1))
                   !if(n.ne.0) my_dvals(indx)=my_dvals(indx) - real(polesi(e_nods(n),a-e1+1))
                end if
                
                if(p==b) then
                   if(m.ne.0) then
                      do j=1,4
                         my_dvals(indx)=my_dvals(indx) - real(polesi(source_nodes(m,j),b-e1+1)*source_currents(m,j))
                      end do
                   end if
                   if(n.ne.0) then
                      do j=1,4
                         my_dvals(indx)=my_dvals(indx) + real(polesi(source_nodes(n,j),b-e1+1)*source_currents(n,j))
                      end do
                   end if
                   !if(m.ne.0) my_dvals(indx)=my_dvals(indx) - real(polesi(e_nods(m),b-e1+1))
                   !if(n.ne.0) my_dvals(indx)=my_dvals(indx) + real(polesi(e_nods(n),b-e1+1))
                   
                end if
             end do
          end if
       end do
    end if


  end subroutine assemble_data
  !___________________________________________________________

end module assemble
