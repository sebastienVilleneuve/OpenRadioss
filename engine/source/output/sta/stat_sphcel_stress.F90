!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2024 Altair Engineering Inc.
!Copyright>
!Copyright>        This program is free software: you can redistribute it and/or modify
!Copyright>        it under the terms of the GNU Affero General Public License as published by
!Copyright>        the Free Software Foundation, either version 3 of the License, or
!Copyright>        (at your option) any later version.
!Copyright>
!Copyright>        This program is distributed in the hope that it will be useful,
!Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
!Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!Copyright>        GNU Affero General Public License for more details.
!Copyright>
!Copyright>        You should have received a copy of the GNU Affero General Public License
!Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
!Copyright>
!Copyright>
!Copyright>        Commercial Alternative: Altair Radioss Software
!Copyright>
!Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
!Copyright>        software under a commercial license.  Contact Altair to discuss further if the
!Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
      module stat_sphcel_stress_mod
        contains
! ======================================================================================================================
!                                                   PROCEDURES
! ======================================================================================================================
!! \sphcel element write stress in state file 
        subroutine stat_sphcel_stress(numsph    ,nisp       ,ngroup     ,nparg         ,sizloc          , npart , &
                                    sizp0     ,nspmd      ,stat_numelsph ,stat_numelsph_g ,  &
                                    kxsp       ,ipartsph   ,ipart_state   ,stat_indxsph    ,  &
                                    iparg     ,elbuf_tab  ,wa         ,wap0          )
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
          use ELBUFDEF_MOD, only: elbuf_struct_
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Included files
! ----------------------------------------------------------------------------------------------------------------------
#include "my_real.inc"
#include "task_c.inc"
#include "units_c.inc"
#include "mvsiz_p.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          integer,                                   intent(in) :: numsph            !< size for array definition number sph cells
          integer,                                   intent(in) :: nisp              !< size for array definition
          integer,                                   intent(in) :: ngroup            !< size for array definition
          integer,                                   intent(in) :: nparg             !< size for array definition
          integer,                                   intent(in) :: sizloc            !< size for array definition
          integer,                                   intent(in) :: sizp0             !< size for array definition
          integer,                                   intent(inout) :: stat_numelsph    !< size for array definition
          integer,                                   intent(inout) :: stat_numelsph_g  !< size for array definition
          integer,                                   intent(in) :: nspmd               !< size for array definition
          integer,                                   intent(in) :: npart               !< size for array definition
          integer,                                   intent(in) :: kxsp(nisp,numsph) 
          integer,                                   intent(in) :: ipartsph(numsph) 
          integer,                                   intent(in) :: ipart_state(npart) 
          integer,                                   intent(inout) :: stat_indxsph(numsph) 
          integer,                                   intent(in) :: iparg(nparg,ngroup) 
          type(elbuf_struct_),                       intent(in) :: elbuf_tab(ngroup)
          double precision,                          intent(inout) :: wa(sizloc)
          double precision,                          intent(inout) :: wap0(sizp0)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer i,j,k,n,jj,len,ioff,ie,ng,nel,nft,lft,llt,ity,id,iprt0,iprt,igtyp,iprop,g_pla
          integer ii(6)
          integer ptwa(stat_numelsph),ptwa_p0(0:max(1,stat_numelsph_g))
          character*100 delimit,line
          data delimit(1:60) &
           /'#---1----|----2----|----3----|----4----|----5----|----6----|'/
          data delimit(61:100) &
           /'----7----|----8----|----9----|----10---|'/
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          jj = 0
          !
          if (stat_numelsph /= 0) then
            !
            ie=0
            do ng=1,ngroup
              ity = iparg(5,ng)
              if (ity == 51) then
                nel  = iparg(2,ng)
                do i=1,6
                  ii(i) = nel*(i-1)
                enddo
                nft  = iparg(3,ng)
                iprop = kxsp(4,nft+1)
                lft=1
                llt=nel
                !
                do i=lft,llt
                  n = i + nft
                  iprt=ipartsph(n)
                  if (ipart_state(iprt) /= 0) then
                    jj = jj + 1
                    wa(jj) = elbuf_tab(ng)%gbuf%off(i)
                    jj = jj + 1
                    wa(jj) = iprt
                    jj = jj + 1
                    wa(jj) = kxsp(nisp,n)
                    jj = jj + 1
                    wa(jj) = elbuf_tab(ng)%gbuf%sig(ii(1)+i)
                    jj = jj + 1
                    wa(jj) = elbuf_tab(ng)%gbuf%sig(ii(2)+i)
                    jj = jj + 1
                    wa(jj) = elbuf_tab(ng)%gbuf%sig(ii(3)+i)
                    !---
                    ie=ie+1
                    !---            pointeur de fin de zone dans wa
                    ptwa(ie)=jj
                  endif ! if (ipart_state(iprt) /= 0)
                enddo  !  do i=lft,llt
                ! end loop over truss elements
              endif ! ity == 51
            enddo ! ng = 1, ngroup
          endif ! if (stat_numelsph == 0) then
!-----------------------------------------------------------------------
!     sphcel - write
!-----------------------------------------------------------------------
            if (nspmd == 1) then
              ptwa_p0(0)=0
              do n=1,stat_numelsph
                ptwa_p0(n)=ptwa(n)
              enddo
              len=jj
              do j=1,len
                wap0(j)=wa(j)
              enddo
            else
              call spmd_stat_pgather(ptwa,stat_numelsph,ptwa_p0,stat_numelsph_g)
              len = 0
              call spmd_rgather9_dp(wa,jj,wap0,sizp0,len)
            endif
!-------------------------------------
            if (ispmd == 0 .and. len > 0) then
              iprt0 = 0
              do n=1,stat_numelsph_g
                k=stat_indxsph(n)
                j=ptwa_p0(k-1)
                !
                ioff = nint(wap0(j + 1))
                if (ioff /= 0) then
                  iprt  = nint(wap0(j + 2)) 
                  id    = nint(wap0(j + 3))
!--------------------------------------
                    if (iprt /= iprt0) then
                      write(iugeo,'(a)') delimit
                      write(iugeo,'(a)')'/INISPHCEL/STRESS'
                      write(iugeo,'(a)')&
                                 '#----------------------------------------------------------'
                      write(iugeo,'(a)')'#SPHCEL_ID'
                      write(iugeo,'(a)')'#format:(1p3e20.13) #(sig1,sig2,sig3)'
                      write(iugeo,'(a)')&
                                 '#----------------------------------------------------------'
                      !
                      iprt0=iprt
                    endif ! if (iprt /= iprt0)
                    !
                    write(iugeo,'(i10)') id
                    write(iugeo,'(1p3e20.13)')wap0(j+4),wap0(j+5),wap0(j+6)
!--------------------------------------
                endif  !  if (ioff >= 1)
              enddo  !  do n=1,stat_numelsph_g
            endif  !  if (ispmd == 0.and.len > 0)

          return
! ----------------------------------------------------------------------------------------------------------------------
        end subroutine stat_sphcel_stress
      end module stat_sphcel_stress_mod
      
            