!copyright>        openradioss
!copyright>        copyright (c) 1986-2024 altair engineering inc.
!copyright>
!copyright>        this program is free software: you can redistribute it and/or modify
!copyright>        it under the terms of the gnu affero general public license as published by
!copyright>        the free software foundation, either version 3 of the license, or
!copyright>        (at your option) any later version.
!copyright>
!copyright>        this program is distributed in the hope that it will be useful,
!copyright>        but without any warranty; without even the implied warranty of
!copyright>        merchantability or fitness for a particular purpose.  see the
!copyright>        gnu affero general public license for more details.
!copyright>
!copyright>        you should have received a copy of the gnu affero general public license
!copyright>        along with this program.  if not, see <https://www.gnu.org/licenses/>.
!copyright>
!copyright>
!copyright>        commercial alternative: altair radioss software
!copyright>
!copyright>        as an alternative to this open-source version, altair also offers altair radioss
!copyright>        software under a commercial license.  contact altair to discuss further if the
!copyright>        commercial version may interest you: https://www.altair.com/radioss/.
!chd|====================================================================
!chd|  hm_read_mat129                 source/materials/mat/mat129/hm_read_mat129.f
!chd|-- called by -----------
!chd|        hm_read_mat                   source/materials/mat/hm_read_mat.f
!chd|-- calls ---------------
!chd|====================================================================
! ======================================================================================================================

      module hm_read_mat129_mod
      contains
  
      subroutine hm_read_mat129( mtag, matparam ,                         &
                   nuvar, unitab   ,lsubmodel,                       &
                   mat_id   ,titr     ,pm         ,                            &
                   iout     ,npropm   )
  ! ----------------------------------------------------------------------------------------------------------------------
  !                                                   modules
  ! ----------------------------------------------------------------------------------------------------------------------
      use elbuftag_mod
      use matparam_def_mod
      use unitab_mod
      use message_mod
      use submodel_mod
      use constant_mod , only : one ,two, zero
! --------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------
!                                                 implicit none
! --------------------------------------------------------------------------------------------------------------------
      implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   included files
! ----------------------------------------------------------------------------------------------------------------------
#include "my_real.inc"
      
!c-----------------------------------------------
!c   d u m m y   a r g u m e n t s
!c-----------------------------------------------
      integer, intent(in)                          :: mat_id
      integer, intent(in)                          :: iout
      integer, intent(in)                          :: npropm
      integer, intent(out)                         :: nuvar 
      type (unit_type_),intent(in) ::unitab 
      type(submodel_data), dimension(nsubmod),intent(in) :: lsubmodel
      character(len=nchartitle) ,intent(in)             :: titr    
      my_real, dimension(npropm)    ,intent(inout)   :: pm     
      type(matparam_struct_) ,intent(inout) :: matparam
      type(mlaw_tag_), intent(inout)  :: mtag
!-----------------------------------------------
!   l o c a l   v a r i a b l e s
!-----------------------------------------------
      logical :: is_available,is_encrypted
      integer :: ilaw,lcid
      my_real :: rho0, kreal,delta1  
      ! -------------------------
!=======================================================================
      is_encrypted = .false.
      is_available = .false.
      ilaw  = 129
!--------------------------------------------------------
!
      call hm_option_is_encrypted(is_encrypted)
!
!--------------------------------------------------------
!     read input fields
!--------------------------------------------------------

      call hm_get_floatv('LSD_RO'           ,rho0      ,is_available, lsubmodel, unitab)
    
      call hm_get_floatv('LSD_K'            ,kreal     ,is_available, lsubmodel, unitab)
      call hm_get_floatv('LSD_DELTA1'       ,delta1    ,is_available, lsubmodel, unitab)
      call hm_get_intv  ('LSD_LCID'         ,lcid      ,is_available, lsubmodel)     

!-------------------------------------
      nuvar = 4
!-------------------------------------
      
      matparam%niparam = 1
      matparam%nuparam = 2
!          
      allocate (matparam%uparam(matparam%nuparam))
      allocate (matparam%iparam(matparam%niparam))
!     
      matparam%iparam(1) =  lcid

      matparam%uparam(1) =   kreal 
      matparam%uparam(2) =   delta1  
!-------------------------------------------------
      pm(1)  = rho0
      pm(89) = rho0
!-------------------------------------------------
      mtag%g_pla  = 1
      mtag%l_pla  = 1
      mtag%l_dmg  = 1
!-------------------------------------------------
      ! properties compatibility  
      call init_mat_keyword(matparam,"LUNG_TISSUE")       
!-------------------------------------------------
      write(iout,1050) trim(titr),mat_id,129
      write(iout,1000)
      if (is_encrypted) then
        write(iout,'(5x,a,//)')'CONFIDENTIAL DATA'
      else
        write(iout,1060) rho0
        write(iout,1100) kreal,delta1, lcid
      endif       

!-----------
      return
!-----------
 1000 format(                                                                 & 
     5x,'MAT LUNG TISSUE                       ',/,                       & 
     5x,'-----------------------------------   ',//)           
 1050 format(/                                                               &
      5x,a,/,                                                                &
      5x,'MATERIAL NUMBER . . . . . . . . . . . . .=',i10/,                  &
      5x,'MATERIAL LAW. . . . . . . . . . . . . . .=',i10/)       
 1060 format(                                                                &
      5x,'INITIAL DENSITY . . . . . . . . . . . . .=',1pg20.13/)  
 1100 format(                                                                & 
      5x,'K . . . . . . . . . . . . . . . . . . . .=',1PG20.13/,  &
      5x,'DELTA 1 . . . . . . . . . . . . . . . . .=',1PG20.13/,  &
      5x,'LCID. . . . . . . . . . . . . . . . . . .=',I10/)
!-----------------
    end subroutine hm_read_mat129
end module hm_read_mat129_mod             