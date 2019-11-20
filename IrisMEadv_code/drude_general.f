! -------------------------------------------------------------------------------------------------------------------------------------
! This file is part of the IrisMEadv implementation. IriMEadv is a fast and accurate program for 
! the calculation of scattering coefficients in arrays of rectangular holes in metals. 
! The following subroutines were implemented by Sergio G Rodrigo sergut@unizar.es
!
! IrisME is licensed under the AGPL and it is free to use. However, if you are using, or plan to use it, 
! specially if it is for research or academic purposes, please send an email with your name, institution and a 
! brief description of your interest for this program. If you use IrisME in a work that leads to a scientific 
! for academic publication, we would appreciate it if you would kindly cite IrisME in your manuscript as:
!      
! > S. Carretero-Palacios, F. J. García-Vidal, L. Martín-Moreno, and S. G. Rodrigo,
! > Effect of film thickness and dielectric environment on optical transmission through subwavelength holes, 
! > Phys. Rev. B 85, 035417 (2012).
!
! Commercial applications may also acquire a commercial license. Please contact Sergio G Rodrigo sergut@unizar.es, for details.
! 
! Copyright (C) 2019 Sergio G Rodrigo sergut@unizar.es
!
! IrisME is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public 
! License as published by the Free Software Foundation, either version 3 of the License,or (at your option) any later version.
!
! IrisME is distributed in the hope that it will be useful for research or/and academic purpouses, but WITHOUT ANY WARRANTY; 
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public
! License for more details. You should have received a copy of the GNU Affero General Public License 
! along with IrisME. If not,see http://www.gnu.org/licenses/.
! -------------------------------------------------------------------------------------------------------------------------------------

      
      subroutine init_metal(w,eps,inFDTD,msg_in)
	implicit none
	! ********************************************************************** !
	! En unidades atómicas, en base a que E = h-barra*w como h-barra = 1 se  !
	! mide en eV y de hecho E = w    
	real*8                     :: w !in eV
	complex*16                 :: eps
	integer                    :: inFDTD

	complex                    :: e_inf   
	real                       :: mu_inv_bbox  	       
	integer                    :: nDrude  
	character(1)              :: msg_in

	character(1)              :: msg
	real,pointer               :: gamma(:)   ! En Hartree´s
	real,pointer               :: wplasma(:)   ! En hartree´s 
	integer                    :: mLorentz
	real,pointer               :: de(:)
	real,pointer               :: fil(:)
	real,pointer               :: gammal(:)
	
	integer::nD,mL
      complex,parameter          :: ci=(0.0,1.0)
	real                       :: er,ei

	integer::i_metal,i_metal_fit,old_access,logfile,i,j

	!**************************************************************************!
    
	rewind(unit=inFDTD)  
	 do while(trim(msg)/=trim(msg_in))
	   read(inFDTD,*) msg
	 end do
	 	
	 read(inFDTD,*) old_access		 
	    if(old_access==1)then
			 read(inFDTD,*) e_inf
			 read(inFDTD,*) mu_inv_bbox
			 mu_inv_bbox= 1.0/mu_inv_bbox

			 !>>> Drude terms
			 read(inFDTD,*) nDrude	 	 			 
			 allocate(wplasma(nDrude),gamma(nDrude))

			 if(nDrude/=0)then
			 do nD=1,nDrude
				 read(inFDTD,*) wplasma(nD)							 
				 read(inFDTD,*) gamma(nD)							 
			 end do
			 end if
			 !    Drude terms

			 !>>> Lorentz terms
			 read(inFDTD,*) mLorentz			 
			 if(mLorentz/=0)then	
			   allocate(fil(mLorentz),gammal(mLorentz),de(mLorentz))			
				 do mL=1,mLorentz
					 read(inFDTD,*) fil(mL)					 
					 read(inFDTD,*) gammal(mL)					 
					 read(inFDTD,*) de(mL)					 
				 end do
			 end if	 
		 else	
		   ! read(inFDTD,*) i_metal !>>> Metal type
		   ! read(inFDTD,*) i_metal_fit   !>>> Old/new tabulated values
		   ! call select_metal	     			 
	   end if	

!      >>> Drude - Lorentz
      er = 0.0
	ei = 0.0
      er = real(e_inf)
	ei = imag(e_inf)
	do i=1,nDrude
      er =er - wplasma(i)**2/(w**2 + gamma(i)**2)	    			
	ei =ei + (gamma(i)/w)*(wplasma(i)**2/(w**2 + gamma(i)**2))	
      end do	
	do j=1,mLorentz
	er =er-real(de(j)*gammal(j)**2/(w**2 - gammal(j)**2+ci*fil(j)*w))		
	ei =ei-imag(de(j)*gammal(j)**2/(w**2 - gammal(j)**2+ci*fil(j)*w))	
	end do	
	eps=er+ci*ei
      return
      end subroutine
