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
      
      subroutine inttoASCII(ino,charno)	
	 implicit none

	 integer,intent(in)::ino
	 character(5),intent(out)::charno
	 integer::sCase,U,D,C,M,iInter  

	if (ino.LE.9) then
	     sCase=0  	 
	  end if

	  if ((ino.GT.9).AND.(ino.LE.99)) then
		   sCase=1
	  end if

	  if (ino.GT.99) then
	     sCase=2 	 
	  end if  

	  if (ino.GT.999) then
          sCase=3 	 
        end if  
   
        if (ino.GT.9999) then
          sCase=4 	 
        end if    

	select case(sCase)
	   case(0)
		    U=ino; charno=ACHAR(INT(U+48))
	   case(1)
		   D=ino/10
	       U=Mod(ino,10)
           charno=ACHAR(INT(D+48))//ACHAR(INT(U+48))
	   case(2)
		   C=ino/100
	       U=Mod(ino,10)
	  	   D=(ino-(100*C+U))/10			
           charno=ACHAR(INT(C+48))
		   charno=trim(charno)//ACHAR(INT(D+48))//ACHAR(INT((U+48)))		
	   case(3)
		   M=ino/1000
		   iInter = ino-1000*M
		   C=iInter/100
		   U=Mod(iInter,10)
		   D=(iInter-(100*C+U))/10				
         charno=ACHAR(INT(M+48))//ACHAR(INT(C+48))
         charno=trim(charno)//ACHAR(INT(D+48))//ACHAR(INT((U+48)))  

	   case default
		  !  write(logfile,*) "...Error - Subroutine:inttoASCII"
	end select   
      return
      end subroutine inttoASCII
