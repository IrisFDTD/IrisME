! -------------------------------------------------------------------------------------------------------------------------------------
! This file is part of the IrisMEadv implementation. IriMEadv is a fast and accurate program for 
! the calculation of scattering coefficients in arrays of rectangular holes in metals. 
! The following subroutines belongs to the "Multiple scattering formalism" on which is based IrisMEadv.
! The subroutines were implemented by Luis Martin-Moreno lmm@unizar.es
!
! IrisME is licensed under the AGPL and it is free to use. However, if you are using, or plan to use it, 
! specially if it is for research or academic purposes, please send an email with your name, institution and a 
! brief description of your interest for this program. If you use IrisME in a work that leads to a scientific 
! for academic publication, we would appreciate it if you would kindly cite IrisME in your manuscript as:
!      
! > L. Martín-Moreno, F.J. García-Vidal, Minimal model for optical transmission through holey metal films, 
! > J. Phys.: Condens. Matter 20, 304214 (2008). 
! > and/or
! > F.J. Garcia-Vidal, L. Martin-Moreno, T.W. Ebbesen, L. Kuipers, Light passing through subwavelength apertures, 
! > Rev. Mod. Phys. 82, 729–787 (2010)
!
! Commercial applications may also acquire a commercial license. Please contact lmm@unizar.es, for details.
! 
! Copyright (C) 2019 Luis Martin-Moreno lmm@unizar.es
!
! IrisME is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public 
! License as published by the Free Software Foundation, either version 3 of the License,or (at your option) any later version.
!
! IrisME is distributed in the hope that it will be useful for research or/and academic purpouses, but WITHOUT ANY WARRANTY; 
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public
! License for more details. You should have received a copy of the GNU Affero General Public License 
! along with IrisME. If not,see http://www.gnu.org/licenses/.
! -------------------------------------------------------------------------------------------------------------------------------------

      
      subroutine index2(ndimk2x, ndimk2y, ndim2, n2x, n2y, width, 
     &                   cutoff, vk2z, ip2, im2, fackz, nummod2, numte)
c
c      indexes modes in region II. 
c      The convention is: first mode TE, later mode TM. 
c      starting with, TE: (i,0), then (0,j), then (1,j), (2,j)...(n2x,j)
c      following with TM: (1,j), (2,j)...(n2x,j), (j:1->n2y)
c
c      ip2(pol,i,j) gives the index for mode (pol,i,j).
c      im2(ind,0) gives polarization for mode indexed by 'ind'
c      im2(ind,1) gives m for mode indexed by 'ind'
c      im2(ind,2) gives n for mode indexed by 'ind'
c 
       implicit real*8(a-h,o-z)
       dimension ip2(0:1, 0:ndimk2x, 0:ndimk2y), im2(ndim2, 0:2) 
       complex*16 ui, vk2z(0:ndimk2x,0:ndimk2y), fackz(ndim2)
c
       ui = dcmplx(0.d0, 1.d0)
c
       ind = 0
       indte = 0
       do ipol = 0, 1
       do i = 0, n2x
       do j = 0, n2y
          ip2(ipol, i, j) = 0
       end do 
       end do 
       end do 
c
       ipol = 1
       j = 0
       do i = 1, n2x
          if( dimag(vk2z(i,j))*width .le. cutoff) then
            ind = ind + 1
            indte = indte + 1
            ip2(ipol, i, j) = ind 
            im2(ind,0) = ipol
            im2(ind,1) = i
            im2(ind,2) = j
            fackz(ind) = cdexp(ui*vk2z(i,j)*width)
          end if
       end do
c
       ipol = 1
       i = 0
       do j = 1, n2y
          if( dimag(vk2z(i,j))*width .le. cutoff) then
            ind = ind + 1
            indte = indte + 1
            ip2(ipol, i, j) = ind
            im2(ind,0) = ipol
            im2(ind,1) = i
            im2(ind,2) = j
            fackz(ind) = cdexp(ui*vk2z(i,j)*width)
          end if
       end do
c
       ipol = 1
       do i = 1, n2x
       do j = 1, n2y
          if( dimag(vk2z(i,j))*width .le. cutoff) then
            ind = ind + 1
            indte = indte + 1
            ip2(ipol, i, j) = ind
            im2(ind,0) = ipol
            im2(ind,1) = i
            im2(ind,2) = j
            fackz(ind) = cdexp(ui*vk2z(i,j)*width)
          end if
       end do
       end do
c
       ipol = 0
       do i = 1, n2x
       do j = 1, n2y
          if( dimag(vk2z(i,j))*width .le. cutoff) then
            ind = ind + 1
            ip2(ipol, i, j) = ind
            im2(ind,0) = ipol
            im2(ind,1) = i
            im2(ind,2) = j
            fackz(ind) = cdexp(ui*vk2z(i,j)*width)
          end if
       end do
       end do
c
       nummod2 = ind
       numte = indte
c
       return
       end
