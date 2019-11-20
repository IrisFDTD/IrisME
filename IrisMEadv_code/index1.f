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

      
      subroutine index1(ndimk1x, ndimk1y, ndim1, n1x, n1y, ip1, 
     &                   im1, nummod1)
c
c      indexes modes in region II. 
c      The convention is: first mode P, later mode S. 
c      starting with, P: (-n1x,j), (0,j)...(n1x,j)
c      following with S: (-n1y,j), (0,j)...(n1y,j), (j:-n1y->n1y)
c
c      ip1(pol,i,j) gives the index for mode (pol,i,j).
c      im1(ind,0) gives polarization for mode indexed by 'ind'
c      im1(ind,1) gives m for mode indexed by 'ind'
c      im1(ind,2) gives n for mode indexed by 'ind'
c 
       implicit real*8(a-h,o-z)
       dimension ip1(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y),
     &           im1(ndim1, 0:2)
c
       do ipol = 0, 1
       do i = -n1x, n1x
       do j = -n1y, n1y
          ip1(ipol, i, j) = 0
       end do 
       end do 
       end do 
c
       ind = 0
c
       ipol = 1
       do i = -n1x, n1x
       do j = -n1y, n1y
            ind = ind + 1
            ip1(ipol, i, j) = ind
            im1(ind,0) = ipol
            im1(ind,1) = i
            im1(ind,2) = j
       end do
       end do
c
       ipol = 0
       do i = -n1x, n1x
       do j = -n1y, n1y
            ind = ind + 1
            ip1(ipol, i, j) = ind
            im1(ind,0) = ipol
            im1(ind,1) = i
            im1(ind,2) = j
       end do
       end do
c
       nummod1 = ind
c
       return
       end
