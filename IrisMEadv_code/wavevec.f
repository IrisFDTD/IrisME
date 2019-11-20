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

      subroutine wavevec(ndimk1x, ndimk1y, n1x, n1y, ndimk2x,
     &                    ndimk2y, n2x, n2y, rlx, rly, ra, rb,
     &                    theta, fi, gvaccuum, eps1, eps2, eps3,  
     &                    vk1x, vk1y, vk1z, vk2x, vk2y, vk2z,
     &                    vk3x, vk3y, vk3z)
c
       implicit real*8 (a-h,o-z)
       dimension vk1x(-ndimk1x:ndimk1x), vk1y(-ndimk1y:ndimk1y),
     &           vk3x(-ndimk1x:ndimk1x), vk3y(-ndimk1y:ndimk1y),
     &           vk2x(0:ndimk2x), vk2y(0:ndimk2y)
       complex*16 vk1z(-ndimk1x:ndimk1x,-ndimk1y:ndimk1y),
     &            vk2z(0:ndimk2x,0:ndimk2y),
     &            vk3z(-ndimk1x:ndimk1x,-ndimk1y:ndimk1y)

c
       pi = 4.d0 * datan(1.d0)
c
       g1 = gvaccuum * dsqrt(eps1)
       g12 = g1**2
       do m = -n1x, n1x
          vk1x(m) = dsin(theta)*dcos(fi)*g1 + 2.d0*m*pi/rlx
       end do
       do n = -n1y, n1y
          vk1y(n) = dsin(theta)*dsin(fi)*g1 + 2.d0*n*pi/rly
       end do
       do m = -n1x, n1x
       do n = -n1y, n1y
        vk1z(m,n) = cdsqrt(dcmplx(g12 - vk1x(m)**2 - vk1y(n)**2,0.d0))
       end do
       end do
c
       g0x = pi/ra
       g0y = pi/rb
       do i = 0 , n2x
          vk2x(i) = i * g0x
       end do
       do j = 0 , n2y
          vk2y(j) = j * g0y
       end do
       do i = 0 , n2x
       do j = 0 , n2y
        vk2z(i,j) = cdsqrt(dcmplx(eps2*gvaccuum**2 - vk2x(i)**2 - 
     &                            vk2y(j)**2,0.d0))
       end do
       end do
c
       g3 = gvaccuum * dsqrt(eps3)
       g32 = g3**2
       do m = -n1x, n1x
          vk3x(m) = vk1x(m)
       end do
       do n = -n1y, n1y
          vk3y(n) = vk1y(n)
       end do
       do m = -n1x, n1x
       do n = -n1y, n1y
        vk3z(m,n) = cdsqrt(dcmplx(g32 - vk3x(m)**2 - vk3y(n)**2,0.d0))
       end do
       end do
c
       return
       end
