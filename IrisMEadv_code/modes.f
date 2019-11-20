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

      
      subroutine modes(ndimk1x, ndimk1y, ndimk2x, ndimk2y, n1x, n1y, 
     &                 n2x, n2y, vk1x, vk1y, vk1z, vk2x, vk2y, vk2z, 
     &                 vk3x, vk3y, vk3z, theta, gvaccuum, eps1, eps2, 
     &                 eps3, f1p, f1m, f2p, f2m, f3p, f3m)
c
      implicit real*8(a-h,o-z)
      dimension vk1x(-ndimk1x:ndimk1x), vk1y(-ndimk1y:ndimk1y),
     &          vk3x(-ndimk1x:ndimk1x), vk3y(-ndimk1y:ndimk1y),
     &          vk2x(0:ndimk2x), vk2y(0:ndimk2y)
      complex*16 cz, ur, vk1z(-ndimk1x:ndimk1x, -ndimk1y:ndimk1y),
     &           vk2z(0:ndimk2x, 0:ndimk2y),
     &           vk3z(-ndimk1x:ndimk1x, -ndimk1y:ndimk1y),
     &           f1p(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f1m(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f2p(0:1, 0:ndimk2x, 0:ndimk2y, 4),
     &           f2m(0:1, 0:ndimk2x, 0:ndimk2y, 4),
     &           f3p(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f3m(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4)
c
       cz = dcmplx(0.d0,0.d0)
       ur = dcmplx(1.d0,0.d0)
       gv2 = gvaccuum**2
       pi = 4.d0*datan(1.d0)
       dospi = 2.d0 * pi
c
c S Polarization
c
       do m = -n1x, n1x
       do n = -n1y, n1y
          f1p(0,m,n,1) = + vk1y(n) 
          f1p(0,m,n,2) = - vk1x(m)
          f1p(0,m,n,3) = + vk1x(m) * vk1z(m,n) / dospi 
          f1p(0,m,n,4) = + vk1y(n) * vk1z(m,n) / dospi
c
          f1m(0,m,n,1) = + vk1y(n)
          f1m(0,m,n,2) = - vk1x(m)
          f1m(0,m,n,3) = - vk1x(m) * vk1z(m,n) / dospi
          f1m(0,m,n,4) = - vk1y(n) * vk1z(m,n) / dospi
c
          f3p(0,m,n,1) = + vk3y(n)
          f3p(0,m,n,2) = - vk3x(m)
          f3p(0,m,n,3) = + vk3x(m) * vk3z(m,n) / dospi
          f3p(0,m,n,4) = + vk3y(n) * vk3z(m,n) / dospi
c
          f3m(0,m,n,1) = + vk3y(n)
          f3m(0,m,n,2) = - vk3x(m)
          f3m(0,m,n,3) = - vk3x(m) * vk3z(m,n) / dospi
          f3m(0,m,n,4) = - vk3y(n) * vk3z(m,n) / dospi
c
       end do
       end do
c
       if (theta.eq.0.d0) then
          f1p(0,0,0,1) = cz
          f1p(0,0,0,2) = ur
          f1p(0,0,0,3) = - dsqrt(eps1)
          f1p(0,0,0,4) = cz
c
          f1m(0,0,0,1) = cz
          f1m(0,0,0,2) = ur
          f1m(0,0,0,3) = dsqrt(eps1)
          f1m(0,0,0,4) = cz
c
          f3p(0,0,0,1) = cz
          f3p(0,0,0,2) = ur
          f3p(0,0,0,3) = - dsqrt(eps3)
          f3p(0,0,0,4) = cz
c
          f3m(0,0,0,1) = cz
          f3m(0,0,0,2) = ur
          f3m(0,0,0,3) = + dsqrt(eps3)
          f3m(0,0,0,4) = cz
       end if
c
c P Polarization
c
       do m = -n1x, n1x
       do n = -n1y, n1y
          f1p(1,m,n,1) = + vk1x(m) * vk1z(m,n) / dospi
          f1p(1,m,n,2) = + vk1y(n) * vk1z(m,n) / dospi
          f1p(1,m,n,3) = - eps1 * vk1y(n)
          f1p(1,m,n,4) = + eps1 * vk1x(m)
c
          f1m(1,m,n,1) = + vk1x(m) * vk1z(m,n) / dospi
          f1m(1,m,n,2) = + vk1y(n) * vk1z(m,n) / dospi
          f1m(1,m,n,3) = + eps1 * vk1y(n)
          f1m(1,m,n,4) = - eps1 * vk1x(m)
c
          f3p(1,m,n,1) = + vk3x(m) * vk3z(m,n) / dospi
          f3p(1,m,n,2) = + vk3y(n) * vk3z(m,n) / dospi
          f3p(1,m,n,3) = - eps3 * vk3y(n)
          f3p(1,m,n,4) = + eps3 * vk3x(m)
c
          f3m(1,m,n,1) = + vk3x(m) * vk3z(m,n) / dospi
          f3m(1,m,n,2) = + vk3y(n) * vk3z(m,n) / dospi
          f3m(1,m,n,3) = + eps3 * vk3y(n)
          f3m(1,m,n,4) = - eps3 * vk3x(m)
c
       end do
       end do
c
       if (theta.eq.0.d0) then
          f1p(1,0,0,1) = ur
          f1p(1,0,0,2) = cz
          f1p(1,0,0,3) = cz
          f1p(1,0,0,4) = dsqrt(eps1)
c
          f1m(1,0,0,1) = ur
          f1m(1,0,0,2) = cz
          f1m(1,0,0,3) = cz
          f1m(1,0,0,4) = - dsqrt(eps1)
c
          f3p(1,0,0,1) = ur
          f3p(1,0,0,2) = cz
          f3p(1,0,0,3) = cz
          f3p(1,0,0,4) = dsqrt(eps3)
c
          f3m(1,0,0,1) = ur
          f3m(1,0,0,2) = cz
          f3m(1,0,0,3) = cz
          f3m(1,0,0,4) = - dsqrt(eps3)
       end if
c
c      REGION II
c
c
c      TM POLARIZATION
c
       do i = 0, n2x
       do j = 0, n2y
          f2p(0,i,j,1) = + vk2x(i) * vk2z(i,j) / dospi
          f2p(0,i,j,2) = + vk2y(j) * vk2z(i,j) / dospi
          f2p(0,i,j,3) = - eps2 * vk2y(j)
          f2p(0,i,j,4) = + eps2 * vk2x(i)
c
          f2m(0,i,j,1) = + vk2x(i) * vk2z(i,j) / dospi
          f2m(0,i,j,2) = + vk2y(j) * vk2z(i,j) / dospi
          f2m(0,i,j,3) = + eps2 * vk2y(j)
          f2m(0,i,j,4) = - eps2 * vk2x(i)
c
       end do
       end do
c
c      TE POLARIZATION
c
       do i = 0, n2x
       do j = 0, n2y
          f2p(1,i,j,1) = + vk2y(j)
          f2p(1,i,j,2) = - vk2x(i)
          f2p(1,i,j,3) = + vk2x(i) * vk2z(i,j) / dospi
          f2p(1,i,j,4) = + vk2y(j) * vk2z(i,j) / dospi
c
          f2m(1,i,j,1) = + vk2y(j)
          f2m(1,i,j,2) = - vk2x(i)
          f2m(1,i,j,3) = - vk2x(i) * vk2z(i,j) / dospi
          f2m(1,i,j,4) = - vk2y(j) * vk2z(i,j) / dospi
c
       end do
       end do
c
       return
       end
