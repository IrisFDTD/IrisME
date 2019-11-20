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

      subroutine redef_field(ndimk1x, ndimk1y, ndimk2x, ndimk2y, n1x, 
     &                        n1y, n2x, n2y, zs, f1p, f1m, f2p, f2m, 
     &                        f3p, f3m)
c
      implicit real*8(a-h,o-z)
      complex*16 f1p(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f1m(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f2p(0:1, 0:ndimk2x, 0:ndimk2y, 4),
     &           f2m(0:1, 0:ndimk2x, 0:ndimk2y, 4),
     &           f3p(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f3m(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4), zs
c
       do ipol = 0, 1
       do m = -n1x, n1x
       do n = -n1y, n1y
          f1p(ipol,m,n,1) = f1p(ipol,m,n,1) - zs * f1p(ipol,m,n,4)
          f1p(ipol,m,n,2) = f1p(ipol,m,n,2) + zs * f1p(ipol,m,n,3)
c
          f1m(ipol,m,n,1) = f1m(ipol,m,n,1) - zs * f1m(ipol,m,n,4)
          f1m(ipol,m,n,2) = f1m(ipol,m,n,2) + zs * f1m(ipol,m,n,3)
c
          f3p(ipol,m,n,1) = f3p(ipol,m,n,1) - zs * f3p(ipol,m,n,4)
          f3p(ipol,m,n,2) = f3p(ipol,m,n,2) + zs * f3p(ipol,m,n,3)
c
          f3m(ipol,m,n,1) = f3m(ipol,m,n,1) - zs * f3m(ipol,m,n,4)
          f3m(ipol,m,n,2) = f3m(ipol,m,n,2) + zs * f3m(ipol,m,n,3)
c
       end do
       end do
       end do
c
c      REGION II
c
       do ipol = 0, 1
       do i = 0, n2x
       do j = 0, n2y
          f2p(ipol,i,j,1) = f2p(ipol,i,j,1) - zs * f2p(ipol,i,j,4)
          f2p(ipol,i,j,2) = f2p(ipol,i,j,2) + zs * f2p(ipol,i,j,3)
c
          f2m(ipol,i,j,1) = f2m(ipol,i,j,1) - zs * f2m(ipol,i,j,4)
          f2m(ipol,i,j,2) = f2m(ipol,i,j,2) + zs * f2m(ipol,i,j,3)
c
       end do
       end do
       end do
c
       return
       end
