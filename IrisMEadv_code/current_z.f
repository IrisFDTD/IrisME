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
        

      
      subroutine current_z(ndimk1x, ndimk1y, ndim1, n1x, n1y, 
     &                     vk1z, vk3z, f1p, f1m, f3p, ip1, ipol0, 
     &                     m0, n0, refle,
     &                     trans, trans2, refcurr, refcurr0s, refcurr0p,
     &                     tracurr, tracurr0s, tracurr0p, tranum,tl,al)
c
       implicit real*8(a-h,o-z)
       dimension ip1(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y)
       complex*16 vk1z(-ndimk1x:ndimk1x,-ndimk1y:ndimk1y),
     &           vk3z(-ndimk1x:ndimk1x,-ndimk1y:ndimk1y),
     &           f1p(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f1m(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f3p(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           trans(ndim1), trans2(ndim1), refle(ndim1),tl
c



       curr0 = dreal(f1p(ipol0,m0,n0,1) * f1p(ipol0,m0,n0,4) -
     &         f1p(ipol0,m0,n0,2) * f1p(ipol0,m0,n0,3))
c
c     reflexion
c
      refcurr = 0.d0
      do ipol =0,1
      do m = -n1x, n1x
      do n = -n1y, n1y
       if( dimag(vk1z(m,n)).eq.0.d0) then
         refcurr = refcurr + dreal(
     &                       f1m(ipol,m,n,1)*f1m(ipol,m,n,4) -
     &                       f1m(ipol,m,n,2) * f1m(ipol,m,n,3))*
     &                       cdabs(refle(ip1(ipol,m,n)))**2
       end if
      end do
      end do
      end do
      refcurr = refcurr/curr0
c
      refcurr0p = dreal( f1m(1,0,0,1)*f1m(1,0,0,4) -
     &                   f1m(1,0,0,2) * f1m(1,0,0,3))*
     &                   cdabs(refle(ip1(1,0,0)))**2 / curr0
      refcurr0s = dreal( f1m(0,0,0,1)*f1m(0,0,0,4) -
     &                   f1m(0,0,0,2) * f1m(0,0,0,3))*
     &                   cdabs(refle(ip1(0,0,0)))**2 / curr0
c
c     transmission
c
      tracurr = 0.d0	

      do ipol =0,1	
      do m = -n1x, n1x
      do n = -n1y, n1y
       if( dimag(vk3z(m,n)).eq.0.d0) then
         tracurr = tracurr + dreal(
     &                       f3p(ipol,m,n,1)*f3p(ipol,m,n,4) -
     &                       f3p(ipol,m,n,2) * f3p(ipol,m,n,3))*
     &                       cdabs(trans(ip1(ipol,m,n)))**2 / curr0
       end if
      end do
      end do
      end do
c
      m=0
	n=0
	ipol=0  !>>> p-polarization

      tracurr0p = dreal(f3p(ipol,m,n,1)*f3p(ipol,m,n,4) -
     &                       f3p(ipol,m,n,2) * f3p(ipol,m,n,3))*
     &                       cdabs(trans(ip1(ipol,m,n)))**2 / curr0

      tranum = dreal( f3p(1,0,0,1)*f3p(1,0,0,4) -
     &                   f3p(1,0,0,2) * f3p(1,0,0,3))*
     &                   cdabs(trans2(ip1(1,0,0)))**2 / curr0

      tracurr0s = dreal( f3p(0,0,0,1)*f3p(0,0,0,4) -
     &                   f3p(0,0,0,2) * f3p(0,0,0,3))*
     &                   cdabs(trans(ip1(0,0,0)))**2 / curr0
c
      return
      end

         
   



