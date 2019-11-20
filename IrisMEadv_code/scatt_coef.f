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

      subroutine scatt_coef(ndimk1x, ndimk1y, ndim1, ndimk2x, ndimk2y,
     &                       ndim2, n1x, n1y, n2x, n2y, nummod1, 
     &                       nummod2, fackz, wunit, wresden, wa, rho11, 
     &                       tau12, rhol22, tau21, rhor22, tau23, 
     &                       trans, trans2, refle)
c
      implicit real*8(a-h,o-z)
      dimension  wa(ndim2)
      complex*16 cz, ur, fackz(ndim2), 
     &           wunit(ndim2, ndim2), wide(ndim2,ndim2), 
     &           wresden(ndim2, ndim2), trans2(ndim1), 
     &           rho11(ndim1), tau12(ndim2), rhol22(ndim2, ndim2), 
     &           tau21(ndim1, ndim2), rhor22(ndim2, ndim2), 
     &           tau23(ndim1, ndim2), trans(ndim1), refle(ndim1)
c
       cz = dcmplx(0.d0, 0.d0)
       ur = dcmplx(1.d0, 0.d0)
c
       do k2 = 1, nummod2
       do l2 = 1, nummod2
         wunit(k2,l2) = cz
         wide(k2,l2) = cz
       end do
       end do
       do k2 = 1, nummod2
         wunit(k2,k2) = ur
         wide(k2,k2) = ur
       end do
c
       do k21 = 1, nummod2
       do k23 = 1, nummod2
         wresden(k23,k21) = cz
         do k22 = 1, nummod2
          wresden(k23,k21) = wresden(k23,k21) + 
     &                     rhol22(k23,k22)*fackz(k22)*rhor22(k22,k21)
         end do
         wresden(k23,k21) = wunit(k23,k21) - wresden(k23,k21)*fackz(k23)
       end do
       end do
c
       call LEQT1C(wresden,nummod2,ndim2,wunit,nummod2,ndim2,0,WA,IER) 
c
c      transmssion coefficients (1->3)
c
        do k1 = 1, nummod1
         trans(k1) = cz
         trans2(k1) = cz
         do k21 = 1, nummod2
         do k22 = 1, nummod2
            trans(k1) = trans(k1) + tau23(k1,k22) * wunit(k22,k21)
     &                              * fackz(k21) * tau12(k21)
            trans2(k1) = trans2(k1) + tau23(k1,k22) * wide(k22,k21)
     &                              * fackz(k21) * tau12(k21)
         end do
         end do
        end do
c
c      reflection coefficients (1->1)
c
        do k1 = 1, nummod1
         refle(k1) = rho11(k1)
         do k21 = 1, nummod2
         do k22 = 1, nummod2
         do k23 = 1, nummod2
            refle(k1) = refle(k1) + tau21(k1,k23) * wunit(k23,k22) *
     &         fackz(k22) * rhor22(k22,k21) * fackz(k21) * tau12(k21)
         end do
         end do
         end do
        end do
c
       return
       end
