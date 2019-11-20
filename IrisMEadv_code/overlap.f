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

      
      subroutine overlap(ndimk1x, ndimk1y, ndimk2x, ndimk2y, ndim2,
     &                   n1x, n1y, n2x, n2y, rlx, rly, ra, rb, nummod2,
     &                   ip2, im2, vk1x, vk1y, vk2x, vk2y, vsx, 
     &                   vcx, vsy, vcy, sol)  
c
       implicit real*8(a-h,o-z)
       dimension vk1x(-ndimk1x:ndimk1x), vk1y(-ndimk1y:ndimk1y),
     &           vk2x(0:ndimk2x), vk2y(0:ndimk2y),
     &           ip2(0:1, 0:ndimk2x, 0:ndimk2y),
     &           im2(ndim2, 0:2)
       complex*16 cz, cosint, sinint,
     &            vsx(-ndimk1x:ndimk1x, 0:ndimk2x), 
     &            vcx(-ndimk1x:ndimk1x, 0:ndimk2x),
     &            vsy(-ndimk1y:ndimk1y, 0:ndimk2y),
     &            vcy(-ndimk1y:ndimk1y, 0:ndimk2y),
     &            sol(-ndimk1x:ndimk1x, -ndimk1y:ndimk1y, ndim2, 2)
c
       cz = dcmplx(0.d0, 0.d0)
c
       cnx = dsqrt(2.d0/ra/rlx)
       cny = dsqrt(2.d0/rb/rly)
       cnx0 = dsqrt(1.d0/ra/rlx)
       cny0 = dsqrt(1.d0/rb/rly)
c
       do m = -n1x, n1x
          vsx(m,0) = 0.d0
          vcx(m,0) = cnx0 * cosint(vk1x(m), vk2x(0), ra)
          do i = 1, n2x
             vsx(m,i) = cnx * sinint(vk1x(m), vk2x(i), ra)
             vcx(m,i) = cnx * cosint(vk1x(m), vk2x(i), ra)
          end do
       end do
       do n = -n1y, n1y
          vsy(n,0) = 0.d0
          vcy(n,0) = cny0 * cosint(vk1y(n), vk2y(0), rb)
          do j = 1, n2y
             vsy(n,j) = cny * sinint(vk1y(n), vk2y(j), rb)
             vcy(n,j) = cny * cosint(vk1y(n), vk2y(j), rb)
          end do
       end do
c
       do m = -n1x, n1x
       do n = -n1y, n1y
       do k = 1, nummod2
          sol(m,n,k,1) = vcx(m,im2(k,1)) * vsy(n,im2(k,2))
          if(cdabs(sol(m,n,k,1)).lt.1.d-15) sol(m,n,k,1)=cz
          sol(m,n,k,2) = vsx(m,im2(k,1)) * vcy(n,im2(k,2))
          if(cdabs(sol(m,n,k,2)).lt.1.d-15) sol(m,n,k,2)=cz
       end do
       end do
       end do
c
       return
       end
c
       function cosint(vk, vg, x)
         implicit real*8(a-h,o-z)
         complex*16 ui, cosint, cx
c
         ui = dcmplx(0.d0, 1.d0)
         cx = ui * x
         if(vk. eq. vg. or. vk. eq. -vg) then
          if(vk.eq.0.d0) then
           cosint = x
          else
           cosint = 0.5*x - (cdexp(-2.d0*vk*cx)-1.d0)/(4*vk*ui)
          end if
         else
           cosint = 0.5d0*ui*( (cdexp(-(vk-vg)*cx)-1)/(vk-vg) +
     &                         (cdexp(-(vk+vg)*cx)-1)/(vk+vg) )
         end if
c
         return
         end
c
       function sinint(vk, vg, x)
         implicit real*8(a-h,o-z)
         complex*16 ui, aid, sinint, cx
c
         ui = dcmplx(0.d0, 1.d0)
         cx = ui * x
         if(vk. ne. vg. and. vk. ne. -vg) then
            sinint = 0.5d0 * ( + (cdexp(-(vk-vg)*cx)-1.d0)/(vk-vg) 
     &                         - (cdexp(-(vk+vg)*cx)-1.d0)/(vk+vg) )
         else
           if(vk.eq.0.d0) then
             sinint = 0.d0
           else
             aid = - 0.5d0*ui*x - (cdexp(-2.d0*cx*vk)-1.d0)/(4.d0*vk)
             if(vk. eq. vg) then
               sinint = aid
             else
               sinint = - aid
             end if
           end if
         end if
c
         return
         end
