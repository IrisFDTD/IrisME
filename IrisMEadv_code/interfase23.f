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

      subroutine interfase23(ndimk1x, ndimk1y, ndim1, ndimk2x, ndimk2y, 
     &                        ndim2, n1x, n1y, n2x, n2y, f3m, f2p, f2m, 
     &                        sol, nummod1, nummod2,
     &                        numte, im1, im2, amat3, conex32, conex23, 
     &                        wbiga, wa, rhor22, tau23)
c
      implicit real*8(a-h,o-z)
      dimension  im1(ndim1, 0:2), im2(ndim2, 0:2), wa(ndim2)
      complex*16 cz, ur,
     &           f3m(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f2p(0:1, 0:ndimk2x, 0:ndimk2y, 4),
     &           f2m(0:1, 0:ndimk2x, 0:ndimk2y, 4),
     &           amat3(ndim2, ndim2, 3:4), conex32(ndim1, ndim2),
     &           conex23(ndim1, ndim2), wbiga(ndim2, ndim2), 
     &           rhor22(ndim2, ndim2), tau23(ndim1, ndim2),
     &           sol(-ndimk1x:ndimk1x, -ndimk1y:ndimk1y, ndim2, 2)
c
      cz = dcmplx(0.d0, 0.d0)
c
      iequation = 0
c
      do k2eq = 1, numte
        ieq = im2(k2eq,1)
        jeq = im2(k2eq,2)
        do 10 neta = 3, 4
           if(jeq.eq.0 .and. neta.eq.4) goto 10
           if(ieq.eq.0 .and. neta.eq.3) goto 10
           iequation = iequation + 1
           do k2 = 1, nummod2
              wbiga(iequation,k2) = amat3(k2eq,k2,neta)
           end do  
10      continue
      end do
c
c   independent term (one for each k2_0)
c
      do k20 = 1, nummod2
        ipol20 = im2(k20,0)
        i0 = im2(k20,1)
        j0 = im2(k20,2)
        iequation = 0
c
        do k2eq = 1, numte
           ieq = im2(k2eq,1)
           jeq = im2(k2eq,2)
           do 20 neta = 3, 4
              netas = 5 - neta
              if(jeq.eq.0 .and. neta.eq.4) goto 20
              if(ieq.eq.0 .and. neta.eq.3) goto 20
              iequation = iequation + 1
              if (jeq.eq.j0 .and. ieq.eq.i0) then
                 rhor22(iequation, k20) = f2m(ipol20,i0,j0,neta)
              else
                 rhor22(iequation, k20) = cz
              end if
              do k3 = 1, nummod1
                 ipol3 = im1(k3,0)
                 m = im1(k3,1)
                 n = im1(k3,2)
                 rhor22(iequation, k20) = rhor22(iequation, k20) -
     &               dconjg(sol(m,n,k2eq,netas)) * f3m(ipol3,m,n,neta) *
     &               conex23(k3,k20)
              end do
20         continue
        end do
      end do
c
c      reflection 2->2
c
       call LEQT1C(wbiga,nummod2,ndim2,rhor22,nummod2,ndim2,0,WA,IER)
c
c      transmission 2->3
c
       do k20 = 1, nummod2
        do k3 = 1, nummod1
         tau23(k3,k20) = conex23(k3,k20)
         do k2 = 1, nummod2
          tau23(k3,k20) = tau23(k3,k20) + conex32(k3,k2)*rhor22(k2,k20)
         end do
        end do
       end do
c
       return
       end
