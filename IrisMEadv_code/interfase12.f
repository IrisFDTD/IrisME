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

      
      subroutine interfase12(ndimk1x, ndimk1y, ndim1, ndimk2x,ndimk2y,
     &                        ndim2, n1x, n1y, n2x, n2y, f1p, f1m, f2p,
     &                        sol, ipol0, m0, n0, v0, nummod1, nummod2,
     &                        numte, ip1, im2, amat1, conex12,
     &                        wbiga, wa, rho11, tau12)
c
      implicit real*8(a-h,o-z)
      dimension  ip1(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y),
     &           im2(ndim2, 0:2), wa(ndim2)
      complex*16 cz, b0(3:4), v0(0:1), 
     &           f1p(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f1m(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f2p(0:1, 0:ndimk2x, 0:ndimk2y, 4),
     &           f2m(0:1, 0:ndimk2x, 0:ndimk2y, 4),
     &           sol(-ndimk1x:ndimk1x, -ndimk1y:ndimk1y, ndim2, 2),
     &           amat1(ndim2, ndim2, 3:4), conex12(ndim1, ndim2),
     &           wbiga(ndim2, ndim2), rho11(ndim1), tau12(ndim2)
c
       cz = dcmplx(0.d0, 0.d0)
c
       do neta = 3, 4
        b0(neta) = - f1p(ipol0,m0,n0,neta) + v0(0)*f1m(0,m0,n0,neta) +
     &             v0(1)*f1m(1,m0,n0,neta)
       end do
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
          netas = 5 - neta
c
c   independent term
c
          tau12(iequation) = dconjg(sol(m0,n0,k2eq,netas)) * b0(neta)
c
c   matrix term
c
          do k2 = 1, nummod2
             wbiga(iequation,k2) = amat1(k2eq,k2,neta)
          end do   
10       continue
       end do
c
       nequation = iequation
       if(nequation .ne. nummod2) then
          write(6,*)'faltan modos y/o ecuaciones'
          write(6,*)' # ecuaciones = ', nequation
          write(6,*)' # modos = ', nummod2
       end if
c
c
c      transmission 1->2
c
       call LEQT1C(wbiga,nummod2,ndim2,tau12,1,ndim2,0,WA,IER)
c
c      reflection 1->1
c
       do k1 = 1, nummod1
          rho11(k1) = cz
       end do
       do ipol1 = 0, 1
          rho11(ip1(ipol1,m0,n0)) = - v0(ipol1)
       end do
       do k1 = 1, nummod1
          do k2 = 1, nummod2
             rho11(k1) = rho11(k1) + conex12(k1, k2) * tau12(k2) 
          end do
       end do
c
       return
       end
