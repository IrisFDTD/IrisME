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

      subroutine matrix(ndimk1x, ndimk1y, ndim1, ndimk2x, ndimk2y, 
     &                 ndim2, n1x, n1y, n2x, n2y, f1p, f1m, f2p, f2m,
     &                 f3p, f3m, ipol0, m0, n0, eps1, eps3, nummod1,
     &                 nummod2, numte, im1, im2, sol, v0, den1, den3,
     &                 conex12, conex21, conex32, conex23, amat1, amat3)
c
      implicit real*8(a-h,o-z)
      dimension  im1(ndim1, 0:2), im2(ndim2, 0:2)
      complex*16 cz, v0(0:1),
     &           f1p(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f1m(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f2p(0:1, 0:ndimk2x, 0:ndimk2y, 4),
     &           f2m(0:1, 0:ndimk2x, 0:ndimk2y, 4),
     &           f3p(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f3m(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           sol(-ndimk1x:ndimk1x, -ndimk1y:ndimk1y, ndim2, 2),
     &           den1(-ndimk1x:ndimk1x, -ndimk1y:ndimk1y),
     &           den3(-ndimk1x:ndimk1x, -ndimk1y:ndimk1y),
     &           conex12(ndim1, ndim2), conex21(ndim1, ndim2),
     &           conex32(ndim1, ndim2), conex23(ndim1, ndim2),
     &           amat1(ndim2, ndim2, 3:4), amat3(ndim2, ndim2, 3:4)
c
       cz = dcmplx(0.d0, 0.d0)
c
       do m = -n1x, n1x
       do n = -n1y, n1y
          den1(m,n) = f1m(0,m,n,1) * f1m(1,m,n,2) -
     &                f1m(0,m,n,2) * f1m(1,m,n,1)
       end do
       end do
c
       v0(0) = (f1p(ipol0,m0,n0,1) * f1m(1,m0,n0,2) -
     &         f1p(ipol0,m0,n0,2) * f1m(1,m0,n0,1)) / den1(m0,n0)
       v0(1) = (f1p(ipol0,m0,n0,1) * f1m(0,m0,n0,2) -
     &         f1p(ipol0,m0,n0,2) * f1m(0,m0,n0,1)) / (-den1(m0,n0))
c
       do k1 = 1, nummod1
        ipol1 = im1(k1,0)
        nipol1 = 1 - ipol1
        m = im1(k1,1)
        n = im1(k1,2)
        do k2 = 1, nummod2
c
         ipol2 = im2(k2,0)
         i = im2(k2,1)
         j = im2(k2,2)
         conex12(k1,k2) = ( f2p(ipol2,i,j,1) * sol(m,n,k2,1) * 
     &                      f1m(nipol1,m,n,2) -
     &         f2p(ipol2,i,j,2) * sol(m,n,k2,2) * f1m(nipol1,m,n,1) )
     &         / ( (-1)**ipol1 * den1(m,n) )
         conex21(k1,k2) = ( f2m(ipol2,i,j,1) * sol(m,n,k2,1) * 
     &                      f1m(nipol1,m,n,2) -
     &         f2m(ipol2,i,j,2) * sol(m,n,k2,2) * f1m(nipol1,m,n,1) )
     &         / ( (-1)**ipol1 * den1(m,n) )
        end do
       end do
c
       do k2 = 1, nummod2
       do k2eq = 1, numte
       do neta = 3, 4
        netas = 5 - neta
        if(im2(k2,1).eq.im2(k2eq,1) .and. im2(k2,2).eq.im2(k2eq,2)) then
         amat1(k2eq,k2,neta)= - f2p(im2(k2,0),im2(k2,1),im2(k2,2),neta)
        else
         amat1(k2eq,k2,neta) = cz
        end if
        do k1 = 1, nummod1
          ipol1 = im1(k1,0)
          m = im1(k1,1)
          n = im1(k1,2)
c
c defino amat1 para (i,j,neta) donde i,j son los correspondientes
c a los del modo TE, aun sabiendo que (i,0) y (0,j) solo generan una
c ecuacion (neta=3 o 4) y no dos. Esto se tiene en cuenta al generar
c la ecuacion matricial que genera rho, tau. 
c
          amat1(k2eq,k2,neta) = amat1(k2eq,k2,neta) + 
     &                conex12(k1,k2) * dconjg(sol(m,n,k2eq,netas)) * 
     &                f1m(ipol1,m,n,neta)
        end do
       end do
       end do
       end do
c
c      matrixes related to interfase 2-3
c
       if (eps3 .eq. eps1) then
        do m = -n1x, n1x
        do n = -n1y, n1y
           den3(m,n) = den1(m,n)
        end do
        end do
        do k3 = 1, nummod1
        do k2 = 1, nummod2
           conex23(k3,k2) = conex21(k3,k2)
           conex32(k3,k2) = conex12(k3,k2)
        end do
        end do
        do k2 = 1, nummod2
        do k2eq = 1, numte
        do neta = 3, 4
           amat3(k2eq,k2,neta) = amat1(k2eq,k2,neta) 
        end do
        end do
        end do
       else
c
        do m = -n1x, n1x
        do n = -n1y, n1y
           den3(m,n) = f3m(0,m,n,1) * f3m(1,m,n,2) -
     &                 f3m(0,m,n,2) * f3m(1,m,n,1)
        end do
        end do
c
        do k3 = 1, nummod1
         ipol3 = im1(k3,0)
         nipol3 = 1 - ipol3
         m = im1(k3,1)
         n = im1(k3,2)
         do k2 = 1, nummod2
          ipol2 = im2(k2,0)
          i = im2(k2,1)
          j = im2(k2,2)
          conex32(k3,k2) = 
     &     ( f2p(ipol2,i,j,1) * sol(m,n,k2,1) * f3m(nipol3,m,n,2) -
     &       f2p(ipol2,i,j,2) * sol(m,n,k2,2) * f3m(nipol3,m,n,1) )
     &       / ( (-1)**ipol3 * den3(m,n) )
          conex23(k3,k2) =  ( f2m(ipol2,i,j,1) * sol(m,n,k2,1) * 
     &                        f3m(nipol3,m,n,2) -
     &       f2m(ipol2,i,j,2) * sol(m,n,k2,2) * f3m(nipol3,m,n,1) )
     &       / ( (-1)**ipol3 * den3(m,n) )
        end do
        end do
c
        do k2 = 1, nummod2
        do k2eq = 1, numte
        do neta = 3, 4
           netas = 5 - neta
        if(im2(k2,1).eq.im2(k2eq,1) .and. im2(k2,2).eq.im2(k2eq,2)) then
         amat3(k2eq,k2,neta)= - f2p(im2(k2,0),im2(k2,1),im2(k2,2),neta)
        else
         amat3(k2eq,k2,neta) = cz
        end if
           do k3 = 1, nummod1
              ipol3 = im1(k3,0)
              m = im1(k3,1)
              n = im1(k3,2)
              amat3(k2eq,k2,neta) = amat3(k2eq,k2,neta) + 
     &                    conex32(k3,k2) * dconjg(sol(m,n,k2eq,netas)) * 
     &                    f3m(ipol3,m,n,neta)
             end do
          end do
        end do
        end do
       end if
c
       return
       end
