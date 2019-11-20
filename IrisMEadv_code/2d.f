! -------------------------------------------------------------------------------------------------------------------------------------
! This is the main file of the IrisMEadv program implementation. IriMEadv is a fast and accurate program for 
! the calculation of scattering coefficients in arrays of rectangular holes in metals. 
! IrisME is based on the Modal Expansion (ME) of the electromagnetic fields. Two approaches are implemented: 
! 1) "Multiple Scattering" formalism.
! 2) "Amplitudes of the Electric Field at the Openings" formalism.
!
! IrisME is licensed under the AGPL and it is free to use. However, if you are using, or plan to use it, 
! specially if it is for research or academic purposes, please send an email with your name, institution and a 
! brief description of your interest for this program. If you use IrisME in a work that leads to a scientific 
! for academic publication, we would appreciate it if you would kindly cite IrisME in your manuscript as:
!      
! > L. Martín-Moreno, F.J. García-Vidal, Minimal model for optical transmission through holey metal films, 
! > J. Phys.: Condens. Matter 20, 304214 (2008). (minimal model that includes a full description of both formalisms).
! > and/or
! > F.J. Garcia-Vidal, L. Martin-Moreno, T.W. Ebbesen, L. Kuipers, Light passing through subwavelength apertures, 
! > Rev. Mod. Phys. 82, 729–787 (2010). (more general approach to ME, including other systems...) 
! > and/or
! > S. Carretero-Palacios, F. J. García-Vidal, L. Martín-Moreno, and S. G. Rodrigo,
! > Effect of film thickness and dielectric environment on optical transmission through subwavelength holes, 
! > Phys. Rev. B 85, 035417 (2012). 
!
! Commercial applications may also acquire a commercial license. Please contact sergut@unizar.es & lmm@unizar.es, for details.
! 
! Copyright (C) 2019 Sergio G Rodrigo sergut@unizar.es & Luis Martin-Moreno lmm@unizar.es
!
! IrisME is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public 
! License as published by the Free Software Foundation, either version 3 of the License,or (at your option) any later version.
!
! IrisME is distributed in the hope that it will be useful for research or/and academic purpouses, but WITHOUT ANY WARRANTY; 
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
! Affero General Public License for more details. You should have received a copy of the GNU Affero General Public License 
! along with IrisME. If not,see http://www.gnu.org/licenses/.
! -------------------------------------------------------------------------------------------------------------------------------------
   
     
      
      parameter(ndimk1x=30, ndimk1y=30, ndimk2x=10, ndimk2y=10,
     &           ndim1= 2*(2*ndimk1x+1)*(2*ndimk1y+1),
     &           ndim2= 2*(ndimk2x+1)*(ndimk2y+1))


c
      implicit real*8 (a-h,o-z)
c
      dimension  ip1(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y),
     &           im1(ndim1, 0:2), im2(ndim2, 0:2),
     &           ip2(0:1, 0:ndimk2x, 0:ndimk2y)
c
      dimension  vk1x(-ndimk1x:ndimk1x), vk1y(-ndimk1y:ndimk1y),
     &           vk3x(-ndimk1x:ndimk1x), vk3y(-ndimk1y:ndimk1y),
     &           vk2x(0:ndimk2x), vk2y(0:ndimk2y), wa(ndim2),
     &           h(1:9)
c
      complex*16 ui, eps, epsmetal,zs,zsp,zsm,  
     &           vk1z(-ndimk1x:ndimk1x,-ndimk1y:ndimk1y),
     &           vk2z(0:ndimk2x,0:ndimk2y), fackz(ndim2),
     &           vk3z(-ndimk1x:ndimk1x,-ndimk1y:ndimk1y),
     &           f1p(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f1m(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f2p(0:1, 0:ndimk2x, 0:ndimk2y, 4),
     &           f2m(0:1, 0:ndimk2x, 0:ndimk2y, 4),
     &           f3p(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           f3m(0:1, -ndimk1x:ndimk1x, -ndimk1y:ndimk1y, 4),
     &           vsx(-ndimk1x:ndimk1x, 0:ndimk2x),
     &           vcx(-ndimk1x:ndimk1x, 0:ndimk2x),
     &           vsy(-ndimk1y:ndimk1y, 0:ndimk2y),
     &           vcy(-ndimk1y:ndimk1y, 0:ndimk2y),
     &           sol(-ndimk1x:ndimk1x, -ndimk1y:ndimk1y, ndim2, 2),
     &           den1(-ndimk1x:ndimk1x, -ndimk1y:ndimk1y),
     &           den3(-ndimk1x:ndimk1x, -ndimk1y:ndimk1y),
     &           conex12(ndim1, ndim2), conex21(ndim1, ndim2),
     &           conex32(ndim1, ndim2), conex23(ndim1, ndim2),
     &           amat1(ndim2, ndim2, 3:4), amat3(ndim2, ndim2, 3:4),
     &           v0(0:1), wbiga(ndim2, ndim2), wunit(ndim2, ndim2),
     &           rho11(ndim1), tau12(ndim2), rhol22(ndim2, ndim2), 
     &           tau21(ndim1, ndim2), rhor22(ndim2, ndim2), 
     &           tau23(ndim1, ndim2), trans(ndim1), 
     &           trans2(ndim1), refle(ndim1),denden
     &        ,g1,g2,g3,greenI,greenIII,phase,sigma,gv,expo,phaserho
     &        ,gva,sigmaa,Io,g,gc,con_max,con_max1,Taprx,Taprx2,Tanl,
     &         So,illu,e0p,e0,e0pa,e0a,e0pa2,e0a2
     &        ,gmx,Tp,con_max_valid,con_max3,phase3,phaseGsigma,con_max4
     &       ,Y2,con_maxA,gva2,sigmaa2,den,dena,dena2,Too,Tanl2,fa,fb,fc
     &       ,Gsum,Gdif,
     &        denWITHminus,denOUTminus,con_maxDIF,con_maxDIFext,
     &        con_luis,dPLUS,dMINUS,aPLUS,aMINUS,bden,rho1,rhoeps
      
 !>>> Cálculo too para transmisión directa a través de la lámina 
	complex*16, dimension(2) :: e0pfac,e0fac
	real*8,    dimension(11) :: transfac

	integer SIBC,sim,nosim,loop

	real   w,wh,
     &       nr,ni,   !n´=nr+ini
     &       nII,nIII,
     &       del,epsI,epsII,epsIII,delsim,transFULL,transaFULL,refFULL
      
	complex*16 cn,cnp,
     &           rl,tl,bt
	character(5)::charno
c
       logical, parameter :: flag_open=.true.
	 parameter(inME=5)
	 parameter(transFile=19)
	 ndimk1xs = ndimk1x
       ndimk1ys = ndimk1y
       ndimk2xs = ndimk2x
       ndimk2ys = ndimk2y
       ndim1s = ndim1
       ndim2s = ndim2	
c
c values for constants
c   
       pi = 4.d0*datan(1.0d0)
       clight = 0.1957d0
       ui = dcmplx(0.d0, 1.d0)
c
c      ! To compare with IrisME.nb Mathematica script
       ! Minimal Model Mathematica No:=-1; Nf:=1; Mo:= 0; Mf:=0;
       ! Minimal Model Fortran n2x,n2y (0,1,1,0)===>Minimal model (n1x,n1y->I y III;n2x,n2y->II)
       
       print*,"...Loading input"
	 open(unit=inME,file='input_IrisMEadv.dat',status='old',action='read')
	 read(inME,*) rlx0, rly0, ax0, ay0, width0
       read(inME,*) eps3, theta, fi, eps1, eps2        
       read(inME,*) vlan0, delvlan, nvlan
       read(inME,*) n1x, n1y, n2x, n2y      
       read(inME,*) cutoff    
	 read(inME,*) SIBC   	
	 read(inME,*) loop
	 read(inME,*) sim,deltal,nosim   
            
      !>>> Loop calculation 	
       !>>> Redefine parameters
	 epsI=eps1; epsII=eps2; epsIII=eps3
	 Lx = rlx0; Ly = rly0	      
	 ax=ax0; ay=ay0; widthaux=width0	 
       !<<< Redefine parameters          
       
	 do w =1,nosim
	 select case (loop)
	 case (1)  !>>> Running the program for several substrate dielectric constants
	   eps3 = epsIII	 
         eps3 = eps3+deltal*(w-1)
	   if(sim) eps1   = eps3
	 case (2)  !>>> Running the program for several hole dielectric constants
	   eps2 = epsII	 
         eps2 = eps2+deltal*(w-1)	   
	 case (3)  !>>> Running the program for several hole sizes (ax changes)
	   ax0 = ax
	   ax0 = ax0+deltal*(w-1)
	 case (4)  !>>> Running the program for several hole sizes (ay changes)
	   ay0 = ay
         ay0 = ay0+deltal*(w-1)
       case (5)  !>>> Running the program for several film widths
	   width0=widthaux
	   width0=width0+deltal*(w-1)
	 case default 
         !....
       end select      

	call inttoASCII(int(w),charno)
	if(flag_open)then
      if(nosim/=1)then
      open(unit=15,file=trim("T-R_IrisMEadv_step"//charno)//".dat",
     & action='write')
c	open(unit=16,file=trim("Tanalytic_step"//charno)//".dat",action='write')
!	    open(unit=17,file=trim("Sigma_step"//charno)//".dat",action='write')
!	    open(unit=18,file=trim("Gv_step"//charno)//".dat",action='write')
!c	open(unit=19,file=trim("T_step"//charno)//".dat",action='write')
!c	open(unit=20,file=trim("den_fact="//charno)//".dat",action='write')	
!	    open(unit=21,file=trim("con_step"//charno)//".dat",action='write')
!c	open(unit=22,file=trim("den_step"//charno)//".dat",action='write')	
!	    open(unit=24,file=trim("GreenI_step"//charno)//".dat",action='write')
!c     open(unit=25,file=trim("GreenIII_step"//charno)//".dat",action='write')
!	    open(unit=26,file=trim("phase_step"//charno)//".dat",action='write')
!	    open(unit=27,file=trim("expo_step"//charno)//".dat",action='write')
!	    open(unit=28,file=trim("Ep_step"//charno)//".dat",action='write')
!	    open(unit=29,file=trim("E_step"//charno)//".dat",action='write')
!	    open(unit=30,file=trim("T_for_step"//charno)//".dat",action='write')
!	    open(unit=31,file=trim("qzh_step"//charno)//".dat",action='write')
!	    open(unit=32,file=trim("T_den_step"//charno)//".dat",action='write')
      else
          open(unit=15,file="T-R_IrisMEadv.dat",action='write') 	
      end if
      end if	
      
c	 open(unit=19,file=trim("trans_h="//charno)//".dat",action='write')
c   
c      incidence in P-polarization, mode m=0, n=0
c
       ipol0 = 0
       m0 = 0 
       n0 = 0 
**       eps1 = 1.d0
c 
c      LOOP IN WAVELENGTHS
c
       ipaso=0
       do ivlan = 0, nvlan
       ipaso=ipaso+1
          vlan = vlan0 + delvlan*ivlan
          gvaccuum = 2.d0 * pi / vlan
c          ew = clight * gvaccuum
          ew = 1240.d0 / vlan
          ewh=ew/27.2d0

       if(SIBC/=0) then   	
     	   !>>> Drude-Lorentz model for metals from S. G. Rodrigo, F. J. García-Vidal, and L. Martín-Moreno, 
         !    “Influence of material properties on extraordinary optical transmission through hole arrays," 
         !    Phys. Rev. B 77, 075401 (2008).
       select case(SIBC) !The metal chosen is coming...
       case(1)
          !Ag 
          wplasmon = 9.0d0/27.2d0
	    gamma    = 0.07d0/27.2d0
	    el       = 1.10d0 
	    wpl      = 4.9d0/27.2d0
		gammal   = 1.2d0/54.4d0
		 
          eps = 4.6d0 - wplasmon**2/(ewh*(ewh+ui*gamma))
	    eps = eps + (el*wpl**2)/(wpl**2-2*ui*ewh*gammal-ewh**2)
       case(2)
        ! Au
	    wplasmon = 8.729/27.2d0
	    gamma    = 0.065d0/27.2d0
	    el       = 1.09d0 
	    wpl      = 2.684d0/27.2d0
		gammal   = 0.433d0/54.4d0
		 
          eps = 5.967d0 - wplasmon**2/(ewh*(ewh+ui*gamma))
	    eps = eps + (el*wpl**2)/(wpl**2-2*ui*ewh*gammal-ewh**2)      
       end select
	   !<<< Drude-Lorentz model for metals   
          
	   	zs = 1.d0/cdsqrt(eps+1.0d0)   !>>> SIBC conditios <<<!
		dr=dreal(eps)
          di=dimag(eps)	
      else 
          !>>> PEC 
          zs = 0.d0                                  
      end if
c
c      ALL LENGTHS IN UNIS OF LANDA
c
       !delta = vlan/pi/dsqrt(-dr)
       delta = 0.d0
       write(6,*)vlan, 'zs ', zs, delta
c
       rlx = rlx0 / vlan
       rly = rly0 / vlan

	!>>> Ensanchamiento efectivo de las dimensiones laterales del agujero
       ra = ( ax0 + 1.2d0*delta ) / vlan
       rb = ( ay0 + 1.2d0*delta ) / vlan
       al = (rlx*rly - ra*rb)/rlx*rly
      !    Ensanchamiento efectivo de las dimensiones laterales del agujero <<<

       width = width0 / vlan
       gvaccuum = 2.d0 * pi

!---------------------------------!
! Multiple Scattering" formalism  !
!---------------------------------!
       
c wavevectors in all three regions
c
      call wavevec(ndimk1xs, ndimk1ys, n1x, n1y, ndimk2xs, ndimk2ys,
     &             n2x, n2y, rlx, rly, ra, rb, theta, fi, gvaccuum,eps1,
     &             eps2, eps3, vk1x, vk1y, vk1z, vk2x, vk2y, vk2z, vk3x,
     &             vk3y, vk3z)
c
c expression for modes in all 3 regions
c
      call modes(ndimk1xs, ndimk1ys, ndimk2xs, ndimk2ys, n1x, n1y,
     &           n2x, n2y, vk1x, vk1y, vk1z, vk2x, vk2y, vk2z, vk3x,
     &           vk3y, vk3z, theta, gvaccuum, eps1, eps2, eps3,
     &           f1p, f1m, f2p, f2m, f3p, f3m)
c
c Redifining F, f, for Surface Impedance Boundary Conditions.
c
      call redef_field(ndimk1xs, ndimk1ys, ndimk2xs, ndimk2ys, n1x, n1y,
     &                  n2x, n2y, zs, f1p, f1m, f2p, f2m, f3p, f3m)
c
c assigning indexes to modes in regions I and II
c
      call index1(ndimk1xs, ndimk1ys, ndim1s, n1x, n1y, ip1,
     &            im1, nummod1)
c
      call index2(ndimk2xs, ndimk2ys, ndim2s, n2x, n2y, width,
     &            cutoff, vk2z, ip2, im2, fackz, nummod2, numte)
c       write(6,*)'nummod2 ', nummod2
c
c
c computation of overlaps between free space and waveguide modes.
c
      call overlap(ndimk1xs, ndimk1ys, ndimk2xs, ndimk2ys, ndim2s,
     &             n1x, n1y, n2x, n2y, rlx, rly, ra, rb, nummod2,
     &             ip2, im2, vk1x, vk1y, vk2x, vk2y, vsx,
     &             vcx, vsy, vcy, sol)
c  
c auxiliary functions
c
      call matrix(ndimk1xs, ndimk1ys, ndim1s, ndimk2xs, ndimk2ys,
     &            ndim2s, n1x, n1y, n2x, n2y, f1p, f1m, f2p, f2m,
     &            f3p, f3m, ipol0, m0, n0, eps1, eps3, nummod1,
     &            nummod2, numte, im1, im2, sol, v0, den1, den3,
     &            conex12, conex21, conex32, conex23, amat1, amat3)
c
      call interfase12(ndimk1xs, ndimk1ys, ndim1s, ndimk2xs, ndimk2ys,
     &                 ndim2s, n1x, n1y, n2x, n2y, f1p, f1m, f2p,
     &                 sol, ipol0, m0, n0, v0, nummod1, nummod2,
     &                 numte, ip1, im2, amat1, conex12, wbiga, wa, 
     &                 rho11, tau12)
c
      call interfase21(ndimk1xs, ndimk1ys, ndim1s, ndimk2xs, ndimk2ys,
     &                ndim2s, n1x, n1y, n2x, n2y, f1m, f2p, f2m, sol,
     &                nummod1, nummod2, numte, im1, im2, amat1, conex12,
     &                conex21, wbiga, wa, rhol22, tau21)
c
c
      call interfase23(ndimk1xs, ndimk1ys, ndim1s, ndimk2xs, ndimk2ys,
     &                ndim2s, n1x, n1y, n2x, n2y, f3m, f2p, f2m, sol,
     &                nummod1, nummod2, numte, im1, im2, amat3, conex32,
     &                conex23, wbiga, wa, rhor22, tau23)
c
      call scatt_coef(ndimk1xs, ndimk1ys, ndim1s, ndimk2xs, ndimk2ys,
     &                ndim2s, n1x, n1y, n2x, n2y, nummod1, nummod2, 
     &                fackz, wunit, wbiga, wa, rho11, tau12, rhol22, 
     &                tau21, rhor22, tau23, trans, trans2, refle)
c
      call modes(ndimk1xs, ndimk1ys, ndimk2xs, ndimk2ys, n1x, n1y,
     &           n2x, n2y, vk1x, vk1y, vk1z, vk2x, vk2y, vk2z, vk3x,
     &           vk3y, vk3z, theta, gvaccuum, eps1, eps2, eps3,
     &           f1p, f1m, f2p, f2m, f3p, f3m)
      call current_z(ndimk1xs, ndimk1ys, ndim1s, n1x, n1y, vk1z, vk3z, 
     &               f1p, f1m, f3p, ip1, ipol0, m0, n0, refle, trans, 
     &               trans2, refcurr, refcurr0s, refcurr0p, tracurr,  
     &               tracurr0s, tracurr0p,tranum,tl,al)

      
      
!------------------------------------------------------------!
! Amplitudes of the Electric Field at the Openings formalism !
!------------------------------------------------------------!
      
c       denden = 1.d0-rhor22(1,1)*rhol22(1,1) *fackz(1)**2
        denden = rhor22(1,1)*rhol22(1,1) *fackz(1)**2

             
	  !Warging: gvaccumm is not go=2pi/vlan
	  !It works for G1 and GIII because vk2z is normalized by gvaccumm
	  vk2z(1,0)= vk2z(1,0)/vlan
	  if(dimag(vk2z(1,0))==0.0)then
	   expo = exp(ui*vk2z(1,0)*width0)
	  else
         expo = exp(-abs(vk2z(1,0))*width0)
	  end if

	  Y2= vk2z(1,0)*vlan/(2*pi)
	  gv  = (vk2z(1,0)*vlan/(2*pi))*2.d0/(expo-1.0d0/expo)
       sigma=(vk2z(1,0)*vlan/(2*pi))*(expo+1.0d0/expo)/(expo-1.0d0/expo)
  

	  phaserho = atan(imag(rhol22(1,1))/real(rhol22(1,1)))
	  phase=2.0*vk2z(1,0)*width0+2*phaserho

	  vk2z(1,0)= vk2z(1,0)*vlan
        
	 
	  greenI=(vk2z(1,0)/gvaccuum)*
     & (1.d0-rhol22(1,1))/(1.d0+rhol22(1,1))
	  
	  greenIII=(vk2z(1,0)/gvaccuum)*
     & (1.d0-rhor22(1,1))/(1.d0+rhor22(1,1))

	  g=2.0*pi/vlan
	  gc=2.0*pi/(2.0*sqrt(eps2)*ax0) 
	  lanc= 2.0*sqrt(eps2)*ax0
	  
c	  gva2 = -ui*(1.0/(width0*gc))+
c     &       ui*(-eps2*width0/3.0+1.0/(width0*gc**2))*(g-gc)    	  
        
c  	  sigmaa2 = -ui*(1.0/(width0*gc))+
c     &       ui*(2.0*eps2*width0/3.0+1.0/(width0*gc**2))*(g-gc)

   	  sigmaa = -ui*lanc/(2.d0*pi*width0)+
     &        -ui*(1.0d0/(2.d0*pi*width0)+
     &      4.d0*pi*width0*eps2/(3.d0*lanc**2))*(vlan-lanc)    	  
        
  	  gva = -ui*lanc/(2.d0*pi*width0)+
     &        -ui*(1.0d0/(2.d0*pi*width0)-
     &         2.d0*pi*width0*eps2/(3.d0*lanc**2))*(vlan-lanc)    	  
        
	  Gsum=(greenI+greenIII)/2.d0
	  Gdif=(greenI-greenIII)/2.d0 

c	  con_max  = g - gc - dimag(greenI)/(eps2*width0)
c	  con_maxA = g - gc - dimag(greenI+greenIII)/(2.d0*eps2*width0)
c	  con_max_valid = -real(greenI)**2/(imag(greenI)-abs(sigma)-abs(gv))
	  
c	  con_max4= Abs(Gsum-sigma)-abs(gv)

        con_max= Abs(greenI-sigma)*Abs(greenIII-sigma)-abs(gv)**2	  
	  con_max1 =vlan-lanc +
     &   (lanc**2/(2.d0*pi))*dimag(greenI+greenIII)/(2.d0*eps2*width0)
        bden=(abs(Gsum-sigma)**2-abs(gv)**2)/(2.d0*abs(gv))
	   
	  phase3 = atan(imag(greenI-sigma -gv)/real(greenI-sigma -gv))
	  phaseGsigma=atan((imag(greenI)-abs(sigma))/real(greenI))

	  So =(2.0*dsqrt(2.d0)/pi)*(sqrt(ay0*ax0/(rlx0*rly0))) ! 06/sep/17: must be changed to work for off-normal incidence
        illu = 2.0*sqrt(eps1)*So                             ! 06/sep/17: must be changed to work for off-normal incidence
         
	   
         dMINUS=Gsum-sigma-gv
	   aMINUS=1.d0/dMINUS
	   dPLUS =Gsum-sigma+gv         
	   aPLUS = 1.d0/dPLUS

	   e0fac(1) = illu/(2.0*dPLUS)  !>> Resonant (G-Sigma+Gv) 
	   e0fac(2) = illu/(2.0*dMINUS)
	   e0pfac(1)=-illu/(2.0*dPLUS)  !>>Resonant (G-Sigma+Gv) 
	   e0pfac(2)= illu/(2.0*dMINUS)

         denWITHminus=(Gsum-sigma)**2-gv**2-Gdif**2
	   denOUTminus=(Gsum-sigma)**2-gv**2
       
	!>>> Exact transmission and reflection  
	   den = denWITHminus
         e0  = illu * (greenIII - sigma) / den
         e0p = illu * gv / den
         !transFULL  = dreal(e0   * dconjg(e0p)   * gv)/sqrt(eps1)! 06/sep/17: must be changed to work for off-normal incidence
         transFULL  = real(greenIII)*(e0p*dconjg(e0p))/sqrt(eps1) ! Changed 13/nov/19
         refFULL  = 1.0d0 -2.d0*real(So*e0)+
     & real(greenI)*(e0*dconjg(e0))/sqrt(eps1)       


      !>>> Exact transmission (different versions)         
         e0 =e0fac(1)+e0fac(2)
	   e0p=e0pfac(1)+e0pfac(2)
	   transfac(1) = dreal(e0 * dconjg(e0p) * gv)/sqrt(eps1)
	   e0 =e0fac(1)
	   e0p=e0pfac(2)
	   transfac(2) = 2.0*dreal(e0 * dconjg(e0p) * gv)/sqrt(eps1)
         e0 =e0fac(2)
	   e0p=e0pfac(1)
	   transfac(3) = 2.0*dreal(e0 * dconjg(e0p) * gv)/sqrt(eps1)
  
	transfac(4)=(illu**2)*dreal(gv*(abs(aMINUS)**2-abs(aPLUS)**2 +
     &                2.d0*ui*dimag(aPLUS*dconjg(aMINUS))))/4.d0
      transfac(4)=transfac(4)/sqrt(eps1)  	   	   
	
	transfac(5)=-(illu**2/2.d0)*
     &  	        dimag(gv)*dimag(aPLUS*dconjg(aMINUS))
      transfac(5)=transfac(5)/sqrt(eps1)

	
	transfac(6)=illu**2*dreal(greenIII)*
     &               	(abs(aMinus)**2+abs(aPLUS)**2 -
     &                2.d0*dreal(aMinus*dconjg(aPlus)))/4.d0
	transfac(6)=transfac(6)/sqrt(eps1)
      
	transfac(7)=illu**2*dreal(greenIII)*abs(gv)**2/
     &   ((abs(Gsum-sigma)**2-abs(gv)**2)**2+
     &    4.d0*(dreal(Gsum)*abs(gv))**2) 
	transfac(7)=transfac(7)/sqrt(eps1)
      !<<< Exact transmission

      !>>> Approximated transmission from: S. Carretero-Palacios et al. Phys. Rev. B 85, 035417 (2012)
         den = denOUTminus
         e0  = illu * (greenIII - sigma) / den
         e0p = illu * gv / den
         transaFULL = dreal(e0  * dconjg(e0p)  * gv)/sqrt(eps1)
      
	!>>> Approximated transmission ---> from 8 to 11
	transfac(8)=(illu/2.d0)**2*dreal(greenIII)/
     &((dimag(Gsum)-dimag(sigma)+sign(abs(gv),dimag(sigma)))**2+
     &dreal(Gsum)**2) 
	transfac(8)=transfac(8)/sqrt(eps1)
      
	
      transfac(9)=(illu/2.d0)**2*dreal(greenIII)/
     &((dimag(Gsum)-dimag(sigmaa)+sign(abs(gva),dimag(sigmaa)))**2+
     &dreal(Gsum)**2) 
	transfac(9)=transfac(9)/sqrt(eps1)

	transfac(10)=(illu/2.d0)**2*dreal(greenIII)/
     & (dimag(Gsum)**2+dreal(Gsum)**2) 
	transfac(10)=transfac(10)/sqrt(eps1)

	transfac(11)=(illu/2.d0)**2*dreal(greenIII)/
     & ((vlan - sqrt(eps1/eps2)*lanc)**2+dreal(Gsum)**2) 
	transfac(11)=transfac(11)/sqrt(eps1)
!<<< Approximated transmission from: S. Carretero-Palacios et al. Phys. Rev. B 85, 035417 (2012)

       !>>> output
	if(flag_open)then
        write(15,*) vlan,transFULL,refFULL !,tracurr0p,transFULL,transaFULL
        !write(15,*) vlan,tracurr !,tracurr0p,transFULL,transaFULL
        !write(15,*) tracurr,tracurr0p,refcurr
!	  write(21,*) dreal(con_max),dreal(con_max1),dreal(bden)
!       write(17,*) dimag(sigma),dimag(sigmaa)
!	  write(18,*) dimag(gv),dimag(gva) 
!	  write(24,*) dreal(greenI),dimag(greenI)       
!	  write(25,*) dreal(greenIII),dimag(greenIII)      
!c	  write(22,*) dreal(den),dreal(dena),dreal(dena2)
!c	  write(20,*) dreal(fa),dreal(fb),dreal(fc)	  
!	  write(26,4) abs(rhol22(1,1)),dreal(phaserho),
!     &              dreal(rhol22(1,1)),dimag(rhol22(1,1)) 
!	  write(27,*) dreal(1.d0/(expo**2)),dimag(1.d0/(expo**2))
!	  write(28,5) abs(e0p),dreal(e0pfac(1)),dimag(e0pfac(1)),
!     &	                    dreal(e0pfac(2)),dimag(e0pfac(2))
!	   write(29,*) abs(e0)
!         write(32,*) transfac(1),transfac(2),transfac(3)
!	   write(30,7) transfac(5),transfac(6),transfac(7),transfac(8),
!     &        	   transfac(9),transfac(10),transfac(11)
	   
!	  write(31,*) abs(vk2z(1,0)*width0/vlan)	
	end if
1001    format(2x, 5(e12.6,1x))
5      format(5F30.12)
4      format(4F30.12)
6      format(6F30.12)
7      format(7F30.12)
c
c      end of loop in wavelengths
c
       end do
c
1000   format(2x,4(1x,e12.6),2x,i5)      
       close(19)
       
	end do !<< End of loop calculations
	close(inME)
	 stop
       end



	