*IrisME* 
============
*IrisME* is a fast and accurate software for calculating scattering coefficients from arrays of rectangular holes perforating a metallic film. *IrisME* is based on the Modal Expansion (ME) of the electromagnetic fields (see below). <img style="float: right;" width="400" src="./HA_geometry.png">

*IrisME* is licensed under the AGPL and it is free to use.  However, if you are using, or plan to use it, specially if it is for research or academic purposes, please send an email with your name, institution and a brief description of your interest for this program. If you use *IrisME* in a work that leads to a scientific or academic publication, we would appreciate it if you would kindly cite *IrisME* in your manuscript as:

> L. Martín-Moreno, F.J. García-Vidal, *Minimal model for
> optical transmission through holey metal films*, J. Phys.:
> Condens. Matter 20, 304214 (2008) *and*
> S.G. Rodrigo, F. de León-Pérez, and L. Martín-Moreno, 
>*Extraordinary Optical Transmission: fundamentals and applications*, 
> Proceedings of the IEEE 104, 2288 (2016) 

Commercial applications may also acquire a commercial license. Please contact <sergut@unizar.es> & <lmm@unizar.es>, for details.

Copyright (C) 2019 Sergio G Rodrigo <sergut@unizar.es> & Luis Martin-Moreno <lmm@unizar.es>

*IrisME* is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License,or (at your option) any later version.
  
*IrisME* is distributed in the hope that it will be useful for research or/and academic purpouses, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details. You should have received a copy of the GNU Affero General Public License along with *IrisME*. If not,see <http://www.gnu.org/licenses/>.

***

*Modal expansion (briefly)* 
--------------------------------

Within the ME method, Maxwells equations are solved by expanding the electromagnetic (EM) fields in the different regions of space, transmission and reflection coefficients are thus calculated by imposing appropriate boundary conditions (see Ref. [1], and references therein). The whole space is divided in three regions: (I) the cover, (II) holes and, (III) the substrate. In regions (I) and (III) the EM fields are expanded in planes waves. In region (II) in terms of TE and TM waveguide eigenmodes. However, good convergence is attained in many situations only by considering the less decaying transversal electric mode. Within this minimal model results can be worked out analytically.  

The ME method (on which *IrisME* is based) was developed in context of Extraordinary optical transmission (EOT) research [1]. Roughly speaking, EOT is a family of EM resonances through subwavelength apertures, in either a flat or a corrugated metallic film. These resonances provide high transmission of light, much more that would be expected for such tiny apertures as compared to the wavelength of light [2]. EOT was discovered in 1998 [3] and since then it has been a very active research field, leading both to the discovery of new ways of enhancing the optical transmission and to its application to sensing, color filters, metamaterials, lenses, optical trapping, enhancement of nonlinear effects, among others [4]. 

 *References:*
 
[1] F.J. Garcia-Vidal, L. Martin-Moreno, T.W. Ebbesen, L. Kuipers, [*Light passing through subwavelength apertures*](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.82.729), Rev. Mod. Phys. 82, 729–787 (2010).

[2] H. A. Bethe, [*Theory of difraction by small holes*](https://journals.aps.org/pr/abstract/10.1103/PhysRev.66.163), Phys. Rev. 66, 163–182 (1944).

[3] T. W. Ebbesen, H. L. Lezec, H. F. Ghaemi, T. Thio, and P. A. Wolff, [*Extraordinary optical transmission through subwavelength
hole arrays*](https://www.nature.com/articles/35570), Nature 391, 667–669 (1998).

[4] S.G. Rodrigo, F. de León-Pérez, and L. Martín-Moreno, [*Extraordinary Optical Transmission: fundamentals and applications*](https://ieeexplore.ieee.org/document/7592449), Proceedings of the IEEE 104, 2288 (2016). 

***

*IrisME (the programs)*
--------------------------------

**Getting and installing *IrisME***

The library can be download from: <https://github.com/IrisFDTD/IrisME>

Two different implementations of the ME method are made available:

1. **IrisME** is a Wolfram Mathematica (C) script for the calculation of scattering coefficients from arrays of rectangular holes in metals. The metal is approximated as a Perfect Electric Conductor, which is a good approximation in the THz regime for geometrical features much larger than the skin depth. A brief description of the method can be found in the Appendix of Ref. [1]. The formalism is based on the amplitudes of the electric field at the openings (see Ref. [2], for further details). This version has been checked to work from Mahtematica 8.0 version.

2. **IrisMEadv** is a program written in Fortran with extended capabilities: additional EM modes inside the holes, realistic optical properties of the metals under the Suface Impedance Boundary Conditions (SIBC), multiple scattering formalism (see Ref. [2], for further details),... If you are interested in these add-ons contact to Sergio G Rodrigo <sergut@unizar.es> & Luis Martin-Moreno <lmm@unizar.es>. While Windows and Linux executable versions of the program are provided (IrisMEadv64_Windows & IrisMEadv64_Linux), you can build it by your own. The code is shared at ./IrisMEadv_code. It is written in a "mix" of old and new versions of Fortran language. See the README file at same directory, for additional details. 

*References:*

[1] S.G. Rodrigo, [*Terahertz gas sensor based on absorption-induced transparency*](https://epjam.edp-open.org/articles/epjam/full_html/2016/01/epjam160014/epjam160014.html), EPJ Applied Metamaterials 3, 11 (2016).

[2] L. Martín-Moreno, F.J. García-Vidal, [*Minimal model for optical transmission through holey metal films*](https://iopscience.iop.org/article/10.1088/0953-8984/20/30/304214/meta), J. Phys.: Condens. Matter 20, 304214 (2008)

***

*Examples*
------------------------------------------------------------

|Available|Example|Published in (see [Publons](https://publons.com/researcher/2827893/sergio-gutierrez-rodrigo/))|
|:-----:|:-----:|:-----:|
|✅  |(1) Extraordinary Optical Transmission| S. G. Rodrigo, [*Optical Properties of Nanostructured Metallic Systems: Studied with the Finite-Difference Time-Domain Method*](https://www.springer.com/gp/book/9783642230844), Springer-Verlag, Berlin, (2012) *(Fig 1.12)*|
|✅  |(2)  Localized Extraordinary Optical Transmission| S. Carretero-Palacios, F. J. García-Vidal, L. Martín-Moreno, and S. G. Rodrigo [*Effect of film thickness and dielectric environment on optical transmission through subwavelength holes*](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.035417), Phys. Rev. B 85, 035417 (2012) *(Fig 7)* |
|✅  |(3)  Absoption Induced Transparency| S.G. Rodrigo, and L. Martín-Moreno, [*Absorption-induced transparency metamaterials in the terahertz regime*](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-41-2-293), Opt.  Lett. 41, 293-296 (2016) *(Fig 2)* |




```python

```
