*IrisME* 
============
This file is part of *IrisME*, a software for the calculation of scattering coefficients in arrays of rectangular holes perforating a metal film. IrisME is based on the Modal Expansion (ME) of the electromagnetic fields (see below). 

*IrisME* is licensed under the AGPL and it is free to use.  However, if you are using, or plan to use it, specially if it is for research or academic purposes, please send an email with your name, institution and a brief description of your interest for this program. Or course, if you use *IrisME* in a work that leads to a scientific or academic publication, we would appreciate it if you would kindly cite *IrisME* in your manuscript as:

> L. Martín-Moreno, F.J. García-Vidal, *Minimal model for
> optical transmission through holey metal films*, J. Phys.:
> Condens. Matter 20, 304214 (2008) *and*
> S.G. Rodrigo, F. de León-Pérez, and L. Martín-Moreno, 
>*Extraordinary Optical Transmission: fundamentals and applications*, 
> Proceedings of the IEEE 104, 2288 (2016) 

Commercial applications may also acquire a commercial license. Please contact <sergut@unizar.es> & <lmm@unizar.es>, for details.

Copyright (C) 2019 Sergio G Rodrigo <sergut@unizar.es> & Luis Martin-Moreno <lmm@unizar.es>

IrisME is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License,or (at your option) any later version.
  
IrisME is distributed in the hope that it will be useful for research or/and academic purpouses, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details. You should have received a copy of the GNU Affero General Public License along with IrisME. If not,see <http://www.gnu.org/licenses/>.

***

*Modal expansion (briefly)* 
--------------------------------

*IrisME* is a good tool to obtain fast and accurate results calculating scattering coefficients in Perfect Electric Conductor rectangular hole arrays, as compared to other numerical models widely used in computational electromagnetism.

Within ME, Maxwells equations are solved by expanding the EM fields in the different regions of space, transmission and reflection coefficients are thus calculated by imposing appropriate boundary conditions (see Ref. [1], and references therein). The whole space is divided in three regions: (I) the substrate, (II) holes and, (III) the superstrate. In region (I) the EM field is expanded in terms of TE and TM waveguide eigenmodes. However, good convergence is attained for many problems only by considering the less decaying transversal electric mode, $TE_{0,1}$ for subwavelength rectangular holes (the fundamental waveguide mode). Within this minimal model, results can be worked out analytically.  

The ME method was firstly applied and developed by our group in the context of Extraordinary optical transmission (EOT) research [1]. Roughly speaking, EOT are electromagnetic resonances through subwavelength apertures in either a flat or a corrugated metal film. These resonances provide a higher transmission of radiation than would be expected from such a small aperture size. EOT was discovered in 1998 [3] and since then it has been a very active research field, leading both to the discovery of new ways of enhancing
optical transmission and to its application to sensing, color filters, metamaterials, lenses, optical trapping, enhancement of nonlinear effects, among others. 

 *References:*
 
[1] F.J. Garcia-Vidal, L. Martin-Moreno, T.W. Ebbesen, L. Kuipers, *Light passing through subwavelength apertures*, Rev. Mod.
Phys. 82, 729–787 (2010).

[2] L. Martín-Moreno, F.J. García-Vidal, *Minimal model for optical transmission through holey metal films*, J. Phys.:
Condens. Matter 20, 304214 (2008).

[3] T. W. Ebbesen, H. L. Lezec, H. F. Ghaemi, T. Thio, and P. A. Wolff, *Extraordinary optical transmission through subwavelength
hole arrays*, Nature 391, 667–669 (1998).

***

IrisME (the programs) 
--------------------------------

**Getting and installing *IrisME***

The library can be download from: <https://github.com/IrisFDTD/IrisME>

Two different implementations of ME are made available:

1. **IrisME** is a Wolfram Mathematica (C) script for the calculation of scattering coefficients from arrays of rectangular holes in metals. The metal is approximated as a Perfect Electric Conductor, which is a good approximation in the THz regime for geometrical features much larger than the skin depth. A brief description of the method can be found in the Appendix of Ref. [3]. The formalism is based on the amplitudes of the electric field at the openings (see Ref. [2], for further details). This version has been checked to work from Mahtematica 8.0 version.

2. **IrisMEadv** is a program written in Fortran with extended capabilities: additional EM modes inside the holes, realistic optical properties of the metals under the Suface Impedance Boundary Conditions (SIBC), multiple scattering formalism (see Ref. [2], for further details),... If you are interested in these add-ons contact to Sergio G Rodrigo <sergut@unizar.es> & Luis Martin-Moreno <lmm@unizar.es>. IrisMEadv32_Windows and IrisMEadv64_Linux are versions for Windows and Linux, respectively.

Several examples of the IrisME program in action are provided. See:

+ **Example 1**: *Localized Extraordinary Optical Transmission*

