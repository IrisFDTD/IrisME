Copyright (C) 2019 Sergio G Rodrigo <sergut@unizar.es> & Luis Martin-Moreno <lmm@unizar.es>

**Important:** 
*IrisME* is licensed under the AGPL and it is free to use. if you are using, or plan to use this example, specially if it is for research or academic purposes, please send an email with your name, institution and a brief description of your interest for this program.  If you use this example of *IrisME* in a work that leads to a scientific or academic publication, we would appreciate it if you would kindly cite Refs. [1,2]  in your manuscript (see below).


***
- - -

> **Example 2: Localized Extraordinary Optical Transmission**

This example is extracted from Ref. [1]. The transmission calculated can be compared with one of those shown in **Figure 7-(top)**. In that work we investigated the so-called Localized Extraordinary Optical Tranmission (LEOT) phenomenom, for which the spectral position of transmission and reflection resonances appearing in subwavelength apertures drilled in metallic films drastically changes when the dielectric environment is modified, for extremely subwavelength fim widths.  These results are quantitatively correct in the microwave and terahertz regimes, but they also have qualitative validity in the optical regime. Additional details about LEOT can be found in Ref. [2], and references therein.

**Description:**
This example shows you how to obtain transmission and reflection of terahertz waves through a given hole array (perforating a Perfect Electric Conductor metal screen). The period chosen is $L_x=L_y=400 \mu$m, and the rectangles are defined by $a_x = 10 \mu$m and  $a_y = 350 \mu$m sides.  The PEC screen is $h = 25 \mu$m and the holes are filled with a dielectric, $\epsilon_{II} = 4$. The cover and substrate dielectric constants are the same ($\epsilon_{I}= \epsilon_{III}= 2$).  Light iluminates the structure at normal incidence, being the electric field polarized along the x direction.

<img style="float: left;" width="400" src="../HA_geometry.png">

**Steps:**
1. Run IrisME.nb, included in this version, at the same directory this script in launched. The file input_IrisME.dat contains the geometrical and rest of parameters for IrisME.nb. The output file T-R_IrisME.dat has 3 data columns: wavelengths ( $\mu$m), transmission and reflection. Mathematica has to be installed (from 8.0 version). If not, you will still use IrisMEadv. 

2. Run the iPython script. The script calls the IrisMEAdv 64-bits executable for Windows/Linux, included in this release.  The file Input_IrisMEadv.dat contais the geometrical and rest of parameters for IrisMEadv. The output file is T-R_IrisMEadv.dat, which contains 3 data columns: wavelengths($\mu$ m), transmission and reflection. The iPython script finally plots all the outputs in a figure. Note that you can alternatively run IrsMEAdv programs in from the console.
    
*References:*

[1] S. Carretero-Palacios, F. J. García-Vidal, L. Martín-Moreno, and S. G. Rodrigo [*Effect of film thickness and dielectric environment on optical transmission through subwavelength holes*](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.035417),
Phys. Rev. B 85, 035417 (2012).

[2] S.G. Rodrigo, F. de León-Pérez, and L. Martín-Moreno, [*Extraordinary Optical Transmission: fundamentals and applications*](https://ieeexplore.ieee.org/document/7592449), Proceedings of the IEEE 104, 2288 (2016). 

    
***
- - -
