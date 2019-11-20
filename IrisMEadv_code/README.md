*FAQ*
============

> While Windows and Linux executable versions of the program are provided (IrisMEadv64_Windows & IrisMEadv64_Linux), you can build it by your own. The code is written in a kind of messy "mix" of old and new versions of Fortran language, and this fact must be taken into account when compiling and linking the files.

**How can we compiled and linked the code?**
- A general rule can not be provided, but here we go with two recipies that work in Windows and Linux
1. Windows with Microsoft Visual Studio - set Fortran properties to: /nologo /O2 /module:"x64\Release\\" /object:"x64\Release\\" /Fd"x64\Release\vc100.pdb" /libs:dll /threads /c
2. Linux with ifort - call  ifort -O3 -o IrisMEAdv64_Linux.exe inttoAsscii.for index2.f redef_field.f interfase23.f index1.f overlap.f drude_general.f modes.f current_z.f matrix.f interfase21.f 2d.f wavevec.f leqt1c.f interfase12.f scatt_coef.f
