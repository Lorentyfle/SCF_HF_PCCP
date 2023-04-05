==========================================================================
=                       How to use the program?                          =
==========================================================================

=========================Inside the input file=============================

M
2           <= You should put the exact size of your matrix under M.
--Begin--
1.45        <= Each Slater coefficients are below one to each other.
2.90        <= No comment and no text should be put in between "--Begin--" and "--End--".
--End--

--SCF_Begin--
p =         <= Any modification isn't recommanded.
1           <= Parameter for the starting coefficient. It says at which point the value of the coefficient vary.
cp1 =
1           <= Value of the coefficient if we are at a rank <= p
0           <= Value of the coefficient if we are at a rank > p
thr_SCF =   <= It's to be noted that if the value exceed 1e-9, the threshold read will not be exact. Please, refer to the threshold outputed to be certain of the precision of the SCF.
1.0E-8      <= Value of the threshold in scientific notation.
--SCF_End--

===========================================================================

When the input is finalised, the "compile.sh" could be runned with the
terminal inside the folder with all the fortran files to compile and 
compute the values.
The "compile.sh" was tested with the compiler gfortran on a linux system and
with the mingw64 application.

If any problem during the launch is to be remarked, verify that you have
gfortran and use the following command:

gfortran main.f90 module_cst.f90 module_sub.f90 module_diag.f90 -o HF_computation_compiled.x

--No copyright--
Thank you for using our HF program.

Timothee Jamin and Alex Delhumeau