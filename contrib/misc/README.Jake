Analysis code for confined water project:

Compile by running: bash complile.sh
Change compiler flag in compile.sh as needed
Code (a.out) was originally compiled on calhoun with ifort


Rename movie file to fort.10
Need fort.10 and confined.inp in the subdirectory you are running code from
Example input/output files are given in the subdirectory TEST/


Variable definitions in confined.inp (input file):

lrhoz: true if you want to calculate density profile

rhoz_width: bin width for density profile

ldipolez: true if you want to calculate the z-profile of the cosine of
the angle between the water dipole vector and the z-axis and its order
parameter  

dipolez_width: bin width for dipole profile

lhbondz: true if you want to calculate number of hydrogen bonds as
function of z 

hbondz_width: bin width for H-bond profile

lQ4z: true if you want to calculate tetrahedral order parameters as
function of z 

Q4z_width: bin width for order profiles



Output files:

rhoz.dat: density profiles (in g/mL) for water, based on COM

cos_dipole: z-profile of the cosine of the dipole vector relative to
the z-axis of the box 

p2_dipole: order parameter for above cosine
(1 = perpendicular to z-axis, 0 = random, -0.5 = parallel to z-axis)
(see one of our chromatography papers)

hbondz.dat: number of hydrogen bonds per water molecule as a function of z

donorz.dat: number of h-bonds donated per water molecule as a function of z

acceptz.dat: number of h-bonds accepted per water molecule as a function of z

tetra_angle.dat: angle based tetrahedral order parameter as function of z 
(0 for random, 1 for perfect order)
(see Nature, vol. 49, p. 318) 

tetra_angle.dat: distance based tetrahedral order parameter as function 
of z 
(from Molecular Physics, vol. 93, p. 511, but adjusted such that the
order parameter approaches 1 for perfect order, and approaches 0 for 
the ideal gas [rather than the inverse])   

All profiles are shifted by half the box length, i.e., z=0 corresponds
to the center of the box. 

NOTE: Written specifically for this project, not for general use.  Analyzes
only box 1, molecule type 1, assumes molecule type 1 is water.


Feel free to modify any of this as you see fit.
