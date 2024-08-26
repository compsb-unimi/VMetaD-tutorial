# Input files preparation

In the zip file available in the GitHub repository you will find all the initial files needed to run this tutorial.

## Preliminary molecular dynamics run
As a first step for the VMetaD run, we have to choose the atoms which will constitute the reference frame on which we will calculate the position of the ligand with respect to the protein. To be sure in choosing residues that remain stable (i.e. do not belong to a moving loop), we run a 100 ns-long plain molecular dynamics (MD) simulation (it took 2 h on 4 CPUs of a M1 Max MacBook Pro).

After the simulation, we calculate the per-residue root mean square fluctuations (RMSF), which highlights the more dynamic residues.

<p align="center">
  <img src="img/rmsf.jpg" alt="Alt text" width="50%">
  <br>
  <em>Per-residue root mean square fluctuation (RMSF) computed from a plain MD simulation of lysozyme-benzene complex. The portion of the plot with light blue background represents the C-terminal domain, where the benzene binds.</em>
</p>

We can see that almost all the C-terminal domain (residues 71-162) does not show large fluctuations, and only few residues are above the (arbitrary) RMSF threshold of 0.15 nm. We thus consider to define the reference frame considering all the residues 71-162 but the five ones with RMSF > 0.15 nm (residues 135, 136, 137, 140, and 141).

## Choice of the restraining potential size

We now need to choose the size of the restraint potential. In the [original paper](https://doi.org/10.1021/acs.jpclett.9b01183) we showed that the reliability of the estimates is not affected by the size of the potential. However, we have to keep in mind that we need a part of the box where the ligand can stay far away from the protein in order to represent the unbound state in a satisfactory way. An important point here is that the sphere constraint __must__ be inside the box, otherwise the entropic correction will not accurately account for the loss of configurational space. To visually inspect how large the potential is, and to get a feel for the possible movements of the ligand, we can visualize both the system and the restraint with VMD (downloadable [here](https://www.ks.uiuc.edu/Research/vmd/)).

We can open VMD and load the `starting.gro` structure file in the GitHub folder. After the structure is loaded and the visualization has been set up at your taste, you can open the Tk console and write
```
pbc box
```
which draws the cubic box in which the system is inserted. 

Now we can generate the atom selection we defined after checking the RMSF:
```
set sel [atomselect top "resid 71 to 134 or resid 138 to 139 or resid 142 to 162"]
```
The console should answer
```
atomselectXX
```
Where `XX` is a number. The selection for the reference frame has been defined and named `$sel`. Now we can compute the position of the center of mass of this set of atoms:
```
measure center $sel weight mass
```
The console should answer with the position of the center of mass (in Ångstrom)
```
35.28876876831055 34.06196594238281 32.041622161865234
```
Knowing this information, we can draw the sphere with a radius of (for example) 2 nm  (20 Å) with the following command
```
draw sphere {35.28876876831055 34.06196594238281 32.041622161865234} radius 20 resolution 100
```
Receiving a number as an answer from the console. Such number is the ID of the 3D object we just draw. The drawn sphere is opaque, not allowing us to see inside it; to make it transparent, we need to specify that we want a transparent material
```
draw material Transparent
```
We can see that the sphere contains the entire domain, but it is probably too small to represent the unbound state in a precise way. Let's delete the sphere using the ID of the 3D object (let's say that it is `14`), and plot a new sphere of radius 2.8 nm
```
draw delete 14
draw sphere {35.28876876831055 34.06196594238281 32.041622161865234} radius 28 resolution 100
```
You can see the expected result below
<p align="center">
  <img src="img/sphere_box.jpg" alt="Alt text" width="50%">
  <br>
  <em>Cartoon representation of the lysozyme-benzene complex, including the restraining potential applied within a 3 nm radius of the reference frame center of mass. The boundaries of the simulation box are also highlighted to show that the sphere is entirely contained by the box. </em>
</p>
