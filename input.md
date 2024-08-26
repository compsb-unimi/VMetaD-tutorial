# Input files preparation

In the zip file available in the GitHub repository you will find all the initial files needed to run this tutorial.

## Preliminary molecular dynamics run
As a first step for the VMetaD run, we have to choose the atoms which will constitute the reference frame on which we will calculate the position of the ligand with respect to the protein. To be sure in choosing residues that remain stable (i.e. do not belong to a moving loop), we run a 100 ns-long plain molecular dynamics (MD) simulation (it took 2 h on 4 CPUs of a M1 Max MacBook Pro).

After the simulation, we calculate the per-residue root mean square fluctuations (RMSF), which highlights the more dynamic residues.

<p align="center">
  <img src="img/rmsf.jpg" alt="Alt text" width="50%">
  <br>
  <em>Per-residue root mean square fluctuation (RMSF) computed from a plain MD simulation of lysozyme-benzene complex.</em>
</p>
