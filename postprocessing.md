# Post-processing and analysis
Here we assume that you checked the convergence of your system and that all the run ended successfully (no protein unfolding, no other issues). Here we will not focus on convergence, error calculation and metadynamics technical details, considering that such topics are already covered by the [Metadynamics tutorial](https://www.plumed-tutorials.org/lessons/21/004/data/NAVIGATION.html). 

### Reweighting 
After the end of the simulation, we need to reweight our free energy landscape on apt collective variables. As anticipated in the previous step, we will use the distance from the origin of the reference frame $\rho$ and the number of contacts $c$. We thus prepare step-by-step a reweighting script (you can find it in the GitHub folder, called `reweight.dat`).

We begin from the reading of the relevant data files (refer to the `plumed.dat` file if you do not remember what is contained in them):
```plumed
rho: READ FILE=results/coord_rho.dat VALUES=rho IGNORE_FORCES IGNORE_TIME
c: READ FILE=results/coord_rho.dat VALUES=c IGNORE_FORCES IGNORE_TIME
metad: READ FILE=results/metad_data.dat VALUES=metad.* IGNORE_FORCES IGNORE_TIME
restr_rmsd: READ FILE=results/rmsd_restraint.dat VALUES=restr_rmsd.* IGNORE_FORCES IGNORE_TIME
```
After loading the data, we have to perform the reweighting of the metadynamics potential via the [Tiwary-Parrinello estimator](https://doi.org/10.1021/jp504920s). Concurrently, we remove the (almost negligible) contribution of the RMSD restraining, obtaining the final weights for the histogram
```plumed
weights: REWEIGHT_BIAS TEMP=300 ARG=metad.rbias,restr_rmsd.bias
```
Having the weights, we can compute the histogram, considering the data from 200 ns to 1 µs to ignore the out-of-equilibrium portion of the simulation:
```plumed
HISTOGRAM ...
  ARG=rho,c
  GRID_MIN=0.,0
  GRID_MAX=3.0,160
  GRID_BIN=300,160
  KERNEL=DISCRETE
  LOGWEIGHTS=weights
  LABEL=histo
  UPDATE_FROM=200000
  UPDATE_UNTIL=1000000
... HISTOGRAM
```
This histogram can be converted to a free energy landscape
```plumed
ff: CONVERT_TO_FES GRID=histo TEMP=300
```
and finally printed
```plumed
DUMPGRID GRID=ff FILE=reweighted_fes.dat
```
From this data we can plot the 2D free energy landscape reweighted on $\rho$ and $c$:

<p align="center">
  <img src="https://github.com/riccardocapelli/VMetaD-tutorial/blob/main/img/fes.jpg?raw=true" alt="Alt text" width="75%">
  <br>
  <em>Free energy surface projected along the new CVs.</em>
</p>

### ∆G calculation
In the GitHub folder you can find a simple python script to get the free energy difference between two basins, `deltaG.py`. You have to input the path of the reweighted FES, the minima and the maxima in x and y for both the basins and it returns the free energy difference. For my run, we defined the two minima and obtained this result:
```
$> python3 deltaG.py reweighted_fes.dat 0.45 0.55 105 115 2. 2.8 0 0.5
The free energy difference between basin A and basin B is -3.53 kcal/mol
```
We now have $\Delta G_{\text{MetaD}}$, the last step is the calculation of the entropic correction.

### Entropic correction calculation
We have to compute the volume of the protein included in the sphere. To do so, we can again use VMD. After loading the structure, we open the Tk console and select all the atom that are within 2.8 nm from the center of mass we computed in the previous step of the tutorial:
```
set sel [atomselect top "sqrt((x-35.28876876831055)^2 + (y-34.06196594238281)^2 + (z-32.041622161865234)^2) < 28.0"]
```
having set up this selection, we can write an output PDB file.
```
$sel writepdb atoms_in_sphere.pdb
```
Having this file, we can use `gmx sasa` to have an estimate of the volume of these atoms:
```
gmx sasa -f atoms_in_sphere.pdb -s atoms_in_sphere.pdb -tv volume.xvg
```
obtaining a volume of 27.294 nm$^3$.

Putting this value and $\rho_{\text{s}}=2.8$ nm in the formula we get the entropic correction

$$
\begin{align*}
RT \log\left( \frac{V^{0}}{\frac{4}{3}\pi \rho_{\text{s}}^3 -V_{\text{host}}}\right)
&= -2.18 \text{ kcal/mol}
\end{align*}
$$

And we have a final binding free energy value $\Delta G^{0}=-(3.53+2.18)$  kcal/mol $= -5.71$ kcal/mol, which is close to the experimental value of $-5.2$ kcal/mol obtained by Morton _et al._ in [this work](https://doi.org/10.1021/bi00027a006).

##### [Back to VMetad home](NAVIGATION.md)
