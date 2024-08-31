# Volume-based Metadynamics intro and theory

In Metadynamics, one of the most important decisions the user must make is selecting an appropriate collective variable (CV) that can comprehensively describe the slow degrees of freedom relevant to the phenomenon of interest. When focusing on protein-ligand (or more generally, host-guest) binding, significant effort has been dedicated to identifying CVs, often based on human intuition and/or data-driven techniques.

In VMetaD, the choice is simplified: the CVs are the 3D coordinates of the ligand's center of mass relative to a reference frame defined by the host. The primary goal of VMetaD is to identify all possible binding and unbinding pathways quickly and reliably, and this approach is tailored to that objective.

The main problem in this approach is given by the time needed for VMetaD to "fill" the unbound state, allowing recrossing events: if the host is immersed in a reasonably large box, the (long) time needed to observe a favorable binding event makes this technique completely useless. For this reason, a restraining potential is applied (usually at the center of mass of the host), which limits the conformational space and allows the guest to enter and exit the host in an efficient manner. The shape of such restraining potential strongly depends in what is the problem or how much information the user knows about the system. In the [original implementation](https://doi.org/10.1021/acs.jpclett.9b01183) (and in this tutorial) we implemented a spherical-shaped restraint and used spherical coordinates as CVs. In a [subsequent work](https://doi.org/10.1021/acs.jctc.1c00649) we used a parabolic potential and parabolic coordinates as CVs.
In the case of the sphere, we have the restraining potential defined as

$$
U_{\text{s}}(\rho(t)) = 
\begin{cases}
  k \left(\rho(t) - \rho_{\text{s}}\right)^{2} ~~~ \text{if }\rho(t) > \rho_{\text{s}} \\
  0 ~~~~~~~~~~~~~~~~~~~~~\text{otherwise}
  \end{cases}
$$

where $k$ is the harmonic constant, $\rho(t)$ is the distance of the guest center of mass from the center of the restraining sphere, and $\rho_{\text{s}}$ is the radius of the restraining sphere.

Given enough simulation time, this setup alone can return all the possible binding and unbinding pathways (if that is the only goal). However, considering the fact that VMetaD builds a history-dependent potential, we can also compute the free energy difference between the bound and unbound states (i.e. the binding free energy) when we are at convergence. The main point here is to identify the unbound state in a comprehensive way (we assume that the bound state is well defined in volume). To do this, we reweight the converged free energy landscape and project it onto some new CVs that allow better identification and discrimination of both the bound and unbound states.

Finally, we have to take into account the non-negligible entropic contribution given by the presence of the restraining potential in the simulation (more generally, this should be done in any simulation where the ligand concentration is not the standard one) to make the _in silico_ ∆G compatible with its experimental evaluation. 
In general, the free energy difference in standard condition is given by

$$
\Delta G^{0} = -RT \log\left( C^{0} K_{b} \right)
$$

where $R$ is the gas constant, $T$ is the system temperature, $K_{b}$ is the binding constant, $C^{0} = (1660 \text{Å}^{3})^{-1}$ is the standard concentration. In the context of our _in silico_ approach, this is equal to

$$
\Delta G^{0} = \Delta G_{\text{MetaD}} + RT \log\left( \frac{V^{0}}{V_{\text{restr}} -V_{\text{host}}}\right)
$$

where $\Delta G_{\text{MetaD}}$ is the free energy difference obtained from the VMetaD calculation, $V^{0}$ is the reciprocal of $C^{0}$, $V_{\text{restr}}$ is the volume included in the restraining potential, and $V_{\text{host}}$ is the portion of the host included in the restraining potential. This, in case of a spherical potential, results in

$$
\Delta G^{0} = \Delta G_{\text{MetaD}} + RT \log\left( \frac{V^{0}}{\frac{4}{3}\pi \rho_{\text{s}}^3 -V_{\text{host}}}\right)
$$

where $\rho_{\text{s}}$ is the radius of the restraining sphere.

##### [Back to VMetad home](NAVIGATION.md)
