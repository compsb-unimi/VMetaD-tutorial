# VMetaD-tutorial
A tutorial to setup and run Volume-based Metadynamics simulations. The original work which explains the theory for Volume-based Metadynamics (VMetaD) is [this](https://doi.org/10.1021/acs.jpclett.9b01183). Some time later, in [this](https://doi.org/10.1021/acs.jctc.1c00649) paper we added some details (the RMSD restraining) which we will use in this tutorial.

This guide is divided in different sections

#### [Volume-based Metadynamics theory](theory.md)
A brief introduction to the rationale behind the use of VMetaD.

#### [Case study: L99A T4 Lysozyme:Benzene](lysozyme_benzene.md)
We will present the case study for this tutorial, a small engineered protein which binds benzene. We will use, to favor the possibility to test the method without that need of HPC resources, a simplified model called multi-eGO (more details [here](https://doi.org/10.26434/chemrxiv-2024-jcmgc)) which can run on laptop/workstation at the timescales of interest. 

#### [Input files preparation](input.md)
We will discuss all the preliminary steps to set up VMetaD, from the restraining potential size to the selection of the atoms that define the reference frame. We will also discuss all the PLUMED instructions to implement the method.

#### [Post-processing, reweighting and entropic correction](postprocessing.md)
Finally, we will discuss how to analyze the VMetaD files to get the binding free energy difference, reweighting the free energy surface on apt CVs. We will also discuss how to compute the entropic correction due to the presence of the restraining potential. 

___NB:___ This tutorial assumes that you know metadynamics theory and practice in PLUMED (if not, a good start can be the [PLUMED tutorial](https://www.plumed-tutorials.org/lessons/21/004/data/NAVIGATION.html) about it).


```mermaid
flowchart TB
  A[PLUMED syntax] -.-> C[Theory];
  B[Metadynamics Tutorial] -.->  C
  C --> D[Case study]
  D --> E[Input files]
  E --> F[Post-processing]
  
  click A "syntax" "This lesson teaches you the basic features of the PLUMED input syntax"
  click B "metad" "This lesson teaches you how to perform and analyze Metadynamics simulations"
  click C "theory.md" "Volume-based Metadynamics theory intro & overview"
  click D "lysozyme_benzene.md" "Case study: L99A T4 Lysozyme:Benzene"
  click E "input.md" "Input files preparation"
  click F "postprocessing.md" "Post-processing, reweighting and entropic correction"
```
