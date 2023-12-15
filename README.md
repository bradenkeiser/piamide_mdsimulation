# piamide_mdsimulation

These scripts analyze the lifetime presence of piamide interaction between an aromatic ring and an amide region. The calculation is accessible in both Python and R environments, utilizing Bio3D and MDAnalysis for trajectory reading, respectively. Both GROMACS and Amber trajectory inputs are suitable; Bio3D offers support for reading CHARMM DCD trajectories only while MDAnalysis can work with both GROMACS XTC and Amber DCD inputs.

The examples folder provides analysis over the _EV_Xyn11TS xylanase at a favorable region between a phenyalanine ring and FF backbone region. The pre-rendered cutoff is for distances of 4.5 angstroms between the ring's center of mass (COM) and the amide COM with an angle less than 40 degrees in the standard configuration. The angular is deduced from a point on the normal plane above the amide's COM, the amide COM, and the ring's COM. 

<img src="https://github.com/bradenkeiser/piamide_mdsimulation/assets/108278411/9652180e-955e-4ffc-bff3-5e9530f7a9ef" width="450" height="450">

A new configuration is being worked with that omits T-shaped conformations from registering as positive counts. This follows the ring's COM to a point on the normal plane facing the amide region and to the amide's COM. Cutoffs are 120 degrees and 40 degrees if the normal point happens to be calculated for the opposite face of the amide ring through the trajectory frames. 

<img src="https://github.com/bradenkeiser/piamide_mdsimulation/assets/108278411/7da90b2c-c099-4a18-acaa-3a1b2ade406f" width="450" height="450">
