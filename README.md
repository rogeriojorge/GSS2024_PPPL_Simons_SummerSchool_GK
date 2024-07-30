# GSS2024_PPPL_Simons_SummerSchool_GK
 Files for the Simons/PPPL Summer School Lectures on Gyrokinetics

###
These are the scripts used to show a stellarator fluxtube geometry and a corresponding linear gyrokinetic simulation.
It uses VMEC output (wout) files for magnetic field equilibria (in `equilibria` folder),
the [GS2 code](https://bitbucket.org/gyrokinetics/gs2/src/master/) for the gyrokinetic simulations,
the [SIMSOPT code](https://github.com/hiddenSymmetries/simsopt) to obtain the geometry coefficients,
which requires the [VMEC2000 code](https://github.com/hiddenSymmetries/VMEC2000).

To run this repository, use
`python main.py --type X`
where `X` is the stellarator configuration.
Type 1 is HSX, type 2 is W7-X, type 3 is QI, type 4 is QH, type 5 is QA.

The code will create figures in a newly created `output` folder.

The GS2 executable location should be added in the `main.py` file, at the following line:
`gs2_executable = f'{home_directory}/local/gs2/bin/gs2'`

A Mathematica notebook is added to show how to easily check the simulation output parameters.