README CHANGES

###############################################################################################################
Modifications from vn1.0.0_monc_belluin:

Rewritting of the IO-server:
- rewritting of timestep callback stage in iobridge.F90
- rewriting of the server function in io_server.F90
Compared to the previous version no federator and no manager are needed. The data are gathered inside one MONC process and send to the IO-server synchronously.

In simplesetup.F90
Flag Casim removed because it deactivates qv,.... if Casim disabled

Rewritting of the system to collect data from configuration files:
- only use of one file global_config.mcf-like which regroup all commands on MONC/CASIM/SOCRATES
- no use of pointer structures which is replaced by a basic array structure with a limited number of cells to scan to find the desired information

Rewritting of the q variables inside MONC, no loops are used as they are explicitely coded

Pdf analysis can be disabled without crashing

Meanprofiles must enabled for profile_diagnostic to work

Removing of old makefile files

Please use the global_config_cloud.mcf file to launch a complete simulation with both warm and cold phases microphysics

###############################################################################################################
Modifications from vn1.0.1_monc_belluin_delbeke:

Add of RH and RHI as variables (available in 3d diagnostics output)

Full modularity of diagnostics output

Casim new physics for ice nucleation (Karcher et al 2021), "iopt_inuc = 11"

Connection of Socrates to MONC

###############################################################################################################
Modifications from vn1.0.2_monc_belluin_delbeke:

Correction of developper errors (lines commented for previous analysis for example)

Completion of modules to interact with all the q species

Casim: autoconversion.F90 and evaporation.F90 Heisenbugs detected during the calculation of dnumber. Resolved by optimizing Casim and Socrates in O1 instead of O3
To compile, please execute the command:  export LC_ALL=C ; fcm make -f ./fcm-make/monc-ubuntu-22.04-gnu-debug.cfg

Application of unitary tests to determine which part of the code execute similarly to vn1.0.0
Determination of the point where the vn1.0.0 et vn1.0.3 could not work the same du to compilation differences

Correction of some conditions to allow some diagnostics to be calculated if not all the diagnostics options are activated

Addition of an option to output diagnostics and checkpoints files in function of time or in function of timestep (default is timestep)

Reactivation of the option restart from a checkpoint. To restart from a checkpoint, a restart configuration file is needed even no changes occured between the different launches. Please use checkpt_run_time_15_restart.nc and global_ref_restart.mcf to test the restart of a simulation using the command:
mpirun -np 5 ./build/bin/monc_driver.exe --reconfig=./global_ref_restart.mcf --checkpoint=./checkpt_run_time_15_restart.nc

Casim: strengthening of conditions of the AbdulRazzakGhan2000_dust scheme in activation.F90, not enough constrained on the minimum mass required

Karcher case: modification of AbdulRazzakGhan2000_dust scheme to avoid coarse_dust mode calculation. Add of a switch option "opt_karcher"

Casim: aerosol parameter for all modes are completely available for the user. A flag "opt_aerosol_param_user" allows to turn it on/off. By default, it is settled to false and the model will use default values.

