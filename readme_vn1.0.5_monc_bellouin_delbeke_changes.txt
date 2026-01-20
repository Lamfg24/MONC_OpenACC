README CHANGES

###############################################################################################################
Modifications from vn1.0.0_monc_belluin:

Rewritting of the IO-server:
- rewritting of timestep callback stage in iobridge.F90
- rewriting of the server function in io_server.F90
Compared to the previous version no federator and no manager are needed. The data are gathered inside one MONC process and send to the IO-server synchronously.

In simplesetup.F90
Flag Casim removed because it deactivates qv,.... if Casim disabled. Not necessary as all q and n species are desired to be activated

Rewritting of the system to collect data from configuration files:
- only use of one file global_config.mcf-like which regroup all commands on MONC/CASIM/SOCRATES
- no use of pointer structures which is replaced by a basic array structure with a limited number of cells to scan to find the desired information

Rewritting of the q variables inside MONC, no loops are used as they are explicitely coded. Some parts are missing due to prudent completion of the code. Following versions will complete even finish these parts of the code.

Pdf analysis can be disabled without crashing

Meanprofiles must enabled for profile_diagnostic to work

Removing of old makefile files

###############################################################################################################
Modifications from vn1.0.1_monc_belluin_delbeke:

Correction of developper errors (for example in pressuresource: reactivation of some source term for wind components)

Modifications of modules to interact with all hydrometeor and aerosol species:
- buoyancy: addition of missing mass mixing ratio and number concentration species
- damping: addition of missing mass mixing ratio and number concentration species
- set_consistent_lowbc: addition of missing number concentration species
- swapsmooth: addition of missing mass mixing ratio and number concentration species
- tvdadvection: addition of missing mass mixing ratio and number concentration species

Heisenbug detected in gridmanager.F90 from vn1.0.0 about the defintion of the vertical variables like prefn. It concerns module which use these variable like forcing.
The bug has been corrected since vn1.0.1 using more constrained compilation option but it will generate some different results it the user want to compare vn1.0.0 and vn1.0.1 (or higher version)

Addition of RH and RHI as variables (available in 3d diagnostics output)

Full modularity of diagnostics output

Casim new physics for ice nucleation (Karcher et al 2021), "iopt_inuc = 11"

Connection of Socrates to MONC

Completion of modules to interact with all the q species:

###############################################################################################################
Modifications from vn1.0.2_monc_belluin_delbeke:

Correction of developper errors (for example in swapsmooth: supression of the double processing of zqCoarseDust)

Completion of modules to interact with all hydrometeor and aerosol species:
- pwadvection: addition of missing number concentration species
- stepfields: addition of missing reset and remove_negative for q and n

Documentation of the Heisenbug detected in gridmanager.F90 about the defintion of the vertical variables like prefn.

Casim: autoconversion.F90 and evaporation.F90 Heisenbugs detected during the calculation of dnumber. Resolved by optimizing Casim and Socrates in O1 instead of O3
To compile, please execute the command:  export LC_ALL=C ; fcm make -f ./fcm-make/monc-ubuntu-22.04-gnu-debug.cfg

Casim: addition of autoconversion parameters in global_config file, DImax and DI2S

Application of unitary tests to determine which part of the code execute similarly to vn1.0.0
Determination of the point where the vn1.0.0 et vn1.0.3 could not work the same du to compilation differences

Correction of some conditions to allow some diagnostics to be calculated if not all the diagnostics options are activated

Addition of an option to output diagnostics and checkpoints files in function of time or in function of timestep (default is timestep)

Reactivation of the option restart from a checkpoint. To restart from a checkpoint, a restart configuration file is needed even no changes occured between the different launches. Please use checkpt_run_time_15_restart.nc and global_ref_restart.mcf to test the restart of a simulation using the command:
mpirun -np 5 ./build/bin/monc_driver.exe --reconfig=./global_ref_restart.mcf --checkpoint=./checkpt_run_time_15_restart.nc

Casim: strengthening of conditions of the AbdulRazzakGhan2000_dust scheme in activation.F90, not enough constrained on the minimum mass required

Karcher case: modification of AbdulRazzakGhan2000_dust scheme to avoid coarse_dust mode calculation. Addition of a switch option "opt_karcher"

Casim: aerosol parameter for all modes are completely available for the user. A flag "opt_aerosol_param_user" allows to turn it on/off. By default, it is settled to false and the model will use default values (all modes can now be personalized).

Diagnostics: addition of options in global_config file to select a specific layer of the domain in order to get the corresponding water path for the different hydrometeors species

Gridmanager: addition of a initial profile of vertical wind speed w, see global_config file

constant_w_profile_enabled: an option to force MONC to reach the intialization profile of w wind component at each timestep, see global_config

Casim: implementation of the homogeneous freezing scheme of Seifert and Beheng (2006), user can choose between Wisener and Seifert scheme in global_config file

Io-server - iobridge: addition of a security to ensure the completion of netcdf writing (synchronous)

Condiational compilation for the different hydrometeor and aerosol species for the modules which deal with them. It allows more flexibility for users who want to allow/cut off specific part of the code concerning the hydrometeor and aerosol species. To operate, conditionner are available in comp-gnu-4.4.7-debug.cfg with the complete reference flags described in comments

###############################################################################################################
Modifications from vn1.0.3_monc_belluin_delbeke:

Clarification of ioserver.F90 and transfer of the netcdf writing functions into ioserver_netcdf_writer.F90

Connection of flux_budget module to MONC and cleaning among available fluxes. Modification/addition to advection schemes to get structured variables (like u_advection)

Connection of conditionnal diagnostics

Connection of pdf diagnostics

Generation of a folder containing checkpoint and diagnostics subfolders allowing a better tracability and management of data results.
Results folder title format is: YYYY_(M)M_(D)D_(H)H_(M)M_(S)S_(re)start

Revision/completion of MONC documentation

Revision and cleaning of the Forcing.F90 time-varying forcing (with and without forcing file). Time varying forcing of theta, q and subsidence where reviewed and geostrophic wind comopenents are added. A prototype python file has been developped to help users to create their own forcing netcdf files following MONC nomenclature

Reorganization of time-varying and constant forcing option for better readability in file forcing.F90

Correction and reactivation of time-varying surface values in file setfluxlook.F90 (with and without surface forcing file). A Python file has been developped to generate such surface forcing file

A tool folder has been created to regroup all Python routine prototype for users to develop their own netcdf forcing files

Generation of a diagnostic routine capable to save at n timestep/time prognostic and diagnostic variables and mean them to obtain a average of the variable on a certain period/number of timesteps

Discovery of a non-consertive issue when using anelastic equations for density and using the w forcing at each timestep. A correcting factor has been applied on tvd scheme (conservative by nature) precisely on the calculation of q scalars transported by w. Not definitive solution

Addition of a feature capable to choose to add random noise at initialization or at a time the user chooses

Addition of a timer for application of the w constant profile applied all along the simulation

###############################################################################################################
Modifications from vn1.0.4_monc_belluin_delbeke:

Correction for density, thref, prefn allocation as previously these variables were defined in grid.F90 and state.F90...only remains allocation in grid.F90 in vertical structure

Addition of a starting time and duration time of forcings allowing more flexibility to apply forcing during a simulation

Correction of the determination of ice crystal diameter by including the ice number concentration in the formula in disbribution.F90

Fix of haloswapper.F90 concerning "copy_fields_to_halo_buffer" function, wrong order of q species corrected
Which order is: qv,ql,qr,qi,qs,qg,qAikt,qAccumSol,qAccumInsol,qCoarse,qDust and same for n species

Optional: addition of WENO5 scheme and improvment of halocommunication.F90 to handle 3 points in order to make WENO5 functionnal
Validated for wind components and temperature but not for q and n species (still conservation issues). All files available in components/weno5/src/

Addition of scenario 12 in ice_nucleation.F90

###############################################################################################################
Modifications from vn1.0.5_monc_belluin_delbeke:

Addition of a duration for application of random noise during the simulation at a certain time


