# MONC_OpenACC

MONC is a highly scalable Large Eddy Simulation (LES) model that has been developed to simulate clouds and turbulent flows at high resolution (~ 10s of metres) on large domains (https://code.metoffice.gov.uk/trac/monc/).
MPI is used for parallelism and the actual configuration allows a great scalability on large HPC machines (https://arxiv.org/pdf/2009.12849). The component based architecture is also adapted to GPU porting.

MONC relies on different modules like buoyancy, viscosity, diffusivity and others. Two main modules called CASIM and Socrates are included in this version.
A version of CASIM parallelized with OpenACC is already available: https://code.metoffice.gov.uk/trac/monc/browser/casim/branches/dev/weizhang/um13.6_copy/src?rev=11536

This project aims to provide a new insight of the use of OpenACC programming standard for parallel computing following the previous work of Brown et al 2020 (https://arxiv.org/pdf/2009.12850).

Source: https://code.metoffice.gov.uk/trac/monc/browser/main/branches/dev/lambertdelbeke/del_new_monc?rev=12051
