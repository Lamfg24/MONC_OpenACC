MONC is a highly scalable Large Eddy Simulation (LES) model that has been developed to simulate clouds and turbulent flows at high resolution (~ 10s of metres) on large domains (https://code.metoffice.gov.uk/trac/monc/). \
\
MPI is used for parallelism and the actual configuration allows a great scalability on large HPC machines (https://arxiv.org/pdf/2009.12849). The component based architecture is also adapted to GPU porting. \
\
MONC relies on different modules like buoyancy, viscosity, diffusivity and others. Two main modules called CASIM and Socrates are included in this version. \
A version of CASIM parallelized with OpenACC is already available: https://code.metoffice.gov.uk/trac/monc/browser/casim/branches/dev/weizhang/um13.6_copy/src?rev=11536 \
\
This project aims to provide a new insight of the use of OpenACC programming standard for parallel computing following the previous work of Brown et al 2020 (https://arxiv.org/pdf/2009.12850). \
\
Source: https://code.metoffice.gov.uk/trac/monc/browser/main/branches/dev/lambertdelbeke/del_new_monc?rev=12051
\
\
Development:
- build of librairies (zlib, curl, libxml2, pnetcdf, netcdf-C, netcdf-Fortran)\
  Issues encountered with the use of nv# compilers:\
  Trouble to build libxml2, contact with NVIDIA support: https://forums.developer.nvidia.com/t/libxml2-manual-installation/360152/3

\
Build librairies: \
First after installing the SDK v13 package from NVIDIA please load the following module: module load cuda/sdk_13/modulefiles/nvhpc-hpcx-cuda13/25.11 \

################################# 

ZLIB

Download v1.2.12 at https://www.zlib.net/fossils/

Then:\
CFLAGS="-fpic" ./configure --prefix=/usr/local/Modules/modulefiles/zlib/1.2.12_nvidia \
make -j 4 \
sudo make install 

################################# 

CURL

Download v8.9.0 at https://curl.se/download/ \
Then:\
./configure --prefix=/usr/local/Modules/modulefiles/curl/8.9.0_nvidia --without-ssl --enable-shared=yes --enable-static=yes \
make -j 4 \
sudo make install 

#################################
