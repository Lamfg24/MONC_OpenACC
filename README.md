# MONC Introduction: 

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
  - Issues encountered with the use of nv# compilers:
  - Trouble to build libxml2, contact with NVIDIA support: https://forums.developer.nvidia.com/t/libxml2-manual-installation/360152 
  - Trouble with zlib symbols which troubles the HDF5 installation, contact with NVIDIA support: https://forums.developer.nvidia.com/t/zlib-manual-installation-conflict/360230 
  - Trouble with NETCDF libraries, contact with NVIDIA support: https://forums.developer.nvidia.com/t/netcdf-fortran-manual-installation/360255 (solved by myself)
- compilation of MONC
 

# Build librairies: 
  
First after installing the SDK package with CUDA v13 from NVIDIA please load the following module: module load cuda/sdk_13/modulefiles/nvhpc-hpcx-cuda13/25.11 \

################################# 

ZLIB

Download v1.2.12 at https://www.zlib.net/fossils/

Then:\
module load cuda/sdk_13/modulefiles/nvhpc-hpcx-cuda13/25.11
export LDSHARED='mpicc -shared -Wl,-soname,libz.so.1,--version-script,zlib.map' \
CC=mpicc CFLAGS="-fpic" ./configure --prefix=/usr/local/Modules/modulefiles/zlib/1.2.12_nvidia \
make -j 4 \
sudo make install 

################################# 

CURL

Download v8.9.0 at https://curl.se/download/ 

Then:\
module load cuda/sdk_13/modulefiles/nvhpc-hpcx-cuda13/25.11
./configure --prefix=/usr/local/Modules/modulefiles/curl/8.9.0_nvidia --without-ssl --enable-shared=yes --enable-static=yes \
make -j 4 \
sudo make install 

#################################

LIBXML2

Download v2.16.0 at https://gitlab.gnome.org/GNOME/libxml2/-/tree/master 

Then:\
module load cuda/sdk_13/modulefiles/nvhpc-hpcx-cuda13/25.11

In folder ./ open configure and at line: \
13859 \
EXTRA_CFLAGS="${EXTRA_CFLAGS} -pedantic -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes -Wnested-externs -Winline -Wredundant-decls" \
replace by \
EXTRA_CFLAGS="${EXTRA_CFLAGS} " 

13861 \
EXTRA_CFLAGS="${EXTRA_CFLAGS} -Wno-long-long -Wno-format-extra-args" \
replace by \
EXTRA_CFLAGS="${EXTRA_CFLAGS} " 

./configure --prefix=/usr/local/Modules/modulefiles/libxml2/2.16.0_nvidia \
make -j 4 \
sudo make install 

#################################

HDF5

Download v1.10.3 at https://support.hdfgroup.org/downloads/ 

Then: \
module load cuda/sdk_13/modulefiles/nvhpc-hpcx-cuda13/25.11 \
module load zlib/1.2.12_nvidia/zlib_1.2.12_nvidia

CC=mpicc FC=mpif90 CFLAGS="-fPIC" CXXFLAGS="-fPIC" FCFLAGS="-fPIC" ./configure --with-zlib=/usr/local/Modules/modulefiles/zlib/1.2.12_nvidia prefix=/usr/local/Modules/modulefiles/hdf5/1.10.3_nvidia --enable-parallel --enable-fortran 

make -j 4 

Go to hdf5-1.10.3 folder and open libtool \
Go to line 316, then pass a line and insert line: \
PATH="$PATH:/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/bin:/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/clusterkit/bin:/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/tests/imb:/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/sharp/bin:/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/hcoll/bin:/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ucc/bin:/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ucx/bin:/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/compilers/extras/qd/bin:/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/compilers/bin:/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/cuda/13.0/bin:/usr/local/Modules/modulefiles/zlib/1.2.12_nvidia" 

sudo make install 

#################################

PNETCDF

Download v1.12.3 at https://parallel-netcdf.github.io/wiki/Download.html

Then: \
module load cuda/sdk_13/modulefiles/nvhpc-hpcx-cuda13/25.11

MPICC=mpicc MPICXX=mpicxx MPIF77=mpif77 MPIF90=mpif90 CFLAGS="-fPIC" CXXFLAGS="-fPIC" FFLAGS="-fPIC" FCFLAGS="-fPIC" ./configure --prefix=/usr/local/Modules/modulefiles/pnetcdf/1.12.3_nvidia --with-mpi=/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi --enable-shared

make -j 4\
sudo make install

Please open in pnetcdf_install/lib the file libpnetcdf.la \
At line: \
dependency_libs='  /libmpi_usempif08.la  /libmpi_usempi_ignore_tkr.la  /libmpi_mpifh.la -lnvf -lnsnvc  /libmpi.la -lrt -lutil -lz' \
replacy by: \
dependency_libs='  /usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/lib/libmpi_usempif08.la  /usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/lib/libmpi_usempi_ignore_tkr.la  /usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/lib/libmpi_mpifh.la -lnvf -lnsnvc  /usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/lib/libmpi.la -lrt -lutil -lz' \
This procedure is important to enable the installation of NETCDF-C using PNETCDF

#################################

NETCDF-C

Download v4.9.2 at https://downloads.unidata.ucar.edu/netcdf/ 

Then: \
module load cuda/sdk_13/modulefiles/nvhpc-hpcx-cuda13/25.11 \
module load curl/8.9.0_nvidia/curl-8.9.0_nvidia \
module load zlib/1.2.12_nvidia/zlib_1.2.12_nvidia \
module load libxml2/2.16.0_nvidia/libxml2_2.16.0_nvidia \
module load hdf5/1.10.3_nvidia/hdf5_1.10.3_nvidia  \

CC=mpicc CXX=mpicxx CFLAGS="-fPIC" CXXFLAGS="-fPIC" CPPFLAGS='-I/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/include -I/usr/local/Modules/modulefiles/curl/8.9.0_nvidia/include -I/usr/local/Modules/modulefiles/libxml2/2.16.0_nvidia/include -I/usr/local/Modules/modulefiles/hdf5/1.10.3_nvidia/include -I/usr/local/Modules/modulefiles/zlib/1.2.12_nvidia/include' LDFLAGS='-L/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/lib -L/usr/local/Modules/modulefiles/curl/8.9.0_nvidia/lib -L/usr/local/Modules/modulefiles/libxml2/2.16.0_nvidia/lib -L/usr/local/Modules/modulefiles/hdf5/1.10.3_nvidia/lib -L/usr/local/Modules/modulefiles/zlib/1.2.12_nvidia/lib ./configure --prefix=/usr/local/Modules/modulefiles/netcdf-c/4.9.2_nvidia --enable-parallel

make -j 4 \
sudo make install


Please open in netcdf-c_install/lib the file libnetcdf.la \
At line: \
dependency_libs=' -L/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/lib -L/usr/local/Modules/modulefiles/curl/8.9.0_nvidia/lib -L/usr/local/Modules/modulefiles/libxml2/2.16.0_nvidia/lib -L/usr/local/Modules/modulefiles/hdf5/1.10.3_nvidia/lib -L/usr/local/Modules/modulefiles/zlib/1.2.12_nvidia/lib -L/usr/local/Modules/modulefiles/pnetcdf/1.12.3_nvidia/lib /usr/local/Modules/modulefiles/pnetcdf/1.12.3_nvidia/lib/libpnetcdf.la  /libmpi_usempif08.la  /libmpi_usempi_ignore_tkr.la  /libmpi_mpifh.la -lnvf -lnsnvc  /libmpi.la -lrt -lutil /usr/local/Modules/modulefiles/hdf5/1.10.3_nvidia/lib/libhdf5_hl.la /usr/local/Modules/modulefiles/hdf5/1.10.3_nvidia/lib/libhdf5.la -ldl /usr/local/Modules/modulefiles/libxml2/2.16.0_nvidia/lib/libxml2.la -lm /usr/local/Modules/modulefiles/curl/8.9.0_nvidia/lib/libcurl.la -lz' \
replace by: \
dependency_libs=' -L/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/lib -L/usr/local/Modules/modulefiles/curl/8.9.0_nvidia/lib -L/usr/local/Modules/modulefiles/libxml2/2.16.0_nvidia/lib -L/usr/local/Modules/modulefiles/hdf5/1.10.3_nvidia/lib -L/usr/local/Modules/modulefiles/zlib/1.2.12_nvidia/lib -L/usr/local/Modules/modulefiles/pnetcdf/1.12.3_nvidia/lib /usr/local/Modules/modulefiles/pnetcdf/1.12.3_nvidia/lib/libpnetcdf.la /usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/lib/libmpi_usempif08.la  /usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/lib/libmpi_usempi_ignore_tkr.la  /usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/lib/libmpi_mpifh.la -lnvf -lnsnvc  /usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/lib/libmpi.la -lrt -lutil /usr/local/Modules/modulefiles/hdf5/1.10.3_nvidia/lib/libhdf5_hl.la /usr/local/Modules/modulefiles/hdf5/1.10.3_nvidia/lib/libhdf5.la -ldl /usr/local/Modules/modulefiles/libxml2/2.16.0_nvidia/lib/libxml2.la -lm /usr/local/Modules/modulefiles/curl/8.9.0_nvidia/lib/libcurl.la -lz' \
This procedure is important to enable the installation of NETCDF-FORTRAN from NETCDF-C

# MONC Compilation: 
#################################

NETCDF-FORTRAN

Download v4.6.1 at https://downloads.unidata.ucar.edu/netcdf/

Then:\
module load cuda/sdk_13/modulefiles/nvhpc-hpcx-cuda13/25.11 \
module load curl/8.9.0_nvidia/curl-8.9.0_nvidia \
module load zlib/1.2.12_nvidia/zlib_1.2.12_nvidia \
module load libxml2/2.16.0_nvidia/libxml2_2.16.0_nvidia \
module load hdf5/1.10.3_nvidia/hdf5_1.10.3_nvidia \
module load netcdf-c/4.9.2_nvidia/netcdf_c_4.9.2_nvidia \

./configure CC=mpicc FC=mpif90 F77=mpif77 CFLAGS="-fPIC" FFLAGS="-fPIC" FCFLAGS="-fPIC" CPPFLAGS='-I/usr/local/Modules/modulefiles/hdf5/1.10.3_nvidia/include -I/usr/local/Modules/modulefiles/zlib/1.2.12/include -I/usr/local/Modules/modulefiles/netcdf-c/4.9.2_nvidia/include -I/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/include' LDFLAGS='-L/usr/local/Modules/modulefiles/hdf5/1.10.3_nvidia/lib -L/usr/local/Modules/modulefiles/zlib/1.2.12/lib -L/usr/local/Modules/modulefiles/netcdf-c/4.9.2_nvidia/lib -L/usr/local/Modules/modulefiles/cuda/sdk_13/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/lib'  --prefix=/usr/local/Modules/modulefiles/netcdf-fortran/4.6.1_nvidia

make -j 4 \
sudo make install
