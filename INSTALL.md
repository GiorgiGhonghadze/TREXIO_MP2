Installation and Build Instructions

1. Prerequisites

To compile and run this program, you need a standard C compiler, the HDF5 development library and the TREXIO C library.

Install HDF5

Ubuntu / Debian:

sudo apt update
sudo apt install build-essential libhdf5-dev


macOS (Homebrew):

brew install hdf5


Install TREXIO

Download and install TREXIO from source into /usr/local:

wget https://github.com/TREX-CoE/trexio/releases/download/v2.6.1/trexio-2.6.1.tar.gz
tar -zxvf trexio-2.6.1.tar.gz
cd trexio-2.6.1
./configure
make
sudo make install


On Linux, update the dynamic linker cache so libtrexio.so can be located at runtime:

sudo ldconfig


2. Compilation

From the root directory of the repository, compile src/mp2.c and link against libtrexio:

gcc -I/usr/local/include src/mp2.c -L/usr/local/lib -Wl,-rpath,/usr/local/lib -ltrexio -lm -o mp2

3. Running the Program

Execute the compiled binary with the path to a TREXIO input file:

./mp2 data/h2o.h5

Expected Reference Output (data/h2o.h5)

Number of MOs: 24
Number of occupied MOs: 5

Hartree-Fock energy:       -76.0267987082 a.u.
MP2 correlation energy:     -0.2039599741 a.u.
Total MP2 energy:          -76.2307586823 a.u.
