# IFEM OpenFrac


## Introduction

This module contains the Fracture Dynamics library and applications built
using the IFEM library.

### Getting all dependencies

1. Install IFEM from https://github.com/OPM/IFEM

### Getting the code

This is done by first navigating to a folder `<App root>` in which you want
the application and typing

    git clone https://github.com/OPM/IFEM-Elasticity
    git clone https://github.com/OPM/IFEM-PoroElasticity
    git clone https://github.com/OPM/IFEM-OpenFrac

The build system uses sibling directory logic to locate the IFEM-Elasticity
and IFEM-PoroElasticity modules.

### Compiling the code

To compile, first navigate to the root catalogue `<App root>`.

1. `cd IFEM-OpenFrac`
2. `mkdir Debug`
3. `cd Debug`
5. `cmake -DCMAKE_BUILD_TYPE=Debug ..`
6. `make`

This will compile the library and applications.
The executables can be found in the 'bin' sub-folder.
Change all instances of `Debug` with `Release` to drop debug-symbols,
and get a faster running code.

### Testing the code

IFEM uses the cmake test system.
To compile and run all regression- and unit-tests, navigate to your build
folder (i.e. `<App root>/IFEM-OpenFrac/Debug`) and type

    make check
