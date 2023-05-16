
Building MSIEVE
===============

Msieve depends on several other packages, which have to be built
before msieve itself can be built. To achieve this the Visual
Studio build files for msieve and these other packages are laid
out in a specific directory structure:

  common_root_directory
    msieve        MSIEVE
    mpir          MPIR
    gmp-ecm       GMP-ECM
    pthreads      PTHREADS

each of which have the following sub-directories:

      build.vc10  - build files for Visual Studio 2010
      build.vc11  - build files for Visual Studio 2012
      bin         - any executable output from the build
      lib         - any static library output from the build
      dll         - any DLL library output from the build

The above directory names are important and this may require 
some renaming after source files are extracted.

MPIR, GMP-ECM and PTHREADS must be built before MSIEVE (although
GMP-ECM is optional).

Building MSIEVE without GMP and GMP-ECM
=======================================

These files build MSIEVE in a way that includes the GMP-ECM library. 
 
If you don't want to use GMP-ECM, you can build MSIEVE without using
this library in the following way:

1. Load the MSIEVE solution - msieve.sln - in the Visual Studio 
   2010 IDE.
2. Open the Property Manager and expand each of the four sub-projects
   in turn and perform the following operation.
3. Open each of the four project configurations in turn and remove
   the 'gmp_ecm_config' item (right click and select 'Remove').

4. The build will then complete without GMP and GMP-ECM support.

Blame Microsoft for the apparent lack of an easier way to do this.

     Brian Gladman, April 2012
