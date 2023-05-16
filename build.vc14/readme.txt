
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

      build.vc14       - build files for Visual Studio 2015
      bin              - any executable output from the build
      lib              - any static library output from the build
      dll              - any DLL library output from the build

In addition MSIEVE has this sub-directory for building the CUDA
version of MSIEVE

      build.cuda.vc14  - build files for Visual Studio 2015

The above directory names are important and this may require 
some renaming after the source files are extracted.

The build fieles for older versions of Visual Studio are still
present but are not being maintained by me.  Morover I am now
only maintaining builds for 64-bit versions of Windows.  

Building MSIEVE (non GPU version)
=================================

MPIR, GMP-ECM and PTHREADS must be built before MSIEVE (GMP-ECM is 
optional - see below).

Open the project file msieve.sln in Visual Studio 2015 and build
the msieve project (this will build he other projects as well if
they are out of date).

The executable will be placed in the msieve/bin directory.

Building MSIEVE without GMP and GMP-ECM
=======================================

These files build MSIEVE in a way that includes the GMP-ECM library. 
 
If you don't want to use GMP-ECM, you can build MSIEVE without using
this library in the following way:

1. Load the MSIEVE solution - msieve.sln - in the Visual Studio 
   2015 IDE.
2. Open the Property Manager and expand each of the four sub-projects
   in turn and perform the following operation.
3. Open each of the four project configurations in turn and remove
   the 'gmp_ecm_config' item (right click and select 'Remove').

4. The build will then complete without GMP and GMP-ECM support.

Building MSIEVE (GPU version)
=============================

The build files for the GPU version of msieve are in the directory

      build.cuda.vc14  - build files for Visual Studio 2015

where the file msieve.sln will load the following projects:

      common
	  gnfs
	  mpqs
	  msieve
	  zlib
	  sort_engine_sm
	  stage1_core_sm

The first five projects are built when the msieve project is
built but the last two have to be built individually as they
have to be configured for the GPU (CUDA) display driver that
is being used.

Each Nvidia disply driver has a compute capability which has
to be set for building both sort_engine_sm and stage1_core_sm.
To do this, open the file gpu_sm.props in the build.cuda.vc14
directory and edit the following two lines to set the digits of 
the compute capability:     

    <CC_major>5</CC_major>
    <CC_minor>0</CC_minor>

for the build.  Then save the updated file and then build 
sort_engine_sm and stage1_core_sm (it may be necessary to
reload these two projects in Visual Studio in order for 
the edited gpu_sm.props file to be reloaded).

The resulting executables wiill be put in the bin directory
with the sort_engine in the sub-directory 'cub'. 

Using GMP-ECM with GMP
======================

GMP-ECM relies on a multiple precision integer arithmetic library,
either GMP or MPIR, a windows friendly fork of GMP.  The Msieve
Visual Studio build is set up to use the MPIR library but GMP
can be used by editing the VC++ gmp_ecm_config.vsprops sheet to
set the user macro to 'gmp' (without the quotes) instead of 'mpir'. 

The layout and naming of the directories holding MSIEVE, MPIR and 
GMP-ECM is assumed to be:

    common_root_directory 
        mpir
        gmp-ecm
        msieve
        
(although the name of the msieve directory doesn't matter). If the 
MPIR and GMP-ECM directories are named and/or located differently
it will then be necessary to rename these directories as above 
or modify the gmp_ecm_config.vsprops properties file so that it
uses the revised names.  The 'Additional Include Directories'
(under Properties|C/C++) and the 'Additional Dependencies' (under
Properties|Linker|Input) will need to be changed to match the
names and locations for MPIR and GMP-ECM.

    Brian Gladman, September 2016
