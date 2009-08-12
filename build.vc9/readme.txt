
Building MSIEVE without GMP and GMP-ECM
=======================================

These files build MSIEVE in a way that includes the GMP-ECM library. 
 
If you don't want to use GMP-ECM, you can build MSIEVE without using
this library in the following way:

1. Load the MSIEVE solution - msieve.sln - in the Visual Studio 
   2008 IDE.
2. Select all the five sub-projects in the Solution Explorer window, 
   right click and select 'Properties'.
3. Select Configuration: 'All Configurations' and Platform:
   'All Platforms' (at the top of the properties dialog).
4. In the 'Inherited Project Property Sheets' property item change:
      
         ..\gmp_ecm_config.vsprops
   to:
         ..\no_gmp_ecm_config.vsprops
   
   to change the build configuration.
   
5. Under the file menu, select 'Save All'.

6. The build will then complete without GMP and GMP-ECM support.

Using GMP-ECM
=============

GMP-ECM relies on a multiple precision integer arithmetic library,
either GMP or MPIR, a windows freiendly fork of GMP.  The Msieve
Visual Studio build is now set up to use the MPIR library as the 
default as thi9s pprovides fast x64 assembler code support on 
Windows that GMP does not provide.  GMP can still be used by 
editing the VC++ gmp_ecm_config.vsprops sheet to set the user
macro to 'gmp' (without the quotes) instead of 'mpir'. 

The layout and naming of the directories holding MSIEVE, MPIR and 
GMP-ECM is assumed to be:

    common_root_directory 
        mpir
        gmp-ecm
        msieve
        
(although the name of the msieve directory doesn't matter). If the 
MPIR and GMP-ECM directories are named and/or located differently
it will then be necessary to rename thyese directories as above 
or modify the gmp_ecm_config.vsprops properties file so that it
uses the revised names.  The 'Additional Include Directories'
(under Properties|C/C++) and the 'Additional Dependencies' (under
Properties|Linker|Input) will need to be changed to match the
names and locations for MPIR and GMP-ECM.

     Brian Gladman, August 2009

