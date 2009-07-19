
Building MSIEVE without GMP and GMP-ECM
=======================================

These files build MSIEVE in a way that includes support for both GMP 
and GMP-ECM.  If you don't have GMP and GMP-ECM available as Visual 
Studio build projects you can build MSIEVE without using these 
libraries in the following way:

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

The layout and naming of the directories holding MSIEVE, GMP and 
GMP-ECM is assumed to be:

    common_root_directory 
        gmp
        gmp-ecm
        msieve
        
(although the name of the msieve directory doesn't matter). If the 
GMP and GMP-ECM directories are named and/or located differently
it will then be necessary to modify the gmp_ecm_config.vsprops
properties file accordingly. The 'Additional Include Directories'
(under Properties|C/C++) and the 'Additional Dependencies' (under
Properties|Linker|Input) will need to be changed to match the
names and locations for GMP and GMP-ECM.

     Brian Gladman, November 2008
