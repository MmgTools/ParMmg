# Basic use of the **parmmg** library

## I/ Implementation
To call the **parmmg** library, you must:  
  1. build mesh and sol at PMMG format;
  2. call the ParMMG library;
  3. get the final mesh and sol.

### sequential_IO  

   The test cases of this folder start and return a sequential mesh (on rank 0).
  
#### automatic_IO  
  We read mesh and solution files using the **PMMG_loadMesh** and **PMMG_loadSol** functions.
  Results are saved using **PMMG_saveMesh** and **PMMG_saveSol** functions.

#### manual_IO  
  The mesh and solution are hard coded.    
  They are build in PMMG format using API functions and are recovered by the same way.  
  We show how to recover the mesh/sol then we write it in a file.

## parallel_IO  
  The test cases of this folder start and return a distributed mesh.

#### automatic_IO  (not yet implemented)
  Each MPI process read the mesh and the solution files using the **PMMG_loadMesh**
  and **PMMG_loadSol** functions. Then we read the communication data using the
  **PMMG_loadComm** function.
  Results and communicators are saved using **PMMG_saveMesh**, **PMMG_saveSol**
  and **PMMG_saveComm** functions.

#### manual_IO  
  The mesh, solution and communicators are hard coded on each processor.    
  They are build in PMMG format using API functions and are recovered by the
  same way.  
  We show how to recover the mesh/sol and the communicators then we write it in
  a file.

## II/ Compilation
  1. Build and install the **mmg3d** shared or static library. We suppose in the following that you have installed the **mmg3d** library in the **_$CMAKE_INSTALL_PREFIX_** directory (see the [installation](https://github.com/MmgTools/Mmg/wiki/Setup-guide#iii-installation) section of the setup guide);
  2. Build and install the **parmmg** shared or static library. We suppose in the following that you have installed the **pmmg** library in the **_$CMAKE_INSTALL_PREFIX_** directory;
  2. compile the main.c file specifying:
    * the **mmg3d** and **parmmg** include directories with the **-I** option;
    * the **mmg3d** and **parmmg** library locations with the **-L** option;
    * the **mmg3d** and **parmmg** library names with the **-l** option;
    * with the shared library, you must add the ***_$CMAKE_INSTALL_PREFIX_** directory to your **LD_LIBRARY_PATH**.

> Example 1  
>  Command line to link the application with the **parmmg** static library
> ```Shell
> gcc -I$CMAKE_INSTALL_PREFIX/include/parmmg -I$CMAKE_INSTALL_PREFIX/include/mmg/mmg3d main.c -L$CMAKE_INSTALL_PREFIX/lib -lparmmg -lmmg3d -lm
> ```

> Example 2  
>  Command line to link the application with the **parmmg** shared library:  
> ```Shell
> gcc -I$CMAKE_INSTALL_PREFIX/include/parmmg -I$CMAKE_INSTALL_PREFIX/include/mmg/mmg3d main.c -L$CMAKE_INSTALL_PREFIX/lib -lparmmg
> export LD_LIBRARY_PATH=$CMAKE_INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
> ```
