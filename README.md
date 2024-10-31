# Parmmg - Parallel 3D volume remesher
[![Build Status](https://ci.inria.fr/parmmg/buildStatus/icon?job=ParMmg-multiConf)](https://ci.inria.fr/parmmg/job/ParMmg-multiConf/)

ParMmg is an open source software for parallel mesh adaptation of 3D volume meshes.
The parallel mesh adaptation algorithm is based on iterative remeshing and repartitioning of the distributed mesh (see [1](https://hal.inria.fr/hal-02386837)). It is inspired by the [Yales2](https://www.coria-cfd.fr/index.php/YALES2) algorithm (see [2](https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.4204)).
ParMmg uses the [Mmg](http://mmgtools.org) software to perform the sequential remeshing steps.

:warning: `ParMmg` is built on `Mmg` but it is not officially patronized by the _Mmg Consortium_, thus it does not benefit from the same professional user support. Support is currently done with love by the developers and the community, on their own initiative and free time, and without related funding.

Feel free to contact the maintainers if you think that an advanced usage of `ParMmg` could be useful for your project, or to discuss opportunities for collaborations on specific developments. :warning:

## Get and compile the ParMmg project
### Needed tools
To get and build ParMmg, you will need:
 * **Git**: to download the code you will have to use a git manager. You can install a git manager from the link below but there are many other git clients that you can use:
    * [Official Git client](https://git-scm.com/download) (command line program)
    * [GitKraken](https://www.gitkraken.com/)
    * [SourceTree](https://www.sourcetreeapp.com/)  

    Note that if you uses Microsoft Visual Studio (Windows OS), you can simply activate the Git Module of the application.

  * **CMake** : ParMmg uses the CMake building system that can be downloaded on the
    following web page:
    [https://cmake.org/download/](https://cmake.org/download/). On Windows OS,
    once CMake is installed, please <span style="color:red"> do not forget to
    mark the option: 
    ```
    "Add CMake to the system PATH for all users"
    ```
    </span>  

### ParMmg download and compilation
#### Unix-like OS (Linux, MacOS...)

  1. Get the repository:  
```Shell
      git clone https://github.com/MmgTools/ParMmg.git
```

  The project sources are available under the **_src/_** directory.

  2. By default ParMmg download and install automatically Mmg and Metis (as CMake external projects):
```Shell
      cd ParMmg  
      mkdir build  
      cd build  
      cmake ..  
      make  
      make install
```
  If the `make install` command fail, try to run the `sudo make install` command.
  If you don't have root access, please refers to the [Installation section](https://github.com/MmgTools/Mmg/wiki/Setup-guide#iii-installation) of the [setup guide](https://github.com/MmgTools/Mmg/wiki/Setup-guide#setup-guide) of Mmg and follow the same steps.

  If you don't have internet access and/or want to use your own installation of
  Mmg (resp. Metis), you can disable the automatic download of Mmg setting the
  `DOWNLOAD_MMG` (resp. `DOWNLOAD_METIS`) CMake variable to `OFF`. In this case:

  - Mmg has to be built with the private header installation enabled (turn `ON` the ` MMG_INSTALL_PRIVATE_HEADERS` CMake variable in Mmg at cmake configuration step);
  - you can help CMake to find Mmg by specifying the installation
  directory of Mmg in the `MMG_DIR` variable;
  - you can help CMake to find Metis by specifying the installation directory of Metis in the
  `METIS_DIR` variable.
  
  Example:
  ```Shell
      cd ParMmg  
      mkdir build  
      cd build  
      cmake -DDOWNLOAD_MMG=OFF -DMMG_DIR=~/mmg -DMMG_BUILDDIR=~/mmg/build ..  
      cmake -DDOWNLOAD_METIS=OFF -DMETIS_DIR=~/metis-install ..  
      make  
      make install
```

  The **parmmg** application is available under the `parmmg_O3` command.

Note that if you use some specific options and want to set it easily, you can use a shell script to execute the previous commands. An example is provided in the Mmg wiki [here](https://github.com/MmgTools/mmg/wiki/Configure-script-for-CMake-(UNIX-like-OS)).

#### Windows OS
For now we do not provide Windows portability.

## Documentation
A short documentation is available [here](https://github.com/MmgTools/ParMmg/wiki#user-guide).
To get ParMmg help, run:
```Shell
  parmmg_O3 -h
```
To see the ParMmg default option:
```Shell
  parmmg_O3 -val
```

### Algorithm
The parallel algorithm used in ParMmg has been developed by the Coria team and is described [here](https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.4204). Its current implementation and differences are detailed [here](https://hal.inria.fr/hal-02386837).

### About the team
ParMmg's current developers and maintainers are [Luca Cirrottola](mailto:luca.cirrottola@inria.fr) and [Algiane Froehly](mailto:algiane.froehly@inria.fr).

## License and copyright
Code is under the [terms of the GNU Lesser General Public License](https://raw.githubusercontent.com/MmgTools/mmg/master/LICENSE).

Copyright © Bx INP/Inria/UBordeaux, 2018- .

## Acknowledgements
This work has received funding from the ExaQUte project (http://exaqute.eu/), founded by the European Union’s Horizon 2020 research and innovation programme under grant agreement No 800898.
This work has also received funding from the FUI22 ICARUS project (http://i-carus.fr/), founded by Région Nouvelle- Aquitaine (RNA) under grant agreement 17001260-031 and the Fonds Unique intérministeriel (FUI) under grant agreement No DOS0050687100.

The authors would like to thank the scientific support of the YALES2 development team at the CORIA research center (http://www.coria.fr/) in the implementation of the parallel remeshing algorithm.

## Reference

[Parallel unstructured mesh adaptation using iterative remeshing and repartitioning - _L. Cirrottola and A. Froehly_ - November 29, 2019 - _Inria research report 9307_](https://hal.inria.fr/hal-02386837)


[Three-dimensional adaptive domain remeshing, implicit domain meshing, and applications to free and moving boundary problems - _C. Dapogny, C. Dobrzynski and P. Frey_ - April 1, 2014 - _JCP_](http://www.sciencedirect.com/science/article/pii/S0021999114000266)


[Mesh adaptation for large-eddy simulations in complex geometries. - _P. Benard, G. Balarac, V. Moureau, C. Dobrzynski, G. Lartigue, et al._ - 2016 - _International Journal for Numerical Methods in Fluids_](https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.4204)
