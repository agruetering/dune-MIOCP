# dune-MIOCP

This package contains the implementation of the outer approximation algorithm in [2], 
using the DUNE-library [3] for the discretization of the PDE, 
and an interface to use the MINLP solver Gurobi 9.1.2 [9] to obtain exact optimal solutions for comparison; see Section 4 in [1]. 

See the file Overview.pdf to get information about the different header and source code files.

Dependencies
------------

The currect development version of dune-MIOCP depends on the following software: 

* PDELab library [5], the needed core libraries (dune-common, dune-geometry, dune-grid, dune-istl, dune-localfunctions) [4] 
* and the extension libraries (dune-functions, dune-typetree) [6]; all in version 2.7.0.  
* The dune-uggrid module (release/2.7.0) [7]. 
* Gurobi 9.1.2 [9]
* A compiler with support for C++17, at least GCC >= 7.1 or clang >= 5.0 
* CMake >= 3.1

Getting started
---------------

Download all needed DUNE modules to your computer and extract them in one common directory. 

To compile the DUNE modules run 
    ./dune-common/bin/dunecontroll all 
The command configures and builds all installed dune modules as well as all dune modules (not installed) 
which sources reside in a subdirectory of the current directory. 

To separately configure and build all modules run the commands
    ./dune-common/bin/dunecontroll cmake 
    ./dune-common/bin/dunecontroll make 

If you want to compile only a specific DUNE module X, then run
    ./dune-common/bin/dunecontroll --only=X all 
    
If you want to compile a specific DUNE module X and the modules it depends on, you must run 
    ./dune-common/bin/dunecontroll --module=X all 

If you'll have to provide additional information to dunecontrol
(e.g., compilers, configure options) and/or make options, the most convenient way is to use options files, 
which specify the options via the variable
    CMAKE_FLAGS=<flags>
    
An example of an option file example.opts is
    CMAKE_FLAGS=" -DCMAKE_CXX_COMPILER=g++-7 -DCMAKE_INSTALL_PREFIX='/install/path' -DCMAKE_CXX_FLAGS='-Wall -pedantic' "

which
  * uses a specific compiler,
  * installs the files to a custom directory, default is /usr/local/bin and 
  * uses specific compiler flags.
  * 
You can pass the opts file to dunecontrol via the --opts option, e.g.,
    ./dune-common/bin/dunecontrol --opts=example.opts all

See [8] for a more comprehensive introduction to the dune build system.

An installation guide for Gurobi 9.1.2 can be found on [10]. 

Links 
-----

 [1]: C. Buchheim, A. Grütering, and C. Meyer, Parabolic optimal control problems with combinatorial switching constraints – Part I: Convex relaxations, arXiV preprint arXiv:2203.07121 (2022)

 [2]: C. Buchheim, A. Grütering, and C. Meyer, Parabolic optimal control problems with combinatorial switching constraints – Part II: Outer approximation algorithm, arXiV preprint arXiv:2204.07008 (2022)
 
 [3]: http://www.dune-project.org
 
 [4]: https://www.dune-project.org/groups/core
 
 [5]: https://www.dune-project.org/modules/dune-pdelab
 
 [6]: https://www.dune-project.org/groups/extension
 
 [7]: https://www.dune-project.org/modules/dune-uggrid
 
 [8]: https://www.dune-project.org/doc/installation
 
 [9]: https://www.gurobi.com
 
 [10]: https://www.gurobi.com/documentation/quickstart.html
