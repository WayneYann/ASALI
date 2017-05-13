# ASALI

1D/2D models for reacting flows with heterogeneous/homogeneous reactions.

The compiled folder contains pre-compiled solvers, which can be used without any external libraries. Instead, the external libraries required for compiling ASALI are:

1/ compulsory libraries

    Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
    RapidXML (http://rapidxml.sourceforge.net/)
    Boost C++ (http://www.boost.org/)
    OpenSMOKE++ (alberto.cuoci@polimi.it)

2/ optional libraries

    Intel MKL (https://software.intel.com/en-us/intel-mkl)
    BzzMath   (http://super.chem.polimi.it/)

The CHEMKIN kinetic schemes required by ASALI can be compiled by using the chemikin-preprocessor distributed with:

    catalyticFOAM (https://www.catalyticfoam.polimi.it)
