/*##############################################################################################
#                                                                                              #
#     #############       #############       #############       ####                ####     #
#    #             #     #             #     #             #     #    #              #    #    #
#    #    #####    #     #    #########      #    #####    #     #    #              #    #    #
#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #
#    #    #####    #     #    #              #    #####    #     #    #              #    #    #
#    #             #     #    #########      #             #     #    #              #    #    #
#    #             #     #             #     #             #     #    #              #    #    #
#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #
#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #
#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #
#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #
#     ####     ####       #############       ####     ####       #############       ####     #
#                                                                                              #
#   Department of Energy                                                                       #
#   Politecnico di Milano                                                                      #
#   Author: Stefano Rebughini <stefano.rebughini@polimi.it>                                    #
#                                                                                              #
################################################################################################
#                                                                                              #
#   License                                                                                    #
#                                                                                              #
#   This file is part of ASALI.                                                                #
#                                                                                              #
#   ASALI is free software: you can redistribute it and/or modify                              #
#   it under the terms of the GNU General Public License as published by                       #
#   the Free Software Foundation, either version 3 of the License, or                          #
#   (at your option) any later version.                                                        #
#                                                                                              #
#   ASALI is distributed in the hope that it will be useful,                                   #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                              #
#   GNU General Public License for more details.                                               #
#                                                                                              #
#   You should have received a copy of the GNU General Public License                          #
#   along with ASALI. If not, see <http://www.gnu.org/licenses/>.                              #
#                                                                                              #
##############################################################################################*/

namespace ASALI
{
class ODESystem
#if ASALI_USE_BZZ == 1
 : public BzzOdeSystemObject
#endif
{

public:

    ODESystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&          thermodynamicsMap, 
              OpenSMOKE::KineticsMap_CHEMKIN<double>&                kineticsMap,
              OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&  thermodynamicsSurfaceMap, 
              OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped<double>& kineticsSurfaceMap,
              OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&     transportMap);

    #include "vector.h"

    void setReactorGeometry(const double D,const double L);

    void setFeedValue(const double p, const double T0,
                      const OpenSMOKE::OpenSMOKEVectorDouble x0bulk);

    void setGrid(const OpenSMOKE::OpenSMOKEVectorDouble z);

    void setInert(const std::string inert);

    void resize(const unsigned int NP);

    unsigned int BlockDimension()    const {return NB_;};

    unsigned int NumberOfEquations()    const {return NE_;};

    void start();
    void end();

    unsigned int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

    unsigned int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

    #if ASALI_USE_BZZ == 1
    virtual void GetSystemFunctions(BzzVector &y, double t, BzzVector &dy);
    virtual void ObjectBzzPrint(void);
    virtual ~ODESystem(){};
    #endif


private:

    double p_;
    double T0_;
    double epsi_;
    double D_;
    double L_;
    double t_;


    unsigned int NC_;
    unsigned int SURF_NP_;
    unsigned int SURF_NC_;
    unsigned int NE_;
    unsigned int NB_;
    unsigned int NP_;

    unsigned int iteration;

    std::string inert_;


    OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&            thermodynamicsMap_;             //!< thermodynamic map
    OpenSMOKE::KineticsMap_CHEMKIN<double>&                  kineticsMap_;                   //!< kinetic map
    OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&    thermodynamicsSurfaceMap_;      //!< thermodynamic map
    OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>&          kineticsSurfaceMap_;            //!< kinetic map
    OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&       transportMap_;                  //!< transport map

    OpenSMOKE::OpenSMOKEVectorDouble *omegaBulk_;

    OpenSMOKE::OpenSMOKEVectorDouble *teta_;


    OpenSMOKE::OpenSMOKEVectorDouble  x0bulk_;

    OpenSMOKE::OpenSMOKEVectorDouble  z_;

    OpenSMOKE::OpenSMOKEVectorDouble dyOS_;
    OpenSMOKE::OpenSMOKEVectorDouble  yOS_;
};


ODESystem::ODESystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&          thermodynamicsMap, 
                     OpenSMOKE::KineticsMap_CHEMKIN<double>&                kineticsMap,
                     OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&  thermodynamicsSurfaceMap, 
                     OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped<double>& kineticsSurfaceMap,
                     OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&     transportMap):
    thermodynamicsMap_(thermodynamicsMap), 
    kineticsMap_(kineticsMap),
    thermodynamicsSurfaceMap_(thermodynamicsSurfaceMap), 
    kineticsSurfaceMap_(kineticsSurfaceMap),
    transportMap_(transportMap)
    {
    }

void ODESystem::resize(const unsigned int NP)
{
    NP_      = NP;
    NC_      = thermodynamicsMap_.NumberOfSpecies();
    SURF_NC_ = thermodynamicsSurfaceMap_.number_of_site_species();
    SURF_NP_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
    NB_      = NC_ + SURF_NC_;
    NE_      = (NC_ + SURF_NC_)*NP_;
    
    omegaBulk_     = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    teta_          = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];

    for (unsigned int i=0;i<NP_;i++)
    {
        ChangeDimensions(NC_, &omegaBulk_[i], true);
        ChangeDimensions(SURF_NC_, &teta_[i], true);
    }
    ChangeDimensions(NE_, &dyOS_, true);
    ChangeDimensions(NE_, &yOS_, true);
}

void ODESystem::setInert(const std::string inert)
{
    inert_        = inert;
}

void ODESystem::setReactorGeometry(const double D,const double L)
{
    D_  = D;
    L_  = L;
}

void ODESystem::setFeedValue(const double p, const double T0,
                             const OpenSMOKE::OpenSMOKEVectorDouble x0bulk)
{
    p_   = p;
    T0_  = T0;
    ChangeDimensions(x0bulk.Size(), &x0bulk_, true);
    for (unsigned int j=1;j<=x0bulk.Size();j++)
    {
        x0bulk_[j] = x0bulk[j];
    }
}

void ODESystem::setGrid(const OpenSMOKE::OpenSMOKEVectorDouble z)
{
    ChangeDimensions(z.Size(), &z_, true);
    for (unsigned int k=1;k<=z_.Size();k++)
    {
        z_[k] = z[k]/L_;
    }
}

#if ASALI_USE_BZZ == 1
void ODESystem::GetSystemFunctions(BzzVector &y, double t, BzzVector &dy)
{
    FromBzzToOS(y,yOS_);
    #include "ODEequations.H"
    ObjectBzzPrint();
    FromOSToBzz(dyOS_,dy);
}

void ODESystem::ObjectBzzPrint(void)
{
    #include "printIntegration.H"
}
#endif

unsigned int ODESystem::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
    ChangeDimensions(NE_, &yOS_,  true);
    ChangeDimensions(NE_, &dy,    true);
    ChangeDimensions(NE_, &dyOS_, true);

    for (unsigned int i=1;i<=NE_;i++)
        yOS_[i] = y[i];

    #include "ODEequations.H"

    for (unsigned int i=1;i<=NE_;i++)
        dy[i] = dyOS_[i];
    
    return 0;
}

unsigned int ODESystem::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
    t_ = t;
    #include "printIntegration.H"
    return 0;
}

void ODESystem::start()
{
    std::cout << "\n######################################" << std::endl;
    std::cout << "# ODE system:                 START  #" << std::endl;
    std::cout << "######################################\n" << std::endl;
}

void ODESystem::end()
{
    delete [] omegaBulk_;
    delete [] teta_;
    std::cout << "\n######################################" << std::endl;
    std::cout << "# ODE system:                 END    #" << std::endl;
    std::cout << "######################################\n" << std::endl;
}

}
