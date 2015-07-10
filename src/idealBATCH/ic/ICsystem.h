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
class ICSystem
#if ASALI_USE_BZZ == 1
 : public BzzOdeSystemObject
#endif
{

public:

    ICSystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&          thermodynamicsMap, 
              OpenSMOKE::KineticsMap_CHEMKIN<double>&                kineticsMap,
              OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&  thermodynamicsSurfaceMap, 
              OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped<double>& kineticsSurfaceMap,
              OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&     transportMap);

    #include "vector.h"

    void setHeterogeneusReactions(const bool flag) { heterogeneusReactions_ = flag; }
    
    void setVolume(const double V);

    void setPressure(const double P);

    void setTemperature(const double T);
    
    void setCatalystLoad(const double alfa);
    
    void setResolutionType(const std::string resolution);

    unsigned int NumberOfEquations()    const {return NE_;};

    void start();
    void end();

    unsigned int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

    unsigned int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

    #if ASALI_USE_BZZ == 1
    virtual void GetSystemFunctions(BzzVector &y, double t, BzzVector &dy);
    virtual void ObjectBzzPrint(void);
    virtual ~ICSystem(){};
    #endif


private:

    double MW_;
    double cTot_;
    double rho_;
    double P_;
    double T_;
    double V_;
    double A_;
    double alfa_;
    double mass_;

    unsigned int NC_;
    unsigned int SURF_NP_;
    unsigned int SURF_NC_;
    unsigned int NE_;

    std::string resolution_;

    bool heterogeneusReactions_;

    OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&            thermodynamicsMap_;             //!< thermodynamic map
    OpenSMOKE::KineticsMap_CHEMKIN<double>&                  kineticsMap_;                   //!< kinetic map
    OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&    thermodynamicsSurfaceMap_;      //!< thermodynamic map
    OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>&          kineticsSurfaceMap_;            //!< kinetic map
    OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&       transportMap_;                  //!< transport map

	OpenSMOKE::OpenSMOKEVectorDouble omega_;
	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;
	OpenSMOKE::OpenSMOKEVectorDouble RfromGas_;
	OpenSMOKE::OpenSMOKEVectorDouble RfromSurface_;
	OpenSMOKE::OpenSMOKEVectorDouble Z_;
	OpenSMOKE::OpenSMOKEVectorDouble Gamma_;
	OpenSMOKE::OpenSMOKEVectorDouble Rsurface_;
	OpenSMOKE::OpenSMOKEVectorDouble RsurfacePhases_;
	OpenSMOKE::OpenSMOKEVectorDouble dummy;

    OpenSMOKE::OpenSMOKEVectorDouble dyOS_;
    OpenSMOKE::OpenSMOKEVectorDouble  yOS_;
};


ICSystem::ICSystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&          thermodynamicsMap, 
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
        MW_        = 0.;
        cTot_      = 0.;
        rho_       = 0.;
        P_         = 0.;
        T_         = 0.;
        V_         = 0.;
        alfa_      = 0.;
        A_         = 0.;

        NC_        = 0;
        SURF_NP_   = 0;
        SURF_NC_   = 0;
        NE_        = 0;

        NC_      = thermodynamicsMap_.NumberOfSpecies();
        SURF_NC_ = thermodynamicsSurfaceMap_.number_of_site_species();
        SURF_NP_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
        NE_      = NC_ + SURF_NC_ + SURF_NP_ + 1 + 1;

		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &RfromGas_, true);
		ChangeDimensions(NC_, &RfromSurface_, true);
		
		ChangeDimensions(SURF_NC_, &Z_, true);
		ChangeDimensions(SURF_NP_, &Gamma_, true);
		ChangeDimensions(SURF_NC_, &Rsurface_, true);
		ChangeDimensions(SURF_NP_, &RsurfacePhases_, true);

        ChangeDimensions(NE_, &dyOS_, true);
        ChangeDimensions(NE_, &yOS_,  true);
    }


void ICSystem::setPressure(const double P)
{
    P_ = P;
}

void ICSystem::setTemperature (const double T)
{
    T_ = T;
}

void ICSystem::setVolume (const double V)
{
    V_ = V;
    A_ = std::pow(V_,(2./3.));

}

void ICSystem::setCatalystLoad(const double alfa)
{
	alfa_ = alfa;
}

void ICSystem::setResolutionType(const std::string resolution)
{
	resolution_ = resolution;
}

#if ASALI_USE_BZZ == 1
void ICSystem::GetSystemFunctions(BzzVector &y, double t, BzzVector &dy)
{
    FromBzzToOS(y,yOS_);
    #include "ICequations.H"
    ObjectBzzPrint();
    FromOSToBzz(dyOS_,dy);
}

void ICSystem::ObjectBzzPrint(void)
{
}
#endif

unsigned int ICSystem::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
    ChangeDimensions(NE_, &yOS_,  true);
    ChangeDimensions(NE_, &dy,    true);
    ChangeDimensions(NE_, &dyOS_, true);

    for (unsigned int i=1;i<=NE_;i++)
        yOS_[i] = y[i];

    #include "ICequations.H"

    for (unsigned int i=1;i<=NE_;i++)
        dy[i] = dyOS_[i];
    
    return 0;
}

unsigned int ICSystem::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
    return 0;
}

void ICSystem::start()
{
    std::cout << "\n######################################" << std::endl;
    std::cout << "# IC system:                  START  #" << std::endl;
    std::cout << "######################################\n" << std::endl;
}

void ICSystem::end()
{
    std::cout << "\n######################################" << std::endl;
    std::cout << "# IC system:                  END    #" << std::endl;
    std::cout << "######################################\n" << std::endl;
}

}
