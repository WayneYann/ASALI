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
class equationSystem
{

public:

    equationSystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&          thermodynamicsMap, 
                   OpenSMOKE::KineticsMap_CHEMKIN<double>&                kineticsMap,
                   OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&  thermodynamicsSurfaceMap, 
                   OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped<double>& kineticsSurfaceMap,
                   OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&     transportMap);

    #include "vector.h"

    void setHomogeneusReactions(const bool flag)    {homogeneusReactions_   = flag;}

    void setHeterogeneusReactions(const bool flag)  {heterogeneusReactions_ = flag;}
    
    void setEnergy(const bool flag)                 {energyEquation_        = flag;}

    void setVolume(const double V);
    
    void setArea(const double A);

    void setPressure(const double P);

    void setTemperature(const double T);
    
    void setCatalystLoad(const double alfa);
    
    void setResolutionType(const std::string resolution);

    std::vector<double> getTime()        const {return Time_;};
    std::vector<double> getMass()        const {return Mass_;};
    std::vector<double> getVolume()      const {return Volume_;};
    std::vector<double> getPressure()    const {return Pressure_;};
    std::vector<double> getTemperature() const {return Temperature_;};
    std::vector<double> getSpecie()      const {return Specie_;};
    std::vector<double> getSite()        const {return Site_;};
    std::vector<double> getPhase()       const {return Phase_;};

    unsigned int NumberOfEquations()     const {return NE_;};

    unsigned int OdeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

    unsigned int OdePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

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
    double QRGas_;
    double QRSurface_;
    double t_;

    unsigned int NC_;
    unsigned int SURF_NP_;
    unsigned int SURF_NC_;
    unsigned int NE_;
    
    unsigned int iterator_;

    std::string resolution_;

    bool homogeneusReactions_ ;
    bool heterogeneusReactions_;
    bool energyEquation_;

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

    std::vector<double> Time_;
    std::vector<double> Mass_;
    std::vector<double> Volume_;
    std::vector<double> Pressure_;
    std::vector<double> Temperature_;
    std::vector<double> Specie_;
    std::vector<double> Site_;
    std::vector<double> Phase_;

};


equationSystem::equationSystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&          thermodynamicsMap, 
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
        QRGas_     = 0.;
        QRSurface_ = 0.;

        iterator_  = 0;
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


void equationSystem::setPressure(const double P)
{
    P_ = P;
}

void equationSystem::setTemperature (const double T)
{
    T_ = T;
}

void equationSystem::setVolume (const double V)
{
    V_ = V;
}

void equationSystem::setArea (const double A)
{
    A_ = A;
}

void equationSystem::setCatalystLoad(const double alfa)
{
    alfa_ = alfa;
}

void equationSystem::setResolutionType(const std::string resolution)
{
    resolution_ = resolution;
}

unsigned int equationSystem::OdeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
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

unsigned int equationSystem::OdePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
    t_ = t;
    #include "printIntegration.H"
    return 0;
}

}
