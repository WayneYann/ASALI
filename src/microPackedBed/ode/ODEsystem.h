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

    void setHomogeneusReactions(const bool flag)   { homogeneusReactions_ = flag; }

    void setHeterogeneusReactions(const bool flag) { heterogeneusReactions_ = flag; }

    void setDiffusion(const bool flag)             { gasDiffusion_ = flag; }
    
    void setCorrelation(const std::string correlation);

    void setDiscretizationScheme(const std::string discretizationScheme);

    void setReactorGeometry( const double alfa, const double G, 
                             const double Lcat, const double Linert,
                             const double av);

    void setPackedBedProperties(const double Dp, const double Dt, const double epsi);

    void setHoneyCombProperties(const double epsi, const double Dh);

    void setFeedValue(const double p, const double T0,
                      const OpenSMOKE::OpenSMOKEVectorDouble x0bulk);

    void setGrid(const OpenSMOKE::OpenSMOKEVectorDouble z);

    void setInert(const std::string inert);

    void resize(const unsigned int NP);

    unsigned int BlockDimension()       const {return NB_;};

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

    double MWbulk_;
    double MWwall_;
    double cTotBulk_;
    double cTotWall_;
    double condWall_;
    double rhoWall_;
    double p_;
    double T0_;
    double epsiP_;
    double epsiH_;
    double Dh_;
    double Dp_;
    double Lcat_;
    double Dt_;
    double Linert_;
    double L_;
    double AsymptoticSh_;
    double alfa_;
    double alfaTemp_;
    double G_;
    double SD_;
    double av_;
    double h_;
    double t_;


    unsigned int NC_;
    unsigned int SURF_NP_;
    unsigned int SURF_NC_;
    unsigned int NE_;
    unsigned int NB_;
    unsigned int NP_;

    unsigned int iteration;

    std::string inert_;
    std::string discretizationScheme_;
    std::string correlation_;

    bool homogeneusReactions_ ;
    bool heterogeneusReactions_;
    bool gasDiffusion_;

    OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&            thermodynamicsMap_;             //!< thermodynamic map
    OpenSMOKE::KineticsMap_CHEMKIN<double>&                  kineticsMap_;                   //!< kinetic map
    OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&    thermodynamicsSurfaceMap_;      //!< thermodynamic map
    OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>&          kineticsSurfaceMap_;            //!< kinetic map
    OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&       transportMap_;                  //!< transport map

    OpenSMOKE::OpenSMOKEVectorDouble *omegaBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble *omegaWall_;
    OpenSMOKE::OpenSMOKEVectorDouble *cGas_;
    OpenSMOKE::OpenSMOKEVectorDouble *diffG_;
    OpenSMOKE::OpenSMOKEVectorDouble *jSolid_;

    OpenSMOKE::OpenSMOKEVectorDouble *teta_;

    OpenSMOKE::OpenSMOKEVectorDouble  Tbulk_;
    OpenSMOKE::OpenSMOKEVectorDouble  Th_;
    OpenSMOKE::OpenSMOKEVectorDouble  Tp_;
    OpenSMOKE::OpenSMOKEVectorDouble  cbulk_;
    OpenSMOKE::OpenSMOKEVectorDouble  rhoBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble  cp_;
    OpenSMOKE::OpenSMOKEVectorDouble  condBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble  etaMix_;
    OpenSMOKE::OpenSMOKEVectorDouble  x0bulk_;

    OpenSMOKE::OpenSMOKEVectorDouble  z_;

    OpenSMOKE::OpenSMOKEVectorDouble RfromGas_;
    OpenSMOKE::OpenSMOKEVectorDouble RfromSurface_;
    OpenSMOKE::OpenSMOKEVectorDouble Rsurface_;

    OpenSMOKE::OpenSMOKEVectorDouble Kmat_;

    OpenSMOKE::OpenSMOKEVectorDouble xBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble xWall_;
    OpenSMOKE::OpenSMOKEVectorDouble cBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble cWall_;
    OpenSMOKE::OpenSMOKEVectorDouble RsurfacePhases_;

    OpenSMOKE::OpenSMOKEVectorDouble dyOS_;
    OpenSMOKE::OpenSMOKEVectorDouble  yOS_;

    void MassTransferCoefficient(const double z, const double rho, const double eta, const OpenSMOKE::OpenSMOKEVectorDouble dG);
    OpenSMOKE::OpenSMOKEVectorDouble FirstOrderDerivate (const OpenSMOKE::OpenSMOKEVectorDouble value);

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
    NB_      =  NC_ + NC_ + SURF_NC_ + 1 + 1 + 1;
    NE_      = (NC_ + NC_ + SURF_NC_ + 1 + 1 + 1)*NP_;
    
    omegaBulk_     = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    jSolid_        = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    omegaWall_     = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    teta_          = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    cGas_          = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    diffG_         = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
 
    for (unsigned int i=0;i<NP_;i++)
    {
        ChangeDimensions(NC_, &omegaBulk_[i], true);
        ChangeDimensions(NC_, &omegaWall_[i], true);
        ChangeDimensions(NC_, &cGas_[i],      true);
        ChangeDimensions(NC_, &jSolid_[i],    true);
        ChangeDimensions(NC_, &diffG_[i],     true);
        ChangeDimensions(SURF_NC_, &teta_[i], true);
    }

    ChangeDimensions(NC_, &RfromGas_,      true);
    ChangeDimensions(NC_, &RfromSurface_,  true);
    ChangeDimensions(SURF_NC_, &Rsurface_, true);

    ChangeDimensions(NP_, &Tbulk_,    true);
    ChangeDimensions(NP_, &Th_,       true);
    ChangeDimensions(NP_, &Tp_,       true);
    ChangeDimensions(NP_, &rhoBulk_,  true);
    ChangeDimensions(NP_, &cp_,       true);
    ChangeDimensions(NP_, &condBulk_, true);
    ChangeDimensions(NP_, &etaMix_,   true);
    ChangeDimensions(NP_, &cbulk_,    true);
    
    ChangeDimensions(NC_, &xBulk_, true);
    ChangeDimensions(NC_, &xWall_, true);
    ChangeDimensions(NC_, &cBulk_, true);
    ChangeDimensions(NC_, &cWall_, true);
    ChangeDimensions(NC_, &Kmat_,  true);

    ChangeDimensions(NE_, &dyOS_, true);
    ChangeDimensions(NE_, &yOS_, true);
}

void ODESystem::setInert(const std::string inert)
{
    inert_        = inert;
}

void ODESystem::setDiscretizationScheme(const std::string discretizationScheme)
{
    discretizationScheme_ = discretizationScheme;
}

void ODESystem::setCorrelation(const std::string correlation)
{
    correlation_ = correlation;
}

void ODESystem::setReactorGeometry( const double alfa, const double G, 
                                    const double Lcat, const double Linert,
                                    const double av)
{
    alfaTemp_        = alfa;
    av_              = av;
    Lcat_            = Lcat;
    L_               = Lcat_ + Linert;
    Linert_          = Linert/L_;
    G_               = G;
}

void ODESystem::setPackedBedProperties(const double Dp, const double Dt, const double epsi)
{
    Dp_    = Dp;
    Dt_    = Dt;
    epsiP_ = epsi;
}

void ODESystem::setHoneyCombProperties(const double epsi, const double Dh)
{
    Dh_    = Dh;
    epsiH_ = epsi;
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

void ODESystem::MassTransferCoefficient(const double z, const double rho, const double eta, const OpenSMOKE::OpenSMOKEVectorDouble dG)
{
    if ( correlation_ == "Yoshida" )
    {
        double Re     = G_*Dp_/(epsiP_*eta*(1. - epsiP_)*6.);
        double ReReal = G_*Dp_/(eta*epsiP_);
        OpenSMOKE::OpenSMOKEVectorDouble Sc(NC_);
        OpenSMOKE::OpenSMOKEVectorDouble Sh(NC_);
        OpenSMOKE::OpenSMOKEVectorDouble jM(NC_);
        for (unsigned int i=1;i<=NC_;i++)
        {
            Sc[i] = eta/(rho*dG[i]);
            if ( Re < 50. )
            {
                jM[i]    = 0.91/(std::pow(Re,0.51));
                Sh[i]    = jM[i]*std::pow(Sc[i],(1./3.))*ReReal;
                Kmat_[i] = Sh[i]*dG[i]/Dp_;
            }
            else
            {
                jM[i]    = 0.61/(std::pow(Re,0.41));
                Sh[i]    = jM[i]*std::pow(Sc[i],(1./3.))*ReReal;
                Kmat_[i] = Sh[i]*dG[i]/Dp_;
            }
        }
    }
    else if ( correlation_ == "Petrovic" )
    {
        double Re = G_*Dp_/(eta*epsiP_);
        OpenSMOKE::OpenSMOKEVectorDouble Sc(NC_);
        OpenSMOKE::OpenSMOKEVectorDouble Sh(NC_);
        OpenSMOKE::OpenSMOKEVectorDouble jM(NC_);
        for (unsigned int i=1;i<=NC_;i++)
        {
            Sc[i]    = eta/(rho*dG[i]);
            jM[i]    = 0.357/(epsiP_*std::pow(Re,0.359));
            Sh[i]    = jM[i]*std::pow(Sc[i],(1./3.))*Re;
            Kmat_[i] = Sh[i]*dG[i]/Dp_;
        }
    }
    else if ( correlation_ == "Wakao" )
    {
        double Re = G_*Dp_/(eta*epsiP_);
        OpenSMOKE::OpenSMOKEVectorDouble Sc(NC_);
        OpenSMOKE::OpenSMOKEVectorDouble Sh(NC_);
        for (unsigned int i=1;i<=NC_;i++)
        {
            Sc[i]     = eta/(rho*dG[i]);
            Sh[i]    = 2. + 1.1*std::pow(Re,0.6)*std::pow(Sc[i],(1./3.));
            Kmat_[i] = Sh[i]*dG[i]/Dp_;
        }
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

OpenSMOKE::OpenSMOKEVectorDouble ODESystem::FirstOrderDerivate (const OpenSMOKE::OpenSMOKEVectorDouble value)
{
    OpenSMOKE::OpenSMOKEVectorDouble dvalue_(NP_);

    if ( discretizationScheme_ == "CDS" )
    {
        for (unsigned int k=1;k<=NP_;k++)
        {
            if ( k == 1 )
            {
                dvalue_[k] = ((value[k+1] - value[k])/(z_[k+1] - z_[k]))/L_;
            }
            else if ( k == NP_ )
            {
                dvalue_[k] = ((value[k] - value[k-1])/(z_[k] - z_[k-1]))/L_;
            }
            else
            {
                dvalue_[k] = ((value[k+1] - value[k-1])/(z_[k+1] - z_[k-1]))/L_;
            }
        }
    }
    else if ( discretizationScheme_ == "BDS" )
    {
        for (unsigned int k=1;k<=NP_;k++)
        {
            if ( k == 1 )
            {
                dvalue_[k] = ((value[k+1] - value[k])/(z_[k+1] - z_[k]))/L_;
            }
            else if ( k == NP_ )
            {
                dvalue_[k] = ((value[k] - value[k-1])/(z_[k] - z_[k-1]))/L_;
            }
            else
            {
                dvalue_[k] = ((value[k] - value[k-1])/(z_[k] - z_[k-1]))/L_;
            }
        }
    }
    return dvalue_;
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
    delete [] omegaWall_;
    delete [] jSolid_;
    delete [] teta_;
    std::cout << "\n######################################" << std::endl;
    std::cout << "# ODE system:                 END    #" << std::endl;
    std::cout << "######################################\n" << std::endl;
}

}
