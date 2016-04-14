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
class BVPSystem
#if ASALI_USE_BZZ == 1
 : public BzzDaeSystemObject
#endif
{

public:

    BVPSystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&          thermodynamicsMap, 
              OpenSMOKE::KineticsMap_CHEMKIN<double>&                kineticsMap,
              OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&  thermodynamicsSurfaceMap, 
              OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped<double>& kineticsSurfaceMap,
              OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&     transportMap);

    #include "vector.h"

    void setHoneycombSolid(const double rhoSolid, const double condSolid, const double cpSolid);
    
    void setPackedBedSolid(const double rhoSolid, const double condSolid, const double cpSolid);

    void setHomogeneusReactions(const bool flag)   { homogeneusReactions_ = flag; }

    void setHeterogeneusReactions(const bool flag) { heterogeneusReactions_ = flag; }

    void setEnergyEquation(const bool flag)        { energy_ = flag; }
    
    void setDiffusion(const bool flag)             { gasDiffusion_ = flag; }
 
    void setCorrelation(const std::string correlation);

    void setDiscretizationScheme(const std::string discretizationScheme);

    void setReactorGeometry( const double alfa, const double G, 
                             const double Lcat, const double Linert,
                             const double av,   const double aex);

    void setPackedBedProperties(const double Dp, const double Dt, const double epsi);

    void setHoneyCombProperties(const double epsi, const double Dh);

    void setFeedValue(const double p, const double T0,
                      const OpenSMOKE::OpenSMOKEVectorDouble x0bulk);

    void setGrid(const OpenSMOKE::OpenSMOKEVectorDouble z);

    void setInert(const std::string inert);

    void resize(const unsigned int NP);

    unsigned int BlockDimension()        const {return NB_;};

    unsigned int NumberOfEquations()    const {return NE_;};

    void start();
    void end();

    unsigned int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

    unsigned int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

    void AlgebraicEquations(OpenSMOKE::OpenSMOKEVectorBool& algebraic);

    OpenSMOKE::OpenSMOKEVectorDouble getVelocity()    const {return v_;};

    #if ASALI_USE_BZZ == 1
    void AlgebraicEquations(BzzVectorInt& algebraic);
    virtual void GetSystemFunctions(BzzVector &y, double t, BzzVector &dy);
    virtual void ObjectBzzPrint(void);
    virtual ~BVPSystem(){};
    #endif

private:

    double MWbulk_;
    double MWwall_;
    double cTotBulk_;
    double cTotWall_;
    double rhoWall_;
    double p_;
    double T0_;
    double epsiH_;
    double epsiP_;
    double Dh_;
    double Dp_;
    double Dt_;
    double Lcat_;
    double Linert_;
    double L_;
    double AsymptoticSh_;
    double alfa_;
    double alfaTemp_;
    double G_;
    double SD_;
    double av_;
    double aex_;
    double h_;
    double U_;
    double t_;

    double rhoH_;
    double cpH_;
    double condH_;

    double rhoP_;
    double cpP_;
    double condP_;

    double QfromGas_;
    double QfromSurface_;

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
    bool energy_;
    bool gasDiffusion_;

    OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&              thermodynamicsMap_;             //!< thermodynamic map
    OpenSMOKE::KineticsMap_CHEMKIN<double>&                    kineticsMap_;                   //!< kinetic map
    OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&      thermodynamicsSurfaceMap_;      //!< thermodynamic map
    OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>&            kineticsSurfaceMap_;            //!< kinetic map
    OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&         transportMap_;                  //!< transport map

    OpenSMOKE::OpenSMOKEVectorDouble *omegaBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble *omegaWall_;

    OpenSMOKE::OpenSMOKEVectorDouble *jGas_;
    OpenSMOKE::OpenSMOKEVectorDouble *cGas_;
    OpenSMOKE::OpenSMOKEVectorDouble *diffG_;
    OpenSMOKE::OpenSMOKEVectorDouble *jSolid_;

    OpenSMOKE::OpenSMOKEVectorDouble *teta_;

    OpenSMOKE::OpenSMOKEVectorDouble  Tbulk_;
    OpenSMOKE::OpenSMOKEVectorDouble  Th_;
    OpenSMOKE::OpenSMOKEVectorDouble  Tp_;
    OpenSMOKE::OpenSMOKEVectorDouble  jbulk_;
    OpenSMOKE::OpenSMOKEVectorDouble  cbulk_;
    OpenSMOKE::OpenSMOKEVectorDouble  condHeff_;
    OpenSMOKE::OpenSMOKEVectorDouble  condPeff_;
    OpenSMOKE::OpenSMOKEVectorDouble  jh_;
    OpenSMOKE::OpenSMOKEVectorDouble  jp_;
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
    OpenSMOKE::OpenSMOKEVectorDouble v_;

    OpenSMOKE::OpenSMOKEVectorDouble xBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble xWall_;
    OpenSMOKE::OpenSMOKEVectorDouble cBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble cWall_;
    OpenSMOKE::OpenSMOKEVectorDouble RsurfacePhases_;
    
    OpenSMOKE::OpenSMOKEVectorDouble dyOS_;
    OpenSMOKE::OpenSMOKEVectorDouble  yOS_;

    void HeatTransferCoefficient(const double z, const double cp,  const double eta, const double cond);
    void OverallHeatTransferCoefficient(const double cp, const double eta, const double cond);
    void MassTransferCoefficient(const double z, const double rho, const double eta, const OpenSMOKE::OpenSMOKEVectorDouble dG);
    OpenSMOKE::OpenSMOKEVectorDouble FirstOrderDerivate (const OpenSMOKE::OpenSMOKEVectorDouble value);
    OpenSMOKE::OpenSMOKEVectorDouble SecondOrderDerivate(const OpenSMOKE::OpenSMOKEVectorDouble value, const OpenSMOKE::OpenSMOKEVectorDouble coeff, const std::string type);
    OpenSMOKE::OpenSMOKEVectorDouble SecondOrderDerivate(const OpenSMOKE::OpenSMOKEVectorDouble value, const std::string type);
    double EffectiveThermalConductivity(const double T, const double condG, const std::string type);
};


BVPSystem::BVPSystem(   OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&          thermodynamicsMap, 
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

void BVPSystem::resize(const unsigned int NP)
{
    NP_      = NP;
    NC_      = thermodynamicsMap_.NumberOfSpecies();
    SURF_NC_ = thermodynamicsSurfaceMap_.number_of_site_species();
    SURF_NP_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
    NB_      =  NC_ + NC_ + SURF_NC_ + 1 + 1 + 1;
    NE_      = (NC_ + NC_ + SURF_NC_ + 1 + 1 + 1)*NP_;
    
    omegaBulk_     = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    jGas_          = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    cGas_          = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    jSolid_        = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    omegaWall_     = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    teta_          = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    diffG_         = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];

    for (unsigned int i=0;i<NP_;i++)
    {
        ChangeDimensions(NC_, &omegaBulk_[i], true);
        ChangeDimensions(NC_, &omegaWall_[i], true);
        ChangeDimensions(NC_, &jGas_[i],      true);
        ChangeDimensions(NC_, &cGas_[i],      true);
        ChangeDimensions(NC_, &jSolid_[i],    true);
        ChangeDimensions(NC_, &diffG_[i],     true);
        ChangeDimensions(SURF_NC_, &teta_[i], true);
    }

    ChangeDimensions(NC_, &RfromGas_, true);
    ChangeDimensions(NC_, &RfromSurface_, true);
    ChangeDimensions(SURF_NC_, &Rsurface_, true);

    ChangeDimensions(NP_, &Tbulk_,        true);
    ChangeDimensions(NP_, &Th_,           true);
    ChangeDimensions(NP_, &Tp_,           true);
    ChangeDimensions(NP_, &jbulk_,        true);
    ChangeDimensions(NP_, &cbulk_,        true);
    ChangeDimensions(NP_, &jh_,           true);
    ChangeDimensions(NP_, &jp_,           true);
    ChangeDimensions(NP_, &v_,            true);
    ChangeDimensions(NP_, &rhoBulk_,      true);
    ChangeDimensions(NP_, &condHeff_,     true);
    ChangeDimensions(NP_, &condPeff_,     true);
    ChangeDimensions(NP_, &cp_,           true);
    ChangeDimensions(NP_, &condBulk_,     true);
    ChangeDimensions(NP_, &etaMix_,       true);

    ChangeDimensions(NC_, &xBulk_, true);
    ChangeDimensions(NC_, &xWall_, true);
    ChangeDimensions(NC_, &cBulk_, true);
    ChangeDimensions(NC_, &cWall_, true);
    ChangeDimensions(NC_, &Kmat_, true);

    ChangeDimensions(NE_, &dyOS_, true);
    ChangeDimensions(NE_, &yOS_, true);
}

void BVPSystem::setHoneycombSolid(const double rhoSolid, const double condSolid, const double cpSolid)
{
    rhoH_   = rhoSolid;
    condH_  = condSolid;
    cpH_    = cpSolid;
}

void BVPSystem::setPackedBedSolid(const double rhoSolid, const double condSolid, const double cpSolid)
{
    rhoP_   = rhoSolid;
    condP_  = condSolid;
    cpP_    = cpSolid;
}

void BVPSystem::setInert(const std::string inert)
{
    inert_        = inert;
}

void BVPSystem::setDiscretizationScheme(const std::string discretizationScheme)
{
    discretizationScheme_ = discretizationScheme;
}

void BVPSystem::setCorrelation(const std::string correlation)
{
    correlation_ = correlation;
}

void BVPSystem::setReactorGeometry( const double alfa, const double G, 
                                    const double Lcat, const double Linert,
                                    const double av,   const double aex)
{
    alfaTemp_        = alfa;
    av_              = av;
    aex_             = aex;
    Lcat_            = Lcat;
    L_               = Lcat_ + Linert;
    Linert_          = Linert/L_;
    G_               = G;
}

void BVPSystem::setPackedBedProperties(const double Dp, const double Dt, const double epsi)
{
    Dp_    = Dp;
    Dt_    = Dt;
    epsiP_ = epsi;
}

void BVPSystem::setHoneyCombProperties(const double epsi, const double Dh)
{
    Dh_    = Dh;
    epsiH_ = epsi;
}

void BVPSystem::setFeedValue(const double p, const double T0,
                             const OpenSMOKE::OpenSMOKEVectorDouble x0bulk)
{
    p_  = p;
    T0_ = T0;
    ChangeDimensions(x0bulk.Size(), &x0bulk_, true);
    for (unsigned int j=1;j<=x0bulk.Size();j++)
    {
        x0bulk_[j] = x0bulk[j];
    }
}

void BVPSystem::setGrid(const OpenSMOKE::OpenSMOKEVectorDouble z)
{
    ChangeDimensions(z.Size(), &z_, true);
    for (unsigned int k=1;k<=z_.Size();k++)
    {
        z_[k] = z[k]/L_;
    }
}

double BVPSystem::EffectiveThermalConductivity(const double T, const double condG, const std::string type)
{
    double newCond = 0.;
    if ( type == "honeyComb" )
    {
        newCond = condH_ + (16./3.)*(5.67*1.e-08)*((1.12*Dh_)/2.)*std::pow(T,3.);
    }
    else if ( type == "packedBed" )
    {
        double B = 1.25*std::pow((1-epsiP_)/epsiP_,(10./9.));
        double M = condG/condP_;
        
        newCond = condG*(2./(1.-M*B))*((1-M)*B*std::log(1./(M*B))/std::pow((1-M*B),2.)
                - 0.5*(B+1)
                - (B-1)/(1-M*B));
    }
    return newCond;
}

void BVPSystem::OverallHeatTransferCoefficient(const double cp, const double eta, const double cond)
{
    double phi;
    double A = condP_/cond;
    if ( epsiP_ <= 0.26 )
    {
        phi = 0.072*std::pow((1. - 1./A),2.)/(std::log(A-0.925*(A-1.))-0.075*(1.-1./A)) - 2./(3.*A);
    }
    else if ( epsiP_ >= 0.476 )
    {
        phi = (1/3)*std::pow((1. - 1./A),2.)/(std::log(A-0.577*(A-1.))-0.423*(1.-1./A)) - 2./(3.*A);
    }
    else
    {
        double phi1 = (1./3.)*std::pow((1. - 1./A),2.)/(std::log(A-0.577*(A-1.))-0.423*(1.-1./A)) - 2./(3.*A);
        double phi2 = 0.072*std::pow((1. - 1./A),2.)/(std::log(A-0.925*(A-1.))-0.075*(1.-1./A)) - 2./(3.*A);
        phi = phi2 + (phi1-phi2)*(epsiP_-0.26)/0.216;
    }
    
    double krs = (1.-epsiP_)*(1.+2.66*std::sqrt(Dp_/Dt_))/(2./(3.*A) + phi);
    double krf = G_*cp*Dp_/(8.6*(1+19.4*std::pow((Dp_/Dt_),2.)));
    double Re  = G_*Dp_/(epsiP_*eta);
    double h   = 0.574*G_*cp*std::pow(Re,-0.407)/epsiP_;
    double Ns  = 1.5*(1.-epsiP_)*h*std::pow(Dt_,2.)/(krs*Dp_);
    double Bi  = 5.73*std::sqrt(Dt_/Dp_)*std::pow(Re,-0.26);
    double kr  = krf + krs*(1.+4./Bi)/(1.+8./Ns);
    double hw  = Bi*2.*kr/Dt_;
    
    U_ = 1./(1./hw + 0.5*Dt_*(Bi+3.)/(3.*kr*(Bi+4.)));
}

void BVPSystem::HeatTransferCoefficient(const double z, const double cp, const double eta, const double cond)
{
    if ( correlation_ == "Yoshida" )
    {
        double Re     = G_*Dp_/(epsiP_*eta*(1. - epsiP_)*6.);
        double ReReal = G_*Dp_/(eta*epsiP_);
        double Pr     = cp*eta/cond;
        double Nu;
        double jH;
        if ( Re < 50. )
        {
            jH = 0.91/(std::pow(Re,0.51));
            Nu = jH*std::pow(Pr,(1./3.))*ReReal;
            h_ = Nu*cond/Dp_;
        }
        else
        {
            jH = 0.61/(std::pow(Re,0.41));
            Nu = jH*std::pow(Pr,(1./3.))*ReReal;
            h_ = Nu*cond/Dp_;
        }
    }
    else if ( correlation_ == "Petrovic" )
    {
        double Re = G_*Dp_/(eta*epsiP_);
        double Pr = cp*eta/cond;
        double jH = 0.357/(epsiP_*std::pow(Re,0.359));
        double Nu = jH*std::pow(Pr,(1./3.))*Re;
               h_ = Nu*cond/Dp_;
    }
    else if ( correlation_ == "Wakao" )
    {
        double Re = G_*Dp_/(eta*epsiP_);
        double Pr = cp*eta/cond;
        double Nu = 2. + 1.1*std::pow(Re,0.6)*std::pow(Pr,(1./3.));
               h_ = Nu*cond/Dp_;
    }
}

void BVPSystem::MassTransferCoefficient(const double z, const double rho, const double eta, const OpenSMOKE::OpenSMOKEVectorDouble dG)
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
            Sc[i]    = eta/(rho*dG[i]);
            Sh[i]    = 2. + 1.1*std::pow(Re,0.6)*std::pow(Sc[i],(1./3.));
            Kmat_[i] = Sh[i]*dG[i]/Dp_;
        }
    }
}

#if ASALI_USE_BZZ == 1
void BVPSystem::GetSystemFunctions(BzzVector &y, double t, BzzVector &dy)
{
    FromBzzToOS(y,yOS_);
    #include "BVPequations.H"
    ObjectBzzPrint();
    FromOSToBzz(dyOS_,dy);
}

void BVPSystem::ObjectBzzPrint(void)
{
    #include "printIntegration.H"
}

void BVPSystem::AlgebraicEquations(BzzVectorInt& algebraic)
{
    ChangeDimensions(NE_, &algebraic, true);
    unsigned int counter = 1;
    for (unsigned int i=0;i<NP_;i++)
    {
        if ( i == 0)
        {
            for (unsigned int j=1;j<=NC_;j++)
                algebraic[counter++] = 0;
            for (unsigned int j=1;j<=NC_;j++)
                algebraic[counter++] = 0;
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
                {
                    algebraic[counter++] = 1;
                }
                else
                {
                    algebraic[counter++] = 0;
                }
            }
            if (energy_ == true)
            {
                algebraic[counter++] = 0;
                algebraic[counter++] = 0;
                algebraic[counter++] = 0;
            }
            else
            {
                algebraic[counter++] = 1;
                algebraic[counter++] = 1;
                algebraic[counter++] = 1;
            }
        }
        else if ( i == (NP_ - 1) )
        {
            for (unsigned int j=1;j<=NC_;j++)
            {
                if ( gasDiffusion_ == true )
                {
                    algebraic[counter++] = 0;
                }
                else
                {
                    if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
                    {
                        algebraic[counter++] = 1;
                    }
                    else
                    {
                        algebraic[counter++] = 0;
                    }
                }
            }
            for (unsigned int j=1;j<=NC_;j++)
                algebraic[counter++] = 0;
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
                {
                    algebraic[counter++] = 1;
                }
                else
                {
                    algebraic[counter++] = 0;
                }
            }
            if (energy_ == true)
            {
                if ( gasDiffusion_ == true)
                {
                    algebraic[counter++] = 0;
                }
                else
                {
                    algebraic[counter++] = 1;
                }
                algebraic[counter++] = 0;
                algebraic[counter++] = 0;
            }
            else
            {
                algebraic[counter++] = 1;
                algebraic[counter++] = 1;
                algebraic[counter++] = 1;
            }
        }
        else
        {
            for (unsigned int j=1;j<=NC_;j++)
            {
                if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
                {
                    algebraic[counter++] = 1;
                }
                else
                {
                    algebraic[counter++] = 0;
                }
            }
            for (unsigned int j=1;j<=NC_;j++)
                algebraic[counter++] = 0;
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
                {
                    algebraic[counter++] = 1;
                }
                else
                {
                    algebraic[counter++] = 0;
                }
            }
            algebraic[counter++] = 1;
            algebraic[counter++] = 1;
            algebraic[counter++] = 1;
        }
    }
}
#endif

unsigned int BVPSystem::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
    ChangeDimensions(NE_, &yOS_,  true);
    ChangeDimensions(NE_, &dy,    true);
    ChangeDimensions(NE_, &dyOS_, true);

    for (unsigned int i=1;i<=NE_;i++)
        yOS_[i] = y[i];

    #include "BVPequations.H"

    for (unsigned int i=1;i<=NE_;i++)
        dy[i] = dyOS_[i];
    
    return 0;
}

unsigned int BVPSystem::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
    t_ = t;
    #include "printIntegration.H"
    return 0;
}

void BVPSystem::AlgebraicEquations(OpenSMOKE::OpenSMOKEVectorBool& algebraic)
{
    ChangeDimensions(NE_, &algebraic, true);
    unsigned int counter = 1;
    for (unsigned int i=0;i<NP_;i++)
    {
        if ( i == 0)
        {
            for (unsigned int j=1;j<=NC_;j++)
                algebraic[counter++] = true;
            for (unsigned int j=1;j<=NC_;j++)
                algebraic[counter++] = true;
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
                {
                    algebraic[counter++] = false;
                }
                else
                {
                    algebraic[counter++] = true;
                }
            }
            if (energy_ == true)
            {
                algebraic[counter++] = true;
                algebraic[counter++] = true;
                algebraic[counter++] = true;
            }
            else
            {
                algebraic[counter++] = false;
                algebraic[counter++] = false;
                algebraic[counter++] = false;
            }
        }
        else if ( i == (NP_ - 1) )
        {
            for (unsigned int j=1;j<=NC_;j++)
            {
                if ( gasDiffusion_ == true )
                {
                    algebraic[counter++] = true;
                }
                else
                {
                    if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
                    {
                        algebraic[counter++] = false;
                    }
                    else
                    {
                        algebraic[counter++] = true;
                    }
                }
            }
            for (unsigned int j=1;j<=NC_;j++)
                algebraic[counter++] = true;
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
                {
                    algebraic[counter++] = false;
                }
                else
                {
                    algebraic[counter++] = true;
                }
            }
            if (energy_ == true)
            {
                if ( gasDiffusion_ == true)
                {
                    algebraic[counter++] = true;
                }
                else
                {
                    algebraic[counter++] = false;
                }
                algebraic[counter++] = true;
                algebraic[counter++] = true;
            }
            else
            {
                algebraic[counter++] = false;
                algebraic[counter++] = false;
                algebraic[counter++] = false;
            }
        }
        else
        {
            for (unsigned int j=1;j<=NC_;j++)
            {
                if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
                {
                    algebraic[counter++] = false;
                }
                else
                {
                    algebraic[counter++] = true;
                }
            }
            for (unsigned int j=1;j<=NC_;j++)
                algebraic[counter++] = true;
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
                {
                    algebraic[counter++] = false;
                }
                else
                {
                    algebraic[counter++] = true;
                }
            }
            algebraic[counter++] = false;
            algebraic[counter++] = false;
            algebraic[counter++] = false;
        }
    }
}

void BVPSystem::start()
{
    std::cout << "\n######################################" << std::endl;
    std::cout << "# DAE system:                 START  #" << std::endl;
    std::cout << "######################################\n" << std::endl;
}

void BVPSystem::end()
{
    delete [] omegaBulk_;
    delete [] omegaWall_;
    delete [] jGas_;
    delete [] cGas_;
    delete [] jSolid_;
    delete [] teta_;
    std::cout << "\n######################################" << std::endl;
    std::cout << "# DAE system:                 END    #" << std::endl;
    std::cout << "######################################\n" << std::endl;
}

OpenSMOKE::OpenSMOKEVectorDouble BVPSystem::FirstOrderDerivate (const OpenSMOKE::OpenSMOKEVectorDouble value)
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

OpenSMOKE::OpenSMOKEVectorDouble BVPSystem::SecondOrderDerivate (const OpenSMOKE::OpenSMOKEVectorDouble value, const OpenSMOKE::OpenSMOKEVectorDouble coeff, const std::string type)
{
    OpenSMOKE::OpenSMOKEVectorDouble dvalue_(NP_);

    if ( type == "gas" )
    {
        for (unsigned int k=1;k<=NP_;k++)
        {
            if ( k == 1 )
            {
                dvalue_[k] = 0.;
            }
            else if ( k == NP_ )
            {
                dvalue_[k] = 0.;
            }
            else
            {
                dvalue_[k] = ((coeff[k+1]*value[k+1]-coeff[k]*value[k])/(z[k+1]-z[k])/L_ - (coeff[k]*value[k]-coeff[k-1]*value[k-1])/(z[k]-z[k-1])/L_)/(0.5*(z[k+1]-z[k-1])/L_);
            }
        }
    }
    else if ( type == "solid" )
    {
        for (unsigned int k=1;k<=NP_;k++)
        {
            if ( k == 1 )
            {
                dvalue_[k] = 0.;
            }
            else if ( k == NP_ )
            {
                dvalue_[k] = 0.;
            }
            else
            {
                dvalue_[k] = ((coeff[k+1]*value[k+1]-coeff[k]*value[k])/(z[k+1]-z[k])/L_ - (coeff[k]*value[k]-coeff[k-1]*value[k-1])/(z[k]-z[k-1])/L_)/(0.5*(z[k+1]-z[k-1])/L_);
            }
        }
    }

    return dvalue_;
}

OpenSMOKE::OpenSMOKEVectorDouble BVPSystem::SecondOrderDerivate (const OpenSMOKE::OpenSMOKEVectorDouble value, const std::string type)
{
    OpenSMOKE::OpenSMOKEVectorDouble dvalue_(NP_);

    if ( type == "gas" )
    {
        for (unsigned int k=1;k<=NP_;k++)
        {
            if ( k == 1 )
            {
                dvalue_[k] = 0.;
            }
            else if ( k == NP_ )
            {
                dvalue_[k] = 0.;
            }
            else
            {
                dvalue_[k] = (value[k+1]*(z[k]-z[k-1])/L_ + value[k-1]*(z[k+1]-z[k])/L_ - value[k]*(z[k+1]-z[k-1])/L_)/(0.5*((z[k+1]-z[k-1])/L_)*((z[k+1]-z[k])/L_)*((z[k]-z[k-1])/L_));
            }
        }
    }
    else if ( type == "solid" )
    {
        for (unsigned int k=1;k<=NP_;k++)
        {
            if ( k == 1 )
            {
                dvalue_[k] = 0.;
            }
            else if ( k == NP_ )
            {
                dvalue_[k] = 0.;
            }
            else
            {
                dvalue_[k] = ((value[k+1]-value[k])/(z[k+1]-z[k])/L_ - (value[k]-value[k-1])/(z[k]-z[k-1])/L_)/(0.5*(z[k+1]-z[k-1])/L_);
            }
        }
    }
    return dvalue_;
}
}
