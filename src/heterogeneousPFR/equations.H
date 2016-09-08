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

    void setSolid(const double rhoSolid, const double condSolid, const double cpSolid);

    void setHomogeneusReactions(const bool flag)   { homogeneusReactions_ = flag; }

    void setHeterogeneusReactions(const bool flag) { heterogeneusReactions_ = flag; }

    void setEnergyEquation(const bool flag)        { energy_ = flag; }
    
    void setDiffusion(const bool flag)             { gasDiffusion_ = flag; }

    void setExternalHeatExchange(const bool flag, const double Tex);

    void setReactorType(const std::string reactorType);
    
    void setCorrelation(const std::string correlation);
    
    void setCatalyst(const std::string name);

    void setDiscretizationScheme(const std::string discretizationScheme);

    void setReactorGeometry(const double alfa, const double epsi, 
                            const double Lcat, const double Linert,
                            const double av,   const double aex,
                            const double G);

    void setPackedBedProperties(const double Dh, const double Dt);

    void setHoneyCombProperties(const double AsymptoticSh, const double Dh);

    void setTubularProperties(const double AsymptoticSh, const double Dh);

    void setFeedValue(const double p, const double T0,
                      const OpenSMOKE::OpenSMOKEVectorDouble x0bulk);

    void setGrid(const OpenSMOKE::OpenSMOKEVectorDouble z);

    void setInert(const std::string inert);

    void resize(const unsigned int NP);
    
    void cleanMemory();
    
    OpenSMOKE::OpenSMOKEVectorDouble getVelocity();

    unsigned int BlockDimension()                     const {return NB_;};

    unsigned int NumberOfEquations()                  const {return NE_;};
    
    unsigned int BandDimension()                      const {return NB_ + NB_+ 1;};

    unsigned int BvpEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

    unsigned int BvpPrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

    unsigned int OdeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

    unsigned int OdePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

    void AlgebraicEquations(std::vector<OpenSMOKE::EquationType>& M);
    
    std::vector<double> atomicBalanceError(const OpenSMOKE::OpenSMOKEVectorDouble yF);
    
    double massBalanceError(const OpenSMOKE::OpenSMOKEVectorDouble yF);
    
    double enthalpyBalanceError(const OpenSMOKE::OpenSMOKEVectorDouble yF);
    
    OpenSMOKE::OpenSMOKEVectorDouble getVelocity(const OpenSMOKE::OpenSMOKEVectorDouble yF);


private:

    double MWbulk_;
    double MWwall_;
    double cTotBulk_;
    double cTotWall_;
    double rhoWall_;
    double p_;
    double T0_;
    double epsi_;
    double Dh_;
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
    double t_;
    double Tex_;
    double U_;

    double rhoSolid_;
    double cpSolid_;
    double condSolid_;

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
    std::string reactorType_;
    std::string correlation_;
    std::string catalyst_;

    bool homogeneusReactions_ ;
    bool heterogeneusReactions_;
    bool energy_;
    bool gasDiffusion_;
    bool hex_;

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
    OpenSMOKE::OpenSMOKEVectorDouble  Twall_;
    OpenSMOKE::OpenSMOKEVectorDouble  jbulk_;
    OpenSMOKE::OpenSMOKEVectorDouble  cbulk_;
    OpenSMOKE::OpenSMOKEVectorDouble  newCondSolid_;
    OpenSMOKE::OpenSMOKEVectorDouble  jwall_;
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
    void MassTransferCoefficient(const double z, const double rho, const double eta, const OpenSMOKE::OpenSMOKEVectorDouble dG);
    
    OpenSMOKE::OpenSMOKEVectorDouble FirstOrderDerivate (const OpenSMOKE::OpenSMOKEVectorDouble value);
    OpenSMOKE::OpenSMOKEVectorDouble SecondOrderDerivate(const OpenSMOKE::OpenSMOKEVectorDouble value, const OpenSMOKE::OpenSMOKEVectorDouble coeff, const std::string type);

    double EffectiveThermalConductivity(const double T, const double cond);
    double OverallHeatTransferCoefficient(const double cp, const double eta, const double cond);
    double TaylorArisDisperionCoefficient(const double Pe);

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
    }

void equationSystem::resize(const unsigned int NP)
{
    NP_      = NP;
    NC_      = thermodynamicsMap_.NumberOfSpecies();
    SURF_NC_ = thermodynamicsSurfaceMap_.number_of_site_species();
    SURF_NP_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
    NB_      =  NC_ + NC_ + SURF_NC_ + 1 + 1;
    NE_      = (NC_ + NC_ + SURF_NC_ + 1 + 1)*NP_;
    
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
    ChangeDimensions(NP_, &Twall_,        true);
    ChangeDimensions(NP_, &jbulk_,        true);
    ChangeDimensions(NP_, &cbulk_,        true);
    ChangeDimensions(NP_, &jwall_,        true);
    ChangeDimensions(NP_, &v_,            true);
    ChangeDimensions(NP_, &rhoBulk_,      true);
    ChangeDimensions(NP_, &newCondSolid_, true);
    ChangeDimensions(NP_, &cp_,           true);
    ChangeDimensions(NP_, &condBulk_,     true);
    ChangeDimensions(NP_, &etaMix_,       true);

    ChangeDimensions(NC_, &xBulk_,  true);
    ChangeDimensions(NC_, &xWall_,  true);
    ChangeDimensions(NC_, &cBulk_,  true);
    ChangeDimensions(NC_, &cWall_,  true);
    ChangeDimensions(NC_, &Kmat_,   true);

    ChangeDimensions(NE_, &dyOS_, true);
    ChangeDimensions(NE_, &yOS_, true);
}

void equationSystem::setSolid(const double rhoSolid, const double condSolid, const double cpSolid)
{
    rhoSolid_   = rhoSolid;
    condSolid_  = condSolid;
    cpSolid_    = cpSolid;
}

void equationSystem::setInert(const std::string inert)
{
    inert_        = inert;
}

void equationSystem::setDiscretizationScheme(const std::string discretizationScheme)
{
    discretizationScheme_ = discretizationScheme;
}

void equationSystem::setCatalyst(const std::string name)
{
    catalyst_ = name;
}

void equationSystem::setReactorType(const std::string reactorType)
{
    reactorType_ = reactorType;
}

void equationSystem::setCorrelation(const std::string correlation)
{
    correlation_ = correlation;
}

void equationSystem::setExternalHeatExchange(const bool flag, const double Tex)
{
    hex_ = flag;
    Tex_ = Tex;
}

void equationSystem::setReactorGeometry(const double alfa, const double epsi, 
                                        const double Lcat, const double Linert,
                                        const double av,   const double aex,
                                        const double G)
{
    alfaTemp_        = alfa;
    av_              = av;
    epsi_            = epsi;
    Lcat_            = Lcat;
    L_               = Lcat_ + Linert;
    Linert_          = Linert;
    aex_             = aex;
    G_               = G;
}

void equationSystem::setPackedBedProperties(const double Dh, const double Dt)
{
    Dh_ = Dh;
    Dt_ = Dt;
}

void equationSystem::setHoneyCombProperties(const double AsymptoticSh, const double Dh)
{
    Dh_           = Dh;
    AsymptoticSh_ = AsymptoticSh;
}

void equationSystem::setTubularProperties(const double AsymptoticSh, const double Dh)
{
    Dh_           = Dh;
    AsymptoticSh_ = AsymptoticSh;
}

void equationSystem::setFeedValue(const double p, const double T0,
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

void equationSystem::setGrid(const OpenSMOKE::OpenSMOKEVectorDouble z)
{
    ChangeDimensions(z.Size(), &z_, true);
    for (unsigned int k=1;k<=z_.Size();k++)
    {
        z_[k] = z[k];
    }
}

double equationSystem::EffectiveThermalConductivity(const double T, const double cond)
{
    double newCond = 0.;
    if ( reactorType_ == "honeyComb" )
    {
        newCond = condSolid_ + (16./3.)*(5.67*1.e-08)*((1.12*Dh_)/2.)*std::pow(T,3.);
    }
    else if ( reactorType_ == "packedBed" )
    {
        double ks  = condSolid_/cond;
        double B   = 1.25*std::pow((1.-epsi_)/epsi_,(10./9.));
        double M   = (ks - B)/ks;
        double krs = std::sqrt(1.-epsi_)*(2/M)*(B*(ks-1.)*std::log(ks/B)/(M*M*ks) - 0.5*(B+1)-(B-1)/M);
        
        newCond = cond*krs;
    }
    else if ( reactorType_ == "tubular" )
    {
        newCond = condSolid_;
    }

    return newCond;
}

double equationSystem::OverallHeatTransferCoefficient(const double cp, const double eta, const double cond)
{
    double U = 0.;
    if ( reactorType_ == "honeyComb" )
    {
        U = 500.;
    }
    else if ( reactorType_ == "packedBed" )         // Dixon et al. 1988
    {
        double Re   = G_*Dh_/(eta*epsi_);
        double Pr   = cp*eta/cond;
        double Pe   = Re*Pr;
        double N    = Dt_/Dh_;
        double Bis  = 2.41 +0.15*std::pow((N-1),2.);
        double Nufs = 0.255*std::pow(Pr,(1./3.))*std::pow(Re,0.67)/epsi_;
        double Perf = 1./(0.74*epsi_/Pe + 1./12.);
        double Nwf  =(1.-1./N)*std::pow(Pr,(1./3.))*std::pow(Re,0.6);
        double hw   = Nwf*cond/Dh_;
        double Bif  = Nwf*(2*N)*Perf/Pe;
        double Bi   = 1/(1/Bis + 1/Bif);
        double B    = 1.25*std::pow((1.-epsi_)/epsi_,(10./9.));
        double ks   = condSolid_/cond;
        double M    = (ks - B)/ks;
        double krs0 = std::sqrt(1.-epsi_)*(2/M)*(B*(ks-1.)*std::log(ks/B)/(M*M*ks) - 0.5*(B+1)-(B-1)/M);
        double krs  = krs0*cond;
        double h    = Nufs*cond/Dh_;
        double Ns   = av_*h*Dt_*Dt_*0.25/krs;
        double Per  = (Bi+4)/(Bi*((1/Perf)*(Bif/(Bif+4)) + (krs0/Pe)/(8/Ns + 1)));
        double kr   = Re*Pr/Per;
               U    = 1/(1/hw + (Dt_/(6*kr))*(Bi+3)/(Bi+4));
    }
    else if ( reactorType_ == "tubular" )
    {
        U = AsymptoticSh_*cond/Dh_;
    }

    return U;
}

double equationSystem::TaylorArisDisperionCoefficient(const double Pe)
{
    double diffTA = 0.;
    if ( reactorType_ == "honeyComb" || reactorType_ == "tubular" )
    {
        diffTA = std::pow(Pe,2.)/192.;
    }
    else if ( reactorType_ == "packedBed" )
    {
        diffTA = 0.;
    }
    
    return diffTA;
}


void equationSystem::HeatTransferCoefficient(const double z, const double cp, const double eta, const double cond)
{
    if ( reactorType_ == "honeyComb" || reactorType_ == "tubular" )
    {
        double Re = G_*Dh_/(eta*epsi_);
        double Pr = cp*eta/cond;
        double Nu;
        double zNew;
        double zStar;
        if ( z > Linert_)
        {
            zNew  = std::max( z - Linert_, 1e-06);
            zStar = zNew/(Dh_*Re*Pr);
            zStar = fabs(zStar);
            Nu    = AsymptoticSh_ + 8.827*pow((1000.*zStar),-0.545)*exp(-48.2*zStar);
        }
        else
        {
            Nu = AsymptoticSh_;
        }
        h_ = Nu*cond/Dh_;
    }
    else if ( reactorType_ == "packedBed" )
    {
        if ( correlation_ == "Yoshida" )
        {
            double Re     = G_*Dh_/(epsi_*eta*(1. - epsi_)*6.);
            double ReReal = G_*Dh_/(eta*epsi_);
            double Pr     = cp*eta/cond;
            double Nu;
            double jH;
            if ( Re < 50. )
            {
                jH = 0.91/(std::pow(Re,0.51));
                Nu = jH*std::pow(Pr,(1./3.))*ReReal;
                h_ = Nu*cond/Dh_;
            }
            else
            {
                jH = 0.61/(std::pow(Re,0.41));
                Nu = jH*std::pow(Pr,(1./3.))*ReReal;
                h_ = Nu*cond/Dh_;
            }
        }
        else if ( correlation_ == "Petrovic" )
        {
            double Re = G_*Dh_/(eta*epsi_);
            double Pr = cp*eta/cond;
            double jH = 0.357/(epsi_*std::pow(Re,0.359));
            double Nu = jH*std::pow(Pr,(1./3.))*Re;
                   h_ = Nu*cond/Dh_;
        }
        else if ( correlation_ == "Wakao" )
        {
            double Re = G_*Dh_/(eta*epsi_);
            double Pr = cp*eta/cond;
            double Nu = 2. + 1.1*std::pow(Re,0.6)*std::pow(Pr,(1./3.));
                   h_ = Nu*cond/Dh_;
        }
    }
}

void equationSystem::MassTransferCoefficient(const double z, const double rho, const double eta, const OpenSMOKE::OpenSMOKEVectorDouble dG)
{
    if ( reactorType_ == "honeyComb" || reactorType_ == "tubular" )
    {
        double Re = G_*Dh_/(eta*epsi_);
        OpenSMOKE::OpenSMOKEVectorDouble Sc(NC_);
        OpenSMOKE::OpenSMOKEVectorDouble Sh(NC_);
        double zNew;
        double zStar;
        for (unsigned int i=1;i<=NC_;i++)
        {
            Sc[i] = eta/(rho*dG[i]);
            if ( z > Linert_)
            {
                zNew  = std::max( z - Linert_, 1e-06);
                zStar = zNew/(Dh_*Re*Sc[i]);
                zStar = fabs(zStar);
                Sh[i] = AsymptoticSh_ + 6.874*pow((1000.*zStar),-0.488)*exp(-57.2*zStar);
            }
            else
            {
                Sh[i] = AsymptoticSh_;
            }
            Kmat_[i] = Sh[i]*dG[i]/Dh_;
        }
    }
    else if ( reactorType_ == "packedBed" )
    {
        if ( correlation_ == "Yoshida" )
        {
            double Re     = G_*Dh_/(epsi_*eta*(1. - epsi_)*6.);
            double ReReal = G_*Dh_/(eta*epsi_);
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
                    Kmat_[i] = Sh[i]*dG[i]/Dh_;
                }
                else
                {
                    jM[i]    = 0.61/(std::pow(Re,0.41));
                    Sh[i]    = jM[i]*std::pow(Sc[i],(1./3.))*ReReal;
                    Kmat_[i] = Sh[i]*dG[i]/Dh_;
                }
            }
        }
        else if ( correlation_ == "Petrovic" )
        {
            double Re = G_*Dh_/(eta*epsi_);
            OpenSMOKE::OpenSMOKEVectorDouble Sc(NC_);
            OpenSMOKE::OpenSMOKEVectorDouble Sh(NC_);
            OpenSMOKE::OpenSMOKEVectorDouble jM(NC_);
            for (unsigned int i=1;i<=NC_;i++)
            {
                Sc[i]    = eta/(rho*dG[i]);
                jM[i]    = 0.357/(epsi_*std::pow(Re,0.359));
                Sh[i]    = jM[i]*std::pow(Sc[i],(1./3.))*Re;
                Kmat_[i] = Sh[i]*dG[i]/Dh_;
            }
        }
        else if ( correlation_ == "Wakao" )
        {
            double Re = G_*Dh_/(eta*epsi_);
            OpenSMOKE::OpenSMOKEVectorDouble Sc(NC_);
            OpenSMOKE::OpenSMOKEVectorDouble Sh(NC_);
            for (unsigned int i=1;i<=NC_;i++)
            {
                Sc[i]    = eta/(rho*dG[i]);
                Sh[i]    = 2. + 1.1*std::pow(Re,0.6)*std::pow(Sc[i],(1./3.));
                Kmat_[i] = Sh[i]*dG[i]/Dh_;
            }
        }
    }
}


OpenSMOKE::OpenSMOKEVectorDouble equationSystem::FirstOrderDerivate (const OpenSMOKE::OpenSMOKEVectorDouble value)
{
    OpenSMOKE::OpenSMOKEVectorDouble dvalue_(NP_);

    if ( discretizationScheme_ == "CDS" )
    {
        for (unsigned int k=1;k<=NP_;k++)
        {
            if ( k == 1 )
            {
                dvalue_[k] = ((value[k+1] - value[k])/(z_[k+1] - z_[k]));
            }
            else if ( k == NP_ )
            {
                dvalue_[k] = ((value[k] - value[k-1])/(z_[k] - z_[k-1]));
            }
            else
            {
                dvalue_[k] = ((value[k+1] - value[k-1])/(z_[k+1] - z_[k-1]));
            }
        }
    }
    else if ( discretizationScheme_ == "BDS" )
    {
        for (unsigned int k=1;k<=NP_;k++)
        {
            if ( k == 1 )
            {
                dvalue_[k] = ((value[k+1] - value[k])/(z_[k+1] - z_[k]));
            }
            else if ( k == NP_ )
            {
                dvalue_[k] = ((value[k] - value[k-1])/(z_[k] - z_[k-1]));
            }
            else
            {
                dvalue_[k] = ((value[k] - value[k-1])/(z_[k] - z_[k-1]));
            }
        }
    }
    return dvalue_;
}

OpenSMOKE::OpenSMOKEVectorDouble equationSystem::SecondOrderDerivate (const OpenSMOKE::OpenSMOKEVectorDouble value, const OpenSMOKE::OpenSMOKEVectorDouble coeff, const std::string type)
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
                dvalue_[k] = ((coeff[k+1]+coeff[k])*0.5*(value[k+1]-value[k])/(z[k+1]-z[k]) - (coeff[k]+coeff[k-1])*0.5*(value[k]-value[k-1])/(z[k]-z[k-1]))/(0.5*(z[k+1]-z[k-1])); //((coeff[k+1]+coeff[k]+coeff[k-1])/3.)*(value[k+1]*(z[k]-z[k-1]) + value[k-1]*(z[k+1]-z[k]) - value[k]*(z[k+1]-z[k-1]))/(0.5*(z[k+1]-z[k-1])*(z[k+1]-z[k])*(z[k]-z[k-1]));
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
                dvalue_[k] = ((coeff[k+1]*value[k+1]-coeff[k]*value[k])/(z[k+1]-z[k]) - (coeff[k]*value[k]-coeff[k-1]*value[k-1])/(z[k]-z[k-1]))/(0.5*(z[k+1]-z[k-1]));
            }
        }
    }

    return dvalue_;
}

unsigned int equationSystem::BvpEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
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

unsigned int equationSystem::BvpPrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
    t_ = t;
    #include "printIntegration.H"
    return 0;
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

void equationSystem::AlgebraicEquations(std::vector<OpenSMOKE::EquationType>& M)
{
    M.resize(NE_);
    unsigned int counter = 0;
    for (unsigned int i=0;i<NP_;i++)
    {
        if ( i == 0)
        {
            for (unsigned int j=1;j<=NC_;j++)
            {
                M[counter++] = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
            }
            for (unsigned int j=1;j<=NC_;j++)
            {
                M[counter++] = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
            }
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != catalyst_ )
                {
                    M[counter++] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
                }
                else
                {
                    M[counter++] = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
                }
            }
            if (energy_ == true)
            {
                M[counter++] = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
                M[counter++] = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
            }
            else
            {
                M[counter++] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
                M[counter++] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
            }
        }
        else if ( i == (NP_ - 1) )
        {
            for (unsigned int j=1;j<=NC_;j++)
            {
                if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
                {
                    M[counter++] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
                }
                else
                {
                    M[counter++] = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
                }
            }
            for (unsigned int j=1;j<=NC_;j++)
            {
                M[counter++] = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
            }
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != catalyst_ )
                {
                    M[counter++] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
                }
                else
                {
                    M[counter++] = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
                }
            }
            if (energy_ == true)
            {
                M[counter++] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
                M[counter++] = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
            }
            else
            {
                M[counter++] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
                M[counter++] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
            }
        }
        else
        {
            for (unsigned int j=1;j<=NC_;j++)
            {
                if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
                {
                    M[counter++] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
                }
                else
                {
                    M[counter++] = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
                }
            }
            for (unsigned int j=1;j<=NC_;j++)
            {
                M[counter++] = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
            }
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != catalyst_ )
                {
                    M[counter++] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
                }
                else
                {
                    M[counter++] = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
                }
            }
            M[counter++] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
            M[counter++] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
        }
    }
}


void equationSystem::cleanMemory()
{
    delete [] omegaBulk_;
    delete [] omegaWall_;
    delete [] jGas_;
    delete [] cGas_;
    delete [] jSolid_;
    delete [] teta_;
}

std::vector<double> equationSystem::atomicBalanceError(const OpenSMOKE::OpenSMOKEVectorDouble yF)
{
    unsigned int counter = 1;
    for (unsigned int i=0;i<NP_;i++)
    {
        for (unsigned int j=1;j<=NC_;j++)
            omegaBulk_[i][j] = yF[counter++];
        for (unsigned int j=1;j<=NC_;j++)
            omegaWall_[i][j] = yF[counter++];
        for (unsigned int j=1;j<=SURF_NC_;j++)
            teta_[i][j] = yF[counter++];
        Tbulk_[i+1] = yF[counter++];
        Twall_[i+1] = yF[counter++];
    }

    Eigen::MatrixXd SpecieAtomicComposition;

    std::vector<std::string> ElementsNames;

    std::vector<double> error;

    Eigen::VectorXd Element(NC_);

    ElementsNames           = thermodynamicsMap_.elements();
    SpecieAtomicComposition = thermodynamicsMap_.atomic_composition();

    error.resize(ElementsNames.size());

    for (unsigned int i=0;i<NP_;i++)
    {
        thermodynamicsMap_.SetTemperature(Tbulk_[i+1]);
        thermodynamicsMap_.SetPressure(p_);
        transportMap_.SetPressure(p_);
        transportMap_.SetTemperature(Tbulk_[i+1]);

        thermodynamicsMap_.MoleFractions_From_MassFractions(xBulk_,MWbulk_,omegaBulk_[i]);

        rhoBulk_[i+1] = p_*MWbulk_/(PhysicalConstants::R_J_kmol*Tbulk_[i+1]);

        transportMap_.MassDiffusionCoefficients(diffG_[i], xBulk_);
    }

    for (unsigned int k=0;k<ElementsNames.size();k++)
    {
        double in  = 0.;
        double out = 0.;
        if ( gasDiffusion_ == true )
        {
            Element = SpecieAtomicComposition.col(k);
            for (unsigned int j=0;j<NC;j++)
            {
                in  += (omegaBulk_[0][j+1]*G_- rhoBulk_[1]*diffG_[0][j+1]*(omegaBulk_[1][j+1]-omegaBulk_[0][j+1])/(z[2] - z[1]))*Element[j]/thermodynamicsMap_.MW()[j+1];
                out +=  omegaBulk_[NP_-1][j+1]*G_*Element[j]/thermodynamicsMap_.MW()[j+1];
            }
        }
        else
        {
            Element = SpecieAtomicComposition.col(k);
            for (unsigned int j=0;j<NC;j++)
            {
                in  +=  omegaBulk_[0][j+1]*G_*Element[j]/thermodynamicsMap_.MW()[j+1];
                out +=  omegaBulk_[NP_-1][j+1]*G_*Element[j]/thermodynamicsMap_.MW()[j+1];
            }
        }

        if (in != 0. )
        {
            error[k] = (in - out)*100./in;
        }
        else
        {
            error[k] = 0.;
        }
    }

    return error;
}

double equationSystem::massBalanceError(const OpenSMOKE::OpenSMOKEVectorDouble yF)
{
    double error;
    unsigned int counter = 1;
    for (unsigned int i=0;i<NP_;i++)
    {
        for (unsigned int j=1;j<=NC_;j++)
            omegaBulk_[i][j] = yF[counter++];
        for (unsigned int j=1;j<=NC_;j++)
            omegaWall_[i][j] = yF[counter++];
        for (unsigned int j=1;j<=SURF_NC_;j++)
            teta_[i][j] = yF[counter++];
        Tbulk_[i+1] = yF[counter++];
        Twall_[i+1] = yF[counter++];
    }

    double in  = 0.;
    double out = 0.;
    for (unsigned int j=0;j<NC;j++)
    {
        in  +=  omegaBulk_[0][j+1]*G_;
        out +=  omegaBulk_[NP_-1][j+1]*G_;
    }

    error = (in - out)*100./in;

    return error;
}

double equationSystem::enthalpyBalanceError(const OpenSMOKE::OpenSMOKEVectorDouble yF)
{
    double error;
    unsigned int counter = 1;
    for (unsigned int i=0;i<NP_;i++)
    {
        for (unsigned int j=1;j<=NC_;j++)
            omegaBulk_[i][j] = yF[counter++];
        for (unsigned int j=1;j<=NC_;j++)
            omegaWall_[i][j] = yF[counter++];
        for (unsigned int j=1;j<=SURF_NC_;j++)
            teta_[i][j] = yF[counter++];
        Tbulk_[i+1] = yF[counter++];
        Twall_[i+1] = yF[counter++];
    }

    double in  = 0.;
    {
        OpenSMOKE::OpenSMOKEVectorDouble h(NC);
        thermodynamicsMap_.SetPressure(p_);
        thermodynamicsMap_.SetTemperature(Tbulk_[1]);
        thermodynamicsMap_.hMolar_Species(h);
        for (unsigned int j=0;j<NC;j++)
        {
            in  +=  omegaBulk_[0][j+1]*G_*h[j+1]/thermodynamicsMapXML->MW()[j+1];
        }
    }

    double out  = 0.;
    {
        OpenSMOKE::OpenSMOKEVectorDouble h(NC);
        thermodynamicsMap_.SetPressure(p_);
        thermodynamicsMap_.SetTemperature(Tbulk_[NP_]);
        thermodynamicsMap_.hMolar_Species(h);
        for (unsigned int j=0;j<NC;j++)
        {
            out  +=  omegaBulk_[NP_-1][j+1]*G_*h[j+1]/thermodynamicsMapXML->MW()[j+1];
        }
    }

    error = (in - out)*100./in;

    return error;
}

OpenSMOKE::OpenSMOKEVectorDouble equationSystem::getVelocity(const OpenSMOKE::OpenSMOKEVectorDouble yF)
{
    unsigned int counter = 1;
    for (unsigned int i=0;i<NP_;i++)
    {
        for (unsigned int j=1;j<=NC_;j++)
            omegaBulk_[i][j] = yF[counter++];
        for (unsigned int j=1;j<=NC_;j++)
            omegaWall_[i][j] = yF[counter++];
        for (unsigned int j=1;j<=SURF_NC_;j++)
            teta_[i][j] = yF[counter++];
        Tbulk_[i+1] = yF[counter++];
        Twall_[i+1] = yF[counter++];
    }

    OpenSMOKE::OpenSMOKEVectorDouble vel(NP_);

    for (unsigned int i=0;i<NP_;i++)
    {
        thermodynamicsMap_.SetTemperature(Tbulk_[i+1]);
        thermodynamicsMap_.SetPressure(p_);
        transportMap_.SetPressure(p_);
        transportMap_.SetTemperature(Tbulk_[i+1]);

        thermodynamicsMap_.MoleFractions_From_MassFractions(xBulk_,MWbulk_,omegaBulk_[i]);

        rhoBulk_[i+1] = p_*MWbulk_/(PhysicalConstants::R_J_kmol*Tbulk_[i+1]);

        vel[i+1]   = G_/(rhoBulk_[i+1]*epsi_);
    }
    
    return vel;
}


}
