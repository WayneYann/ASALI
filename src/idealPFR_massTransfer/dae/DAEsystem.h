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
class DAESystem
#if ASALI_USE_BZZ == 1
 : public BzzDaeSystemObject
#endif
{
public:

    DAESystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&          thermodynamicsMap, 
              OpenSMOKE::KineticsMap_CHEMKIN<double>&                kineticsMap,
              OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&  thermodynamicsSurfaceMap, 
              OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped<double>& kineticsSurfaceMap,
              OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&     transportMap);

    #include "vector.h"

    void SetReactor(const double G, const double Dp, const double Dt, const double epsi);

    void SetReactorType(const std::string type);
    
    void SetTemperature(const double T);
    
    void SetInertSpecie(const std::string inert);

    void SetCorrelation(const std::string fmName, const std::string NuName);

    void GetProfile(std::vector<double> &z, std::vector<double> &p, std::vector<double> &comp);

    unsigned int BlockDimension()        const {return NB_;};

    unsigned int NumberOfEquations()     const {return NE_;};

    unsigned int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

    unsigned int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

    void AlgebraicEquations(OpenSMOKE::OpenSMOKEVectorBool& algebraic);

    #if ASALI_USE_BZZ == 1
    void AlgebraicEquations(BzzVectorInt& algebraic);
    virtual void GetSystemFunctions(BzzVector &y, double t, BzzVector &dy);
    virtual void ObjectBzzPrint(void);
    virtual ~DAESystem(){};
    #endif

private:

    std::string fmName_;
    std::string NuName_;
    std::string type_;

    double p_;
    double T_;
    double t_;
    double G_;
    double Dp_;
    double Dt_;
    double rho_;
    double v_;
    double MW_;
    double mu_;
    double av_;
    double fm_;
    double Re_;
    double epsi_;
    double cTot_;

    unsigned int NE_;
    unsigned int NB_;
    unsigned int NC_;
    unsigned int SURF_NC_;
    unsigned int SURF_NP_;
    unsigned int iterator_;

    OpenSMOKE::OpenSMOKEVectorDouble x_;
    OpenSMOKE::OpenSMOKEVectorDouble xW_;
    OpenSMOKE::OpenSMOKEVectorDouble omega_;
    OpenSMOKE::OpenSMOKEVectorDouble omegaW_;
    OpenSMOKE::OpenSMOKEVectorDouble diffG_;
    OpenSMOKE::OpenSMOKEVectorDouble Sc_;
    OpenSMOKE::OpenSMOKEVectorDouble c_;
    OpenSMOKE::OpenSMOKEVectorDouble R_;
    OpenSMOKE::OpenSMOKEVectorDouble Rsurface_;
    OpenSMOKE::OpenSMOKEVectorDouble RsurfacePhases_;
    OpenSMOKE::OpenSMOKEVectorDouble Sh_;
    OpenSMOKE::OpenSMOKEVectorDouble kMat_;

    OpenSMOKE::OpenSMOKEVectorDouble dyOS_;
    OpenSMOKE::OpenSMOKEVectorDouble yOS_;

    OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&              thermodynamicsMap_;             //!< thermodynamic map
    OpenSMOKE::KineticsMap_CHEMKIN<double>&                    kineticsMap_;                   //!< kinetic map
    OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&      thermodynamicsSurfaceMap_;      //!< thermodynamic map
    OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>&            kineticsSurfaceMap_;            //!< kinetic map
    OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&         transportMap_;                  //!< transport map

    std::vector<double> pprofile_;
    std::vector<double> Compprofile_;
    std::vector<double> z_;

    std::string inert_;

    void FrictionFactor(const double Reynolds);
    void Sherwood(const double Reynolds);
    void av();
    void Sc();
    void kMat();
    void error();

    double ReynoldsForFrictionFactor();
    double ReynoldsForMassTransfer();

};

DAESystem::DAESystem(   OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&          thermodynamicsMap, 
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
        p_    = 0.;
        T_    = 0.;
        t_    = 0.;
        G_    = 0.; 
        Dp_   = 0.;
        Dt_   = 0.;
        rho_  = 0.;
        v_    = 0.;
        MW_   = 0.;
        mu_   = 0.;
        av_   = 0.;
        fm_   = 0.;
        Re_   = 0.;
        epsi_ = 0.;
        cTot_ = 0.;

        NE_        = 0.;
        NB_        = 0.;
        NC_        = 0.;
        SURF_NC_   = 0.;
        SURF_NP_   = 0.;
        iterator_  = 0.; 

        NC_      = thermodynamicsMap_.NumberOfSpecies();
        SURF_NC_ = thermodynamicsSurfaceMap_.number_of_site_species();
        SURF_NP_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
        NE_      = NC_ + NC_ + 1;

        ChangeDimensions(NC_, &x_,     true);
        ChangeDimensions(NC_, &xW_,    true);
        ChangeDimensions(NC_, &c_,     true);
        ChangeDimensions(NC_, &omega_, true);
        ChangeDimensions(NC_, &omegaW_,true);
        ChangeDimensions(NC_, &diffG_, true);
        ChangeDimensions(NC_, &Sc_,    true);
        ChangeDimensions(NC_, &Sh_,    true);
        ChangeDimensions(NC_, &kMat_,  true);
        ChangeDimensions(NC_, &R_,     true);

        ChangeDimensions(SURF_NC_, &Rsurface_, true);

        ChangeDimensions(NE_, &yOS_,  true);
        ChangeDimensions(NE_, &dyOS_, true);
        
        Compprofile_.resize(NC_);
    }

void DAESystem::SetReactor(const double G, const double Dp, const double Dt, const double epsi)
{
    G_         = G;
    Dp_        = Dp;
    epsi_      = epsi;
    Dt_        = Dt;
}

void DAESystem::SetCorrelation(const std::string fmName, const std::string NuName)
{
    fmName_    = fmName;
    NuName_    = NuName;
}

void DAESystem::SetReactorType(const std::string type)
{
    type_ = type;
}


void DAESystem::SetInertSpecie(const std::string inert)
{
    inert_ = inert;
}

void DAESystem::SetTemperature(const double T)
{
    T_ = T;
}

#if ASALI_USE_BZZ == 1
void DAESystem::GetSystemFunctions(BzzVector &y, double t, BzzVector &dy)
{
    FromBzzToOS(y,yOS_);
    #include "DAEequations.H"
    #include "DAEprint.H"
    FromOSToBzz(dyOS_,dy);
}

void DAESystem::ObjectBzzPrint(void)
{
}

void DAESystem::AlgebraicEquations(BzzVectorInt& algebraic)
{
    ChangeDimensions(NE_, &algebraic, true);
    for (unsigned int j=1;j<=NC_;j++)
    {
        if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
        {
            algebraic[j]     = 1;
        }
        else
        {
            algebraic[j]     = 0;
        }
    }
    for (unsigned int j=1;j<=NC_;j++)
    {
        algebraic[j+NC_] = 0;
    }
    algebraic[NE_]       = 1;
}
#endif

unsigned int DAESystem::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
    ChangeDimensions(NE_, &dy,    true);
    for (unsigned int i=1;i<=NE_;i++)
        yOS_[i] = y[i];
    #include "DAEequations.H"
    #include "DAEprint.H"
    for (unsigned int i=1;i<=NE_;i++)
        dy[i] = dyOS_[i];
    return 0;
}

unsigned int DAESystem::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
    return 0;
}

void DAESystem::AlgebraicEquations(OpenSMOKE::OpenSMOKEVectorBool& algebraic)
{
    ChangeDimensions(NE_, &algebraic, true);
    for (unsigned int j=1;j<=NC_;j++)
    {
        if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
        {
            algebraic[j]     = false;
        }
        else
        {
            algebraic[j]     = true;
        }
    }
    for (unsigned int j=1;j<=NC_;j++)
    {
        algebraic[j+NC_] = true;
    }
    algebraic[NE_]       = false;
}

void DAESystem::GetProfile(std::vector<double> &z, std::vector<double> &p, std::vector<double> &comp)
{
    p.resize(pprofile_.size());
    for (unsigned int k=0;k<p.size();k++)
        p[k] = pprofile_[k];

    comp.resize(Compprofile_.size());
    for (unsigned int k=0;k<Compprofile_.size();k++)
         comp[k] = Compprofile_[k];

    z.resize(z_.size());
    for (unsigned int k=0;k<z.size();k++)
         z[k] = z_[k];
}

void DAESystem::FrictionFactor(const double Reynolds)
{
    if ( type_ == "PackedBed" )
    {
        if ( fmName_ == "Foumeny" )
        {
            double newReynolds = Reynolds/(1.- epsi_);
            if ( newReynolds > 85000. || newReynolds < 5.)
            {
                error();
                std::cout << fmName_ << " correlation out of range!\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                fm_ = ((1.-epsi_)/std::pow(epsi_,3.))*(((Dt_/Dp_)/(0.335*Dt_/Dp_ + 2.28)) + 130/newReynolds);
            }
        }
        else if ( fmName_ == "Ergun" )
        {
            if ( Reynolds > 1000 || Reynolds < 0.1)
            {
                error();
                std::cout << fmName_ << " correlation out of range!\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                fm_ = (1.-epsi_)*(1.75 + 150*(1.-epsi_)/Reynolds)/std::pow(epsi_,3.);
            }
        }
        else if ( fmName_ == "Lee" )
        {
            if ( Reynolds > 1e05)
            {
                error();
                std::cout << fmName_ << " correlation out of range!\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                double n=0.352 + 0.1*epsi_ + 0.275*epsi_*epsi_;
                fm_ = 6.25*(1.-epsi_)*(1.-epsi_)*(29.32/Reynolds + 1.56/std::pow(Reynolds,n) + 0.1)/std::pow(epsi_,3.);
            }
        }
        else if ( fmName_ == "Hicks" )
        {
            if ( Reynolds/(1.-epsi_) > 60000 || Reynolds/(1.-epsi_) < 300)
            {
                error();
                std::cout << fmName_ << " correlation out of range!\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                fm_ = 6.8*std::pow((1.-epsi_),1.2)*std::pow(Reynolds,-0.2)/std::pow(epsi_,3.);
            }
        }
        else if ( fmName_ == "Eisfeld" )
        {
            if ( Reynolds > 17635 || Reynolds < 0.01)
            {
                error();
                std::cout << fmName_ << " correlation out of range!\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                double K1 = 154.;
                double k1 = 1.15;
                double k2 = 0.87;
                double Aw = 1. + 2./(3.*Dt_*(1.-epsi_)/Dp_);
                double Bw = std::pow((k1*std::pow((Dp_/Dt_),2.) + k2),2.);
                fm_ = K1*std::pow(Aw,2.)*std::pow((1.-epsi_),2.)/(Reynolds*std::pow(epsi_,3.)) + (Aw/Bw)*(1.-epsi_)/std::pow(epsi_,3.);
            }
        }
        else if ( fmName_ == "Achenbach" )
        {
            if ( Reynolds/(1.-epsi_) > 5.*1.e04)
            {
                error();
                std::cout << fmName_ << " correlation out of range!\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                fm_ = (1.-epsi_)*(160./(Reynolds/(1.-epsi_)) + 3./std::pow((Reynolds/(1.-epsi_)),0.1))/std::pow(epsi_,3.);
            }
        }
    }
    else if ( type_ == "Monolith" )
    {
        fm_ = 16./Reynolds;
    }
}

void DAESystem::Sherwood(const double Reynolds)
{
    if ( type_ == "PackedBed" )
    {
        if ( NuName_ == "Wakao" )
        {
            if ( Reynolds > 10000 || Reynolds < 3)
            {
                error();
                std::cout << NuName_ << " correlation out of range!\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                for (unsigned int i=1;i<=NC_;i++)
                {
                    Sh_[i] = 2. + 1.1*std::pow(Sc_[i],(1./3.))*std::pow(Reynolds,0.6);
                }
            }
        }
        else if ( NuName_ == "Gupta" )
        {
            if ( Reynolds > 2140 || Reynolds < 1)
            {
                error();
                std::cout << NuName_ << " correlation out of range!\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                double jM = (0.01 + 0.863/(std::pow(Reynolds,0.58) - 0.483))/epsi_;
                for (unsigned int i=1;i<=NC_;i++)
                {
                    Sh_[i] = jM*Reynolds*std::pow(Sc_[i],(1./3.));
                }
            }
        }
        else if ( NuName_ == "Gamson" )
        {
            double newRe = Reynolds/((1-epsi_));
            if ( newRe <= 100)
            {
                double jM = 17.*std::pow((1.-epsi_),0.2)/newRe;
                for (unsigned int i=1;i<=NC_;i++)
                {
                    Sh_[i] = jM*Reynolds*std::pow(Sc_[i],(1./3.));
                }
            }
            else
            {
                double jM = 1.46*std::pow((1.-epsi_),0.2)/std::pow(newRe,0.41);
                for (unsigned int i=1;i<=NC_;i++)
                {
                    Sh_[i] = jM*Reynolds*std::pow(Sc_[i],(1./3.));
                }
            }
        }
        else if ( NuName_ == "Yoshida" )
        {
            double newRe = Reynolds/(6.*(1.-epsi_));
            if ( newRe <= 50)
            {
                double jM = 0.91/(std::pow(newRe,0.51));
                for (unsigned int i=1;i<=NC_;i++)
                {
                    Sh_[i] = jM*Reynolds*std::pow(Sc_[i],(1./3.));
                }
            }
            else if ( newRe > 50 && newRe <= 1000)
            {
                double jM = 0.61/(std::pow(newRe,0.41));
                for (unsigned int i=1;i<=NC_;i++)
                {
                    Sh_[i] = jM*Reynolds*std::pow(Sc_[i],(1./3.));
                }
            }
            else
            {
                error();
                std::cout << NuName_ << " correlation out of range!\n" << std::endl;
                exit (EXIT_FAILURE);
            }
        }
        else if ( NuName_ == "Petrovic" )
        {
            if ( Reynolds > 3 && Reynolds < 230)
            {
                double jM = 0.357/(std::pow(Reynolds,0.359));
                for (unsigned int i=1;i<=NC_;i++)
                {
                    Sh_[i] = jM*Reynolds*std::pow(Sc_[i],(1./3.))/epsi_;
                }
            }
            else
            {
                error();
                std::cout << NuName_ << " correlation out of range!\n" << std::endl;
                exit (EXIT_FAILURE);
            }
        }
        else if ( NuName_ == "Rebughini" )
        {
            double jM = 0.393/(std::pow(Reynolds,0.384));
            for (unsigned int i=1;i<=NC_;i++)
            {
                Sh_[i] = jM*Reynolds*std::pow(Sc_[i],(1./3.))/epsi_;
            }
        }
    }
    else if ( type_ == "Monolith" )
    {
        for (unsigned int i=1;i<=NC_;i++)
        {
            double zNew  = std::max( 0., 1e-06);
            double zStar = zNew/(Dt_*Reynolds*Sc_[i]);
                   zStar = fabs(zStar);
                   Sh_[i] = 3.659 + 6.874*pow((1000.*zStar),-0.488)*exp(-57.2*zStar);
        }
    }
}

void DAESystem::av()
{
    if ( type_ == "Monolith" )
    {
        av_ = 4.*epsi_/Dt_;
    }
    else if ( type_ == "PackedBed" )
    {
        av_ = 6.*(1-epsi_)/Dp_;
    }
}

double DAESystem::ReynoldsForMassTransfer()
{
    double Reynolds = 0.;
    if ( type_ == "Monolith" )
    {
        Reynolds = G_*Dt_/mu_;
    }
    else if ( type_ == "PackedBed" )
    {
        Reynolds = G_*Dp_/(epsi_*mu_);
    }

    return Reynolds;
}

double DAESystem::ReynoldsForFrictionFactor()
{
    double Reynolds = 0.;
    if ( type_ == "Monolith" )
    {
        Reynolds = G_*Dt_/mu_;
    }
    else if ( type_ == "PackedBed" )
    {
        Reynolds = G_*Dp_/mu_;
    }
    
    return Reynolds;
}

void DAESystem::Sc()
{
    for(unsigned int i=1;i<=NC_;i++)
    {
        Sc_[i] = mu_/(rho_*diffG_[i]);
    }
}

void DAESystem::kMat()
{
    if ( type_ == "Monolith" )
    {
        for (unsigned int i=1;i<=NC_;i++)
        {
            kMat_[i] = Sh_[i]*diffG_[i]/Dt_;
        }
    }
    else if ( type_ == "PackedBed" )
    {
        for (unsigned int i=1;i<=NC_;i++)
        {
            kMat_[i] = Sh_[i]*diffG_[i]/Dp_;
        }
    }
}

void  DAESystem::error()
{
    std::cout << "\nASALI::HEATtransfer::ERROR\n" << std::endl;
}
}
