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

    ODESystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&       thermodynamicsMap, 
              OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&  transportMap);

    #include "vector.h"

    unsigned int NumberOfEquations() const { return NE_; }

    unsigned int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

    unsigned int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

    void SetReactor(const double G, const double Tw, const double Dp, const double Dt, const double epsi);

    void SetReactorType(const std::string type);

    void SetFluid(const std::string name);
    
    void SetCorrelation(const std::string fmName, const std::string NuName);

    void GetProfile(std::vector<double> &z, std::vector<double> &p, std::vector<double> &T);

    void SetTransportLaw(const std::string name);

    #if ASALI_USE_BZZ == 1
    virtual void GetSystemFunctions(BzzVector &y, double t, BzzVector &dy);
    virtual void ObjectBzzPrint(void);
    virtual ~ODESystem(){};
    #endif

private:

    std::string law_;
    std::string name_;
    std::string fmName_;
    std::string NuName_;
    std::string type_;

    double p_;
    double T_;
    double t_;
    double G_;
    double Tw_;
    double Dp_;
    double Dt_;
    double rho_;
    double v_;
    double MW_;
    double n_;
    double cp_;
    double cv_;
    double k_;
    double mu_;
    double av_;
    double fm_;
    double Nu_;
    double Pr_;
    double Re_;
    double epsi_;
    double h_;

    unsigned int NE_;
    unsigned int NC_;
    unsigned int iterator_;

    OpenSMOKE::OpenSMOKEVectorDouble x_;
    OpenSMOKE::OpenSMOKEVectorDouble dyOS_;
    OpenSMOKE::OpenSMOKEVectorDouble yOS_;

    std::vector<double> pprofile_;
    std::vector<double> Tprofile_;
    std::vector<double> z_;

    OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&            thermodynamicsMap_;
    OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&         transportMap_;

    void FrictionFactor(const double Reynolds);
    void Nusselt(const double Reynolds);
    void av();
    void h();
    void error();

    double ReynoldsForFrictionFactor();
    double ReynoldsForHeatTransfer();
    double specificHeatConstantPressure(const double T);
    double specificHeatConstantVolume(const double T);
    double viscosity(const double T);
    double termalConducitivity(const double mu, const double cv, const double MW);
};

ODESystem::ODESystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, 
                     OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>& transportMap):
    thermodynamicsMap_(thermodynamicsMap), 
    transportMap_(transportMap)
    {
        p_    = 0.;
        T_    = 0.;
        t_    = 0.;
        G_    = 0.; 
        Tw_   = 0.; 
        Dp_   = 0.;
        Dt_   = 0.;
        rho_  = 0.;
        v_    = 0.;
        MW_   = 0.;
        n_    = 0.;
        cp_   = 0.;
        cv_   = 0.;
        k_    = 0.;
        mu_   = 0.;
        av_   = 0.;
        fm_   = 0.;
        Nu_   = 0.;
        Pr_   = 0.;
        Re_   = 0.;
        epsi_ = 0.;
        h_    = 0.;

        NE_        = 0.;
        NC_        = 0.;
        iterator_  = 0.; 
        
        
        NC_ = thermodynamicsMap_.NumberOfSpecies();
        NE_ = 2;

        ChangeDimensions(NC_, &x_, true);
        ChangeDimensions(NE_, &yOS_,  true);
        ChangeDimensions(NE_, &dyOS_, true);
    }

void ODESystem::SetReactor(const double G, const double Tw, const double Dp, const double Dt, const double epsi)
{
    G_        = G;
    Tw_        = Tw;
    Dp_        = Dp;
    epsi_    = epsi;
    Dt_        = Dt;
}

void ODESystem::SetFluid(const std::string name)
{
    name_    = name;
}

void ODESystem::SetTransportLaw(const std::string law)
{
    law_ = law;
}

void ODESystem::SetCorrelation(const std::string fmName, const std::string NuName)
{
    fmName_ = fmName;
    NuName_    = NuName;
}

void ODESystem::SetReactorType(const std::string type)
{
    type_ = type;
}

#if ASALI_USE_BZZ == 1
void ODESystem::GetSystemFunctions(BzzVector &y, double t, BzzVector &dy)
{
    FromBzzToOS(y,yOS_);
    #include "ODEequations.H"
    #include "ODEprint.H"
    FromOSToBzz(dyOS_,dy);
}

void ODESystem::ObjectBzzPrint(void)
{
}
#endif

unsigned int ODESystem::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
    ChangeDimensions(NE_, &dy,    true);
    for (unsigned int i=1;i<=NE_;i++)
        yOS_[i] = y[i];
    #include "ODEequations.H"
    #include "ODEprint.H"
    for (unsigned int i=1;i<=NE_;i++)
        dy[i] = dyOS_[i];
    return 0;
}

unsigned int ODESystem::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
    return 0;
}

void ODESystem::GetProfile(std::vector<double> &z, std::vector<double> &p, std::vector<double> &T)
{
    p.resize(pprofile_.size());
    for (unsigned int k=0;k<p.size();k++)
        p[k] = pprofile_[k];

    T.resize(Tprofile_.size());
    for (unsigned int k=0;k<T.size();k++)
         T[k] = Tprofile_[k];

    z.resize(z_.size());
    for (unsigned int k=0;k<z.size();k++)
         z[k] = z_[k];
}

void ODESystem::FrictionFactor(const double Reynolds)
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

void ODESystem::Nusselt(const double Reynolds)
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
                Nu_ = 2. + 1.1*std::pow(Pr_,(1./3.))*std::pow(Reynolds,0.6);
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
                Nu_ = jM*Reynolds*std::pow(Pr_,(1./3.));
            }
        }
        else if ( NuName_ == "Gamson" )
        {
            double newRe = Reynolds/((1-epsi_));
            if ( newRe <= 100)
            {
                double jM = 17.*std::pow((1.-epsi_),0.2)/newRe;
                Nu_ = jM*Reynolds*std::pow(Pr_,(1./3.));
            }
            else
            {
                double jM = 1.46*std::pow((1.-epsi_),0.2)/std::pow(newRe,0.41);
                Nu_ = jM*Reynolds*std::pow(Pr_,(1./3.));
            }
        }
        else if ( NuName_ == "Yoshida" )
        {
            double newRe = Reynolds/(6.*(1.-epsi_));
            if ( newRe <= 50)
            {
                double jM = 0.91/(std::pow(newRe,0.51));
                Nu_ = jM*Reynolds*std::pow(Pr_,(1./3.));
            }
            else if ( newRe > 50 && newRe <= 1000)
            {
                double jM = 0.61/(std::pow(newRe,0.41));
                Nu_ = jM*Reynolds*std::pow(Pr_,(1./3.));
            }
            else
            {
                error();
                std::cout << NuName_ << " correlation out of range!\n" << std::endl;
                exit (EXIT_FAILURE);
            }
        }
    }
    else if ( type_ == "Monolith" )
    {
        double zStar = std::max((t_/Dt_)/(Reynolds*Pr_),1e-06);
        Nu_ = 3.657 + 6.874*std::pow(1e03*zStar,-0.488)*std::exp(-57.2*zStar);
    }
}

void ODESystem::av()
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

double ODESystem::ReynoldsForHeatTransfer()
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

double ODESystem::ReynoldsForFrictionFactor()
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

void ODESystem::h()
{
    if ( type_ == "Monolith" )
    {
        h_ = Nu_*k_/Dt_;
    }
    else if ( type_ == "PackedBed" )
    {
        h_ = Nu_*k_/Dp_;
    }
}

double ODESystem::specificHeatConstantPressure(const double T)
{
    double Tlow     = 200.;
    double Thigh    = 5000.;
    double Tcommon  = 1000.;
    std::vector<double> highCpCoeffs(5);
    highCpCoeffs[0] =  4.45362;
    highCpCoeffs[1] =  0.00314017;
    highCpCoeffs[2] = -1.27841e-06;
    highCpCoeffs[3] =  2.394e-10;
    highCpCoeffs[4] = -1.66903e-14;
    std::vector<double> lowCpCoeffs(5);
    lowCpCoeffs[0] =  2.27572;
    lowCpCoeffs[1] =  0.00992207;
    lowCpCoeffs[2] = -1.04091e-05;
    lowCpCoeffs[3] =  6.86669e-09;
    lowCpCoeffs[4] = -2.11728e-12;
    
    std::vector<double> a;
    a.resize(lowCpCoeffs.size());

    if ( T < Tcommon )
    {
        for (unsigned int i=0;i<a.size();i++)
            a[i] = lowCpCoeffs[i];
    }
    else
    {
        for (unsigned int i=0;i<a.size();i++)
            a[i] = highCpCoeffs[i];
    }

    double cp;
    cp = 8314.*((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0]);
    return cp;
}


double ODESystem::specificHeatConstantVolume(const double T)
{
    double Tlow     = 200.;
    double Thigh    = 5000.;
    double Tcommon  = 1000.;
    std::vector<double> highCpCoeffs(5);
    highCpCoeffs[0] =  4.45362;
    highCpCoeffs[1] =  0.00314017;
    highCpCoeffs[2] = -1.27841e-06;
    highCpCoeffs[3] =  2.394e-10;
    highCpCoeffs[4] = -1.66903e-14;
    std::vector<double> lowCpCoeffs(5);
    lowCpCoeffs[0] =  2.27572;
    lowCpCoeffs[1] =  0.00992207;
    lowCpCoeffs[2] = -1.04091e-05;
    lowCpCoeffs[3] =  6.86669e-09;
    lowCpCoeffs[4] = -2.11728e-12;
    
    std::vector<double> a;
    a.resize(lowCpCoeffs.size());

    if ( T < Tcommon )
    {
        for (unsigned int i=0;i<a.size();i++)
            a[i] = lowCpCoeffs[i];
    }
    else
    {
        for (unsigned int i=0;i<a.size();i++)
            a[i] = highCpCoeffs[i];
    }
    
    double cv;
    cv = 8314.*((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0]);
    cv = cv + 8314;
    return cv;
}

double ODESystem::viscosity(const double T)
{
    double As = 1.67212e-06;
    double Ts = 170.672;
    
    double mu;
    mu = As*std::sqrt(T)/(1 + Ts/T);
    return mu;
}

double ODESystem::termalConducitivity(const double mu, const double cv, const double MW)
{
    double k;
    k = mu*cv*(1.32 + 1.77*8314/(MW*cv));
    return k;
}

void  ODESystem::error()
{
    std::cout << "\nASALI::HEATtransfer::ERROR\n" << std::endl;
}

}
