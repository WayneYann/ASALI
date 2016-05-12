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

    ODESystem();

    unsigned int BlockDimension()       const {return NB_;};

    unsigned int NumberOfEquations()    const {return NE_;};
    
    void setGasProperties(const double cpG,
                          const double kG,
                          const double mu,
                          const double diff);
                          
    void setCatalystProperties(const double cpC,
                               const double kC,
                               const double rhoC);

    void setSupportProperties(const double cpS,
                              const double kS,
                              const double rhoS);
    
    void setGrid(const unsigned int Na, const unsigned int Nr);
    
    void setReactionType(const std::string reactionType);
    
    void setPressure(const double p);
    
    void setFlowRate(const double G);
    
    void setCoolantTemperature(const double Tw);
    
    void setFeedTemperature(const double Tin);
    
    void setMassFraction(const std::vector<double> IN);
    
    void setThermodinamicsData(const std::vector<double> MW, const unsigned int NC);
    
    void setPackedBed(const double Dt,
                      const double Dp,
                      const double L);

    void setHoneyComb(const double Dt,
                      const double CPSI,
                      const double w,
                      const double Sw,
                      const double L,
                      const std::string type);

    void setMicroBed(const double Dt,
                     const double CPSI,
                     const double Dp,
                     const double w,
                     const double L);

    void resize(const std::string type);

    void start(const std::string type);
    void end(const std::string type);

    unsigned int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

    unsigned int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

    #if ASALI_USE_BZZ == 1
    virtual void GetSystemFunctions(BzzVector &y, double t, BzzVector &dy);
    virtual void ObjectBzzPrint(void);
    virtual ~ODESystem(){};
    #endif


private:

    unsigned int NE_;
    unsigned int NB_;
    unsigned int NC_;
    unsigned int Na_;
    unsigned int Nr_;
    unsigned int iteration;


    double kG_;
    double mu_;
    double diff_;
    double cpG_;
    double Pr_;
    double kS_;
    double rhoS_;
    double cpS_;
    double kC_;
    double rhoC_;
    double cpC_;
    double DpP_;
    double DtP_;
    double LP_;
    double epsiP_;
    double aexP_;
    double DtH_;
    double epsiH_;
    double epsiHC_;
    double epsiHS_;
    double DcH_;
    double avH_;
    double LH_;
    double DtM_;
    double LM_;
    double epsiMp_;
    double epsiMh_;
    double DpM_;
    double DcM_;
    double aexM_;
    double p_;
    double Tw_;
    double G_;
    double t_;
    double kr_;
    double ka_;
    double hw_;
    double U_;
    double Tin_;

    std::string reactionType_;
    std::string reactorType_;
    std::string typeH_;
    
    OpenSMOKE::OpenSMOKEVectorDouble  zP_;
    OpenSMOKE::OpenSMOKEVectorDouble  rP_;
    OpenSMOKE::OpenSMOKEVectorDouble  zH_;
    OpenSMOKE::OpenSMOKEVectorDouble  rH_;
    OpenSMOKE::OpenSMOKEVectorDouble  zM_;
    OpenSMOKE::OpenSMOKEVectorDouble  rM_;
    
    OpenSMOKE::OpenSMOKEVectorDouble  MW_;
    
    OpenSMOKE::OpenSMOKEVectorDouble  omegaIN_;

    OpenSMOKE::OpenSMOKEVectorDouble dyOS_;
    OpenSMOKE::OpenSMOKEVectorDouble  yOS_;

    OpenSMOKE::OpenSMOKEVectorDouble **omegaWall_;
    OpenSMOKE::OpenSMOKEVectorDouble **omegaBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble **d1stomegaBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble **domegaBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble **domegaWall_;
    OpenSMOKE::OpenSMOKEVectorDouble **R_;
    
    OpenSMOKE::OpenSMOKEVectorDouble *TBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble *d1stTBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble *d2ndTWallaxial_;
    OpenSMOKE::OpenSMOKEVectorDouble *d2ndTWallradial_;
    OpenSMOKE::OpenSMOKEVectorDouble *d2ndTBulkaxial_;
    OpenSMOKE::OpenSMOKEVectorDouble *d2ndTBulkradial_;
    OpenSMOKE::OpenSMOKEVectorDouble *TWall_;
    OpenSMOKE::OpenSMOKEVectorDouble *dTBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble *dTWall_;
    OpenSMOKE::OpenSMOKEVectorDouble *Q_;

    OpenSMOKE::OpenSMOKEVectorDouble *rho_;
    OpenSMOKE::OpenSMOKEVectorDouble *v_;

    OpenSMOKE::OpenSMOKEVectorDouble *kMat_;
    OpenSMOKE::OpenSMOKEVectorDouble *h_;


    #include "mathematics.h"
    #include "transport.h"
    #include "chemistry.h"

};


ODESystem::ODESystem()
    {
        NE_ = 0;
        NB_ = 0;
        Na_ = 0;
        Nr_ = 0;
        NC_ = 0;

        kG_     = 0.;
        mu_     = 0.;
        diff_   = 0.;
        cpG_    = 0.;
        cpS_    = 0.;
        kS_     = 0.;
        rhoS_   = 0.;
        cpS_    = 0.;
        kC_     = 0.;
        rhoC_   = 0.;
        cpC_    = 0.;
        DpP_    = 0.;
        DtP_    = 0.;
        LP_     = 0.;
        epsiP_  = 0.;
        aexP_   = 0.;
        DtH_    = 0.;
        epsiH_  = 0.;
        epsiHS_ = 0.;
        epsiHC_ = 0.;
        DcH_    = 0.;
        avH_    = 0.;
        LH_     = 0.;
        DtM_    = 0.;
        LM_     = 0.;
        epsiMh_ = 0.;
        epsiMp_ = 0.;
        DpM_    = 0.;
        DcM_    = 0.;
        aexM_   = 0.;
        p_      = 0.;
        Tw_     = 0.;
        G_      = 0.;
        hw_     = 0.;
        U_      = 0.;
        Tin_    = 0.;
    }

void ODESystem::setFlowRate(const double G)
{
    G_ = G;
}

void ODESystem::setFeedTemperature(const double Tin)
{
    Tin_ = Tin;
}

void ODESystem::setMassFraction(const std::vector<double> IN)
{
    ChangeDimensions(IN.size(), &omegaIN_,true);
    for (unsigned int i=0;i<IN.size();i++)
    {
        omegaIN_[i+1] = IN[i];
    }
}

void ODESystem::setCoolantTemperature(const double Tw)
{
    Tw_ = Tw;
}

void ODESystem::setPressure(const double p)
{
    p_ = p;
}

void ODESystem::setThermodinamicsData(const std::vector<double> MW, const unsigned int NC)
{
    NC_ = NC;
    ChangeDimensions(NC_, &MW_, true);
    for (unsigned int i=0;i<NC_;i++)
    {
        MW_[i+1] = MW[i];
    }
}

void ODESystem::setGasProperties(const double cpG,
                                 const double kG,
                                 const double mu,
                                 const double diff)
{
    cpG_  = cpG;
    kG_   = kG;
    mu_   = mu;
    diff_ = diff;
    Pr_   = cpG_*mu_/kG_;
}

void ODESystem::setCatalystProperties(const double cpC,
                                      const double kC,
                                      const double rhoC)
{
    cpC_  = cpC;
    kC_   = kC;
    rhoC_ = rhoC;
}

void ODESystem::setSupportProperties(const double cpS,
                                     const double kS,
                                     const double rhoS)
{
    cpS_  = cpS;
    kS_   = kS;
    rhoS_ = rhoS;
}

void ODESystem::setGrid(const unsigned int Na, const unsigned int Nr)
{
    Na_ = Na;
    Nr_ = Nr;
}

void ODESystem::setReactionType(const std::string reactionType)
{
    reactionType_ = reactionType;
}

void ODESystem::setPackedBed(const double Dt, const double Dp, const double L)
{
    double NP = Dt/Dp;

    DtP_   = Dt;
    DpP_   = Dp;
    LP_    = L;

    epsiP_ = 0.4 + 0.05/NP + 0.412/(NP*NP);
    aexP_  = 4./Dt;

    ChangeDimensions(Na_, &zP_, true);
    zP_[1] = 0;
    for (unsigned int i=2;i<=Na_;i++)
    {
        zP_[i] = zP_[i-1] + LP_/double(Na_ - 1.);
    }

    ChangeDimensions(Nr_, &rP_, true);
    rP_[1] = 0;
    for (unsigned int i=2;i<=Nr_;i++)
    {
        rP_[i] = rP_[i-1] + DtP_*0.5/double(Nr_ - 1.);
    }
}

void ODESystem::setHoneyComb(const double Dt,
                             const double CPSI,
                             const double w,
                             const double Sw,
                             const double L,
                             const std::string type)
{
    DtH_   = Dt;
    LH_    = L;
    typeH_ = type;

    if ( typeH_ == "washcoated" )
    {
        double m = (25.4/std::sqrt(CPSI))*1e-03;
        DcH_     = m - (w*2.54/1000)*0.01 - 2.*Sw;
        epsiH_   = std::pow((DcH_/m),2.);
        epsiHC_  = (std::pow((DcH_ + 2.*Sw),2.) - std::pow(DcH_,2.))/std::pow(m,2.);
        epsiHS_  = 1. - epsiH_ - epsiHC_;
    }
    else if ( typeH_ == "extruded" )
    {
        double m = (25.4/std::sqrt(CPSI))*1e-03;
        DcH_     = m - (w*2.54/1000)*0.01;
        epsiH_   = std::pow((DcH_/m),2.);
        epsiHC_  = 1. - epsiH_;
        epsiHS_  = 1. - epsiH_;
    }

    avH_ = 4.*epsiH_/DcH_;

    ChangeDimensions(Na_, &zH_, true);
    zH_[1] = 0;
    for (unsigned int i=2;i<=Na_;i++)
    {
        zH_[i] = zH_[i-1] + LH_/double(Na_ - 1.);
    }

    ChangeDimensions(Nr_, &rH_, true);
    rH_[1] = 0;
    for (unsigned int i=2;i<=Nr_;i++)
    {
        rH_[i] = rH_[i-1] + DtH_*0.5/double(Nr_ - 1.);
    }
}

void ODESystem::setMicroBed(const double Dt,
                            const double CPSI,
                            const double Dp,
                            const double w,
                            const double L)
{
    DtM_   = Dt;
    DpM_   = Dp;
    LM_    = L;

    double m = (25.4/std::sqrt(CPSI))*1e-03;
    DcM_     = m - (w*2.54/1000)*0.01;
    epsiMh_  = std::pow((DcM_/m),2.);

    double NM = DcM_/DpM_;
    epsiMp_ = 0.4 + 0.05/NM + 0.412/(NM*NM);
    aexM_  = 4./DcM_;

    ChangeDimensions(Na_, &zM_, true);
    zM_[1] = 0;
    for (unsigned int i=2;i<=Na_;i++)
    {
        zM_[i] = zM_[i-1] + LM_/double(Na_ - 1.);
    }

    ChangeDimensions(Nr_, &rM_, true);
    rM_[1] = 0;
    for (unsigned int i=2;i<=Nr_;i++)
    {
        rM_[i] = rM_[i-1] + DtM_*0.5/double(Nr_ - 1.);
    }
}

void ODESystem::resize(const std::string type)
{
    reactorType_ = type;
    iteration    = 0.;
    if ( reactorType_ == "honeyComb" )
    {
        NB_ = NC_ + NC_ + 1 + 1;
    }
    else if ( reactorType_ == "packedBed" )
    {
        NB_ = NC_ + 1;
    }
    else if ( reactorType_ == "microBed" )
    {
        NB_ = NC_ + 1 + 1;
    }

    NE_      = NB_*Na_*Nr_;
    
    omegaBulk_        = new OpenSMOKE::OpenSMOKEVectorDouble*[Na_];
    omegaWall_        = new OpenSMOKE::OpenSMOKEVectorDouble*[Na_];
    domegaBulk_       = new OpenSMOKE::OpenSMOKEVectorDouble*[Na_];
    domegaWall_       = new OpenSMOKE::OpenSMOKEVectorDouble*[Na_];
    d1stomegaBulk_    = new OpenSMOKE::OpenSMOKEVectorDouble*[Na_];
    R_                = new OpenSMOKE::OpenSMOKEVectorDouble*[Na_];
    TWall_            = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    TBulk_            = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    Q_                = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    dTWall_           = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    dTBulk_           = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    d1stTBulk_        = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    d2ndTWallaxial_   = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    d2ndTWallradial_  = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    d2ndTBulkaxial_   = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    d2ndTBulkradial_  = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    rho_              = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    v_                = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    kMat_             = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    h_                = new OpenSMOKE::OpenSMOKEVectorDouble[Na_];
    
    for (unsigned int a=0;a<Na_;a++)
    {
        omegaBulk_[a]     = new OpenSMOKE::OpenSMOKEVectorDouble[Nr_];
        omegaWall_[a]     = new OpenSMOKE::OpenSMOKEVectorDouble[Nr_];
        domegaBulk_[a]    = new OpenSMOKE::OpenSMOKEVectorDouble[Nr_];
        domegaWall_[a]    = new OpenSMOKE::OpenSMOKEVectorDouble[Nr_];
        d1stomegaBulk_[a] = new OpenSMOKE::OpenSMOKEVectorDouble[Nr_];
        R_[a]             = new OpenSMOKE::OpenSMOKEVectorDouble[Nr_];

        for (unsigned int r=0;r<Nr_;r++)
        {
            ChangeDimensions(NC_, &omegaBulk_[a][r],     true);
            ChangeDimensions(NC_, &omegaWall_[a][r],     true);
            ChangeDimensions(NC_, &domegaBulk_[a][r],    true);
            ChangeDimensions(NC_, &domegaWall_[a][r],    true);
            ChangeDimensions(NC_, &d1stomegaBulk_[a][r], true);
            ChangeDimensions(NC_, &R_[a][r],             true);
        }
        
        ChangeDimensions(Nr_, &TWall_[a],           true);
        ChangeDimensions(Nr_, &TBulk_[a],           true);
        ChangeDimensions(Nr_, &Q_[a],               true);
        ChangeDimensions(Nr_, &dTWall_[a],          true);
        ChangeDimensions(Nr_, &dTBulk_[a],          true);
        ChangeDimensions(Nr_, &d1stTBulk_[a],       true);
        ChangeDimensions(Nr_, &d2ndTWallaxial_[a],  true);
        ChangeDimensions(Nr_, &d2ndTWallradial_[a], true);
        ChangeDimensions(Nr_, &d2ndTBulkaxial_[a],  true);
        ChangeDimensions(Nr_, &d2ndTBulkradial_[a], true);
        ChangeDimensions(Nr_, &rho_[a],             true);
        ChangeDimensions(Nr_, &v_[a],               true);
        ChangeDimensions(Nr_, &kMat_[a],            true);
        ChangeDimensions(Nr_, &h_[a],               true);
    }

    ChangeDimensions(NE_, &dyOS_, true);
    ChangeDimensions(NE_, &yOS_, true);
}

#if ASALI_USE_BZZ == 1
void ODESystem::GetSystemFunctions(BzzVector &y, double t, BzzVector &dy)
{
    FromBzzToOS(y,yOS_);
    
    if ( reactorType_ == "honeyComb" )
    {
        #include "ODEhoneycomb.H"
    }
    else if ( reactorType_ == "packedBed" )
    {
        #include "ODEpackedbed.H"
    }
    else if ( reactorType_ == "microBed" )
    {
        #include "ODEmicrobed.H"
    }

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

    if ( reactorType_ == "honeyComb" )
    {
        #include "ODEhoneycomb.H"
    }
    else if ( reactorType_ == "packedBed" )
    {
        #include "ODEpackedbed.H"
    }
    else if ( reactorType_ == "microBed" )
    {
        #include "ODEmicrobed.H"
    }

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

void ODESystem::start(const std::string type)
{
    std::cout << "\n######################################" << std::endl;
    std::cout << "# ODE system for " << type << ":     START  " << std::endl;
    std::cout << "######################################\n" << std::endl;
}

void ODESystem::end(const std::string type)
{
    delete [] omegaWall_;
    delete [] omegaBulk_;
    delete [] d1stomegaBulk_;
    delete [] domegaBulk_;
    delete [] domegaWall_;
    delete [] R_;
    delete [] TBulk_;
    delete [] d1stTBulk_;
    delete [] d2ndTWallaxial_;
    delete [] d2ndTWallradial_;
    delete [] d2ndTBulkaxial_;
    delete [] d2ndTBulkradial_;
    delete [] TWall_;
    delete [] dTBulk_;
    delete [] dTWall_;
    delete [] Q_;
    delete [] rho_;
    delete [] v_;
    delete [] kMat_;
    delete [] h_;

    std::cout << "\n######################################" << std::endl;
    std::cout << "# ODE system for " << type << ":     END  " << std::endl;
    std::cout << "######################################\n" << std::endl;
}

}
