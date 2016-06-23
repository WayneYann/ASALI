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

    equationSystem();

    unsigned int BlockDimension()       const {return NB_;};

    unsigned int NumberOfEquations()    const {return NE_;};
    
    void setGasProperties(const double cpG,
                          const double kG,
                          const double mu,
                          const double diff);
                          
    void setCatalystProperties(const double cpC,
                               const double kC,
                               const double rhoC,
                               const double epsiC,
                               const double tauC);

    void setReactionType(const std::string reactionType);
    
    void setPressure(const double p);

    void setFlowRate(const double G);
    
    void setCoolantTemperature(const double Tw);

    void setThermodinamicsData(const std::vector<double> MW, const unsigned int NC);
    
    void setGeometry(const double Dt,
                     const double Dp);

    void setHeatExchange(const bool flag);

    void resize(const std::string type);

    void GetDaeProfile(std::vector<double> &z, 
                       std::vector<double> &Tb,
                       std::vector<double> &Tw,
                       std::vector<double> &omegab,
                       std::vector<double> &omegaw);


    void GetOdeProfile(std::vector<double> &z,
                       std::vector<double> &Tb,
                       std::vector<double> &omegab);

    unsigned int DaeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

    unsigned int DaePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

    unsigned int OdeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

    unsigned int OdePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

    void AlgebraicEquations(std::vector<OpenSMOKE::EquationType>& M);


private:

    unsigned int NE_;
    unsigned int NB_;
    unsigned int NC_;
    unsigned int dae_iterator_;
    unsigned int ode_iterator_;

    bool extHeat_;

    double kG_;
    double mu_;
    double diff_;
    double cpG_;
    double Pr_;
    double kC_;
    double rhoC_;
    double cpC_;
    double epsiC_;
    double tauC_;
    double Dp_;
    double Dt_;
    double epsi_;
    double aex_;
    double p_;
    double Tw_;
    double G_;
    double t_;
    double kr_;
    double ka_;
    double hw_;
    double U_;
    double Lgeo_;
    double TBulk_;
    double TWall_;
    double dTBulk_;
    double dTWall_;
    double Q_;
    double rho_;
    double v_;
    double kMat_;
    double h_;
    double av_;

    std::string reactionType_;
    std::string reactorType_;

    std::vector<double> Tb_;
    std::vector<double> omegab_;
    std::vector<double> TW_;
    std::vector<double> omegaw_;
    std::vector<double> z_;

    OpenSMOKE::OpenSMOKEVectorDouble  MW_;

    OpenSMOKE::OpenSMOKEVectorDouble dyOS_;
    OpenSMOKE::OpenSMOKEVectorDouble  yOS_;

    OpenSMOKE::OpenSMOKEVectorDouble omegaWall_;
    OpenSMOKE::OpenSMOKEVectorDouble omegaBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble domegaWall_;
    OpenSMOKE::OpenSMOKEVectorDouble domegaBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble R_;
    
    #include "mathematics.h"
    #include "transport.h"
    #include "chemistry.h"

};


equationSystem::equationSystem()
    {
        NE_ = 0;
        NB_ = 0;
        NC_ = 0;

        kG_     = 0.;
        mu_     = 0.;
        diff_   = 0.;
        cpG_    = 0.;
        kC_     = 0.;
        rhoC_   = 0.;
        tauC_   = 0.;
        epsiC_  = 0.;
        cpC_    = 0.;
        Dp_     = 0.;
        Dt_     = 0.;
        epsi_   = 0.;
        aex_    = 0.;
        p_      = 0.;
        Tw_     = 0.;
        G_      = 0.;
        hw_     = 0.;
        U_      = 0.;
        Lgeo_   = 0.;
        av_     = 0.;
    }

void equationSystem::setFlowRate(const double G)
{
    G_ = G;
}

void equationSystem::setCoolantTemperature(const double Tw)
{
    Tw_ = Tw;
}

void equationSystem::setPressure(const double p)
{
    p_ = p;
}

void equationSystem::setHeatExchange(const bool flag)
{
    extHeat_ = flag;
}

void equationSystem::setThermodinamicsData(const std::vector<double> MW, const unsigned int NC)
{
    NC_ = NC;
    ChangeDimensions(NC_, &MW_, true);
    for (unsigned int i=0;i<NC_;i++)
    {
        MW_[i+1] = MW[i];
    }
}

void equationSystem::setGasProperties(const double cpG,
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

void equationSystem::setCatalystProperties(const double cpC,
                                      const double kC,
                                      const double rhoC,
                                      const double epsiC,
                                      const double tauC)
{
    cpC_   = cpC;
    kC_    = kC;
    rhoC_  = rhoC;
    epsiC_ = epsiC;
    tauC_  = tauC;
}

void equationSystem::setReactionType(const std::string reactionType)
{
    reactionType_ = reactionType;
}

void equationSystem::setGeometry(const double Dt, const double Dp)
{
    double NP = Dt/Dp;

    Dt_   = Dt;
    Dp_   = Dp;
    Lgeo_ = Dp_/6.;

    epsi_ = 0.4 + 0.05/NP + 0.412/(NP*NP);
    aex_  = 4./Dt;
    av_   = 6.*(1. - epsi_)/Dp_;
}

void equationSystem::resize(const std::string type)
{
    reactorType_  = type;
    dae_iterator_ = 0;
    ode_iterator_ = 0;

    Tb_.clear();
    omegab_.clear();
    TW_.clear();
    omegaw_.clear();
    z_.clear();

    if ( reactorType_ == "heterogeneous" )
    {
        NE_   = (NC_ + NC_ + 1 + 1);
        NB_   = NC_ + 1;
    }
    else if ( reactorType_ == "pseudoHomogeneous" )
    {
        NE_   = NC_ + 1;
        NB_   = NC_ + 1;
    }

    ChangeDimensions(NC_, &omegaBulk_,  true);
    ChangeDimensions(NC_, &omegaWall_,  true);
    ChangeDimensions(NC_, &domegaBulk_, true);
    ChangeDimensions(NC_, &domegaWall_, true);
    ChangeDimensions(NC_, &R_,          true);

    ChangeDimensions(NE_, &dyOS_, true);
    ChangeDimensions(NE_, &yOS_,  true);
}

void equationSystem::GetOdeProfile(std::vector<double> &z,
                                   std::vector<double> &Tb,
                                   std::vector<double> &omegab)
{
    Tb.resize(Tb_.size());
    for (unsigned int k=0;k<Tb.size();k++)
        Tb[k] = Tb_[k];

    omegab.resize(omegab_.size());
    for (unsigned int k=0;k<omegab_.size();k++)
         omegab[k] = omegab_[k];

    z.resize(z_.size());
    for (unsigned int k=0;k<z.size();k++)
         z[k] = z_[k];
}

void equationSystem::GetDaeProfile(std::vector<double> &z, 
                                   std::vector<double> &Tb,
                                   std::vector<double> &Tw,
                                   std::vector<double> &omegab,
                                   std::vector<double> &omegaw)
{
    Tb.resize(Tb_.size());
    for (unsigned int k=0;k<Tb.size();k++)
        Tb[k] = Tb_[k];

    Tw.resize(TW_.size());
    for (unsigned int k=0;k<Tw.size();k++)
        Tw[k] = TW_[k];

    omegab.resize(omegab_.size());
    for (unsigned int k=0;k<omegab_.size();k++)
         omegab[k] = omegab_[k];

    omegaw.resize(omegaw_.size());
    for (unsigned int k=0;k<omegaw_.size();k++)
         omegaw[k] = omegaw_[k];

    z.resize(z_.size());
    for (unsigned int k=0;k<z.size();k++)
         z[k] = z_[k];
}

unsigned int equationSystem::DaeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
    ChangeDimensions(NE_, &yOS_,  true);
    ChangeDimensions(NE_, &dy,    true);
    ChangeDimensions(NE_, &dyOS_, true);

    for (unsigned int i=1;i<=NE_;i++)
        yOS_[i] = y[i];

    if ( reactorType_ == "heterogeneous" )
    {
        #include "DAEheterogeneous.H"
        #include "DAEprint.H"
    }

    for (unsigned int i=1;i<=NE_;i++)
        dy[i] = dyOS_[i];
    
    return 0;
}

unsigned int equationSystem::DaePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
    return 0;
}

unsigned int equationSystem::OdeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
    ChangeDimensions(NE_, &yOS_,  true);
    ChangeDimensions(NE_, &dy,    true);
    ChangeDimensions(NE_, &dyOS_, true);

    for (unsigned int i=1;i<=NE_;i++)
        yOS_[i] = y[i];

    if ( reactorType_ == "heterogeneous" )
    {
        #include "ODEheterogeneous.H"
    }
    else if ( reactorType_ == "pseudoHomogeneous" )
    {
        #include "ODEpseudoHomogeneous.H"
        #include "ODEprint.H"
    }

    for (unsigned int i=1;i<=NE_;i++)
        dy[i] = dyOS_[i];
    
    return 0;
}

unsigned int equationSystem::OdePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
    return 0;
}

void equationSystem::AlgebraicEquations(std::vector<OpenSMOKE::EquationType>& M)
{
    M.resize(NE_);
    unsigned int counter = 0;
    for (unsigned int s=1;s<=NC_;s++)
    {
        M[counter++] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
    }

    for (unsigned int s=1;s<=NC_;s++)
    {
        M[counter++]  = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
    }

    M[counter++]  = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
    M[counter++]  = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
}
}
