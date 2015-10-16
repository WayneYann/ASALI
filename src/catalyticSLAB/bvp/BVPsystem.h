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
              OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&     transportMap,
              dustyGas&                                              dgm);

    #include "vector.h"

    void setReactions(const bool flag)             { reactions_ = flag; }

    void setReactorGeometry( const double alfa, const double L, const double D);
    
    void setDiffusionModel(const std::string diffusionModel);

    void setCatalystProperties(const double epsi, const double tau, const double pore);

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

    #if ASALI_USE_BZZ == 1
    void AlgebraicEquations(BzzVectorInt& algebraic);
    virtual void GetSystemFunctions(BzzVector &y, double t, BzzVector &dy);
    virtual void ObjectBzzPrint(void);
    virtual ~BVPSystem(){};
    #endif

private:

    double MWbulk_;
    double cTotBulk_;
    double rhoBulk_;
    double p_;
    double T0_;
    double D_;
    double L_;
    double alfa_;
    double SD_;
    double t_;
    double A_;
    double tau_;
    double epsi_;
    double pore_;


    unsigned int NC_;
    unsigned int SURF_NP_;
    unsigned int SURF_NC_;
    unsigned int NE_;
    unsigned int NB_;
    unsigned int NP_;

    unsigned int iteration;

    std::string inert_;
    std::string diffusionModel_;

    bool reactions_;

    OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&              thermodynamicsMap_;             //!< thermodynamic map
    OpenSMOKE::KineticsMap_CHEMKIN<double>&                    kineticsMap_;                   //!< kinetic map
    OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&      thermodynamicsSurfaceMap_;      //!< thermodynamic map
    OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>&            kineticsSurfaceMap_;            //!< kinetic map
    OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&         transportMap_;                  //!< transport map
    dustyGas&                                                  dgm_;

    OpenSMOKE::OpenSMOKEVectorDouble *omegaBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble *jBulk_;
 
    OpenSMOKE::OpenSMOKEVectorDouble *jGas_;
    OpenSMOKE::OpenSMOKEVectorDouble *cGas_;

    OpenSMOKE::OpenSMOKEVectorDouble *teta_;

    OpenSMOKE::OpenSMOKEVectorDouble  x0bulk_;

    OpenSMOKE::OpenSMOKEVectorDouble  z_;

    OpenSMOKE::OpenSMOKEVectorDouble RfromGas_;
    OpenSMOKE::OpenSMOKEVectorDouble RfromSurface_;
    OpenSMOKE::OpenSMOKEVectorDouble Rsurface_;

    OpenSMOKE::OpenSMOKEVectorDouble diffG_;

    OpenSMOKE::OpenSMOKEVectorDouble xBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble cBulk_;
    OpenSMOKE::OpenSMOKEVectorDouble RsurfacePhases_;

    OpenSMOKE::OpenSMOKEVectorDouble dyOS_;
    OpenSMOKE::OpenSMOKEVectorDouble  yOS_;

    OpenSMOKE::OpenSMOKEVectorDouble FirstOrderDerivate (const OpenSMOKE::OpenSMOKEVectorDouble value, std::string type);

};


BVPSystem::BVPSystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&          thermodynamicsMap, 
                     OpenSMOKE::KineticsMap_CHEMKIN<double>&                kineticsMap,
                     OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&  thermodynamicsSurfaceMap, 
                     OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped<double>& kineticsSurfaceMap,
                     OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&     transportMap,
                     dustyGas&                                              dgm):
    thermodynamicsMap_(thermodynamicsMap), 
    kineticsMap_(kineticsMap),
    thermodynamicsSurfaceMap_(thermodynamicsSurfaceMap), 
    kineticsSurfaceMap_(kineticsSurfaceMap),
    transportMap_(transportMap),
    dgm_(dgm)
    {
    }

void BVPSystem::resize(const unsigned int NP)
{
    NP_      = NP;
    NC_      = thermodynamicsMap_.NumberOfSpecies();
    SURF_NC_ = thermodynamicsSurfaceMap_.number_of_site_species();
    SURF_NP_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
    NB_      =  NC_ + SURF_NC_;
    NE_      = (NC_ + SURF_NC_)*NP_;
    
    omegaBulk_     = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    jBulk_         = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    jGas_          = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    cGas_          = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
    teta_          = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];

    for (unsigned int i=0;i<NP_;i++)
    {
        ChangeDimensions(NC_, &omegaBulk_[i], true);
        ChangeDimensions(NC_, &jBulk_[i],     true);
        ChangeDimensions(NC_, &jGas_[i],      true);
        ChangeDimensions(NC_, &cGas_[i],      true);
        ChangeDimensions(SURF_NC_, &teta_[i], true);
    }

    ChangeDimensions(NC_, &RfromGas_, true);
    ChangeDimensions(NC_, &RfromSurface_, true);
    ChangeDimensions(SURF_NC_, &Rsurface_, true);


    ChangeDimensions(NC_, &xBulk_, true);
    ChangeDimensions(NC_, &cBulk_, true);
    ChangeDimensions(NC_, &diffG_, true);

    ChangeDimensions(NE_, &dyOS_, true);
    ChangeDimensions(NE_, &yOS_, true);
}

void BVPSystem::setInert(const std::string inert)
{
    inert_        = inert;
}

void BVPSystem::setDiffusionModel(const std::string diffusionModel)
{
    diffusionModel_ = diffusionModel;
}

void BVPSystem::setReactorGeometry( const double alfa, const double L, const double D)
{
    alfa_ = alfa;
    D_    = D;
    L_    = L;
    A_    = 0.25*3.14*std::pow(D_,2.);
}

void BVPSystem::setCatalystProperties(const double epsi, const double tau, const double pore)
{
    epsi_ = epsi;
    tau_  = tau;
    pore_ = pore;
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
        }
        else if ( i == (NP_ - 1) )
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
            for (unsigned int j=1;j<=SURF_NC_;j++)
                algebraic[counter++] = false;
        }
        else if ( i == (NP_ - 1) )
        {
            for (unsigned int j=1;j<=NC_;j++)
                algebraic[counter++] = true;
            for (unsigned int j=1;j<=SURF_NC_;j++)
                algebraic[counter++] = false;
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
            for (unsigned int j=1;j<=SURF_NC_;j++)
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
    delete [] jBulk_;
    delete [] jGas_;
    delete [] cGas_;
    delete [] teta_;
    std::cout << "\n######################################" << std::endl;
    std::cout << "# DAE system:                 END    #" << std::endl;
    std::cout << "######################################\n" << std::endl;
}

OpenSMOKE::OpenSMOKEVectorDouble BVPSystem::FirstOrderDerivate (const OpenSMOKE::OpenSMOKEVectorDouble value, std::string type)
{
    OpenSMOKE::OpenSMOKEVectorDouble dvalue_(NP_);

    if ( type == "CDS" )
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
    else if ( type == "BDS" )
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
    else if ( type == "FDS" )
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
                dvalue_[k] = ((value[k+1] - value[k])/(z_[k+1] - z_[k]))/L_;
            }
        }
    }

    return dvalue_;
}
}
