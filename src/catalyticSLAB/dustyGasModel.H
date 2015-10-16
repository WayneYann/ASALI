#ifndef DUSTYGAS_H
#define DUSTYGAS_H

class dustyGas
{
    public:

        dustyGas(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&           thermodynamicsMap,
                 OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&      transportMap);

        void resize();

        void setComposition(const OpenSMOKE::OpenSMOKEVectorDouble omega);

        void setCompositionGradient(const Eigen::MatrixXd omegaGrad);
        
        void setPressureGradient(const Eigen::VectorXd pGrad);

        void setCatalystProperties(const double epsi,const double tau,const double poreRadius);

        void setPermeabilityLaw(const std::string law);

        void setTemperature(const double T);

        void setPressure(const double P);

        void solve();

        inline Eigen::MatrixXd getFluxes() const {return fluxes_;}

        virtual ~dustyGas();


    private:

        OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&      thermodynamicsMap_;     //!< thermodynamic map
        OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>& transportMap_;          //!< transport map

        Eigen::MatrixXd A_;
        Eigen::MatrixXd b_;
        Eigen::MatrixXd fluxes_;
        Eigen::MatrixXd Diff_;
        Eigen::MatrixXd xGrad_;

        Eigen::VectorXd pGrad_;
        Eigen::VectorXd knudsenDiffCoeff_;

        unsigned int NC_;

        double T_;
        double P_;
        double epsi_;
        double tau_;
        double correctionFactor_;
        double poreRadius_;
        double K0_;
        double B0_;
        double etaMix_;
        double MW_;

        OpenSMOKE::OpenSMOKEVectorDouble binaryDiffCoeff_;
        OpenSMOKE::OpenSMOKEVectorDouble x_;
        OpenSMOKE::OpenSMOKEVectorDouble omega_;
};

dustyGas::dustyGas(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&        thermodynamicsMap,
                   OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&   transportMap) :
            thermodynamicsMap_(thermodynamicsMap),
            transportMap_(transportMap)
            {
                NC_ = thermodynamicsMap_.NumberOfSpecies();
            }

void dustyGas::resize()
{
    Diff_.resize(NC_,NC_);
    A_.resize(NC_,NC_);
    b_.resize(NC_,3);
    xGrad_.resize(NC_,3);
    fluxes_.resize(NC_,3);

    pGrad_.resize(3);
    knudsenDiffCoeff_.resize(NC_);
    
    ChangeDimensions(NC_, &x_,     true);
    ChangeDimensions(NC_, &omega_, true);

}

void dustyGas::setComposition(const OpenSMOKE::OpenSMOKEVectorDouble omega)
{
    thermodynamicsMap_.MoleFractions_From_MassFractions(x_,MW_,omega);

    for (unsigned int i=1;i<=NC_;i++)
    {
        omega_[i] = omega[i];
    }
}

void dustyGas::setTemperature(const double T)
{
    T_= T;
}

void dustyGas::setPressure(const double P)
{
    P_= P;
}

void dustyGas::setCatalystProperties(const double epsi,const double tau,const double poreRadius)
{
    epsi_             = epsi;
    tau_              = tau;
    correctionFactor_ = epsi_/tau_;
    poreRadius_       = poreRadius;
}

void dustyGas::setPermeabilityLaw(const std::string law)
{
    if (law == "Veldsink")
    {
        B0_ = std::pow(epsi_*(poreRadius_),2.)/(8.*tau_); //Veldsink correlation for B0
        //(ref. J.W. Veldsink, "A catalytically active, non-permselective, membrane reactor for kinetically fast, strongly exothermic, heterogeneous reactions")
    }
    else if (law == "Poiseuille")
    {
        B0_ = std::pow(poreRadius_,2.)/(32.); //Poiseuille correlation for cylindrical pore
        //(ref. R. Krishna, J.A. Wesseling, "The Maxwell-Stefan approach to mass transfer")
    }
    else if (law == "Carman-Kozeny")
    {
        B0_ = std::pow((poreRadius_),2.)/(180.)*std::pow(epsi_/(1.-epsi_),2.); //Carman-Kozeny correlation for aggregated bed of spheres
        //(ref. R. Krishna, J.A. Wesseling, "The Maxwell-Stefan approach to mass transfer")
    }
    else if (law == "none")
    {
        B0_ = 0.;
    }
    else
    {
        std::cout << "DUSTY GAS CLASS: wrong permabilityLaw option: Veldsink || Poiseuille || Carman-Kozeny || none " << std::endl;
        exit(-1);
    }
}


void dustyGas::setCompositionGradient(const Eigen::MatrixXd omegaGrad)
{
    Eigen::VectorXd sum(3);

    for ( unsigned int j=0;j<3;j++)
    {
        sum(j) = 0;
    }

    for ( unsigned int i=0;i<NC_;i++ )
    {
        for ( unsigned int j=0;j<3;j++)
        {
            sum(j) = sum(j) + omegaGrad(i,j)/thermodynamicsMap_.MW()[i+1];
        }
    }

    for ( unsigned int i=0;i<NC_;i++ )
    {
        for ( unsigned int j=0;j<3;j++)
        {
            xGrad_(i,j) = omegaGrad(i,j)*MW_/thermodynamicsMap_.MW()[i+1] + omega_[i+1]*sum(j)/thermodynamicsMap_.MW()[i+1];
        }
    }
}

void dustyGas::setPressureGradient(const Eigen::VectorXd pGrad)
{
    for ( unsigned int i=0;i<3;i++ )
    {
        pGrad_(i) = pGrad(i);
    }
}

void dustyGas::solve()
{
    transportMap_.SetTemperature(T_);
    transportMap_.SetPressure(P_);
    thermodynamicsMap_.SetTemperature(T_);
    thermodynamicsMap_.SetPressure(P_);

    //- Calculate b
    {
        transportMap_.DynamicViscosity(etaMix_,x_);

        {
            double K0_ = (epsi_*poreRadius_)/(2.*tau_); //(ref. J.W. Veldsink, "A catalytically active, non-permselective, membrane reactor for kinetically fast, strongly exothermic, heterogeneous reactions")
            for(unsigned int i=0;i<NC_;i++)
            {
                knudsenDiffCoeff_(i) = (4./3.)*K0_*std::sqrt((8.*PhysicalConstants::R_J_mol*T_)/(PhysicalConstants::pi*thermodynamicsMap_.MW()[i+1]));
            }
        }


        for(unsigned int i=0;i<NC_;i++)
        {
            for(unsigned int j=0;j<3;j++)
            {
                b_(i,j) = -(1./(PhysicalConstants::R_J_mol*T_))*(xGrad_(i,j)+(x_[i+1]/P_)*(1.+(B0_*P_)/(etaMix_*knudsenDiffCoeff_(i)))*pGrad_(j));
            }
        }
    }

    //- Calculate A
    {
        transportMap_.getGammaBinary(binaryDiffCoeff_);
        for(unsigned int i=1;i<=binaryDiffCoeff_.Size();i++)
        {
            binaryDiffCoeff_[i] = (P_/100000.)/binaryDiffCoeff_[i];
        }

        unsigned int counter = 1;

        for(unsigned int i=0;i<NC_;i++)
        {
            Diff_(i,i) = 0.;
            for(unsigned int j=(i+1);j<NC_;j++)
            {
                Diff_(i,j) = binaryDiffCoeff_(counter++);
                Diff_(j,i) = Diff_(i,j);
            }
        }

        Eigen::VectorXd sum(NC_);

        for(unsigned int i=0;i<NC_;i++)
        {
            sum(i) = 0;
        }

        for(unsigned int i=0;i<NC_;i++)
        {
            for(unsigned int j=0;j<NC_;j++)
            {
                if( j != i )
                {
                    sum(i) += x_[j+1]/(Diff_(i,j)*correctionFactor_);
                }
            }
        }

        for(unsigned int i=0;i<NC_;i++)
        {
            for(unsigned int j=0;j<NC_;j++)
            {
                if( j != i )
                {
                    A_(i,j) = -x_[i+1]/(Diff_(i,j)*correctionFactor_);
                }
                else
                {
                    A_(i,j) = 1./(knudsenDiffCoeff_(i)) + sum(i);
                }
            }
        }
    }

    //- Linear System Resolution
    {
        Eigen::PartialPivLU<Eigen::MatrixXd> luDec(A_);
        fluxes_ = luDec.solve(b_);
    }

    //- [mol/m2/s] --> [Kg/m2/s]
    for(unsigned int i=0;i<NC_;i++)
    {
        for(unsigned int j=0;j<3;j++)
        {
            fluxes_(i,j) = fluxes_(i,j)*thermodynamicsMap_.MW()[i+1]*1000.;
        }
    }

}


dustyGas::~dustyGas()
{

}

#endif // DUSTYGAS_H
