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

{
    //- Assigned variable
    {
        for (unsigned int i=1;i<=NC_;i++)
        {
            omega_[i] = yOS_[i];
        }
        for (unsigned int i=1;i<=NC_;i++)
        {
            omegaW_[i] = yOS_[i+NC_];
        }

        p_ = yOS_[NE_];
        t_ = t;
    }

    //-Setting conditions for GAS phase
    {
        thermodynamicsMap_.SetTemperature(T_);
        thermodynamicsMap_.SetPressure(p_);
        transportMap_.SetPressure(p_);
        transportMap_.SetTemperature(T_);
    }

    //- Evaluting properties
    {
        thermodynamicsMap_.MoleFractions_From_MassFractions(x_,MW_,omega_);

        //- Density [Kg/m3]
        rho_ = p_*MW_/(PhysicalConstants::R_J_kmol*T_);

        //- Velocity [m/s]
        v_   = G_/rho_;

        {
            thermodynamicsMap_.MoleFractions_From_MassFractions(xW_,MW_,omegaW_);
            cTot_    = p_/(PhysicalConstants::R_J_kmol*T_);
            Product(cTot_, xW_, &c_); 

            OpenSMOKE::OpenSMOKEVectorDouble fakephase(SURF_NP_);
            OpenSMOKE::OpenSMOKEVectorDouble fakebulk(NC_);
            OpenSMOKE::OpenSMOKEVectorDouble faketeta(SURF_NC_);
            OpenSMOKE::OpenSMOKEVectorDouble dummy;

            for(unsigned int j=1;j<=SURF_NP_;j++)
                fakephase[j] = thermodynamicsSurfaceMap_.matrix_densities_site_phases()[0][j-1];

            for(unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
                { 
                    faketeta[j] = 0.;
                }
                else
                {
                    faketeta[j] = 1.;
                }
            }

            kineticsSurfaceMap_.ReactionEnthalpiesAndEntropies();
            kineticsSurfaceMap_.ArrheniusKineticConstants();
            kineticsSurfaceMap_.ReactionRates(c_, faketeta, fakebulk, fakephase);
            kineticsSurfaceMap_.FormationRates(&R_, &Rsurface_, &fakebulk, &RsurfacePhases_);
        }

        // Diffusity [m2/s]
        transportMap_.MassDiffusionCoefficients(diffG_, x_);

        // Viscosity [kg/m/s]
        transportMap_.DynamicViscosity(mu_, x_);

        //- Specific area [1/m]
        av();

        //- Schimdt number
        Sc();
        
        //- Fanning friction factor
        FrictionFactor(ReynoldsForFrictionFactor());

        //- Sherwood number
        Sherwood(ReynoldsForMassTransfer());

        //- Heat transfer coefficient [W/m2/K]
        kMat();
    }
    
    for (unsigned int i=1;i<=NC_;i++)
    {
        if ( thermodynamicsMapXML->NamesOfSpecies()[i-1] != inert_ )
        {
            dyOS_[i] = -kMat_[i]*av_*(omega_[i] - omegaW_[i])/(v_/epsi_);
        }
        else
        {
            dyOS_[i] = 1. - omega_.SumElements();
        }
    }
    
    for (unsigned int i=1;i<=NC_;i++)
    {
        if ( thermodynamicsMapXML->NamesOfSpecies()[i-1] != inert_ )
        {
            dyOS_[i+NC_] = kMat_[i]*av_*(omega_[i] - omegaW_[i]) + R_[i]*thermodynamicsMap_.MW()[i]/rho_;
        }
        else
        {
            dyOS_[i+NC_] = 1. - omegaW_.SumElements();
        }
    }

    if ( type_ == "Monolith" )
    {
        dyOS_[NE_] = -rho_*v_*v_*fm_/Dt_;
    }
    else if ( type_ == "PackedBed" )
    {
        dyOS_[NE_] = -rho_*v_*v_*fm_/Dp_;
    }
}
