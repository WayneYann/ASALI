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
    //-Assigned variables
    {
        unsigned int counter = 1;
        for (unsigned int j=1;j<=NC_;j++)
            omegaBulk_[j] = y[counter++];
        for (unsigned int j=1;j<=NC_;j++)
            omegaWall_[j] = y[counter++];
        for (unsigned int j=1;j<=SURF_NC_;j++)
            teta_[j] = y[counter++];
        t_ = t;
    }

    //-Setting conditions for GAS phase
    {
        thermodynamicsMap_.SetTemperature(T_);
        thermodynamicsMap_.SetPressure(p_);
        kineticsMap_.SetTemperature(T_);
        kineticsMap_.SetPressure(p_);
        thermodynamicsSurfaceMap_.SetPressure(p_);
        thermodynamicsSurfaceMap_.SetTemperature(T_);
        kineticsSurfaceMap_.SetPressure(p_);
        kineticsSurfaceMap_.SetTemperature(T_);
        transportMap_.SetPressure(p_);
        transportMap_.SetTemperature(T_);
    }

    //-Evaluating GAS phase properties
    {
        thermodynamicsMap_.MoleFractions_From_MassFractions(xBulk_,MWbulk_,omegaBulk_);
        cTotBulk_ = p_/(PhysicalConstants::R_J_kmol*T_);
        rhoBulk_ = cTotBulk_*MWbulk_;
        Product(cTotBulk_, xBulk_, &cBulk_);

        if ( homogeneusReactions_ == true )
        {
            kineticsMap_.ReactionRates(cBulk_);
            kineticsMap_.FormationRates(&RfromGas_);
            ElementByElementProduct(RfromGas_, thermodynamicsMap_.MW(), &RfromGas_);
        }
        else
        {
            for (unsigned int j=1;j<=NC_;j++)
                RfromGas_[j] = 0.;
        }

        //- Diffusity [m2/s]
        transportMap_.MassDiffusionCoefficients(diffG_, xBulk_);

        // Viscosity [kg/m/s]
        transportMap_.DynamicViscosity(etaMix_, xBulk_);

        //Mass Transfer coefficient
        MassTransferCoefficient(t_);

        // Fluxes from gas to solid [m/s]
        for (unsigned int j=1;j<=NC_;j++)
            jSolid_[j] = Kmat_[j]*(omegaBulk_[j] - omegaWall_[j]);
    }
    
    //-Evaluating properties for SOLID phase
    {
        thermodynamicsMap_.MoleFractions_From_MassFractions(xWall_, MWwall_, omegaWall_);
        cTotWall_ = p_/(PhysicalConstants::R_J_kmol*T_);
        rhoWall_ = cTotWall_*MWwall_;
        Product(cTotWall_, xWall_, &cWall_); 

        if ( heterogeneusReactions_ == true )
        {
            OpenSMOKE::OpenSMOKEVectorDouble fakephase(SURF_NP_);
            OpenSMOKE::OpenSMOKEVectorDouble fakebulk(NC_);

            for(unsigned int j=1;j<=SURF_NP_;j++)
                fakephase[j] = thermodynamicsSurfaceMap_.matrix_densities_site_phases()[0][j-1];

            SD_ = fakephase[1];

            kineticsSurfaceMap_.ReactionEnthalpiesAndEntropies();
            kineticsSurfaceMap_.ReactionRates(cWall_, teta_, fakebulk, fakephase);
            kineticsSurfaceMap_.FormationRates(&RfromSurface_, &Rsurface_, &fakebulk, &RsurfacePhases_);
            ElementByElementProduct(RfromSurface_, thermodynamicsMap_.MW(), &RfromSurface_);
        }
        else
        {
            OpenSMOKE::OpenSMOKEVectorDouble fakephase(SURF_NP_);

            for(unsigned int j=1;j<=SURF_NP_;j++)
                fakephase[j] = thermodynamicsSurfaceMap_.matrix_densities_site_phases()[0][j-1];

            SD_ = fakephase[1];

            for (unsigned int j=1;j<=NC_;j++)
                RfromSurface_[j] = 0.;
            for (unsigned int j=1;j<=SURF_NC_;j++)
                Rsurface_[j] = 0.;
        }
    }

    //- Equations
    {
        unsigned int counter = 1;
        for (unsigned int j=1;j<=NC_;j++)
        {
            if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
            {
                dyOS_[counter++] = - av_*rhoBulk_*jSolid_[j]/G_
                                   + RfromGas_[j]/G_; //DIFFERENTIAL
            }
            else
            {
                dyOS_[counter++] = 1. - omegaBulk_.SumElements();// ALGEBRAIC
            }
        }
        for (unsigned int j=1;j<=NC_;j++)
        {
            if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
            {
                dyOS_[counter++] = (+ av_*rhoBulk_*jSolid_[j] + RfromSurface_[j]*alfa_); // ALGEBRAIC
            }
            else
            {
                dyOS_[counter++] = 1. - omegaWall_.SumElements();// ALGEBRAIC
            }
        }
        for (unsigned int j=1;j<=SURF_NC_;j++)
        {
            if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != catalyst_ )
            {
                dyOS_[counter++] = Rsurface_[j]; //ALGEBRAIC
            }
            else
            {
                dyOS_[counter++] = 1. - teta_.SumElements();
            }
        }
    }
}
