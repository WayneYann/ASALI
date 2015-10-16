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
        for (unsigned int i=0;i<NP_;i++)
        {
            for (unsigned int j=1;j<=NC_;j++)
                omegaBulk_[i][j] = yOS_[counter++];
            for (unsigned int j=1;j<=SURF_NC_;j++)
                teta_[i][j] = yOS_[counter++];
        }
        t_ = t;
    }

    //- Set operating conditions
    {
        //-Setting conditions for GAS phase
        thermodynamicsMap_.SetTemperature(T0_);
        thermodynamicsMap_.SetPressure(p_);
        kineticsMap_.SetTemperature(T0_);
        kineticsMap_.SetPressure(p_);
        thermodynamicsSurfaceMap_.SetPressure(p_);
        thermodynamicsSurfaceMap_.SetTemperature(T0_);
        kineticsSurfaceMap_.SetPressure(p_);
        kineticsSurfaceMap_.SetTemperature(T0_);
        transportMap_.SetPressure(p_);
        transportMap_.SetTemperature(T0_);
    }

    //- Create derivate of species
    {
        OpenSMOKE::OpenSMOKEVectorDouble omegabulk_(NP_);
        OpenSMOKE::OpenSMOKEVectorDouble omegagas_(NP_);
        for (unsigned int j=1;j<=NC_;j++)
        {
            for (unsigned int i=0;i<NP_;i++)
            {
                omegabulk_[i+1]   = omegaBulk_[i][j];
            }

            omegagas_ = FirstOrderDerivate(omegabulk_,"FDS");

            for (unsigned int i=0;i<NP_;i++)
            {
                cGas_[i][j] = omegagas_[i+1];
            }
        }
    }

    //- Estimate fluxes
    {
        for (unsigned int i=0;i<NP_;i++)
        {
            thermodynamicsMap_.MoleFractions_From_MassFractions(xBulk_,MWbulk_,omegaBulk_[i]);

            // Diffusity [m2/s]
            transportMap_.MassDiffusionCoefficients(diffG_, xBulk_);

            if ( diffusionModel_ == "Fick" )
            {
                cTotBulk_ = p_/(PhysicalConstants::R_J_kmol*T0_);
                rhoBulk_  = cTotBulk_*MWbulk_;

                for (unsigned int j=1;j<=NC;j++)
                {
                    jBulk_[i][j] = diffG_[j]*cGas_[i][j]*rhoBulk_;
                }
            }
            else if ( diffusionModel_ == "DustyGas" )
            {
                dgm_.resize();
                dgm_.setCatalystProperties(epsi_,tau_,pore_);
                dgm_.setPermeabilityLaw("none");

                Eigen::MatrixXd omegaGrad(NC,3);
                Eigen::VectorXd pGrad(3);

                for (unsigned int j=0;j<3;j++)
                {
                    pGrad(j) = 0.;
                }

                for (unsigned int j=0;j<NC;j++)
                {
                    omegaGrad(j,0) = cGas_[i][j+1];
                    omegaGrad(j,1) = 0.;
                    omegaGrad(j,2) = 0.;
                }

                dgm_.setTemperature(T0_);
                dgm_.setPressure(p_);
                dgm_.setComposition(omegaBulk_[i]);
                dgm_.setCompositionGradient(omegaGrad);
                dgm_.setPressureGradient(pGrad);
                dgm_.solve();

                for (unsigned int j=0;j<NC;j++)
                {
                    jBulk_[i][j+1] = -dgm_.getFluxes()(j,0);
                }
            }
            else if ( diffusionModel_ == "Fick-DustyGas" )
            {
                dgm_.resize();
                dgm_.setCatalystProperties(epsi_,tau_,pore_);
                dgm_.setPermeabilityLaw("none");

                Eigen::MatrixXd omegaGrad(NC,3);
                Eigen::VectorXd pGrad(3);

                for (unsigned int j=0;j<3;j++)
                {
                    pGrad(j) = 0.;
                }

                for (unsigned int j=0;j<NC;j++)
                {
                    omegaGrad(j,0) = cGas_[i][j+1];
                    omegaGrad(j,1) = 0.;
                    omegaGrad(j,2) = 0.;
                }

                dgm_.setTemperature(T0_);
                dgm_.setPressure(p_);
                dgm_.setComposition(omegaBulk_[i]);
                dgm_.setCompositionGradient(omegaGrad);
                dgm_.setPressureGradient(pGrad);
                dgm_.solve();

                cTotBulk_ = p_/(PhysicalConstants::R_J_kmol*T0_);
                rhoBulk_  = cTotBulk_*MWbulk_;

                for (unsigned int j=0;j<NC;j++)
                {
                    if ( cGas_[i][j+1] != 0. )
                    {
                        diffG_[j+1] = -dgm_.getFluxes()(j,0)/(cGas_[i][j+1]*rhoBulk_);
                    }
                }

                for (unsigned int j=1;j<=NC;j++)
                {
                    jBulk_[i][j] = diffG_[j]*cGas_[i][j]*rhoBulk_;
                }

            }
        }
    }

    //- Create derivate of fluxes
    {
        //- Species
        OpenSMOKE::OpenSMOKEVectorDouble jbulk_(NP_);
        OpenSMOKE::OpenSMOKEVectorDouble jgas_(NP_);
        for (unsigned int j=1;j<=NC_;j++)
        {
            for (unsigned int i=0;i<NP_;i++)
            {
                jbulk_[i+1]   = jBulk_[i][j];
            }

            jgas_ = FirstOrderDerivate(jbulk_,"BDS");

            for (unsigned int i=0;i<NP_;i++)
            {
                jGas_[i][j] = jgas_[i+1];
            }
        }
    }

    unsigned int counter = 1;
    for (unsigned int i=0;i<NP_;i++)
    {
        thermodynamicsMap_.MoleFractions_From_MassFractions(xBulk_,MWbulk_,omegaBulk_[i]);
        cTotBulk_ = p_/(PhysicalConstants::R_J_kmol*T0_);
        rhoBulk_  = cTotBulk_*MWbulk_;
        Product(cTotBulk_, xBulk_, &cBulk_);

        if ( reactions_ == true )
        {
            OpenSMOKE::OpenSMOKEVectorDouble fakephase(SURF_NP_);
            OpenSMOKE::OpenSMOKEVectorDouble fakebulk(NC_);
            OpenSMOKE::OpenSMOKEVectorDouble dummy;
            for(unsigned int j=1;j<=SURF_NP_;j++)
                fakephase[j] = thermodynamicsSurfaceMap_.matrix_densities_site_phases()[0][j-1];

            SD_ = fakephase[1];

            kineticsSurfaceMap_.ReactionEnthalpiesAndEntropies();
            kineticsSurfaceMap_.ArrheniusKineticConstants();
            kineticsSurfaceMap_.ReactionRates(cBulk_, teta_[i], fakebulk, fakephase);
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


        //- Equations
        if(i==0)
        {
            //-INLET conditions
            
            for (unsigned int j=1;j<=NC_;j++)
            {
                if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
                {
                    dyOS_[counter++] = x0bulk_[j] - omegaBulk_[i][j];// ALGEBRAIC
                }
                else
                {
                    dyOS_[counter++] = 1. - omegaBulk_[i].SumElements();
                }
            }
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
                {
                    dyOS_[counter++] = Rsurface_[j]/SD_; //DIFFERENTIAL
                }
                else
                {
                    dyOS_[counter++] = 1. - teta_[i].SumElements();
                }
            }
        }
        else if (i==(NP_-1))
        {
            //-OUTLET conditions
            for (unsigned int j=1;j<=NC_;j++)
            {
                if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
                {
                    dyOS_[counter++] = (A_*jGas_[i][j] + RfromSurface_[j]*alfa_)/rhoBulk_; //DIFFERENTIAL
                }
                else
                {
                    dyOS_[counter++] = 1. - omegaBulk_[i].SumElements();
                }
            }
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
                {
                    dyOS_[counter++] = Rsurface_[j]/SD_; //DIFFERENTIAL
                }
                else
                {
                    dyOS_[counter++] = 1. - teta_[i].SumElements();
                }
            }
        }
        else
        {
            for (unsigned int j=1;j<=NC_;j++)
            {
                if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
                {
                    dyOS_[counter++] = (A_*jGas_[i][j] + RfromSurface_[j]*alfa_)/rhoBulk_; //DIFFERENTIAL
                }
                else
                {
                    dyOS_[counter++] = 1. - omegaBulk_[i].SumElements();// ALGEBRAIC
                }
            }
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
                {
                    dyOS_[counter++] = Rsurface_[j]/SD_; //DIFFERENTIAL
                }
                else
                {
                    dyOS_[counter++] = 1. - teta_[i].SumElements();
                }
            }
        }
    }
}
