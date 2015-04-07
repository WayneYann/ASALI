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
			for (unsigned int j=1;j<=NC_;j++)
				omegaWall_[i][j] = yOS_[counter++];
			for (unsigned int j=1;j<=SURF_NC_;j++)
				teta_[i][j] = yOS_[counter++];
			Tbulk_[i+1] = yOS_[counter++];
			Twall_[i+1] = yOS_[counter++];
		}
		t_ = t;
	}

	unsigned int counter = 1;
	for (unsigned int i=0;i<NP_;i++)
	{
		//-Setting conditions for GAS phase
		thermodynamicsMap_.SetTemperature(Tbulk_[i+1]);
		thermodynamicsMap_.SetPressure(p_);
		kineticsMap_.SetTemperature(Tbulk_[i+1]);
		kineticsMap_.SetPressure(p_);
		thermodynamicsSurfaceMap_.SetPressure(p_);
		thermodynamicsSurfaceMap_.SetTemperature(Tbulk_[i+1]);
		kineticsSurfaceMap_.SetPressure(p_);
		kineticsSurfaceMap_.SetTemperature(Tbulk_[i+1]);
		transportMap_.SetPressure(p_);
		transportMap_.SetTemperature(Tbulk_[i+1]);

		//-Evaluating GAS phase properties
		thermodynamicsMap_.MoleFractions_From_MassFractions(xBulk_,MWbulk_,omegaBulk_[i]);

		cTotBulk_ = p_/(PhysicalConstants::R_J_kmol*Tbulk_[i+1]); 	// [Kmol/m3]
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

		// Diffusity [m2/s]
		transportMap_.MassDiffusionCoefficients(diffG_, xBulk_);

		// Fluxes in gas[m/s]
		for (unsigned int j=1;j<=NC_;j++)
		{
			if ( i == 0)
			{
				jGas_[i][j] = 0.;
			}
			else if ( i == (NP_ - 1) )
			{
				jGas_[i][j] = 0.;
			}
			else
			{
				jGas_[i][j] = (omegaBulk_[i+1][j]*Dz_[i+1-1] + omegaBulk_[i-1][j]*Dz_[i+1+1-1]-omegaBulk_[i][j]*(Dz_[i+1+1-1]+Dz_[i+1-1]))/(Dz_[i+1+1-1]*Dz_[i+1+1-1]*Dz_[i+1-1]);
			}
		}

		// Viscosity [kg/m/s]
		transportMap_.DynamicViscosity(etaMix_, xBulk_);

		//- Specific heat [J/Kg/K]
		thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(cp_, xBulk_);
		cp_ = cp_/MWbulk_;

		//- Thermal conductivity [W/m/K]
		transportMap_.ThermalConductivity(condBulk_,xBulk_);

		MassTransferCoefficient(z_[i]);

		// Fluxes from gas to solid [m/s]
		for (unsigned int j=1;j<=NC_;j++)
			jSolid_[i][j] = Kmat_[j]*(omegaBulk_[i][j] - omegaWall_[i][j]);

		//-Setting conditions for SOLID phase
		thermodynamicsMap_.SetTemperature(Twall_[i+1]);
		thermodynamicsMap_.SetPressure(p_);
		kineticsMap_.SetTemperature(Twall_[i+1]);
		kineticsMap_.SetPressure(p_);
		thermodynamicsSurfaceMap_.SetPressure(p_);
		thermodynamicsSurfaceMap_.SetTemperature(Twall_[i+1]);
		kineticsSurfaceMap_.SetPressure(p_);
		kineticsSurfaceMap_.SetTemperature(Twall_[i+1]);
		transportMap_.SetPressure(p_);
		transportMap_.SetTemperature(Twall_[i+1]);
		
		//-Evaluating properties for SOLID phase
		thermodynamicsMap_.MoleFractions_From_MassFractions(xWall_, MWwall_, omegaWall_[i]);
		cTotWall_ = p_/(PhysicalConstants::R_J_kmol*Twall_[i+1]);
		rhoWall_ = cTotWall_*MWwall_;
		Product(cTotWall_, xWall_, &cWall_); 

		if ( heterogeneusReactions_ == true )
		{
			OpenSMOKE::OpenSMOKEVectorDouble fakephase(SURF_NP_);
			OpenSMOKE::OpenSMOKEVectorDouble fakebulk(NC_);
			OpenSMOKE::OpenSMOKEVectorDouble dummy;
			for(unsigned int j=1;j<=SURF_NP_;j++)
				fakephase[j] = thermodynamicsSurfaceMap_.matrix_densities_site_phases()[0][j-1];

			SD_ = fakephase[1];

			kineticsSurfaceMap_.ReactionEnthalpiesAndEntropies();
			kineticsSurfaceMap_.ArrheniusKineticConstants();
			kineticsSurfaceMap_.ReactionRates(cWall_, teta_[i], fakebulk, fakephase);
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

		if ( z_[i+1] <= Linert_)
		{
			for (unsigned int j=1;j<=NC_;j++)
				RfromSurface_[j] = 0.;
			for (unsigned int j=1;j<=SURF_NC_;j++)
				Rsurface_[j] = 0.;
			alfa_ = 0;
		}
		else
		{
			alfa_ = alfaTemp_;
		}

		//Gas thermal flux
		if (i==0)
		{
			jbulk_[i+1] = 0.;
		}
		else if( i == (NP_ - 1) )
		{
			jbulk_[i+1] = 0.;
		}
		else
		{
			jbulk_[i+1] = (Tbulk_[i+1+1]*Dz_[i+1-1] + Tbulk_[i-1+1]*Dz_[i+1+1-1]-Tbulk_[i+1]*(Dz_[i+1+1-1]+Dz_[i+1-1]))/(Dz_[i+1+1-1]*Dz_[i+1+1-1]*Dz_[i+1-1]);
		}
		
		//Solid thermal flux
		if (i==0)
		{
			jwall_[i+1] = 0.;
		}
		else if( i == (NP_ - 1) )
		{
			jwall_[i+1] = 0.;
		}
		else
		{
			jwall_[i+1] = (Twall_[i+1+1]*Dz_[i+1-1] + Twall_[i-1+1]*Dz_[i+1+1-1]-Twall_[i+1]*(Dz_[i+1+1-1]+Dz_[i+1-1]))/(Dz_[i+1+1-1]*Dz_[i+1+1-1]*Dz_[i+1-1]);
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
					dyOS_[counter++] = 1e06*(1. - omegaBulk_[i].SumElements());
				}
			}
			for (unsigned int j=1;j<=NC_;j++)
			{
				if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
				{
					dyOS_[counter++] = (+ av_*rhoBulk_*jSolid_[i][j]*epsi_ + RfromSurface_[j]*alfa_*epsi_); // ALGEBRAIC
				}
				else
				{
					dyOS_[counter++] = 1e06*(1. - omegaWall_[i].SumElements());// ALGEBRAIC
				}
			}
			for (unsigned int j=1;j<=SURF_NC_;j++)
			{
				if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
				{
					dyOS_[counter++] = 0.; //DIFFERENTIAL
				}
				else
				{
					dyOS_[counter++] = 1e06*(1. - teta_[i].SumElements());
				}
			}
			dyOS_[counter++] = T0_ - Tbulk_[i+1]; // + (condBulk_/(G_*cp_))*((Tbulk_[i+1+1]-Tbulk_[i+1])/Dz_[i+1]); // ALGEBRAIC
			dyOS_[counter++] = Twall_[i+1+1] - Twall_[i+1]; // ALGEBRAIC
		}
		else if (i==(NP_-1))
		{
			//-OUTLET conditions
			for (unsigned int j=1;j<=NC_;j++)
			{
				if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
				{
					dyOS_[counter++] = omegaBulk_[i][j]-omegaBulk_[i-1][j]; // ALGEBRAIC
				}
				else
				{
					dyOS_[counter++] = 1e06*(1. - omegaBulk_[i].SumElements());
				}
			}
			for (unsigned int j=1;j<=NC_;j++)
			{
				if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
				{
					dyOS_[counter++] = (+ av_*rhoBulk_*jSolid_[i][j]*epsi_ + RfromSurface_[j]*alfa_*epsi_); // ALGEBRAIC
				}
				else
				{
					dyOS_[counter++] = 1e06*(1. - omegaWall_[i].SumElements());// ALGEBRAIC
				}
			}
			for (unsigned int j=1;j<=SURF_NC_;j++)
			{
				if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
				{
					dyOS_[counter++] = 0.; //DIFFERENTIAL
				}
				else
				{
					dyOS_[counter++] = 1e06*(1. - teta_[i].SumElements());
				}
			}
			dyOS_[counter++] = Tbulk_[i+1] - Tbulk_[i-1+1]; // ALGEBRAIC
			dyOS_[counter++] = Twall_[i+1] - Twall_[i-1+1]; // ALGEBRAIC
		}
		else
		{
			for (unsigned int j=1;j<=NC_;j++)
			{
				if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
				{
					dyOS_[counter++] = 0.; //DIFFERENTIAL
				}
				else
				{
					dyOS_[counter++] = 1e06*(1. - omegaBulk_[i].SumElements());// ALGEBRAIC
				}
			}
			for (unsigned int j=1;j<=NC_;j++)
			{
				if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
				{
					dyOS_[counter++] = (+ av_*rhoBulk_*jSolid_[i][j]*epsi_ + RfromSurface_[j]*alfa_*epsi_); // ALGEBRAIC
				}
				else
				{
					dyOS_[counter++] = 1e06*(1. - omegaWall_[i].SumElements());// ALGEBRAIC
				}
			}
			for (unsigned int j=1;j<=SURF_NC_;j++)
			{
				if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
				{
					dyOS_[counter++] = 0.; //DIFFERENTIAL
				}
				else
				{
					dyOS_[counter++] = 1e06*(1. - teta_[i].SumElements());
				}
			}
			if (energy_ == true)
			{
				dyOS_[counter++]=0.; //DIFFERENTIAL
				dyOS_[counter++]=0.; //DIFFERENTIAL
			}
			else
			{
				dyOS_[counter++]=0.; //DIFFERENTIAL
				dyOS_[counter++]=0.; //DIFFERENTIAL
			}
		}
	}
}
