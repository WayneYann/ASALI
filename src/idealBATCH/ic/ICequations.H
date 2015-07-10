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

if ( resolution_ == "volume" )
{
	// Recover unknowns
	{
		unsigned int k=1;

		for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = yOS_[k++];

		mass_ = yOS_[k++];

		for(unsigned int i=1;i<=SURF_NP_;++i)
			Gamma_[i] = yOS_[k++];

		for(unsigned int i=1;i<=SURF_NC_;++i)
			Z_[i] = yOS_[k++];

		T_ = yOS_[k++];
	}

	// Calculates the pressure and the concentrations of species
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);
	const double rho_  = mass_/V_;
	const double cTot_ = rho_/MW_;
	Product(cTot_, x_, &c_);
	P_ = cTot_ * PhysicalConstants::R_J_kmol * T_;

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T_);
	thermodynamicsMap_.SetPressure(P_);

	// Calculates heterogeneous kinetics
	{
		thermodynamicsSurfaceMap_.SetPressure(P_);
		thermodynamicsSurfaceMap_.SetTemperature(T_);
		kineticsSurfaceMap_.SetPressure(P_);
		kineticsSurfaceMap_.SetTemperature(T_);

		if ( heterogeneusReactions_ == true)
		{
			kineticsSurfaceMap_.ReactionEnthalpiesAndEntropies();
			kineticsSurfaceMap_.ArrheniusKineticConstants();
			kineticsSurfaceMap_.ReactionRates(c_, Z_, dummy, Gamma_);
			kineticsSurfaceMap_.FormationRates(&RfromSurface_, &Rsurface_, &dummy, &RsurfacePhases_);
		}
		else
		{
			for(unsigned int i=1;i<=NC_;++i)
				RfromSurface_[i] = 0.;

			for(unsigned int i=1;i<=SURF_NP_;++i)
				RsurfacePhases_[i] = 0.;

			for(unsigned int i=1;i<=SURF_NC_;++i)
				Rsurface_[i] = 0.;
		}
	}

	// Recovering residuals
	{
		unsigned int k=1;
		
		// Gas phase species
		for (unsigned int i=1;i<=NC_;++i)
		{
			dyOS_[k++] = 0.;
		}

		// Total mass
		dyOS_[k++] = 0.;

		// Phases
		for (unsigned int i=1;i<=SURF_NP_;++i)
		{
			dyOS_[k++] = 0.;
		}
		
		// Surface site species
		for (unsigned int i=1;i<=SURF_NC_;++i)	
		{
			const unsigned int index_phase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[i-1]+1;
			dyOS_[k++] = (thermodynamicsSurfaceMap_.vector_occupancies_site_species()[i-1]*Rsurface_[i] - 
									Z_[i]*dy[NC_+1+index_phase])/Gamma_[index_phase];
		}
		
		dyOS_[k++] = 0.;
	}
}
else if ( resolution_ == "pressure" )
{
	// Recover unknowns
	{
		unsigned int k=1;

		for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = yOS_[k++];

		mass_ = yOS_[k++];

		for(unsigned int i=1;i<=SURF_NP_;++i)
			Gamma_[i] = yOS_[k++];

		for(unsigned int i=1;i<=SURF_NC_;++i)
			Z_[i] = yOS_[k++];

		T_ = yOS_[k++];

	}

	// Calculates the volume and the concentrations of species
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);
	const double cTot_ = P_/(PhysicalConstants::R_J_kmol * T_);
	Product(cTot_, x_, &c_);

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T_);
	thermodynamicsMap_.SetPressure(P_);

	// Calculates heterogeneous kinetics
	{
		thermodynamicsSurfaceMap_.SetPressure(P_);
		thermodynamicsSurfaceMap_.SetTemperature(T_);
		kineticsSurfaceMap_.SetPressure(P_);
		kineticsSurfaceMap_.SetTemperature(T_);

		if ( heterogeneusReactions_ == true)
		{
			kineticsSurfaceMap_.ReactionEnthalpiesAndEntropies();
			kineticsSurfaceMap_.ArrheniusKineticConstants();
			kineticsSurfaceMap_.ReactionRates(c_, Z_, dummy, Gamma_);
			kineticsSurfaceMap_.FormationRates(&RfromSurface_, &Rsurface_, &dummy, &RsurfacePhases_);
		}
		else
		{
			for(unsigned int i=1;i<=NC_;++i)
				RfromSurface_[i] = 0.;

			for(unsigned int i=1;i<=SURF_NP_;++i)
				RsurfacePhases_[i] = 0.;

			for(unsigned int i=1;i<=SURF_NC_;++i)
				Rsurface_[i] = 0.;
		}
	}

	// Recovering residuals
	{
		unsigned int k=1;
		
		// Gas phase species
		for (unsigned int i=1;i<=NC_;++i)
		{
			dyOS_[k++] = 0.;
		}

		// Total mass
		dyOS_[k++] = 0.;

		// Phases
		for (unsigned int i=1;i<=SURF_NP_;++i)
		{
			dyOS_[k++] = 0.;
		}
		
		// Surface site species
		for (unsigned int i=1;i<=SURF_NC_;++i)	
		{
			const unsigned int index_phase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[i-1]+1;
			dyOS_[k++] = (thermodynamicsSurfaceMap_.vector_occupancies_site_species()[i-1]*Rsurface_[i] - 
									Z_[i]*dy[NC_+1+index_phase])/Gamma_[index_phase];
		}
		
		dyOS_[k++] = 0.;
	}
}
