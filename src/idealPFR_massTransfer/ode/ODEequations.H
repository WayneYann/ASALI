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


		// Diffusity [m2/s]
		transportMap_.MassDiffusionCoefficients(diffG_, x_);

		// Viscosity [kg/m/s]
		transportMap_.DynamicViscosity(mu_, x_);

		//- Specific area [1/m]
		av();

		//- Reynolds number
		Re();
		
		//- Schimdt number
		Sc();
		
		//- Fanning friction factor
		FrictionFactor();
		
		//- Nusselt number
		Sherwood();

		//- Heat transfer coefficient [W/m2/K]
		kMat();
		
	}
	
	for (unsigned int i=1;i<=NC_;i++)
	{
		dyOS_[i] = -kMat_[i]*av_*(omega_[i] - omegaW_[i]);
	}
	dyOS_[NE_] = -rho_*v_*v_*fm_/Dp_;
}
