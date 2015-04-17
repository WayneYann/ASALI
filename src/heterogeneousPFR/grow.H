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
	growGrid->resize(N);
	growGrid->setGrid(z);
	growGrid->setVariables(yF);

	growGrid->calculateErrors();
	growGrid->addition();
	growGrid->interpolation();

	unsigned int Nnew    = growGrid->getInterpolatedPoints();
				 NC      = thermodynamicsMapXML->NumberOfSpecies();
				 SURF_NC = thermodynamicsSurfaceMapXML->number_of_site_species();
				 NE      = (NC + NC + SURF_NC + 1 + 1)*Nnew;

	ChangeDimensions(Nnew,    &z,      true);
	for (unsigned int k=1;k<=Nnew;k++)
		z[k] = growGrid->getInterpolatedGrid()[k];

	ChangeDimensions(NE, &x0, true);
	unsigned int counter = 1;
	for (unsigned int i=1;i<=Nnew;i++)
	{
		for (unsigned int j=1;j<=NC;j++)
		{
			x0[counter] = growGrid->getInterpolatedVariables()[counter];
			counter++;
		}
		for (unsigned int j=1;j<=NC;j++)
		{
			x0[counter] = growGrid->getInterpolatedVariables()[counter];
			counter++;
		}
		for (unsigned int j=1;j<=SURF_NC;j++)
		{
			x0[counter] = growGrid->getInterpolatedVariables()[counter];
			counter++;
		}
		{
			x0[counter] = growGrid->getInterpolatedVariables()[counter];
			counter++;
		}
		{
			x0[counter] = growGrid->getInterpolatedVariables()[counter];
			counter++;
		}
	}

	ended = growGrid->getExitCondition();

	growGrid->recapOnScreen();

	N = Nnew;
	
	ASALI::discretizationSchemeOnScreen(discretizationScheme);
}
