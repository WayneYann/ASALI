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
	ASALI::POSTprocessing  post(NC, SURF_NC);
	post.resize(N,input.getSamplingGridDimension());
	post.setGrid(z);
	post.setVariables(yF);
	post.setSampledGrid(input.getSamplingGrid());
	
	if ( input.getSamplingGridDimension() > N )
	{
		growGrid = new ASALI::GROWgrid(NC, SURF_NC);

		growGrid->setAddPoints(input.getAddPointsNumber());
		growGrid->setMaxPoints(input.getMaxPointsNumber());
		growGrid->setTol(input.getGridError());
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
				yF[counter] = growGrid->getInterpolatedVariables()[counter];
				counter++;
			}
			for (unsigned int j=1;j<=NC;j++)
			{
				yF[counter] = growGrid->getInterpolatedVariables()[counter];
				counter++;
			}
			for (unsigned int j=1;j<=SURF_NC;j++)
			{
				yF[counter] = growGrid->getInterpolatedVariables()[counter];
				counter++;
			}
			{
				yF[counter] = growGrid->getInterpolatedVariables()[counter];
				counter++;
			}
			{
				yF[counter] = growGrid->getInterpolatedVariables()[counter];
				counter++;
			}
		}

		growGrid->recapOnScreen();

		N = Nnew;
	}
	else if ( input.getSamplingGridDimension() < N )
	{
		post.sampling();
		unsigned int Nnew    = input.getSamplingGridDimension();
					 NC      = thermodynamicsMapXML->NumberOfSpecies();
					 SURF_NC = thermodynamicsSurfaceMapXML->number_of_site_species();
					 NE      = (NC + NC + SURF_NC + 1 + 1)*Nnew;

		ChangeDimensions(Nnew,    &z,      true);
		for (unsigned int k=1;k<=Nnew;k++)
			z[k] = input.getSamplingGrid()[k-1];

		ChangeDimensions(NE, &yF, true);
		unsigned int counter = 1;
		for (unsigned int i=1;i<=Nnew;i++)
		{
			for (unsigned int j=1;j<=NC;j++)
			{
				yF[counter] = post.getSampledVariables()[counter];
				counter++;
			}
			for (unsigned int j=1;j<=NC;j++)
			{
				yF[counter] = post.getSampledVariables()[counter];
				counter++;
			}
			for (unsigned int j=1;j<=SURF_NC;j++)
			{
				yF[counter] = post.getSampledVariables()[counter];
				counter++;
			}
			{
				yF[counter] = post.getSampledVariables()[counter];
				counter++;
			}
			{
				yF[counter] = post.getSampledVariables()[counter];
				counter++;
			}
		}

		N = Nnew;
	}
}
