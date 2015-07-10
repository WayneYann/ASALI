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
    iterator_++;

    Time_.resize(iterator_);
    Mass_.resize(iterator_);
    Volume_.resize(iterator_);
    Pressure_.resize(iterator_);
    Temperature_.resize(iterator_);

    Specie_.resize(iterator_*NC_);

    Site_.resize(iterator_*SURF_NC_);

    Phase_.resize(iterator_*SURF_NP_);

	{
		unsigned int k=1;

		for (unsigned int i=(iterator_-1)*NC_;i<(NC_*iterator_);i++)
			Specie_[i] = yOS_[k++];

		Mass_[iterator_-1] = yOS_[k++];

		for (unsigned int i=(iterator_-1)*SURF_NP_;i<(SURF_NP_*iterator_);i++)
			Phase_[i] = yOS_[k++];

		for (unsigned int i=(iterator_-1)*SURF_NC_;i<(SURF_NC_*iterator_);i++)
			Site_[i] = yOS_[k++];

		Temperature_[iterator_-1] = yOS_[k++];
	}

	Pressure_[iterator_-1] = P_;
	Volume_[iterator_-1]   = V_;
	Time_[iterator_-1]     = t_;

    std::cout.precision(6);
    std::cout.setf(std::ios::scientific);
    std::cout << "Integration time: " << t_ << "\tResidual: " << dyOS_.Norm2()/dyOS_.Size() << "\n" << std::endl;

}
