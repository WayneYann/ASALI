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
	std::vector<double> p;
	std::vector<double> comp;
	std::vector<double> z;
	ode->GetProfile(z,p,comp);
	
	OpenSMOKE::OpenSMOKEVectorDouble *omega;

	omega = new OpenSMOKE::OpenSMOKEVectorDouble[z.size()];

	for (unsigned int i=0;i<z.size();i++)
	{
		ChangeDimensions(NC,      &omega[i], true);
	}

	{
		unsigned int conter = 0;
		for (unsigned int k=0;k<z.size();k++)
		{
			for (unsigned int i=1;i<=NC;i++)
			{
				omega[k][i] = comp[conter++];
			}
		}
	}


	boost::filesystem::remove_all("results");
	boost::filesystem::create_directory("results");
	std::ofstream output;
	output.setf(std::ios::scientific);
	output.precision(6);
	output.open("results/pressure.txt",std::ios::out);
	output << "#Length" << "\t" << "#Pressure" << std::endl;
	for (unsigned int k=0;k<p.size();k++)
		output << z[k] << "\t" << p[k] << std::endl;
	output.close();

	std::ofstream specie;
	specie.setf(std::ios::scientific);
	specie.precision(6);
	specie.open("results/massFraction.txt",std::ios::out);
	specie << "#Length" << "\t";

	for (unsigned int k=1;k<=NC;k++)
	{
		specie << thermodynamicsMapXML->NamesOfSpecies()[k-1] << "\t";
	}
	specie << std::endl;

	for (unsigned int k=0;k<p.size();k++)
	{
		specie << z[k] << "\t";
		for (unsigned int j=1;j<=NC;j++)
		{
			specie << omega[k][j] << "\t";
		}
		specie << std::endl;
	}
	specie.close();

	std::ofstream recap;
	recap.setf(std::ios::scientific);
	recap.precision(6);
	recap.open("results/information.txt",std::ios::out);
	recap << "Input:\n" << std::endl;
	recap << " 1/ pressure correlation:      " << input.getPressureCorrelation() << std::endl;
	recap << " 2/ mass transfer correlation: " << input.getMassTransferCorrelation() << std::endl;
	recap << " 3/ feed temperature:          " << input.getFeedTemperature() << " [K]" << std::endl;
	recap << " 4/ feed pressure:             " << input.getPressureCorrelation() << " [Pa]" << std::endl;
	recap << " 5/ feed velocity:             " << input.getFeedVelocity() << " [m/s]" << std::endl;

	if ( input.getReactorType() == "Monolith" )
	{
		recap << " 6/ void fraction:             " << input.getVoidFraction() << " [-]" << std::endl;
	}
	else if ( input.getReactorType() == "PackedBed" )
	{
		recap << " 6/ particle diameter:         " << input.getParticleDiameter() << " [m]" << std::endl;
	}

	recap << " 7/ tube diameter:             " << input.getShellDiameter() << " [m]" << std::endl;
	recap << " 8/ tube length:               " << input.getReactorLength()  << " [m]" << std::endl;
	recap << std::endl;
	recap << "Output:\n" << std::endl;
	recap << " 1/ pressure drops:            " << input.getFeedPressure() - yF[NC+1] << " [Pa]" << std::endl;
	recap << " 2/ pressure drops:            " << (input.getFeedPressure() - yF[NC+1])*1e-05/input.getReactorLength() << " [bar/m]" << std::endl;
	recap.close();
}

