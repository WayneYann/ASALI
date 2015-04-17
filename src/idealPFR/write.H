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
	std::string folderName;

	folderName = "results/";

	std::ofstream mass;
	std::string fileNameMass = folderName + "/massFraction.txt";
	const char *pathMass = fileNameMass.c_str();
	remove(pathMass);
	mass.open(pathMass,std::ios::out);

	std::ofstream mole;
	std::string fileNameMole = folderName + "/moleFraction.txt";
	const char *pathMole = fileNameMole.c_str();
	remove(pathMole);
	mole.open(pathMole,std::ios::out);

	std::ofstream wall;
	std::string fileNameWall = folderName + "/solid.txt";
	const char *pathWall = fileNameWall.c_str();
	remove(pathWall);
	wall.open(pathWall,std::ios::out);

	std::ofstream site;
	std::string fileNameSite = folderName + "/adsorbed.txt";
	const char *pathSite = fileNameSite.c_str();
	remove(pathSite);
	site.open(pathSite,std::ios::out);

	std::ofstream extra;
	std::string fileNameExtra = folderName + "/information.txt";
	const char *pathExtra = fileNameExtra.c_str();
	remove(pathExtra);
	extra.open(pathExtra,std::ios::out);

	{
		unsigned int WP = 6;
		mass.setf(std::ios::scientific);
		mole.setf(std::ios::scientific);
		wall.setf(std::ios::scientific);
		site.setf(std::ios::scientific);
		extra.setf(std::ios::scientific);
		mass.precision(WP);
		mole.precision(WP);
		wall.precision(WP);
		site.precision(WP);
		extra.precision(WP);
	}

	OpenSMOKE::OpenSMOKEVectorDouble omegaBulk(NC);
	OpenSMOKE::OpenSMOKEVectorDouble omegaWall(NC);
	OpenSMOKE::OpenSMOKEVectorDouble conv(NC);
	OpenSMOKE::OpenSMOKEVectorDouble teta(SURF_NC);

	unsigned int counter = 1;
	{
		for (unsigned int j=1;j<=NC;j++)
		{
			omegaBulk[j] = std::max(0.,yF[counter++]);
			if ( omegaBulk[j] == 0. )
			{
				conv[j] = 0.;
			}
			else
			{
				conv[j] = (x0bulk[j] - omegaBulk[j])/x0bulk[j];
			}
		}
		for (unsigned int j=1;j<=NC;j++)
		{
			omegaWall[j] = std::max(0.,yF[counter++]);
		}
		for (unsigned int j=1;j<=SURF_NC;j++)
		{
			teta[j] = std::max(0.,yF[counter++]);
		}
	}
	

	//- GAS mass fraction
	{
		mass << "Specie" << "\t" << "Fraction IN [-]" << "\t" << "Fraction OUT [-]" << "\t" << "Conversion [-]" << std::endl;
		for (unsigned int j=1;j<=NC;j++)
		{
			mass << thermodynamicsMapXML->NamesOfSpecies()[j-1] << "\t" << x0bulk[j] << "\t" << omegaBulk[j] << "\t";
			if ( x0bulk[j] != 0.)
			{
				mass << conv[j] << std::endl;
			}
			else
			{
				mass << " " << std::endl;
			}
		}
		mass.close();
	}

	//- GAS mole fraction
	{
		OpenSMOKE::OpenSMOKEVectorDouble x(NC);
		OpenSMOKE::OpenSMOKEVectorDouble y(NC);
		double MW;
		thermodynamicsMapXML->MoleFractions_From_MassFractions(x, MW, omegaBulk);
		thermodynamicsMapXML->MoleFractions_From_MassFractions(y, MW, x0bulk);
		mole << "Specie" << "\t" << "Fraction IN [-]" << "\t" << "Fraction OUT [-]" << "\t" << "Conversion [-]" << std::endl;
		for (unsigned int j=1;j<=NC;j++)
		{
			mole << thermodynamicsMapXML->NamesOfSpecies()[j-1] << "\t" << y[j] << "\t" << x[j] << "\t";
			if ( x0bulk[j] != 0.)
			{
				mole << conv[j] << std::endl;
			}
			else
			{
				mole << " " << std::endl;
			}
		}
		mole.close();
	}

	//- SOLID fraction
	{
		OpenSMOKE::OpenSMOKEVectorDouble x(NC);
		double MW;
		thermodynamicsMapXML->MoleFractions_From_MassFractions(x, MW, omegaWall);
		wall << "Specie" << "\t" << "MASS fraction [-]" << "\t" << "MOLE fraction [-]" << std::endl;
		for (unsigned int j=1;j<=NC;j++)
		{
			wall << thermodynamicsMapXML->NamesOfSpecies()[j-1] << "\t" << omegaWall[j] << "\t" << x[j] << std::endl;
		}
		wall.close();
	}

	//- SITE fraction
	{
		site << "Specie" << "\t" << "SITE fraction [-]" << std::endl;
		for (unsigned int j=1;j<=SURF_NC;j++)
		{
			site << thermodynamicsSurfaceMapXML->NamesOfSpecies()[j-1+NC]<< "\t" << teta[j] << std::endl;
		}
		site.close();
	}

	//- EXTRA information
	{
		double MassFluxIn;
		double MassFluxOut;

		std::vector<double> IN;
		std::vector<double> OUT;

		{
			thermodynamicsMapXML->SetPressure(input.getPressure());
			thermodynamicsMapXML->SetTemperature(input.getTemperature());
			MassFluxIn = 0.;
			for (unsigned int j=1;j<=NC;j++)
			{
				MassFluxIn += x0bulk[j]*G;
			}

			Eigen::MatrixXd SpecieAtomicComposition;
			Eigen::VectorXd Element;
			std::vector<std::string> ElementsNames;
			ElementsNames = thermodynamicsMapXML->elements();
			SpecieAtomicComposition = thermodynamicsMapXML->atomic_composition();
			Element.resize(NC);
			IN.resize(ElementsNames.size());
			for (unsigned int k=0;k<ElementsNames.size();k++)
			{
				IN[k] = 0.;
				Element = SpecieAtomicComposition.col(k);
				for (unsigned int j=0;j<NC;j++)
					IN[k] += x0bulk[j+1]*Element[j]/thermodynamicsMapXML->MW()[j+1];
			}
		}

		{
			thermodynamicsMapXML->SetPressure(input.getPressure());
			thermodynamicsMapXML->SetTemperature(input.getTemperature());
				MassFluxOut = 0.;
			for (unsigned int j=1;j<=NC;j++)
			{
				MassFluxOut += omegaBulk[j]*G;
			}

			Eigen::MatrixXd SpecieAtomicComposition;
			Eigen::VectorXd Element;
			std::vector<std::string> ElementsNames;
			ElementsNames = thermodynamicsMapXML->elements();
			SpecieAtomicComposition = thermodynamicsMapXML->atomic_composition();
			Element.resize(NC);
			OUT.resize(ElementsNames.size());
			for (unsigned int k=0;k<ElementsNames.size();k++)
			{
				OUT[k] = 0.;
				Element = SpecieAtomicComposition.col(k);
				for (unsigned int j=0;j<NC;j++)
					OUT[k] += omegaBulk[j+1]*Element[j]/thermodynamicsMapXML->MW()[j+1];;
			}
		}

		extra << "\nMass flux IN             = " << MassFluxIn << "	[Kg/m2/s]" << std::endl;
		extra << "Mass flux OUT            = " << MassFluxOut << "	[Kg/m2/s]\n" << std::endl;
		extra << "ERROR on MASS BALANCE    = " << (MassFluxIn - MassFluxOut)*100./MassFluxIn << "	[%]\n" << std::endl;
		
		extra << "ATOMIC BALANCE               " << std::endl;
		for (unsigned int k=0;k<IN.size();k++)
		{
			if ( IN[k] != 0. )
				extra << thermodynamicsMapXML->elements()[k] << " :\t" << (IN[k] - OUT[k])*100./IN[k] << "\t[%]" << std::endl;
		}
		extra << "\n" << std::endl;
		
		extra.close();

		std::cout << "\n#################################################" << std::endl;
		std::cout << "\nMass flux IN             = " << MassFluxIn << "	[Kg/m2/s]" << std::endl;
		std::cout << "Mass flux OUT            = " << MassFluxOut << "	[Kg/m2/s]\n" << std::endl;
		std::cout << "ERROR on MASS BALANCE    = " << fabs((MassFluxIn - MassFluxOut)*100./MassFluxIn) << "	[%]\n" << std::endl;
		std::cout << "#################################################" << std::endl;
		std::cout << "\n############################" << std::endl;
		std::cout << "ATOMIC BALANCE               " << std::endl;
		for (unsigned int k=0;k<IN.size();k++)
		{
			if ( IN[k] != 0. )
				std::cout << thermodynamicsMapXML->elements()[k] << " :\t" << (IN[k] - OUT[k])*100./IN[k] << "\t[%]" << std::endl;
		}
		std::cout << "############################\n" << std::endl;
	}
}

