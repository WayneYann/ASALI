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

ASALI::check(argc);
std::string* file;
std::string* options;

{
	if ( argc == 2 )
	{
		file    = new std::string(argv[1]);
		options = new std::string("none");
	}
	else if ( argc == 3 )
	{
		file    = new std::string(argv[1]);
		options = new std::string(argv[2]);
	}
}

ASALI::READinput  input(*file,*options);
ASALI::ODESystem *ode;
ASALI::BVPSystem *bvp;
ASALI::GROWgrid  *growGrid;

{
	std::string folderPath = input.getKineticsPath() + "/";
	boost::filesystem::path path_kinetics(folderPath);
	{
		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics / "kinetics.xml");

		thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>(doc); 
		transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>(doc); 
		kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN<double>(*thermodynamicsMapXML, doc);
	}
	{
		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics / "kinetics.surface.xml");

		thermodynamicsSurfaceMapXML = new OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>(doc); 
		kineticsSurfaceMapXML = new OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped<double>(*thermodynamicsSurfaceMapXML, doc);
	}

	ASALI::logo(input.getStart());

	if ( input.getStart() == "new" || input.getStart() == "latest" )
	{
		input.recapOnScreen();
	}

	{
		unsigned int WP = 6;
		std::cout.setf(std::ios::scientific);
		std::cout.precision(WP);
	}

	grow = input.getGridType();
}

{
	if ( input.getStart() == "new")
	{
		std::string folderName = "results/";
		boost::filesystem::remove_all(folderName);
		boost::filesystem::create_directory(folderName);

		N       = input.getMinPointsNumber();
		NC      = thermodynamicsMapXML->NumberOfSpecies();
		SURF_NC = thermodynamicsSurfaceMapXML->number_of_site_species();
		NE      = (NC + NC + SURF_NC + 1 + 1)*N;

		ChangeDimensions(N,       &z,      true);
		ChangeDimensions(NC,      &x0bulk, true);
		ChangeDimensions(NC,      &x0wall, true);
		ChangeDimensions(SURF_NC, &x0site, true);

		if ( input.getEnergy() == false )
			T0 = input.getSolidTemperature();
		else
			T0 = input.getGasTemperature();

		thermodynamicsMapXML->SetTemperature(T0);
		thermodynamicsMapXML->SetPressure(input.getPressure());

		x0site[1] = 1.;
		for (unsigned int k=1;k<=SURF_NC;k++)
			x0site[k]=std::max(x0site[k],0.);

		std::vector<double> InletCOMPOSITION(input.getFraction().size());
		double sum = 0.;
		for(unsigned int k=0; k<input.getFraction().size(); k++)
			sum = sum + input.getFraction()[k];

		if ( sum != 1. )
		{
			std::cout << "\nASALI::WARNING Sum of inlet composition is not 1!\n" << std::endl;
			std::cout << "\nPress Enter to normalized them..." << std::endl;
			getchar();
			for(unsigned int k=0; k<input.getFraction().size(); k++)
				InletCOMPOSITION[k]=input.getFraction()[k]/sum;
		}
		else
		{
			for(unsigned int k=0; k<input.getFraction().size(); k++)
				InletCOMPOSITION[k]=input.getFraction()[k];
		}

		double MW;
		if ( input.getFeed() == "Mole" )
		{
			OpenSMOKE::OpenSMOKEVectorDouble x(NC);
			for(unsigned int k=0;k<input.getFractionName().size();k++)
				x[thermodynamicsMapXML->IndexOfSpecies(input.getFractionName()[k])] = InletCOMPOSITION[k];

			thermodynamicsMapXML->MassFractions_From_MoleFractions(x0bulk,MW,x);
		}
		else if ( input.getFeed() == "Mass" )
		{
			for(unsigned int k=0;k<input.getFractionName().size();k++)
				x0bulk[thermodynamicsMapXML->IndexOfSpecies(input.getFractionName()[k])] = InletCOMPOSITION[k];

			thermodynamicsMapXML->MolecularWeight_From_MassFractions(MW,x0bulk);
		}

		for (unsigned int k=1;k<=NC;k++)
			x0wall[k] = 0.;

		x0wall[thermodynamicsMapXML->IndexOfSpecies(input.getInert())] = 1.;

		G = input.getVelocity()*input.getPressure()*MW/(8314.*T0);

		double step = input.getReactorLength()/(N-1);
		z[1] = 0.;
		for (unsigned int k=2;k<=z.Size();k++)
			z[k] = z[k-1] + step;

		ChangeDimensions(NE, &x0, true);
		unsigned int counter = 1;
		for (unsigned int k=1;k<=N;k++)
		{
			for (unsigned int j=1;j<=NC;j++)
			{
				if (k==1)
					x0[counter++] = x0bulk[j];
				else
					x0[counter++] = x0bulk[j];
			}
			for (unsigned int j=1;j<=NC;j++)
				x0[counter++] = x0bulk[j];
			for (unsigned int j=1;j<=SURF_NC;j++)
				x0[counter++] = x0site[j];
			x0[counter++] = T0;
			x0[counter++] = input.getSolidTemperature();
		}
	}
	else if ( input.getStart() == "latest" || input.getStart() == "converter" || input.getStart() == "sampling")
	{
		N       = input.getRestartPoints();
		NC      = thermodynamicsMapXML->NumberOfSpecies();
		SURF_NC = thermodynamicsSurfaceMapXML->number_of_site_species();
		NE      = (NC + NC + SURF_NC + 1 + 1)*N;

		ChangeDimensions(N,       &z,      true);
		ChangeDimensions(NC,      &x0bulk, true);
		ChangeDimensions(NC,      &x0wall, true);
		ChangeDimensions(SURF_NC, &x0site, true);

		OpenSMOKE::OpenSMOKEVectorDouble Tbulk(N);
		OpenSMOKE::OpenSMOKEVectorDouble Twall(N);

		OpenSMOKE::OpenSMOKEVectorDouble* omegaBulk;
		OpenSMOKE::OpenSMOKEVectorDouble* omegaWall;
		OpenSMOKE::OpenSMOKEVectorDouble* teta;

		omegaBulk = new OpenSMOKE::OpenSMOKEVectorDouble[N];
		omegaWall = new OpenSMOKE::OpenSMOKEVectorDouble[N];
		teta      = new OpenSMOKE::OpenSMOKEVectorDouble[N];

		for (unsigned int i=0;i<N;i++)
		{
			ChangeDimensions(NC,      &omegaBulk[i], true);
			ChangeDimensions(NC,      &omegaWall[i], true);
			ChangeDimensions(SURF_NC, &teta[i],      true);
		}

		unsigned int counter = 0;
		for (unsigned int k=0;k<N;k++)
		{
			for (unsigned int j=1;j<=NC;j++)
				omegaBulk[k][j] = input.getRestartBulk()[counter++];
		}
		counter = 0;
		for (unsigned int k=0;k<N;k++)
		{
			for (unsigned int j=1;j<=NC;j++)
				omegaWall[k][j] = input.getRestartWall()[counter++];
		}
		counter = 0;
		for (unsigned int k=0;k<N;k++)
		{
			for (unsigned int j=1;j<=SURF_NC;j++)
				teta[k][j] = input.getRestartSite()[counter++];
		}

		if ( input.getEnergy() == false )
			T0 = input.getSolidTemperature();
		else
			T0 = input.getGasTemperature();

		thermodynamicsMapXML->SetTemperature(T0);
		thermodynamicsMapXML->SetPressure(input.getPressure());

		if ( input.getRestartResults() == "mole" )
		{
			double dummy;
			OpenSMOKE::OpenSMOKEVectorDouble x(NC);
			OpenSMOKE::OpenSMOKEVectorDouble y(NC);
			for (unsigned int i=0;i<N;i++)
			{
				thermodynamicsMapXML->MassFractions_From_MoleFractions(x,dummy,omegaBulk[i]);
				for (unsigned int j=1;j<=NC;j++)
					omegaBulk[i][j] = x[j];

				thermodynamicsMapXML->MassFractions_From_MoleFractions(x,dummy,omegaWall[i]);
				for (unsigned int j=1;j<=NC;j++)
					omegaWall[i][j] = x[j];
			}
		}

		for (unsigned int k=1;k<=N;k++)
		{
			Tbulk[k] = input.getRestartGasTemperature()[k-1];
			Twall[k] = input.getRestartSolidTemperature()[k-1];
		}

		ChangeDimensions(NE, &yF, true);
		counter = 1;
		for (unsigned int k=0;k<N;k++)
		{
			for (unsigned int j=1;j<=NC;j++)
				yF[counter++] = omegaBulk[k][j];
			for (unsigned int j=1;j<=NC;j++)
				yF[counter++] = omegaWall[k][j];
			for (unsigned int j=1;j<=SURF_NC;j++)
				yF[counter++] = teta[k][j];
			yF[counter++] = Tbulk[k+1];
			yF[counter++] = Twall[k+1];
		}

		x0site[1] = 1.;
		for (unsigned int k=1;k<=SURF_NC;k++)
			x0site[k]=std::max(x0site[k],0.);

		std::vector<double> InletCOMPOSITION(input.getFraction().size());
		double sum = 0.;
		for(unsigned int k=0; k<input.getFraction().size(); k++)
			sum = sum + input.getFraction()[k];

		if ( sum != 1. )
		{
			std::cout << "\nASALI::WARNING Sum of inlet composition is not 1!\n" << std::endl;
			std::cout << "\nPress Enter to normalized them..." << std::endl;
			getchar();
			for(unsigned int k=0; k<input.getFraction().size(); k++)
				InletCOMPOSITION[k]=input.getFraction()[k]/sum;
		}
		else
		{
			for(unsigned int k=0; k<input.getFraction().size(); k++)
				InletCOMPOSITION[k]=input.getFraction()[k];
		}

		double MW;
		if ( input.getFeed() == "Mole" )
		{
			OpenSMOKE::OpenSMOKEVectorDouble x(NC);
			for(unsigned int k=0;k<input.getFractionName().size();k++)
				x[thermodynamicsMapXML->IndexOfSpecies(input.getFractionName()[k])] = InletCOMPOSITION[k];

			thermodynamicsMapXML->MassFractions_From_MoleFractions(x0bulk,MW,x);
		}
		else if ( input.getFeed() == "Mass" )
		{
			for(unsigned int k=0;k<input.getFractionName().size();k++)
				x0bulk[thermodynamicsMapXML->IndexOfSpecies(input.getFractionName()[k])] = InletCOMPOSITION[k];

			thermodynamicsMapXML->MolecularWeight_From_MassFractions(MW,x0bulk);
		}

		for (unsigned int k=1;k<=NC;k++)
			x0wall[k] = 0.;

		x0wall[thermodynamicsMapXML->IndexOfSpecies(input.getInert())] = 1.;

		G = input.getVelocity()*input.getPressure()*MW/(8314.*T0);

		for (unsigned int k=1;k<=z.Size();k++)
			z[k] = input.getRestartGrid()[k-1];

		delete [] omegaBulk;
		delete [] omegaWall;
		delete [] teta;
	}
}

if ( input.getStart() == "new" || input.getStart() == "latest" )
{
	ode = new ASALI::ODESystem(*thermodynamicsMapXML,        *kineticsMapXML,
							   *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML,
							   *transportMapXML);

	ode->setHomogeneusReactions(input.getHomogeneousReactions());
	ode->setHeterogeneusReactions(input.getHeterogenousReactions());
	ode->setEnergyEquation(input.getEnergy());
	ode->setReactorGeometry(input.getAlfa(),             input.getVoidFraction(), input.getSherwood(),
							 input.getCatalyticLength(),  input.getInertLength(),
							 input.getChannelDiameter(),  G);
	ode->setFeedValue(input.getPressure(),T0,x0bulk);
	ode->setInert(input.getInert());

	bvp = new ASALI::BVPSystem(*thermodynamicsMapXML,        *kineticsMapXML,
							   *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML,
							   *transportMapXML);

	bvp->setSolid(input.getSolidDensity(),input.getSolidConductivity(),input.getSolidSpecificHeat());
	bvp->setHomogeneusReactions(input.getHomogeneousReactions());
	bvp->setHeterogeneusReactions(input.getHeterogenousReactions());
	bvp->setEnergyEquation(input.getEnergy());
	bvp->setReactorGeometry(input.getAlfa(),             input.getVoidFraction(), input.getSherwood(),
							input.getCatalyticLength(),  input.getInertLength(),
							input.getChannelDiameter(),  G);
	bvp->setFeedValue(input.getPressure(),T0,x0bulk);
	bvp->setInert(input.getInert());

	growGrid = new ASALI::GROWgrid(NC, SURF_NC);

	growGrid->setAddPoints(input.getAddPointsNumber());
	growGrid->setMaxPoints(input.getMaxPointsNumber());
	growGrid->setTol(input.getGridError());
}

if ( *options == "none" && input.getStart() == "sampling" )
{
	ASALI::error();
	std::cout << " || sampling || options required a file named || samplingGrid.txt ||" << std::endl;
	std::cout << "\n " << std::endl;
	std::cout << "INPUT ORDER SHOULD BE:  " << std::endl;
	std::cout << "\n " << std::endl;
	std::cout << " 1/ input file " << std::endl;
	std::cout << " 2/ samplingGrid.txt " << std::endl;
	std::cout << "\n " << std::endl;
	exit(EXIT_FAILURE);
}
