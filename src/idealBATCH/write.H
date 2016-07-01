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

//- On Screen
{
    unsigned int N = eq->getTime().size();
    OpenSMOKE::OpenSMOKEVectorDouble omega0(NC);
    OpenSMOKE::OpenSMOKEVectorDouble omegaf(NC);

    for (unsigned int k=0;k<NC;k++)
    {
        omega0[k+1] = eq->getSpecie()[k];
    }

    for (unsigned int k=0;k<NC;k++)
    {
        omegaf[k+1] = eq->getSpecie()[NC*(N-1) + k];
    }

    double h0 = 0.;
    {
        OpenSMOKE::OpenSMOKEVectorDouble h(NC);
        thermodynamicsMapXML->SetPressure(eq->getPressure()[0]);
        thermodynamicsMapXML->SetTemperature(eq->getTemperature()[0]);
        thermodynamicsMapXML->hMolar_Species(h);
        for (unsigned int j=1;j<=NC;j++)
        {
            h0 += omega0[j]*h[j]/thermodynamicsMapXML->MW()[j];
        }
    }

    double hf = 0.;
    {
        OpenSMOKE::OpenSMOKEVectorDouble h(NC);
        thermodynamicsMapXML->SetPressure(eq->getPressure()[N-1]);
        thermodynamicsMapXML->SetTemperature(eq->getTemperature()[N-1]);
        thermodynamicsMapXML->hMolar_Species(h);
        for (unsigned int j=1;j<=NC;j++)
        {
            hf += omegaf[j]*h[j]/thermodynamicsMapXML->MW()[j];
        }
    }

    std::cout << "\n#################################################" << std::endl;
    std::cout << "Enthalpy flux IN         = " << fabs(h0)                << "    [W/m2]" << std::endl;
    std::cout << "Enthalpy flux OUT        = " << fabs(hf)                << "    [W/m2]\n" << std::endl;
    std::cout << "ERROR on ENERGY BALANCE  = " << fabs((h0 - hf)*100./h0) << "    [%]" << std::endl;
    std::cout << "#################################################" << std::endl;

}


{
    boost::filesystem::remove_all("results");
    boost::filesystem::create_directory("results");

    //- Mass
    {
        std::ofstream output;
        output.setf(std::ios::scientific);
        output.precision(6);
        output.open("results/mass.txt",std::ios::out);
        output << "time\tmass" << std::endl;
        
        for (unsigned int i=0;i<eq->getTime().size();i++)
        {
            output << eq->getTime()[i] << "\t" << eq->getMass()[i] << std::endl;
        }
        output.close();
    }

    //- Volume
    {
        std::ofstream output;
        output.setf(std::ios::scientific);
        output.precision(6);
        output.open("results/volume.txt",std::ios::out);
        output << "time\tvolume" << std::endl;

        for (unsigned int i=0;i<eq->getTime().size();i++)
        {
            output << eq->getTime()[i] << "\t" << eq->getVolume()[i] << std::endl;
        }
        output.close();
    }

    //- Pressure
    {
        std::ofstream output;
        output.setf(std::ios::scientific);
        output.precision(6);
        output.open("results/pressure.txt",std::ios::out);
        output << "time\tpressure" << std::endl;
        
        for (unsigned int i=0;i<eq->getTime().size();i++)
        {
            output << eq->getTime()[i] << "\t" << eq->getPressure()[i] << std::endl;
        }
        output.close();
    }

    //- Temperature
    {
        std::ofstream output;
        output.setf(std::ios::scientific);
        output.precision(6);
        output.open("results/temperature.txt",std::ios::out);
        output << "time\ttemperature" << std::endl;
        
        for (unsigned int i=0;i<eq->getTime().size();i++)
        {
            output << eq->getTime()[i] << "\t" << eq->getTemperature()[i] << std::endl;
        }
        output.close();
    }

    //- Specie
    {
        std::ofstream output;
        output.setf(std::ios::scientific);
        output.precision(6);
        output.open("results/specie.txt",std::ios::out);
        output << "time";

        for (unsigned int i=0;i<NC;i++)
        {
            output << "\t" << thermodynamicsMapXML->NamesOfSpecies()[i];
        }
        output << std::endl;

        unsigned int h = 0;
        for (unsigned int i=0;i<eq->getTime().size();i++)
        {
            output << eq->getTime()[i];
            for (unsigned int k=0;k<NC;k++)
            {
                output << "\t" << eq->getSpecie()[h++];
            }
            output << std::endl;
        }
        output.close();
    }

    //- Site
    {
        std::ofstream output;
        output.setf(std::ios::scientific);
        output.precision(6);
        output.open("results/site.txt",std::ios::out);
        output << "time";

        for (unsigned int i=0;i<SURF_NC;i++)
        {
            output << "\t" << thermodynamicsSurfaceMapXML->NamesOfSpecies()[i+NC];
        }
        output << std::endl;

        unsigned int h = 0;
        for (unsigned int i=0;i<eq->getTime().size();i++)
        {
            output << eq->getTime()[i];
            for (unsigned int k=0;k<SURF_NC;k++)
            {
                output << "\t" << eq->getSite()[h++];
            }
            output << std::endl;
        }
        output.close();
    }
}
