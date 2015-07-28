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
    std::vector<double> T;
    std::vector<double> z;
    ode->GetProfile(z,p,T);
    boost::filesystem::remove_all("results");
    boost::filesystem::create_directory("results");
    std::ofstream output;
    output.setf(std::ios::scientific);
    output.precision(6);
    output.open("results/profile.txt",std::ios::out);
    output << "#Length" << "\t" << "#Pressure" << "\t" << "#Temperature" << std::endl;
    for (unsigned int k=0;k<p.size();k++)
        output << z[k] << "\t" << p[k] << "\t" << T[k] << std::endl;
    output.close();

    std::ofstream recap;
    recap.setf(std::ios::scientific);
    recap.precision(6);
    recap.open("results/information.txt",std::ios::out);
    recap << "Input:\n" << std::endl;
    recap << " 1/ specie name:               " << input.getSpecie() << std::endl;
    recap << " 2/ pressure correlation:      " << input.getPressureCorrelation() << std::endl;
    recap << " 3/ heat transfer correlation: " << input.getHeatTransferCorrelation() << std::endl;
    recap << " 4/ feed temperature:          " << input.getFeedTemperature() << " [K]" << std::endl;
    recap << " 5/ wall temperature:          " << input.getSolidTemperature() << " [K]" << std::endl;
    recap << " 6/ feed pressure:             " << input.getPressureCorrelation() << " [Pa]" << std::endl;
    recap << " 7/ feed velocity:             " << input.getFeedVelocity() << " [m/s]" << std::endl;

    if ( input.getReactorType() == "Monolith" )
    {
        recap << " 8/ void fraction:             " << input.getVoidFraction() << " [-]" << std::endl;
    }
    else if ( input.getReactorType() == "PackedBed" )
    {
        recap << " 8/ particle diameter:         " << input.getParticleDiameter() << " [m]" << std::endl;
    }

    recap << " 9/ tube diameter:             " << input.getShellDiameter() << " [m]" << std::endl;
    recap << "10/ tube length:               " << input.getReactorLength()  << " [m]" << std::endl;
    recap << std::endl;
    recap << "Output:\n" << std::endl;
    recap << " 1/ out temperature:           " << yF[2] << " [K]" << std::endl;
    recap << " 2/ out pressure:              " << yF[1] << " [Pa]" << std::endl;
    recap << " 3/ pressure drops:            " << input.getFeedPressure() - yF[1] << " [Pa]" << std::endl;
    recap << " 4/ pressure drops:            " << (input.getFeedPressure() - yF[1])*1e-05/input.getReactorLength() << " [bar/m]" << std::endl;
    recap.close();
}

