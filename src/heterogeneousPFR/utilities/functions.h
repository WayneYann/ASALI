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

namespace ASALI
{
    void logo(const std::string type)
    {
        std::cout << "\033[2J\033[1;1H" << std::endl;
        std::cout << "################################################################################################" << std::endl;
        std::cout << "#                                                                                              #" << std::endl;
        std::cout << "#     #############       #############       #############       ####                ####     #" << std::endl;
        std::cout << "#    #             #     #             #     #             #     #    #              #    #    #" << std::endl;
        std::cout << "#    #     ###     #     #    #########      #     ###     #     #    #              #    #    #" << std::endl;
        std::cout << "#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #" << std::endl;
        std::cout << "#    #     ###     #     #    #              #     ###     #     #    #              #    #    #" << std::endl;
        std::cout << "#    #             #     #    #########      #             #     #    #              #    #    #" << std::endl;
        std::cout << "#    #             #     #             #     #             #     #    #              #    #    #" << std::endl;
        std::cout << "#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #" << std::endl;
        std::cout << "#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #" << std::endl;
        std::cout << "#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #" << std::endl;
        std::cout << "#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #" << std::endl;
        std::cout << "#     ####     ####       #############       ####     ####       #############       ####     #" << std::endl;
        std::cout << "#                                                                                              #" << std::endl;
        std::cout << "#   Department of Energy                                                                       #" << std::endl;
        std::cout << "#   Politecnico di Milano                                                                      #" << std::endl;
        std::cout << "#   Author: Stefano Rebughini <stefano.rebughini@polimi.it>                                    #" << std::endl;
        std::cout << "#                                                                                              #" << std::endl;
        std::cout << "################################################################################################" << std::endl;
        
        if ( type == "converter" )
        {
            std::cout << "\n############################### " << std::endl;
            std::cout << "#          CONVERTER          # " << std::endl;
            std::cout << "#       mole <---> mass       # " << std::endl;
            std::cout << "###############################\n " << std::endl;
        }
        else if ( type == "sampling" )
        {
            std::cout << "\n################################" << std::endl;
            std::cout << "#           SAMPLING           # " << std::endl;
            std::cout << "################################\n " << std::endl;
        }
        else if ( type == "kinetic" )
        {
            std::cout << "\n################################" << std::endl;
            std::cout << "#       KINETIC analysis       # " << std::endl;
            std::cout << "################################\n " << std::endl;
        }
        else if ( type == "heat" )
        {
            std::cout << "\n#################################" << std::endl;
            std::cout << "#  HEAT of REACTION evaluation  # " << std::endl;
            std::cout << "#################################\n " << std::endl;
        }
        else if ( type == "RPA" )
        {
            std::cout << "\n################################" << std::endl;
            std::cout << "#    REACTION PATH ANALYSIS    # " << std::endl;
            std::cout << "#                              # " << std::endl;
            std::cout << "#                              # " << std::endl;
            std::cout << "#            WARNING           # " << std::endl;
            std::cout << "#     This utility has been    #" << std::endl;
            std::cout << "#   created only for kinetic   #" << std::endl;
            std::cout << "# schemes where every reaction #" << std::endl;
            std::cout << "#  is followed by its inverse  #" << std::endl;
            std::cout << "#                              # " << std::endl;
            std::cout << "#           A ---> B           #" << std::endl;
            std::cout << "#           B ---> A           #" << std::endl;
            std::cout << "#                              # " << std::endl;
            std::cout << "################################\n " << std::endl;
        }
        else if ( type == "conversion" )
        {
            std::cout << "\n################################" << std::endl;
            std::cout << "#          CONVERSION          # " << std::endl;
            std::cout << "################################\n " << std::endl;
        }
    }

    void check( const int argc )
    {
        if ( argc != 3 && argc != 2)
        {
            std::cout << "\n " << std::endl;
            std::cout << "INPUT ORDER SHOULD BE:  " << std::endl;
            std::cout << "\n " << std::endl;
            std::cout << " 1/ input file " << std::endl;
            std::cout << " 2/ [options] " << std::endl;
            std::cout << "\n " << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    
    void discretizationSchemeOnScreen(const std::string scheme)
    {
        std::cout << "\n#################################" << std::endl;
        std::cout << "# Discretization scheme:    " << scheme << " #" << std::endl;
        std::cout << "#################################" << std::endl;
    }

    template < typename T > std::string to_string( const T& v )
    {
        std::ostringstream stm ;
        return ( stm << v ) ? stm.str() : "{*** error ***}" ;
    }

    void CPUtime(const double tStart, const double tEnd)
    {
        std::cout.setf( std::ios::fixed, std:: ios::floatfield );
        std::cout.precision(6);
        std::cout << "\n###################### " << std::endl;
        std::cout << "#  Simulation time: " << std::endl;
        std::cout << "#  " << (tEnd - tStart)       << " [s] " << std::endl;
        std::cout << "#  " << (tEnd - tStart)/60.   << " [min] " << std::endl;
        std::cout << "#  " << (tEnd - tStart)/3600. << " [h] " << std::endl;
        std::cout << "######################\n " << std::endl;
        
        std::ofstream timeFile;
        std::string timeName = "results/time.txt";
        const char *pathTime = timeName.c_str();
        timeFile.open(pathTime,std::ios::out);
        timeFile.setf( std::ios::fixed, std:: ios::floatfield );
        timeFile.precision(6);
        timeFile << "Simulation time: " << std::endl;
        timeFile << (tEnd - tStart)       << " [s] " << std::endl;
        timeFile << (tEnd - tStart)/60.   << " [min] " << std::endl;
        timeFile << (tEnd - tStart)/3600. << " [h] " << std::endl;
        timeFile.close();
    }
    
    void error()
    {
        std::cout << "\nASALI::READinput::ERROR\n" << std::endl;
    }

    void ODEstart()
    {
        std::cout << "\n######################################" << std::endl;
        std::cout << "# ODE system:                 START  #" << std::endl;
        std::cout << "######################################\n" << std::endl;
    }
    
    void ODEend()
    {
        std::cout << "\n######################################" << std::endl;
        std::cout << "# ODE system:                 END    #" << std::endl;
        std::cout << "######################################\n" << std::endl;
    }

    void DAEstart()
    {
        std::cout << "\n######################################" << std::endl;
        std::cout << "# DAE system:                 START  #" << std::endl;
        std::cout << "######################################\n" << std::endl;
    }
    
    void DAEend()
    {
        std::cout << "\n######################################" << std::endl;
        std::cout << "# DAE system:                 END    #" << std::endl;
        std::cout << "######################################\n" << std::endl;
    }
}
