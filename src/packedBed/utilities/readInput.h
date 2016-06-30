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
    class readInput 
    {
        public:

            readInput(std::string& file);

            inline std::vector<bool>        getModels()                     const {return models_;};

            inline std::vector<double>      getMW()                         const {return MW_;};
            inline std::vector<double>      getMassFraction()               const {return omega_;};

            inline std::vector<std::string> getSpecieName()                 const {return name_;};
            inline std::vector<std::string> getModelName()                  const {return modelsName_;};
            
            inline unsigned int             numberOfSpecies()               const {return NC_;};

            inline bool                     externalExchange()              const {return extHeat_;};
            inline bool                     energy()                        const {return energy_;};

            inline std::string              getReactionType()               const {return reactionType_;};
            inline std::string              getSolver()                     const {return solver_;};

            inline double                   getCoolantTemperature()         const {return Tw_;};
            inline double                   getFeedTemperature()            const {return Tin_;};
            inline double                   getPressure()                   const {return p_;};
            inline double                   getSpecificMassFlowRate()       const {return G_;};
            inline double                   getTubeDiameter()               const {return Dt_;};
            inline double                   getParticleDiameter()           const {return Dp_;};
            inline double                   getReactorLength()              const {return L_;};
            inline double                   getGasSpecificHeat()            const {return cpG_;};
            inline double                   getGasViscosity()               const {return mu_;};
            inline double                   getGasDiffusivity()             const {return diff_;};
            inline double                   getGasConductivity()            const {return kG_;};
            inline double                   getCatalystSpecificHeat()       const {return cpC_;};
            inline double                   getCatalystDensity()            const {return rhoC_;};
            inline double                   getCatalystVoidFraction()       const {return epsiC_;};
            inline double                   getCatalystTortuosity()         const {return tauC_;};
            inline double                   getCatalystConductivity()       const {return kC_;};

            void recapOnScreen();

        private:

            std::string file_;
            std::string reactionType_;
            std::string solver_;

            std::vector<bool> models_;

            std::vector<double> MW_;
            std::vector<double> omega_;

            std::vector<std::string> name_;
            std::vector<std::string> modelsName_;

            unsigned int NC_;

            double Tw_;
            double Tin_;
            double p_;
            double G_;
            double Dp_;
            double Dt_;
            double L_;
            double cpG_;
            double mu_;
            double kG_;
            double diff_;
            double kC_;
            double cpC_;
            double rhoC_;
            double kS_;
            double tauC_;
            double epsiC_;

            bool   extHeat_;
            bool   energy_;

            void error() { remove("BzzFile.txt"); std::cout << "\nASALI::readInput::ERROR\n" << std::endl;};
            
            bool fileExist(const std::string &file);
            
            std::string boolOnScreen(bool value);
    };
    
    readInput::readInput(std::string& file):
    file_(file)
    {

        if ( !fileExist(file_) )
        {
            error();
            std::cout << file_ << " doesn't exist!" << std::endl;
            exit(EXIT_FAILURE);
        }

        boost::property_tree::ptree tree;
        
        try 
        {
            boost::property_tree::read_xml(file, tree);
        }
        catch (std::exception &e)
        {
            error();
            std::cout << file_ << " isn't a correct XML file!" << std::endl;
            exit(EXIT_FAILURE);
        }
        
        try
        {
            // models
            {
                models_.resize(2);
                models_[0] = tree.get<bool>("models.heterogeneous");
                models_[1] = tree.get<bool>("models.pseudoHomogeneous");
                
                modelsName_.resize(2);
                modelsName_[0] = "heterogeneous";
                modelsName_[1] = "pseudoHomogeneous";
            }

            //reaction
            {
                reactionType_ = tree.get<std::string>("reaction");

                if ( reactionType_ == "O-xylene-to-phthalic" )
                {
                    NC_    = 5;

                    MW_.resize(NC_);
                    MW_[0] = 32.;
                    MW_[1] = 106.16;
                    MW_[2] = 148.12;
                    MW_[3] = 18.;
                    MW_[4] = 28.;

                    name_.resize(NC_);
                    name_[0] = "O2";
                    name_[1] = "XYLENE";
                    name_[2] = "PHTHALIC";
                    name_[3] = "H2O";
                    name_[4] = "N2";
                }
                else if ( reactionType_ == "O-xylene-to-phthalic-complex" )
                {
                    NC_    = 6;

                    MW_.resize(NC_);
                    MW_[0] = 32.;
                    MW_[1] = 106.16;
                    MW_[2] = 148.12;
                    MW_[3] = 44.;
                    MW_[4] = 28.;
                    MW_[5] = 18.;

                    name_.resize(NC_);
                    name_[0] = "O2";
                    name_[1] = "XYLENE";
                    name_[2] = "PHTHALIC";
                    name_[3] = "CO2";
                    name_[4] = "N2";
                    name_[5] = "H2O";
                }
                else
                {
                    error();
                    std::cout << "node || reaction || could be only || O-xylene-to-phthalic || O-xylene-to-phthalic-complex ||\n" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }

            // mole
            {
                omega_.resize(NC_);

                std::vector<double> x(NC_);
                double sum = 0;
                if ( reactionType_ == "O-xylene-to-phthalic" )
                {
                     for (unsigned int i=0;i<NC_;i++)
                    {
                       x[i] = 0.;
                    }

                    x[0]     = tree.get<double>("mole.O2");
                    x[1]     = tree.get<double>("mole.XYLENE");
                    x[4]     = tree.get<double>("mole.N2");
                }
                else if ( reactionType_ == "O-xylene-to-phthalic-complex" )
                {
                    for (unsigned int i=0;i<NC_;i++)
                    {
                       x[i] = 0.;
                    }

                    x[0]     = tree.get<double>("mole.O2");
                    x[1]     = tree.get<double>("mole.XYLENE");
                    x[4]     = tree.get<double>("mole.N2");
                }

                for (unsigned int i=0;i<NC_;i++)
                {
                   sum      = sum + x[i];
                }

                if ( sum != 1. )
                {
                    error();
                    std::cout << "sum of inlet mole fraction MUST BE 1\n" << std::endl;
                    exit(EXIT_FAILURE);
                }

                double MWmix = 0.;
                for (unsigned int i=0;i<NC_;i++)
                {
                    MWmix = MWmix + MW_[i]*x[i];
                }
                
                for (unsigned int i=0;i<NC_;i++)
                {
                    omega_[i] = MW_[i]*x[i]/MWmix;
                }
            }

            // operating conditions
            {
                Tw_       = tree.get<double>("operatingConditions.temperauture.coolant");
                Tin_      = tree.get<double>("operatingConditions.temperauture.feed");
                p_        = tree.get<double>("operatingConditions.pressure");
                G_        = tree.get<double>("operatingConditions.specificMassFlowRate");
                extHeat_  = tree.get<bool>("operatingConditions.externalExchange");
                energy_   = tree.get<bool>("operatingConditions.energy");  
            }
            
            // reactor dimensions
            {
                Dp_     = tree.get<double>("reactorDimensions.particleDiameter"); 
                Dt_     = tree.get<double>("reactorDimensions.tubeDiameter"); 
                L_      = tree.get<double>("reactorDimensions.length");
            }

            // gas properties
            {
                cpG_  = tree.get<double>("gasProperties.specificHeat");  
                mu_   = tree.get<double>("gasProperties.viscosity");
                diff_ = tree.get<double>("gasProperties.diffusivity");
                kG_   = tree.get<double>("gasProperties.conductivity");
            }

            // catalyst properties
            {
                cpC_   = tree.get<double>("catalystProperties.specificHeat");  
                rhoC_  = tree.get<double>("catalystProperties.density");
                kC_    = tree.get<double>("catalystProperties.conductivity");
                epsiC_ = tree.get<double>("catalystProperties.voidFraction");
                tauC_  = tree.get<double>("catalystProperties.tortuosity");
            }


            // solver
            {
                solver_ = tree.get<std::string>("solver");

                if (solver_ == "BzzMath" )
                {
                    #if ASALI_USE_BZZ == 0
                        error();
                        std::cout << "node || solver || cannot be || BzzMath || \n" << std::endl;
                        exit(EXIT_FAILURE);
                    #endif
                }
                else if ( solver_ == "OpenSMOKE" )
                {}
                else
                {
                    error();
                    std::cout << "node || solver || could be || OpenSMOKE ||";
                    #if ASALI_USE_BZZ == 1
                        std::cout << " BzzMath ||";
                    #endif
                    std::cout << "\n" << std::endl;
                    exit(EXIT_FAILURE);
                }

            }

        }
        catch (std::exception &e)
        {
            error();
            std::cout << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    bool readInput::fileExist(const std::string& file)
    {
        if ( !boost::filesystem::exists(file_) || boost::filesystem::is_directory(file_) )
        {
            return false;
        }
        return true;
    }

    void readInput::recapOnScreen()
    {
        std::cout.precision(6);
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                       PACKED BED REACTOR                                       \n" << std::endl;
        std::cout << "Tube     diameter                         = " << Dt_ << "\t[m]" << std::endl;
        std::cout << "Reactor  length                           = " << L_ << "\t[m]" << std::endl;
        std::cout << "Particle diameter                         = " << Dp_ << "\t[m]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                      CATALYST PROPERTIES                                       \n" << std::endl;
        std::cout << "Density                                  = " << rhoC_ << "\t[Kg/m3]" << std::endl;
        std::cout << "Specific heat                            = " << cpC_ << "\t[J/Kg/K]" << std::endl;
        std::cout << "Conductivity                             = " << kC_ << "\t[W/m/K]" << std::endl;
        std::cout << "Tortuosity                               = " << tauC_ << "\t[-]" << std::endl;
        std::cout << "Void fraction                            = " << epsiC_ << "\t[-]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                         GAS PROPERTIES                                         \n" << std::endl;
        std::cout << "Viscosity                                = " << mu_ << "\t[Pas]" << std::endl;
        std::cout << "Specific heat                            = " << cpG_ << "\t[J/Kg/K]" << std::endl;
        std::cout << "Conductivity                             = " << kG_ << "\t[W/m/K]" << std::endl;
        std::cout << "Diffusivity                              = " << diff_ << "\t[m2/s]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                      OPERATING CONDITIONS                                      \n" << std::endl;
        std::cout << "Specific mass flow rate                  = " << G_ << "\t[Kg/m2/s]" << std::endl;
        std::cout << "Pressure                                 = " << p_ << "\t[Pa]" << std::endl;
        std::cout << "Feed  temperature                        = " << Tin_ << "\t[K]" << std::endl;
        std::cout << "Shell temperature                        = " << Tw_ << "\t[K]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                         SOLVER OPTIONS                                         \n" << std::endl;
        std::cout << "External heat exchange:                    " << boolOnScreen(extHeat_) << std::endl;
        std::cout << "Chemistry scheme:                          " << reactionType_ << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                        NUMERICAL SOLVER                                        \n" << std::endl;
        std::cout << "Solver compiled with                        ";
        #if ASALI_USE_BZZ == 1
        std::cout << "|| BzzMath ";
        #endif
        std::cout << "|| OpenSMOKE ";
        std::cout << "||" << std::endl;
        std::cout << "Chosen solver:                              || " << solver_ << " ||" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
    }

    std::string readInput::boolOnScreen(bool value)
    {
        std::string value_;
        if ( value == true)
            value_ = "on";
        else
            value_ = "off";
        
        return value_;
    }

}
