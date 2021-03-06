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

            inline bool                     getEnergy()                     const {return energy_;};

            inline std::vector<bool>        getModels()                     const {return models_;};

            inline std::vector<double>      getMW()                         const {return MW_;};
            inline std::vector<double>      getMassFraction()               const {return omega_;};

            inline std::vector<std::string> getSpecieName()                 const {return name_;};
            inline std::vector<std::string> getModelName()                  const {return modelsName_;};
            
            inline unsigned int             numberOfSpecies()               const {return NC_;};
            inline unsigned int             axialPoints()                   const {return Na_;};
            inline unsigned int             radialPoints()                  const {return Nr_;};

            inline std::string              getReactionType()               const {return reactionType_;};
            inline std::string              getHoneyCombType()              const {return typeH_;};
            inline std::string              getSolver()                     const {return solver_;};

            inline double                   getCoolantTemperature()         const {return Tw_;};
            inline double                   getFeedTemperature()            const {return Tin_;};
            inline double                   getPressure()                   const {return p_;};
            inline double                   getSpecificMassFlowRate()       const {return G_;};
            inline double                   getPackedBedTubeDiameter()      const {return DtP_;};
            inline double                   getPackedBedParticleDiameter()  const {return DpP_;};
            inline double                   getPackedBedLength()            const {return LP_;};
            inline double                   getHoneyCombTubeDiameter()      const {return DtH_;};
            inline double                   getHoneyCombCPSI()              const {return CPSIH_;};
            inline double                   getHoneyCombLength()            const {return LH_;};
            inline double                   getHoneyCombWashCoat()          const {return SwH_;};
            inline double                   getHoneyCombWall()              const {return wH_;};
            inline double                   getMicroBedTubeDiameter()       const {return DtM_;};
            inline double                   getMicroBedCPSI()               const {return CPSIM_;};
            inline double                   getMicroBedLength()             const {return LM_;};
            inline double                   getMicroBedParticleDiameter()   const {return DpM_;};
            inline double                   getMicroBedWall()               const {return wM_;};
            inline double                   getGasSpecificHeat()            const {return cpG_;};
            inline double                   getGasViscosity()               const {return mu_;};
            inline double                   getGasDiffusivity()             const {return diff_;};
            inline double                   getGasConductivity()            const {return kG_;};
            inline double                   getCatalystSpecificHeat()       const {return cpC_;};
            inline double                   getCatalystDensity()            const {return rhoC_;};
            inline double                   getCatalystVoidFraction()       const {return epsiC_;};
            inline double                   getCatalystTortuosity()         const {return tauC_;};
            inline double                   getCatalystConductivity()       const {return kC_;};
            inline double                   getSupportSpecificHeat()        const {return cpS_;};
            inline double                   getSupportDensity()             const {return rhoS_;};
            inline double                   getSupportConductivity()        const {return kS_;};
            inline double                   getInertLength()                const {return Linert_;};
            inline double                   getHeat()                       const {return Q_;};

            void recapOnScreen();

        private:

            std::string file_;
            std::string reactionType_;
            std::string typeH_;
            std::string solver_;

            std::vector<bool> models_;

            std::vector<double> MW_;
            std::vector<double> omega_;

            std::vector<std::string> name_;
            std::vector<std::string> modelsName_;

            unsigned int NC_;
            unsigned int Na_;
            unsigned int Nr_;

            bool energy_;

            double Tw_;
            double Tin_;
            double p_;
            double G_;
            double DpP_;
            double DtP_;
            double LP_;
            double DtH_;
            double CPSIH_;
            double SwH_;
            double LH_;
            double wH_;
            double DtM_;
            double DpM_;
            double CPSIM_;
            double wM_;
            double LM_;
            double cpG_;
            double mu_;
            double kG_;
            double diff_;
            double kC_;
            double cpC_;
            double rhoC_;
            double kS_;
            double cpS_;
            double rhoS_;
            double tauC_;
            double epsiC_;
            double Linert_;
            double Q_;

            void error() { std::cout << "\nASALI::readInput::ERROR\n" << std::endl;};
            
            bool fileExist(const std::string &file);
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
                models_.resize(3);
                models_[0] = tree.get<bool>("models.honeyComb");
                models_[1] = tree.get<bool>("models.packedBed");
                models_[2] = tree.get<bool>("models.microPackedbed");
                
                modelsName_.resize(3);
                modelsName_[0] = "honeyComb";
                modelsName_[1] = "packedBed";
                modelsName_[2] = "microBed";
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
                else if ( reactionType_ == "ethylene-partial-oxidation" )
                {
                    NC_    = 6;

                    MW_.resize(NC_);
                    MW_[0] = 12.01*2. + 4.;
                    MW_[1] = 32.;
                    MW_[2] = 12.01*2. + 4. + 16.;
                    MW_[3] = 44.;
                    MW_[4] = 12.01 + 4.;
                    MW_[5] = 18.;

                    name_.resize(NC_);
                    name_[0] = "C2H4";
                    name_[1] = "O2";
                    name_[2] = "C2H4O";
                    name_[3] = "CO2";
                    name_[4] = "CH4";
                    name_[5] = "H2O";
                }
                else if ( reactionType_ == "heat-generation")
                {
                    NC_    = 1;
                    MW_.resize(NC_);
                    MW_[0] = 44.;

                    name_.resize(NC_);
                    name_[0] = "N2";
                }
                else if ( reactionType_ == "O-xylene-to-phthalic-with-tolualdehyde" )
                {
                    NC_    = 8;

                    MW_.resize(NC_);
                    MW_[0] = 32.;      //O2
                    MW_[1] = 106.16;   //C8H10
                    MW_[2] = 148.12;   //C8H4O3
                    MW_[3] = 12.01*8. + 8. + 16.; //C8H8O
                    MW_[4] = 12.01*2. + 6. + 32.; //C8H6O2
                    MW_[5] = 44.;      //CO2
                    MW_[6] = 28.;      //N2
                    MW_[7] = 18.;      //H2O

                    name_.resize(NC_);
                    name_[0] = "O2";
                    name_[1] = "XYLENE";
                    name_[2] = "PHTHALIC";
                    name_[3] = "TOLULDEHYDE";
                    name_[4] = "PHTHALIDE";
                    name_[5] = "CO2";
                    name_[6] = "N2";
                    name_[7] = "H2O";
                }
                else
                {
                    error();
                    std::cout << "node || reaction || could be only || O-xylene-to-phthalic || O-xylene-to-phthalic-complex || O-xylene-to-phthalic-with-tolualdehyde || ethylene-partial-oxidation || heat-generation ||\n" << std::endl;
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
                else if ( reactionType_ == "ethylene-partial-oxidation" )
                {
                    for (unsigned int i=0;i<NC_;i++)
                    {
                       x[i] = 0.;
                    }

                    x[1] = tree.get<double>("mole.O2");
                    x[0] = tree.get<double>("mole.C2H4");
                    x[3] = tree.get<double>("mole.CO2");
                    x[4] = tree.get<double>("mole.CH4");
                }
                else if ( reactionType_ == "O-xylene-to-phthalic-with-tolualdehyde" )
                {
                    for (unsigned int i=0;i<NC_;i++)
                    {
                       x[i] = 0.;
                    }

                    x[0]     = tree.get<double>("mole.O2");
                    x[1]     = tree.get<double>("mole.XYLENE");
                    x[6]     = tree.get<double>("mole.N2");
                }
                else if ( reactionType_ == "heat-generation" )
                {
                    x[0] = 1.;
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
                Tw_     = tree.get<double>("operatingConditions.temperauture.coolant");
                Tin_    = tree.get<double>("operatingConditions.temperauture.feed");
                p_      = tree.get<double>("operatingConditions.pressure");
                G_      = tree.get<double>("operatingConditions.specificMassFlowRate");
                Linert_ = tree.get<double>("operatingConditions.inertLength");
                energy_ = tree.get<bool>("operatingConditions.energy");

                if ( reactionType_ == "heat-generation")
                {
                    Q_  = tree.get<double>("operatingConditions.heatGenerated");
                }
                else
                {
                    Q_  = 0.;
                }
            }
            
            // packed bed
            {
                if (models_[1] == true )
                {
                    DpP_ = tree.get<double>("packedBed.particleDiameter"); 
                    DtP_ = tree.get<double>("packedBed.tubeDiameter"); 
                    LP_  = tree.get<double>("packedBed.length"); 
                }
            }

            // honeycomb
            {
                if ( models_[0] == true)
                {
                    typeH_ = tree.get<std::string>("honeyComb.type"); 
                    CPSIH_ = tree.get<double>("honeyComb.CPSI"); 
                    DtH_   = tree.get<double>("honeyComb.tubeDiameter");
                    LH_    = tree.get<double>("honeyComb.length");
                    wH_    = tree.get<double>("honeyComb.wallThickness");

                    if ( typeH_ == "washcoated" )
                    {
                        SwH_ = tree.get<double>("honeyComb.washcoatThickness"); 
                    }
                    else if ( typeH_ == "extruded" )
                    {
                        SwH_ = 0.;
                    }
                    else
                    {
                        error();
                        std::cout << "node || honeyComb.type || could be || washcoated || extruded ||\n" << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }

            // micro bed
            {
                if ( models_[2] == true)
                {
                    CPSIM_ = tree.get<double>("microBed.CPSI"); 
                    DtM_   = tree.get<double>("microBed.tubeDiameter");
                    LM_    = tree.get<double>("microBed.length");
                    wM_    = tree.get<double>("microBed.wallThickness");
                    DpM_   = tree.get<double>("microBed.particleDiameter");
                }
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

            // support properties
            {
                cpS_  = tree.get<double>("supportProperties.specificHeat");  
                rhoS_ = tree.get<double>("supportProperties.density");
                kS_   = tree.get<double>("supportProperties.conductivity");
            }
            
            // grid
            {
                Na_  = tree.get<unsigned int>("grid.axial");  
                Nr_  = tree.get<unsigned int>("grid.radial");  
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
                    std::cout << "node || solver || could be || OpenSMOKE ||" << std::endl;
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
        if ( models_[0] == true )
        {
            std::cout << "                                        HONEYCOMB REACTOR                                       \n" << std::endl;
            std::cout << "Type:                                       " << typeH_ << std::endl;
            std::cout << "CPSI                                      = " << CPSIH_ << std::endl;
            std::cout << "Tube     diameter                         = " << DtH_ << "\t[m]" << std::endl;
            std::cout << "Reactor  length                           = " << LH_ << "\t[m]" << std::endl;
            std::cout << "Wall     thickness                        = " << wH_ << "\t[mills]" << std::endl;
            if ( typeH_ == "washcoated" )
            {
                std::cout << "Washcoat thickness                        = " << SwH_ << "\t[m]" << std::endl;
            }
            std::cout << "\n################################################################################################" << std::endl;
        }
        else if ( models_[1] == true )
        {
            std::cout << "                                       PACKED BED REACTOR                                       \n" << std::endl;
            std::cout << "Tube     diameter                         = " << DtP_ << "\t[m]" << std::endl;
            std::cout << "Reactor  length                           = " << LP_ << "\t[m]" << std::endl;
            std::cout << "Particle diameter                         = " << DpP_ << "\t[m]" << std::endl;
            std::cout << "\n################################################################################################" << std::endl;
        }
        else if ( models_[2] == true )
        {
            std::cout << "                                    MICRO PACKED BED REACTOR                                    \n" << std::endl;
            std::cout << "CPSI                                      = " << CPSIM_ << std::endl;
            std::cout << "Tube     diameter                         = " << DtM_ << "\t[m]" << std::endl;
            std::cout << "Particle diameter                         = " << DpM_ << "\t[m]" << std::endl;
            std::cout << "Reactor  length                           = " << LM_ << "\t[m]" << std::endl;
            std::cout << "Wall     thickness                        = " << wM_ << "\t[mills]" << std::endl;
            std::cout << "\n################################################################################################" << std::endl;
        }
        else
        {
            std::cout << "\nNone of the reactor types has been selected" << std::endl;
            std::cout << "\n################################################################################################" << std::endl;
        }
        std::cout << "                                      CATALYST PROPERTIES                                       \n" << std::endl;
        std::cout << "Density                                  = " << rhoC_ << "\t[Kg/m3]" << std::endl;
        std::cout << "Specific heat                            = " << cpC_ << "\t[J/Kg/K]" << std::endl;
        std::cout << "Conductivity                             = " << kC_ << "\t[W/m/K]" << std::endl;
        std::cout << "Tortuosity                               = " << tauC_ << "\t[-]" << std::endl;
        std::cout << "Void fraction                            = " << epsiC_ << "\t[-]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                      SUPPORT  PROPERTIES                                       \n" << std::endl;
        std::cout << "Density                                  = " << rhoS_ << "\t[Kg/m3]" << std::endl;
        std::cout << "Specific heat                            = " << cpS_ << "\t[J/Kg/K]" << std::endl;
        std::cout << "Conductivity                             = " << kS_ << "\t[W/m/K]" << std::endl;
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
        if (energy_ == true)
        {
            std::cout << "Feed  temperature                        = " << Tin_ << "\t[K]" << std::endl;
            std::cout << "Shell temperature                        = " << Tw_ << "\t[K]" << std::endl;
        }
        else
        {
            std::cout << "Energy:                                    " << "off" << std::endl;
        }
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                         SOLVER OPTIONS                                         \n" << std::endl;
        std::cout << "Axial  number of points:                   " << Na_ << std::endl;
        std::cout << "Radial number of points:                   " << Nr_ << std::endl;
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

}
