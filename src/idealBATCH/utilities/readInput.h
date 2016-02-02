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
    class READinput 
    {
        public:

            READinput(std::string& file);

            #include "convert.h"

            inline double getSpecieAbsTol()        const  { return absSpecieTol_;};
            inline double getSpecieRelTol()        const  { return relSpecieTol_;};
            inline double getTemperatureAbsTol()   const  { return absTemperatureTol_;};
            inline double getTemperatureRelTol()   const  { return relTemperatureTol_;};
            inline double getIntegrationTime()     const  { return tMax_;};
            inline double getAlfa()                const  { return alfa_;};
            inline double getPressure()            const  { return p_;};
            inline double getTemperature()         const  { return T_;};
            inline double getVolume()              const  { return V_;};
            inline double getCatalyticArea()       const  { return A_;};

            inline std::string getFeed()           const  { return feed_;};
            inline std::string getOdeSolver()      const  { return odeSolver_;};
            inline std::string getKineticsPath()   const  { return kineticsPath_;};
            inline std::string getResolutionType() const  { return resolution_;};

            inline bool getHomogeneousReactions()  const  { return homo_;};
            inline bool getHeterogenousReactions() const  { return het_;};
            inline bool getConstraints()           const  { return constraints_;};
            inline bool getEnergy()                const  { return energy_;};

            inline std::vector<double> getFraction()                const { return inletValue_;};
            inline std::vector<double> getCoverage()                const { return coverageValue_;};

            inline std::vector<std::string> getFractionName()       const { return inletName_;};
            inline std::vector<std::string> getCoverageName()       const { return coverageName_;};

            void recapOnScreen();

        private:

            double absSpecieTol_;
            double relSpecieTol_;
            double absTemperatureTol_;
            double relTemperatureTol_;
            double alfa_;
            double p_;
            double T_;
            double V_;
            double A_;
            double tMax_;

            std::string feed_;
            std::string type_;
            std::string inert_;
            std::string odeSolver_;
            std::string kineticsPath_;
            std::string resolution_;

            bool constraints_;
            bool homo_;
            bool het_;
            bool energy_;

            const std::string& file_;

            std::vector<double>      inletValue_;
            std::vector<double>      coverageValue_;

            std::vector<std::string> inputVector_;
            std::vector<std::string> inletName_;
            std::vector<std::string> coverageName_;

            unsigned int solverIndex_;
            unsigned int catalystIndex_;
            unsigned int fractionIndex_;
            unsigned int numericalIndex_;
            unsigned int kineticsIndex_;
            unsigned int coverageIndex_;

            void error() { std::cout << "\nASALI::READinput::ERROR\n" << std::endl;};

            void save();
            void check();
            void solver();
            void kinetics();
            void catalyst();
            void coverage();
            void fractions();
            void pressure();
            void temperature();
            void volume();
            void numerical();

            std::string boolOnScreen(bool value);

            template < typename T > std::string to_string( const T& v )
            {
                std::ostringstream stm ;
                return ( stm << v ) ? stm.str() : "{*** error ***}" ;
            }
    };
    
    READinput::READinput(std::string& file):
    file_(file)
    {
        absSpecieTol_      = 0;
        relSpecieTol_      = 0;
        absTemperatureTol_ = 0;
        relTemperatureTol_ = 0;
        alfa_              = 0;
        p_                 = 0;
        T_                 = 0;
        V_                 = 0;
        A_                 = 0;
        tMax_              = 0;

        constraints_       = false;
        homo_              = false;
        het_               = false;
        energy_            = false;

        solverIndex_       = 0;
        catalystIndex_     = 0;
        fractionIndex_     = 0;
        numericalIndex_    = 0;
        kineticsIndex_     = 0;

        //- Check input
        {
            if ( boost::filesystem::exists(file_) && !boost::filesystem::is_directory(file_))
            {}
            else
            {
                error();
                std::cout << file_ << " MUST BE a file" << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        //- Save all 
        save();
        check();
        kinetics();
        solver();
        volume();
        catalyst();
        fractions();
        coverage();
        pressure();
        temperature();
        numerical();
    }

    void READinput::save()
    {
        std::string dummyString;
        const char *path = file_.c_str();
        std::ifstream is(path);
        unsigned int k = 0;
        while (getline(is,dummyString))
        {
            while (getline(is,dummyString))
            {
                std::istringstream iss(dummyString);
                while (iss >> dummyString)
                {
                    k++;
                    inputVector_.resize(k);
                    inputVector_[k-1] = dummyString;
                }
            }
        }
    }

    void READinput::check()
    {
        std::vector<bool>         checkWord(9);
        std::vector<std::string>  words(9);

        for (unsigned int i=0;i<checkWord.size();i++)
            checkWord[i] = false;

        words[0] = "Temperature";
        words[1] = "Pressure";
        words[2] = "Volume";
        words[3] = "Catalyst";
        words[4] = "Mole fractions || Mass fractions";
        words[5] = "Solver options";
        words[6] = "Numerical solvers";
        words[7] = "Kinetics path";
        words[8] = "Coverage";

        for (unsigned int i=0;i<inputVector_.size();i++)
        {
            if      (inputVector_[i]   == "Temperature")        {checkWord[0]  = true;}
            else if (inputVector_[i]   == "Pressure")           {checkWord[1]  = true;}
            else if (inputVector_[i]   == "Volume")             {checkWord[2]  = true;}
            else if (inputVector_[i]   == "Catalyst")           {checkWord[3]  = true; catalystIndex_     = i;}
            else if (inputVector_[i]   == "Mole" &&
                     inputVector_[i+1] == "fractions" )         {checkWord[4]  = true; fractionIndex_     = i;}
            else if (inputVector_[i]   == "Mass" &&
                     inputVector_[i+1] == "fractions" )         {checkWord[4]  = true; fractionIndex_     = i;}
            else if (inputVector_[i]   == "Solver" &&
                     inputVector_[i+1] == "options" )           {checkWord[5]  = true; solverIndex_       = i;}
            else if (inputVector_[i]   == "Numerical" &&
                     inputVector_[i+1] == "solvers" )           {checkWord[6]  = true; numericalIndex_    = i;}
            else if (inputVector_[i]   == "Kinetics" &&
                     inputVector_[i+1] == "path" )              {checkWord[7]  = true; kineticsIndex_     = i;}
            else if (inputVector_[i]   == "Coverage")           {checkWord[8]  = true; coverageIndex_     = i;}
        }
        
        for (unsigned int i=0;i<checkWord.size();i++)
        {
            if ( checkWord[i] == false )
            {
                error();
                std::cout << "key word || " << words[i] << " || " << "is MISSING!\n" << std::endl;
                exit (EXIT_FAILURE);
            }
        }
    }

    void READinput::solver()
    {
        if ( inputVector_[solverIndex_+1+1] != "{" )
        {
            error();
            std::cout << "The Solver options sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=solverIndex_;i<inputVector_.size();i++)
            {
                if (inputVector_[i] == "}")
                {
                    finalCount = i;
                    break;
                }
            }
            
            if ( finalCount == 1e05 )
            {
                error();
                std::cout << "The Solver options sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            std::vector<std::string> dummyVector;
            unsigned int k=0;
            for (unsigned int i=solverIndex_+1+1+1;i<=finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Solver options sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    dummyVector.resize(k+1);
                    dummyVector[k] = inputVector_[i];
                    k++;
                }
            }

            std::vector<bool>        checkWord(7);
            std::vector<std::string> words(7);

            double absIndex;
            double relIndex;
            double reactionIndex;
            double constraintIndex;
            double constantIndex;
            double energyIndex;
            double tIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "Absolute tollerance";
            words[1] = "Relative tollerance";
            words[2] = "Reactions";
            words[3] = "Constraints";
            words[4] = "Constant";
            words[5] = "Energy equation";
            words[6] = "Integration time";

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if      (dummyVector[i]   == "Absolute" &&
                         dummyVector[i+1] == "tollerance")         {checkWord[0] = true; absIndex        = i;}
                else if (dummyVector[i]   == "Relative" &&
                         dummyVector[i+1] == "tollerance")         {checkWord[1] = true; relIndex        = i;}
                else if (dummyVector[i]   == "Reactions")          {checkWord[2] = true; reactionIndex   = i;}
                else if (dummyVector[i]   == "Constraints")        {checkWord[3] = true; constraintIndex = i;}
                else if (dummyVector[i]   == "Constant")           {checkWord[4] = true; constantIndex   = i;}
                else if (dummyVector[i]   == "Energy" &&
                         dummyVector[i+1] == "equation")           {checkWord[5] = true; energyIndex     = i;}
                else if (dummyVector[i] == "Integration" &&
                         dummyVector[i+1] == "time")               {checkWord[6] = true; tIndex           = i;}
            }

            for (unsigned int i=0;i<checkWord.size();i++)
            {
                if ( checkWord[i] == false)
                {
                    error();
                    std::cout << "key word || " << words[i] << " || is MISSING in Solver options sub-dictionary!\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }

            if ( dummyVector[constraintIndex+1] == "on" )
                constraints_ = true;
            else if ( dummyVector[constraintIndex+1] == "true" )
                constraints_ = true;
            else if ( dummyVector[constraintIndex+1] == "yes" )
                constraints_ = true;
            else if ( dummyVector[constraintIndex+1] == "1" )
                constraints_ = true;
            else
                constraints_ = false;

            if ( dummyVector[energyIndex+2] == "on" )
                energy_ = true;
            else if ( dummyVector[energyIndex+2] == "true" )
                energy_ = true;
            else if ( dummyVector[energyIndex+2] == "yes" )
                energy_ = true;
            else if ( dummyVector[energyIndex+2] == "1" )
                energy_ = true;
            else
                energy_ = false;


            resolution_ = dummyVector[constantIndex+1];
            if ( resolution_ == "pressure" )
            {}
            else if ( resolution_ == "volume" )
            {}
            else
            {
                error();
                std::cout << "key word ||" << " Constant " << "|| MUST be || pressure || volume ||\n" << std::endl;
                exit (EXIT_FAILURE);
            }


            tMax_ = boost::lexical_cast<double>(dummyVector[tIndex+2]);
            std::string tdim = dummyVector[tIndex+3];
            ConvertsToSecond(tMax_, tdim);

            if ( dummyVector[reactionIndex+1] != "(" )
            {
                error();
                std::cout << "The Solver options/Reactions sub-dictionary must start with (\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                unsigned int finalCount = 1e05;
                for (unsigned int i=reactionIndex;i<dummyVector.size();i++)
                {
                    if (dummyVector[i] == ")")
                    {
                        finalCount = i;
                        break;
                    }
                }
                
                if ( finalCount == 1e05 )
                {
                    error();
                    std::cout << "The Solver options/Reactions sub-dictionary must finish with )\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    unsigned int k=0;
                    for (unsigned int i=reactionIndex+1+1;i<finalCount;i++)
                    {
                        if (dummyVector[i] == "(")
                        {
                            error();
                            std::cout << "The Solver options/Reactions sub-dictionary must finish with )\n" << std::endl;
                            exit (EXIT_FAILURE);
                        }
                        else
                        {
                            if ( dummyVector[i] == "Homogeneous" )
                            {
                                if ( dummyVector[i+1] == "on" )
                                    homo_ = true;
                                else if ( dummyVector[i+1] == "true" )
                                    homo_ = true;
                                else if ( dummyVector[i+1] == "yes" )
                                    homo_ = true;
                                else if ( dummyVector[i+1] == "1" )
                                    homo_ = true;
                                else
                                    homo_ = false;
                                
                                i++;
                            }
                            else if ( dummyVector[i] == "Heterogeneous" )
                            {
                                if ( dummyVector[i+1] == "on" )
                                    het_ = true;
                                else if ( dummyVector[i+1] == "true" )
                                    het_ = true;
                                else if ( dummyVector[i+1] == "yes" )
                                    het_ = true;
                                else if ( dummyVector[i+1] == "1" )
                                    het_ = true;
                                else
                                    het_ = false;
                                
                                i++;
                            }
                            else
                            {
                                error();
                                std::cout << "key word || Homogeneous || Heterogeneous || is MISSING in Solver options/Reactions sub-dictionary!\n" << std::endl;
                                exit (EXIT_FAILURE);
                            }
                        }
                    }
                }
            }

            if ( dummyVector[absIndex+1+1] != "(" )
            {
                error();
                std::cout << "The Solver options/Absolute tollerance sub-dictionary must start with (\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                unsigned int finalCount = 1e05;
                for (unsigned int i=absIndex;i<dummyVector.size();i++)
                {
                    if (dummyVector[i] == ")")
                    {
                        finalCount = i;
                        break;
                    }
                }
                
                if ( finalCount == 1e05 )
                {
                    error();
                    std::cout << "The Solver options/Absolute tollerance sub-dictionary must finish with )\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    unsigned int k=0;
                    for (unsigned int i=absIndex+1+1+1;i<finalCount;i++)
                    {
                        if (dummyVector[i] == "(")
                        {
                            error();
                            std::cout << "The Solver options/Absolute tollerance sub-dictionary must finish with )\n" << std::endl;
                            exit (EXIT_FAILURE);
                        }
                        else
                        {
                            if ( dummyVector[i] == "specie" )
                            {
                                absSpecieTol_ = boost::lexical_cast<double>(dummyVector[i+1]);
                                i++;
                            }
                            else if ( dummyVector[i] == "temperature" )
                            {
                                absTemperatureTol_ = boost::lexical_cast<double>(dummyVector[i+1]);
                                i++;
                            }
                            else
                            {
                                error();
                                std::cout << "key word || specie || temperature || is MISSING in Solver options/Absolute tollerance sub-dictionary!\n" << std::endl;
                                exit (EXIT_FAILURE);
                            }
                        }
                    }
                }
            }


            if ( dummyVector[relIndex+1+1] != "(" )
            {
                error();
                std::cout << "The Solver options/Relative tollerance sub-dictionary must start with (\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                unsigned int finalCount = 1e05;
                for (unsigned int i=relIndex;i<dummyVector.size();i++)
                {
                    if (dummyVector[i] == ")")
                    {
                        finalCount = i;
                        break;
                    }
                }
                
                if ( finalCount == 1e05 )
                {
                    error();
                    std::cout << "The Solver options/Relative tollerance sub-dictionary must finish with )\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    unsigned int k=0;
                    for (unsigned int i=relIndex+1+1+1;i<finalCount;i++)
                    {
                        if (dummyVector[i] == "(")
                        {
                            error();
                            std::cout << "The Solver options/Relative tollerance sub-dictionary must finish with )\n" << std::endl;
                            exit (EXIT_FAILURE);
                        }
                        else
                        {
                            if ( dummyVector[i] == "specie" )
                            {
                                relSpecieTol_ = boost::lexical_cast<double>(dummyVector[i+1]);
                                i++;
                            }
                            else if ( dummyVector[i] == "temperature" )
                            {
                                relTemperatureTol_ = boost::lexical_cast<double>(dummyVector[i+1]);
                                i++;
                            }
                            else
                            {
                                error();
                                std::cout << "key word || specie || temperature || is MISSING in Solver options/Relative tollerance sub-dictionary!\n" << std::endl;
                                exit (EXIT_FAILURE);
                            }
                        }
                    }
                }
            }
        }
    }

    void READinput::catalyst()
    {
        if ( inputVector_[catalystIndex_+1] != "{" )
        {
            error();
            std::cout << "The Catalyst sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=catalystIndex_;i<inputVector_.size();i++)
            {
                if (inputVector_[i] == "}")
                {
                    finalCount = i;
                    break;
                }
            }
            
            if ( finalCount == 1e05 )
            {
                error();
                std::cout << "The Catalyst sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            std::vector<std::string> dummyVector;
            unsigned int k=0;
            for (unsigned int i=catalystIndex_+1+1;i<=finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Catalyst sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    dummyVector.resize(k+1);
                    dummyVector[k] = inputVector_[i];
                    k++;
                }
            }

            std::vector<bool>        checkWord(2);
            std::vector<std::string> words(2);

            double alfaIndex;
            double areaIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "alfa";
            words[1] = "area";


            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if      (dummyVector[i] == "alfa")                 {checkWord[0] = true; alfaIndex     = i;}
                else if (dummyVector[i] == "area")                 {checkWord[1] = true; areaIndex     = i;}
            }

            for (unsigned int i=0;i<checkWord.size();i++)
            {
                if ( checkWord[i] == false)
                {
                    error();
                    std::cout << "key word || " << words[i] << " || alfa || is MISSING in Reactor sub-dictionary!\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }
            
            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if(dummyVector[i] == "alfa")
                {
                    alfa_ = boost::lexical_cast<double>(dummyVector[alfaIndex+1]);
                    std::string dim = dummyVector[alfaIndex+2];
                    ConvertsToOneOnMeter(alfa_,dim);
                    i++;
                    i++;
                }
                else if(dummyVector[i] == "area")
                {
                    A_ = boost::lexical_cast<double>(dummyVector[areaIndex+1]);
                    std::string dim = dummyVector[areaIndex+2];
                    ConvertsToMeterSquare(A_,dim);
                    i++;
                    i++;
                }
            }
        }
    }

    void READinput::fractions()
    {
        if ( inputVector_[fractionIndex_+1+1] != "{" )
        {
            error();
            std::cout << "The Mole/Mass fractions sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=fractionIndex_;i<inputVector_.size();i++)
            {
                if (inputVector_[i] == "}")
                {
                    finalCount = i;
                    break;
                }
            }
            
            if ( finalCount == 1e05 )
            {
                error();
                std::cout << "The Mole/Mass fractions sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            unsigned int k=0;
            for (unsigned int i=fractionIndex_+1+1+1;i<finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Mole/Mass fractions sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else if ( inputVector_[i] == "inert" )
                {
                    inert_ = inputVector_[i+1];
                    i++;
                }
                else
                {
                    inletValue_.resize(k+1);
                    inletName_.resize(k+1);
                    inletValue_[k] = boost::lexical_cast<double>(inputVector_[i+1]);
                    inletName_[k]  = inputVector_[i];
                    k++;
                    i++;
                }
            }

            bool test = false;
            for (unsigned int i=0;i<inletName_.size();i++)
            {
                if ( inletName_[i] == inert_)
                {
                    test = true;
                }
            }
            
            if ( test == false )
            {
                error();
                std::cout << "key word || inert || MUST BE one of the feed species!\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            feed_ = inputVector_[fractionIndex_];
        }
    }

    void READinput::pressure()
    {
        for (unsigned int i=0;i<inputVector_.size();i++)
        {
            if ( inputVector_[i] == "Pressure")
            {
                p_ = boost::lexical_cast<double>(inputVector_[i+1]);
                std::string dim = inputVector_[i+2];
                ConvertsToPascal(p_, dim);
            }
        }
    }

    void READinput::temperature()
    {
        for (unsigned int i=0;i<inputVector_.size();i++)
        {
            if ( inputVector_[i] == "Temperature")
            {
                T_ = boost::lexical_cast<double>(inputVector_[i+1]);
                std::string dim = inputVector_[i+2];
                if ( dim == "°C")
                    FromCelsiusToKelvin(T_,dim);
            }
        }
    }

    void READinput::volume()
    {
        for (unsigned int i=0;i<inputVector_.size();i++)
        {
            if ( inputVector_[i] == "Volume")
            {
                V_ = boost::lexical_cast<double>(inputVector_[i+1]);
                std::string dim = inputVector_[i+2];
                ConvertsToMeterCube(V_, dim);
            }
        }
    }

    void READinput::kinetics()
    {
        kineticsPath_ = inputVector_[kineticsIndex_+2];
    }

    void READinput::recapOnScreen()
    {
        std::cout.precision(6);
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                          GENERAL INPUT                                         \n" << std::endl;
        std::cout << "Volume                                   = " << V_ << "\t[m3]" << std::endl;
        std::cout << "Catalytic area                           = " << A_ << "\t[m2]" << std::endl;
        std::cout << "Catalytic load                           = " << alfa_ << "\t[1/m]" << std::endl;
        std::cout << "Pressure                                 = " << p_ << "\t[Pa]" << std::endl;
        std::cout << "Temperature                              = " << T_ << "\t[K]"  << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                         SOLVER OPTIONS                                         \n" << std::endl;
        std::cout << "Homogeneous  reactions:                    " << boolOnScreen(homo_) << std::endl;
        std::cout << "Heterogeneus reactions:                    " << boolOnScreen(het_) << std::endl;
        std::cout << "Energy:                                    " << boolOnScreen(het_) << std::endl;
        std::cout << "Feed in:                                   " << feed_ << " fractions" << std::endl;
        std::cout << "Solved with constant:                      " << resolution_ << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                        NUMERICAL SOLVER                                        \n" << std::endl;
        std::cout << "Solver compiled with                        ";
        #if ASALI_USE_BZZ == 1
        std::cout << "|| BzzMath || OpenSMOKE ";
        #endif
        #if ASALI_USE_SUNDIALS == 1
        std::cout << "|| Sundials || OpenSMOKE ";
        #endif
        std::cout << "||" << std::endl;
        std::cout << "ODE:                                        " << odeSolver_ << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
    }
    
    std::string READinput::boolOnScreen(bool value)
    {
        std::string value_;
        if ( value == true)
            value_ = "on";
        else
            value_ = "off";
        
        return value_;
    }

    void READinput::numerical()
    {
        if ( inputVector_[numericalIndex_+1+1] != "{" )
        {
            error();
            std::cout << "The Numerical solvers sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=numericalIndex_;i<inputVector_.size();i++)
            {
                if (inputVector_[i] == "}")
                {
                    finalCount = i;
                    break;
                }
            }
            
            if ( finalCount == 1e05 )
            {
                error();
                std::cout << "The Numerical solvers sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            std::vector<std::string> dummyVector;
            unsigned int k=0;
            for (unsigned int i=numericalIndex_+1+1+1;i<=finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Numerical solvers sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    dummyVector.resize(k+1);
                    dummyVector[k] = inputVector_[i];
                    k++;
                }
            }

            std::vector<bool>        checkWord(1);
            std::vector<std::string> words(1);

            double odeIndex;
            double daeIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "ODE";

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if(dummyVector[i] == "ODE")                 {checkWord[0] = true; odeIndex      = i;};
            }

            for (unsigned int i=0;i<checkWord.size();i++)
            {
                if ( checkWord[i] == false)
                {
                    error();
                    std::cout << "key word || " << words[i] << " || is MISSING in Numerical solvers sub-dictionary!\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }

            odeSolver_ = dummyVector[odeIndex + 1];

            if ( odeSolver_ != "BzzMath" && odeSolver_ != "Sundials" && odeSolver_ != "OpenSMOKE" )
            {
                error();
                std::cout << "key word || " << "ODE" << " || MUST be || BzzMath || Sundials || OpenSMOKE ||\n" << std::endl;
                exit (EXIT_FAILURE);
            }
        }
    }

    void READinput::coverage()
    {
        if ( inputVector_[coverageIndex_+1] != "{" )
        {
            error();
            std::cout << "The Mole/Mass fractions sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=coverageIndex_;i<inputVector_.size();i++)
            {
                if (inputVector_[i] == "}")
                {
                    finalCount = i;
                    break;
                }
            }

            if ( finalCount == 1e05 )
            {
                error();
                std::cout << "The Coverage sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            unsigned int k=0;
            for (unsigned int i=coverageIndex_+1+1;i<finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Coverage sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    coverageValue_.resize(k+1);
                    coverageName_.resize(k+1);
                    coverageValue_[k] = boost::lexical_cast<double>(inputVector_[i+1]);
                    coverageName_[k]  = inputVector_[i];
                    k++;
                    i++;
                }
            }
        }
    }
}
