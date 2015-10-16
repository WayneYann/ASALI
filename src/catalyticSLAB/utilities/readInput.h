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

            inline double getAbsTol()              const  { return absTol_;};
            inline double getRelTol()              const  { return relTol_;};
            inline double getIntegrationTime()     const  { return tMax_;};
            inline double getReactorLength()       const  { return L_;};
            inline double getReactorDiameter()     const  { return D_;};
            inline double getAlfa()                const  { return alfa_;};
            inline double getPressure()            const  { return p_;};
            inline double getTemperature()         const  { return T_;};
            inline double getVoidFraction()        const  { return epsi_;};
            inline double getTourtuosity()         const  { return tau_;};
            inline double getPoreRadius()          const  { return pore_;};

            inline std::string getResults()        const  { return results_;};
            inline std::string getFeed()           const  { return feed_;};
            inline std::string getOdeSolver()      const  { return odeSolver_;};
            inline std::string getDaeSolver()      const  { return daeSolver_;};
            inline std::string getKineticsPath()   const  { return kineticsPath_;};
            inline std::string getDiffusionType()  const  { return gasDiffusionType_;};
            inline std::string getInert()           const { return inert_;};

            inline bool getReactions()             const  { return reactions_;};
            inline bool getConstraints()           const  { return constraints_;};

            inline unsigned int getNumberOfPoints() const  { return N_;};

            inline std::vector<double> getFraction()          const { return inletValue_;};

            inline std::vector<std::string> getFractionName() const { return inletName_;};

            void recapOnScreen();

        private:

            double absTol_;
            double relTol_;
            double tMax_;
            double L_;
            double D_;
            double alfa_;
            double p_;
            double T_;
            double epsi_;
            double tau_;
            double pore_;

            std::string results_;
            std::string feed_;
            std::string start_;
            std::string odeSolver_;
            std::string daeSolver_;
            std::string kineticsPath_;
            std::string gasDiffusionType_;
            std::string inert_;

            bool constraints_;
            bool reactions_;

            const std::string& file_;

            std::vector<double>      inletValue_;

            std::vector<std::string> inputVector_;
            std::vector<std::string> inletName_;

            unsigned int reactorIndex_;
            unsigned int solverIndex_;
            unsigned int catalystIndex_;
            unsigned int temperatureIndex_;
            unsigned int fractionIndex_;
            unsigned int solidIndex_;
            unsigned int numericalIndex_;
            unsigned int kineticsIndex_;
            unsigned int N_;

            void error() { std::cout << "\nASALI::READinput::ERROR\n" << std::endl;};

            void save();
            void check();
            void solver();
            void kinetics();
            void reactor();
            void catalyst();
            void fractions();
            void pressure();
            void temperature();
            void solid();
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
        N_                 = 0;
        absTol_            = 0;
        relTol_            = 0;
        tMax_              = 0;
        L_                 = 0;
        D_                 = 0;
        alfa_              = 0;
        p_                 = 0;
        T_                 = 0;
        epsi_              = 0;
        tau_               = 0;
        pore_              = 0;

        constraints_       = false;
        reactions_         = false;

        reactorIndex_      = 0;
        solverIndex_       = 0;
        catalystIndex_     = 0;
        temperatureIndex_  = 0;
        fractionIndex_     = 0;
        solidIndex_        = 0;
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
        solver();
        kinetics();
        reactor();
        catalyst();
        fractions();
        pressure();
        temperature();
        solid();
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
        words[2] = "Reactor";
        words[3] = "Catalyst";
        words[4] = "Mole fractions || Mass fractions";
        words[5] = "Solver options";
        words[6] = "Solid";
        words[7] = "Numerical solvers";
        words[8] = "Kinetics path";

        for (unsigned int i=0;i<inputVector_.size();i++)
        {
            if      (inputVector_[i]   == "Temperature")         {checkWord[0]  = true; temperatureIndex_  = i;}
            else if (inputVector_[i]   == "Pressure")            {checkWord[1]  = true;}
            else if (inputVector_[i]   == "Reactor")             {checkWord[2]  = true; reactorIndex_      = i;}
            else if (inputVector_[i]   == "Catalyst")            {checkWord[3]  = true; catalystIndex_     = i;}
            else if (inputVector_[i]   == "Mole" &&
                     inputVector_[i+1] == "fractions" )          {checkWord[4]  = true; fractionIndex_     = i;}
            else if (inputVector_[i]   == "Mass" &&
                     inputVector_[i+1] == "fractions" )          {checkWord[4]  = true; fractionIndex_     = i;}
            else if (inputVector_[i]   == "Solver" &&
                     inputVector_[i+1] == "options" )            {checkWord[5]  = true; solverIndex_       = i;}
            else if (inputVector_[i]   == "Solid")               {checkWord[6]  = true; solidIndex_        = i;}
            else if (inputVector_[i]   == "Numerical" &&
                     inputVector_[i+1] == "solvers" )            {checkWord[7]  = true; numericalIndex_       = i;}
            else if (inputVector_[i]   == "Kinetics" &&
                     inputVector_[i+1] == "path" )               {checkWord[8]  = true; kineticsIndex_       = i;}
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

            std::vector<bool>        checkWord(8);
            std::vector<std::string> words(8);

            unsigned int absIndex;
            unsigned int relIndex;
            unsigned int tIndex;
            unsigned int pointsIndex;
            unsigned int resultsIndex;
            unsigned int reactionIndex;
            unsigned int constraintIndex;
            unsigned int diffIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "Absolute tollerance";
            words[1] = "Relative tollerance";
            words[2] = "Number of points";
            words[3] = "Reactions";
            words[4] = "Integration time";
            words[5] = "Results";
            words[6] = "Constraints";
            words[7] = "Diffusion model";

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if      (dummyVector[i] == "Absolute" &&
                         dummyVector[i+1] == "tollerance")         {checkWord[0] = true; absIndex         = i;}
                else if (dummyVector[i] == "Relative" &&
                         dummyVector[i+1] == "tollerance")         {checkWord[1] = true; relIndex         = i;}
                else if (dummyVector[i] == "Number" &&
                         dummyVector[i+1] == "of" &&
                         dummyVector[i+2] == "points")             {checkWord[2] = true; pointsIndex      = i;}
                else if (dummyVector[i] == "Reactions")            {checkWord[3] = true; reactionIndex    = i;}
                else if (dummyVector[i] == "Integration" &&
                         dummyVector[i+1] == "time")               {checkWord[4] = true; tIndex           = i;}
                else if (dummyVector[i] == "Results")              {checkWord[5] = true; resultsIndex     = i;}
                else if (dummyVector[i] == "Constraints")          {checkWord[6] = true; constraintIndex  = i;}
                else if (dummyVector[i] == "Diffusion" &&
                         dummyVector[i+1] == "model")              {checkWord[7] = true; diffIndex         = i;}
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

            if ( dummyVector[reactionIndex+1] == "on" )
                reactions_ = true;
            else if ( dummyVector[reactionIndex+1] == "true" )
                reactions_ = true;
            else if ( dummyVector[reactionIndex+1] == "yes" )
                reactions_ = true;
            else if ( dummyVector[reactionIndex+1] == "1" )
                reactions_ = true;
            else
                reactions_ = false;

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

            results_ = dummyVector[resultsIndex+1];
            if ( results_ != "mass" && results_ != "mole" )
            {
                error();
                std::cout << "key word || " << "Results" << " || MUST be || mass || mole ||\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            
            tMax_ = boost::lexical_cast<double>(dummyVector[tIndex+2]);
            std::string tdim = dummyVector[tIndex+3];
            ConvertsToSecond(tMax_, tdim);

            absTol_ = boost::lexical_cast<double>(dummyVector[absIndex+2]);
            relTol_ = boost::lexical_cast<double>(dummyVector[relIndex+2]);
            
            N_      = boost::lexical_cast<unsigned int>(dummyVector[pointsIndex+3]);

            gasDiffusionType_ = dummyVector[diffIndex+2];
            if ( gasDiffusionType_ != "Fick" && gasDiffusionType_ != "DustyGas" && gasDiffusionType_ != "Fick-DustyGas")
            {
                error();
                std::cout << "key word ||" << " Diffusion model " << "|| MUST be || Fick || DustyGas || Fick-DustyGas ||\n" << std::endl;
                exit (EXIT_FAILURE);
            }


        }
    }

    void READinput::reactor()
    {
        if ( inputVector_[reactorIndex_+1] != "{" )
        {
            error();
            std::cout << "The Reactor sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=reactorIndex_;i<inputVector_.size();i++)
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
                std::cout << "The Reactor sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            std::vector<std::string> dummyVector;
            unsigned int k=0;
            for (unsigned int i=reactorIndex_+1+1;i<=finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Reactor sub-dictionary must finish with }\n" << std::endl;
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

            unsigned int lengthIndex;
            unsigned int diameterIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "length";
            words[1] = "diameter";

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if      (dummyVector[i] == "length")             {checkWord[0] = true; lengthIndex     = i;}
                else if (dummyVector[i] == "diameter"   )        {checkWord[1] = true; diameterIndex   = i;}
            }

            for (unsigned int i=0;i<checkWord.size();i++)
            {
                if ( checkWord[i] == false)
                {
                    error();
                    std::cout << "key word || " << words[i] << " || is MISSING in Reactor sub-dictionary!\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }
            
            L_ = boost::lexical_cast<double>(dummyVector[lengthIndex+1]);
            std::string dimL = dummyVector[lengthIndex+2];
            ConvertsToMeter(L_, dimL);
            
            D_ = boost::lexical_cast<double>(dummyVector[diameterIndex+1]);
            std::string dimD = dummyVector[diameterIndex+2];
            ConvertsToMeter(D_, dimD);
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

            std::vector<bool>        checkWord(5);
            std::vector<std::string> words(5);

            double massIndex;
            double dispersionIndex;
            double RhIndex;
            double deactivationIndex;
            double alfaIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "mass";
            words[1] = "dispersion";
            words[2] = "Rh fraction";
            words[3] = "deactivation factor";
            words[4] = "alfa";


            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if         (dummyVector[i] == "mass")                 {checkWord[0] = true; massIndex           = i;}
                else if (dummyVector[i] == "dispersion")         {checkWord[1] = true; dispersionIndex     = i;}
                else if (dummyVector[i] == "Rh" &&
                         dummyVector[i+1] == "fraction")        {checkWord[2] = true; RhIndex             = i;}
                else if (dummyVector[i] == "deactivation" &&
                         dummyVector[i+1] == "factor")            {checkWord[3] = true; deactivationIndex    = i;}
            }

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if(dummyVector[i] == "alfa")
                {
                    for (unsigned int k=0;k<checkWord.size();k++)
                        checkWord[k] = true;

                    alfaIndex = i;
                    break;
                }
                else if (dummyVector[i] == "mass")                 {checkWord[4] = true; massIndex           = i;}
                else if (dummyVector[i] == "dispersion")         {checkWord[4] = true; dispersionIndex     = i;}
                else if (dummyVector[i] == "Rh" &&
                         dummyVector[i+1] == "fraction")        {checkWord[4] = true; RhIndex             = i;}
                else if (dummyVector[i] == "deactivation" &&
                         dummyVector[i+1] == "factor")            {checkWord[4] = true; deactivationIndex    = i;}
            }

            for (unsigned int i=0;i<checkWord.size();i++)
            {
                if ( checkWord[i] == false)
                {
                    error();
                    std::cout << "key word || " << words[i] << " || alfa || is MISSING in Catalyst sub-dictionary!\n" << std::endl;
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
                    break;
                }
                else
                {
                    double Wcat = boost::lexical_cast<double>(dummyVector[massIndex+1]);
                    std::string dim = dummyVector[massIndex+2];
                    ConvertsToKg(Wcat, dim);
            
                    double RhDispersion      = boost::lexical_cast<double>(dummyVector[dispersionIndex+1]);
                    double RhMassFraction    = boost::lexical_cast<double>(dummyVector[RhIndex+2]);
                    double dectivationFactor = boost::lexical_cast<double>(dummyVector[deactivationIndex+2]);

                    double RhPM = 102.9;
                    double SiteD = 2.49e-08;
                    double ActiveArea = RhMassFraction*RhDispersion*Wcat/(RhPM*SiteD);
                    
                    double Ain = 0.25*3.14*D_*D_;
                    double ReactorVolume = Ain*L_;
                    
                    alfa_ = dectivationFactor*ActiveArea/ReactorVolume;
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

    void READinput::kinetics()
    {
        kineticsPath_ = inputVector_[kineticsIndex_+2];
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

    void READinput::solid()
    {
        if ( inputVector_[solidIndex_+1] != "{" )
        {
            error();
            std::cout << "The Solid sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=solidIndex_;i<inputVector_.size();i++)
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
                std::cout << "The Solid sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            std::vector<std::string> dummyVector;
            unsigned int k=0;
            for (unsigned int i=solidIndex_+1+1;i<=finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Solid sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    dummyVector.resize(k+1);
                    dummyVector[k] = inputVector_[i];
                    k++;
                }
            }

            std::vector<bool>        checkWord(3);
            std::vector<std::string> words(3);

            unsigned int tauIndex;
            unsigned int epsiIndex;
            unsigned int poreIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "void fraction";
            words[1] = "tortuosity";
            words[2] = "pore radius";


            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if      (dummyVector[i] == "void" &&
                         dummyVector[i+1] == "fraction")            {checkWord[0] = true; epsiIndex       = i;}
                else if (dummyVector[i] == "tortuosity")            {checkWord[1] = true; tauIndex     = i;}
                else if (dummyVector[i] == "pore" &&
                         dummyVector[i+1] == "radius")              {checkWord[2] = true; poreIndex       = i;}
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

            epsi_ = boost::lexical_cast<double>(dummyVector[epsiIndex + 2]);
            tau_  = boost::lexical_cast<double>(dummyVector[tauIndex  + 1]);
            pore_ = boost::lexical_cast<double>(dummyVector[poreIndex + 2]);
            std::string dim = dummyVector[poreIndex+3];
            ConvertsToMeter(pore_, dim);
        }
    }


    void READinput::recapOnScreen()
    {
        std::cout.precision(6);
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                          GENERAL INPUT                                         \n" << std::endl;
        std::cout << "Reactor lenght                           = " << L_ << "\t[m]" << std::endl;
        std::cout << "Reactor diameter                         = " << D_ << "\t[m]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                        SOLID PROPERTIES                                        \n" << std::endl;
        std::cout << "Specific catalytic area (alfa)           = " << alfa_ << "\t[1/m]" << std::endl;
        std::cout << "Void fraction                            = " << epsi_ << "\t[-]" << std::endl;
        std::cout << "Tortuosity                               = " << tau_ << "\t[-]" << std::endl;
        std::cout << "Pore radius                              = " << pore_ << "\t[m]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                      OPERATING CONDITIONS                                      \n" << std::endl;
        std::cout << "Pressure                                 = " << p_ << "\t[Pa]" << std::endl;
        std::cout << "Temperature                              = " << T_ << "\t[K]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                         SOLVER OPTIONS                                         \n" << std::endl;
        std::cout << "Reactions:                                 " << boolOnScreen(reactions_) << std::endl;
        std::cout << "Feed in:                                   " << feed_ << " fractions" << std::endl;
        std::cout << "Results in:                                " << results_ << " fractions" << std::endl;
        std::cout << "Diffusion model:                           " << gasDiffusionType_ << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                       AXIAL GRID OPTIONS                                       \n" << std::endl;
        std::cout << "Number of points                          = " << N_ << std::endl;
        std::cout << "Integration time                          = " << tMax_ << "\t[s]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                        NUMERICAL SOLVER                                        \n" << std::endl;
        std::cout << "Solver compiled with                        ";
        #if ASALI_USE_BZZ == 1
        std::cout << "|| BzzMath ";
        #endif
        #if ASALI_USE_SUNDIALS == 1
        std::cout << "|| Sundials ";
        #endif
        std::cout << "||" << std::endl;
        std::cout << "ODE:                                        " << odeSolver_ << std::endl;
        std::cout << "DAE:                                        " << daeSolver_ << std::endl;
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

            std::vector<bool>        checkWord(2);
            std::vector<std::string> words(2);

            double odeIndex;
            double daeIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "ODE";
            words[1] = "DAE";

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if         (dummyVector[i] == "ODE")                 {checkWord[0] = true; odeIndex      = i;}
                else if (dummyVector[i] == "DAE")                 {checkWord[1] = true; daeIndex      = i;}
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
            daeSolver_ = dummyVector[daeIndex + 1];

            if ( odeSolver_ != "BzzMath" && odeSolver_ != "Sundials")
            {
                error();
                std::cout << "key word || " << "ODE" << " || MUST be || BzzMath || Sundials || \n" << std::endl;
                exit (EXIT_FAILURE);
            }

            if ( daeSolver_ != "BzzMath" && daeSolver_ != "Sundials")
            {
                error();
                std::cout << "key word || " << "DAE" << " || MUST be || BzzMath || Sundials || \n" << std::endl;
                exit (EXIT_FAILURE);
            }

        }
    }
}
