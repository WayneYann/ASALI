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

            double getParticleDiameter()    const  { return Dp_;};
            double getShellDiameter()       const  { return Dt_;};
            double getReactorLength()       const  { return L_;};
            double getVoidFraction()        const  { return epsi_;};
            double getFeedPressure()        const  { return p_;};
            double getFeedVelocity()        const  { return v_;};
            double getFeedTemperature()     const  { return Tin_;};

            std::string getReactorType()               const  { return reactorType_;};
            std::string getKineticsPath()              const  { return kineticsPath_;};
            std::string getPressureCorrelation()       const  { return pCorr_;};
            std::string getMassTransferCorrelation()   const  { return hCorr_;};
            std::string getOdeSolver()                 const  { return odeSolver_;};
            std::string getDaeSolver()                 const  { return daeSolver_;};
            std::string getInletType()                 const  { return inlet_;};
            std::string getResultsType()               const  { return results_;};
            std::string getInertSpecie()               const  { return inert_;};

            std::vector<double> getFeedComposition()              const  { return inletValue_; };

            std::vector<std::string> getFeedName()                const  { return inletName_; };

            void recapOnScreen();

        private:

            double Dp_;
            double Dt_;
            double L_;
            double epsi_;
            double p_;
            double v_;
            double Tin_;

            std::string reactorType_;
            std::string kineticsPath_;
            std::string hCorr_;
            std::string pCorr_;
            std::string odeSolver_;
            std::string daeSolver_;
            std::string inert_;
            std::string results_;
            std::string inlet_;

            unsigned int typeIndex_;
            unsigned int operatingIndex_;
            unsigned int propertiesIndex_;
            unsigned int solverIndex_;
            unsigned int numericalIndex_;
            unsigned int kineticsIndex_;
            unsigned int reactionIndex_;

            const std::string& file_;

            std::vector<std::string> inletName_;
            std::vector<std::string> inputVector_;

            std::vector<double>      inletValue_;

            void error() { std::cout << "\nASALI::READinput::ERROR\n" << std::endl;};

            void save();
            void check();
            void type();
            void reactor();
            void solver();
            void operating();
            void numerical();
            void transport();
            void reaction();

            template < typename T > std::string to_string( const T& v )
            {
                std::ostringstream stm ;
                return ( stm << v ) ? stm.str() : "{*** error ***}" ;
            }
    };
    
    READinput::READinput(std::string& file):
    file_(file)
    {

        Dp_      = 0.;
        Dt_      = 0.;
        L_       = 0.;
        epsi_    = 0.;
        Tin_     = 0.;
        v_       = 0.;
        p_       = 0.;


        typeIndex_         = 0;
        operatingIndex_    = 0;
        propertiesIndex_   = 0;
        solverIndex_       = 0;
        numericalIndex_    = 0;
        kineticsIndex_    = 0;
        reactionIndex_     = 0;

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
        type();
        transport();
        reactor();
        operating();
        solver();
        numerical();
        reaction();

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
        if (inputVector_[0] != "#")
        {
            error();
            exit(EXIT_FAILURE);
        }
        else
        {
            std::vector<bool>         checkWord(7);
            std::vector<std::string>  words(7);

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "Reactor type";
            words[1] = "Operating conditions";
            words[2] = "Reactor properties";
            words[3] = "Solver options";
            words[4] = "Numerical solvers";
            words[5] = "Kinetics path";
            words[6] = "Reaction properties";

            for (unsigned int i=0;i<inputVector_.size();i++)
            {
                if      (inputVector_[i]   == "Reactor" &&
                         inputVector_[i+1] == "type" )               {checkWord[0]  = true; typeIndex_            = i;}
                else if (inputVector_[i]   == "Operating" &&
                         inputVector_[i+1] == "conditions" )         {checkWord[1]  = true; operatingIndex_       = i;}
                else if (inputVector_[i]   == "Reactor" &&
                         inputVector_[i+1] == "properties" )         {checkWord[2]  = true; propertiesIndex_      = i;}
                else if (inputVector_[i]   == "Solver" &&
                         inputVector_[i+1] == "options" )            {checkWord[3]  = true; solverIndex_          = i;}
                else if (inputVector_[i]   == "Numerical" &&
                         inputVector_[i+1] == "solvers" )            {checkWord[4]  = true; numericalIndex_       = i;}
                else if (inputVector_[i]   == "Kinetics" &&
                         inputVector_[i+1] == "path" )               {checkWord[5]  = true; kineticsIndex_       = i;}
                else if (inputVector_[i]   == "Reaction" &&
                         inputVector_[i+1] == "properties" )         {checkWord[6]  = true; reactionIndex_        = i;}
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

    }

    void READinput::type()
    {
        if ( inputVector_[typeIndex_+2] == "Monolith" )
        {
            reactorType_ = inputVector_[typeIndex_+2];
        }
        else if ( inputVector_[typeIndex_+2] == "PackedBed" )
        {
            reactorType_ = inputVector_[typeIndex_+2];
        }
        else
        {
            error();
            std::cout << "key word || " << "Reactor type" << " || MUST be || Monolith || PackedBed ||\n" << std::endl;
            exit (EXIT_FAILURE);
        }
    }

    void READinput::reactor()
    {
        if ( inputVector_[propertiesIndex_+2] != "{" )
        {
            error();
            std::cout << "The Reactor properties sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=propertiesIndex_;i<inputVector_.size();i++)
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
                std::cout << "The Reactor properties sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            std::vector<std::string> dummyVector;
            unsigned int k=0;
            for (unsigned int i=propertiesIndex_+3;i<=finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Reactor properties sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    dummyVector.resize(k+1);
                    dummyVector[k] = inputVector_[i];
                    k++;
                }
            }
            
            
            unsigned int void_particleIndex = 0;
            unsigned int lengthIndex        = 0;
            unsigned int shellIndex         = 0;
            
            std::vector<bool>        checkWord(3);
            std::vector<std::string> words(3);

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            if ( reactorType_ == "Monolith" )
            {
                words[0] = "void fraction";
            }
            else if ( reactorType_ == "PackedBed" )
            {
                words[0] = "particle diameter";
            }
            words[1] = "reactor length";
            words[2] = "shell diameter";

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if ( reactorType_ == "Monolith" )
                {
                    if      (dummyVector[i] == "void" &&
                             dummyVector[i+1] == "fraction")         {checkWord[0] = true; void_particleIndex  = i;}
                    else if (dummyVector[i] == "reactor" &&
                             dummyVector[i+1] == "length")           {checkWord[1] = true; lengthIndex         = i;}
                    else if (dummyVector[i] == "shell" &&
                             dummyVector[i+1] == "diameter")          {checkWord[2] = true; shellIndex          = i;}
                }
                else if ( reactorType_ == "PackedBed" )
                {
                    if      (dummyVector[i] == "particle" &&
                             dummyVector[i+1] == "diameter")         {checkWord[0] = true; void_particleIndex  = i;}
                    else if (dummyVector[i] == "reactor" &&
                             dummyVector[i+1] == "length")           {checkWord[1] = true; lengthIndex         = i;}
                    else if (dummyVector[i] == "shell" &&
                             dummyVector[i+1] == "diameter")          {checkWord[2] = true; shellIndex          = i;}
                }
            }

            for (unsigned int i=0;i<checkWord.size();i++)
            {
                if ( checkWord[i] == false)
                {
                    error();
                    std::cout << "key word || " << words[i] << " || is MISSING in Reactor properties sub-dictionary!\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }

            if ( reactorType_ == "Monolith" )
            {
                epsi_ = boost::lexical_cast<double>(dummyVector[void_particleIndex+2]);
            }
            else if ( reactorType_ == "PackedBed" )
            {
                Dp_ = boost::lexical_cast<double>(dummyVector[void_particleIndex+2]);
                std::string dimDp = dummyVector[void_particleIndex+3];
                ConvertsToMeter(Dp_, dimDp);
            }
            
            {
                L_ = boost::lexical_cast<double>(dummyVector[lengthIndex+2]);
                std::string dimL = dummyVector[lengthIndex+3];
                ConvertsToMeter(L_, dimL);
            }
            
            {
                Dt_ = boost::lexical_cast<double>(dummyVector[shellIndex+2]);
                std::string dimDt = dummyVector[shellIndex+3];
                ConvertsToMeter(Dt_, dimDt);
            }

        }
    }

    void READinput::operating()
    {
        if ( inputVector_[operatingIndex_+2] != "{" )
        {
            error();
            std::cout << "The Operating conditions sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=operatingIndex_;i<inputVector_.size();i++)
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
                std::cout << "The Operating conditions sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            std::vector<std::string> dummyVector;
            unsigned int k=0;
            for (unsigned int i=operatingIndex_+3;i<=finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Operating conditions sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    dummyVector.resize(k+1);
                    dummyVector[k] = inputVector_[i];
                    k++;
                }
            }
            
            
            unsigned int feedTindex   = 0;
            unsigned int pIndex       = 0;
            unsigned int vIndex       = 0;
            
            std::vector<bool>        checkWord(3);
            std::vector<std::string> words(3);

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "feed temperature";
            words[1] = "feed pressure";
            words[2] = "velocity";

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if         (dummyVector[i] == "feed" &&
                         dummyVector[i+1] == "temperature")         {checkWord[0] = true; feedTindex      = i;}
                else if (dummyVector[i] == "feed" &&
                         dummyVector[i+1] == "pressure")            {checkWord[1] = true; pIndex          = i;}
                else if (dummyVector[i] == "velocity")                {checkWord[2] = true; vIndex          = i;}
            }

            for (unsigned int i=0;i<checkWord.size();i++)
            {
                if ( checkWord[i] == false)
                {
                    error();
                    std::cout << "key word || " << words[i] << " || is MISSING in Operating conditions sub-dictionary!\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }

            {
                Tin_ = boost::lexical_cast<double>(dummyVector[feedTindex+2]);
                std::string dim = dummyVector[feedTindex+3];
                if ( dim == "°C")
                    FromCelsiusToKelvin(Tin_,dim);
            }

            {
                p_ = boost::lexical_cast<double>(dummyVector[pIndex+2]);
                std::string dim = dummyVector[pIndex+3];
                ConvertsToPascal(p_, dim);
            }

            if ( dummyVector[vIndex + 1 + 1] == "m/s" )
            {
                v_ = boost::lexical_cast<double>(dummyVector[vIndex + 1]);
            }
            else
            {
                error();
                std::cout << "key word || " << "velocity" << " || accepted unit is m/s in Operating conditions sub-dictionary!\n" << std::endl;
                exit (EXIT_FAILURE);
            }

        }
    }

    void READinput::transport()
    {
        kineticsPath_ = inputVector_[kineticsIndex_+2];
    }

    void READinput::solver()
    {
        if ( inputVector_[solverIndex_+2] != "{" )
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
            for (unsigned int i=solverIndex_+3;i<=finalCount;i++)
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
            
            
            unsigned int pIndex      = 0;
            unsigned int hIndex      = 0;
            
            std::vector<bool>        checkWord(2);
            std::vector<std::string> words(2);

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "pressure drops";
            words[1] = "mass transfer";

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if        (dummyVector[i] == "pressure" &&
                         dummyVector[i+1] == "drops")             {checkWord[0] = true; pIndex       = i;}
                else if (dummyVector[i] == "mass" &&
                         dummyVector[i+1] == "transfer")        {checkWord[1] = true; hIndex       = i;}
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

            if ( dummyVector[pIndex + 2] == "Foumeny")
            {
                pCorr_ = dummyVector[pIndex + 2];
            }
            else if ( dummyVector[pIndex + 2] == "Ergun")
            {
                pCorr_ = dummyVector[pIndex + 2];
            }
            else if ( dummyVector[pIndex + 2] == "Lee")
            {
                pCorr_ = dummyVector[pIndex + 2];
            }
            else if ( dummyVector[pIndex + 2] == "Hicks")
            {
                pCorr_ = dummyVector[pIndex + 2];
            }
            else if ( dummyVector[pIndex + 2] == "Eisfeld")
            {
                pCorr_ = dummyVector[pIndex + 2];
            }
            else if ( dummyVector[pIndex + 2] == "Achenbach")
            {
                pCorr_ = dummyVector[pIndex + 2];
            }
            else
            {
                error();
                std::cout << "key word || pressure drops || MUST BE || Foumeny || Ergun || Lee || Hicks || Eisfeld || Achenbach ||\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            if ( dummyVector[hIndex + 2] == "Wakao")
            {
                hCorr_ = dummyVector[hIndex + 2];
            }
            else if ( dummyVector[hIndex + 2] == "Gupta")
            {
                hCorr_ = dummyVector[hIndex + 2];
            }
            else if ( dummyVector[hIndex + 2] == "Gamson")
            {
                hCorr_ = dummyVector[hIndex + 2];
            }
            else if ( dummyVector[hIndex + 2] == "Yoshida")
            {
                hCorr_ = dummyVector[hIndex + 2];
            }
            else if ( dummyVector[hIndex + 2] == "Monolith")
            {
                hCorr_ = dummyVector[hIndex + 2];
            }
            else if ( dummyVector[hIndex + 2] == "Rebughini")
            {
                hCorr_ = dummyVector[hIndex + 2];
            }
             else if ( dummyVector[hIndex + 2] == "Petrovic")
            {
                hCorr_ = dummyVector[hIndex + 2];
            }
            else
            {
                error();
                std::cout << "key word || mass transfer || MUST BE || Wakao || Gupta || Gamson || Yoshida || Rebughini || Petrovic ||\n" << std::endl;
                exit (EXIT_FAILURE);
            }
        }
    }

    void READinput::recapOnScreen()
    {
        std::cout.precision(6);
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                       REACTOR PROPERTIES                                       \n" << std::endl;
        std::cout << "Reactor type                               " << reactorType_ << std::endl;
        std::cout << "Reactor lenght                           = " << L_ << "\t[m]" << std::endl;
        std::cout << "Shell   diameter                         = " << Dt_ << "\t[m]" << std::endl;
        if ( reactorType_ == "Monolith" )
        {
            std::cout << "Void fraction                            = " << epsi_ << std::endl;
        }
        else if ( reactorType_ == "PackedBed" )
        {
            std::cout << "Particle diameter                        = " << Dp_ << "\t[m]" << std::endl;
        }
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                      OPERATING CONDITIONS                                      \n" << std::endl;
        std::cout << "Feed  velocity                            = " << v_ << "\t[m/s]" << std::endl;
        std::cout << "Feed  pressure                            = " << p_ << "\t[Pa]" << std::endl;
        std::cout << "Feed  temperature                         = " << Tin_ << "\t[K]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                         SOLVER OPTIONS                                         \n" << std::endl;
        std::cout << "Pressure drops correlation:                 " << pCorr_ << std::endl;
        std::cout << "Mass  transfer correlation:                 " << hCorr_ << std::endl;
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
                else if    (dummyVector[i] == "DAE")                 {checkWord[1] = true; daeIndex      = i;}
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

            if ( odeSolver_ != "BzzMath" && odeSolver_ != "Sundials")
            {
                error();
                std::cout << "key word || " << "ODE" << " || MUST be || BzzMath || Sundials || \n" << std::endl;
                exit (EXIT_FAILURE);
            }

            daeSolver_ = dummyVector[daeIndex + 1];

            if ( daeSolver_ != "BzzMath" && daeSolver_ != "Sundials")
            {
                error();
                std::cout << "key word || " << "DAE" << " || MUST be || BzzMath || Sundials || \n" << std::endl;
                exit (EXIT_FAILURE);
            }

        }
    }

    void READinput::reaction()
    {
        if ( inputVector_[reactionIndex_+2] != "{" )
        {
            error();
            std::cout << "The Reaction properties sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=reactionIndex_;i<inputVector_.size();i++)
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
                std::cout << "The Reaction properties sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            std::vector<std::string> dummyVector;
            unsigned int k=0;
            for (unsigned int i=reactionIndex_+3;i<=finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Reaction properties sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    dummyVector.resize(k+1);
                    dummyVector[k] = inputVector_[i];
                    k++;
                }
            }
            
            unsigned int fractionIndex   = 0;
            unsigned int inertIndex      = 0;
            unsigned int resultsIndex    = 0;
            
            std::vector<bool>        checkWord(3);
            std::vector<std::string> words(3);

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "mole fraction || mass fraction";
            words[1] = "inert specie";
            words[2] = "results";

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if        (dummyVector[i] == "mole" &&
                         dummyVector[i+1] == "fraction")             {checkWord[0] = true; fractionIndex    = i;}
                else if    (dummyVector[i] == "mass" &&
                         dummyVector[i+1] == "fraction")             {checkWord[0] = true; fractionIndex    = i;}
                else if (dummyVector[i] == "inert" &&
                         dummyVector[i+1] == "specie")               {checkWord[1] = true; inertIndex       = i;}
                else if (dummyVector[i] == "results")                {checkWord[2] = true; resultsIndex     = i;}
            }

            for (unsigned int i=0;i<checkWord.size();i++)
            {
                if ( checkWord[i] == false)
                {
                    error();
                    std::cout << "key word || " << words[i] << " || is MISSING in Reaction properties sub-dictionary!\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }
            
            inert_ = dummyVector[inertIndex + 2];

            results_ = dummyVector[resultsIndex + 1];
            if ( results_ != "mass" && results_ != "mole" )
            {
                error();
                std::cout << "key word || " << "Results" << " || MUST be || mass || mole ||\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            inlet_ = dummyVector[fractionIndex];
            if ( inlet_ != "mass" && inlet_ != "mole" )
            {
                error();
                std::cout << "key word || " << "fraction" << " || MUST be || mass fraction || mole fraction ||\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            if ( dummyVector[fractionIndex+2] != "(" )
            {
                error();
                std::cout << "The Reaction properties/mole fraction || Reaction properties/mass fraction sub-dictionary must start with (\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                unsigned int finalCount = 1e05;
                for (unsigned int i=fractionIndex;i<dummyVector.size();i++)
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
                    std::cout << "The Reaction properties/mole fraction || Reaction properties/mass fraction sub-dictionary must finish with )\n" << std::endl;
                    exit (EXIT_FAILURE);
                }

                unsigned int k=0;
                for (unsigned int i=fractionIndex+1+1+1;i<finalCount;i++)
                {
                    if (dummyVector[i] == "(")
                    {
                        error();
                        std::cout << "The Reaction properties/mole fraction || Reaction properties/mass fraction sub-dictionary must finish with )\n" << std::endl;
                        exit (EXIT_FAILURE);
                    }
                    else
                    {
                        inletValue_.resize(k+1);
                        inletName_.resize(k+1);
                        inletValue_[k] = boost::lexical_cast<double>(dummyVector[i+1]);
                        inletName_[k]  = dummyVector[i];
                        k++;
                        i++;
                    }
                }
            }
        }
    }
}
