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
            inline double getReactorLength()       const  { return L_;};
            inline double getSpecificArea()        const  { return av_;};
            inline double getInnerDiameter()       const  { return Dint_;};
            inline double getOuterDiameter()       const  { return Dext_;};
            inline double getHydraulicDiameter()   const  { return Dh_;};
            inline double getAlfa()                const  { return alfa_;};
            inline double getPressure()            const  { return p_;};
            inline double getTemperature()         const  { return T_;};
            inline double getVelocity()            const  { return v_;};

            inline std::string getFeed()           const  { return feed_;};
            inline std::string getReactorType()    const  { return type_;};
            inline std::string getInert()          const  { return inert_;}
            inline std::string getOdeSolver()      const  { return odeSolver_;};
            inline std::string getDaeSolver()      const  { return daeSolver_;};
            inline std::string getKineticsPath()   const  { return kineticsPath_;};

            inline bool getHomogeneousReactions()  const  { return homo_;};
            inline bool getHeterogenousReactions() const  { return het_;};
            inline bool getConstraints()           const  { return constraints_;};

            inline std::vector<double> getFraction()                const { return inletValue_;};

            inline std::vector<std::string> getFractionName()       const { return inletName_;};

            void recapOnScreen();

        private:

            double absSpecieTol_;
            double relSpecieTol_;
            double L_;
            double av_;
            double Dint_;
            double Dext_;
            double Dh_;
            double alfa_;
            double p_;
            double T_;
            double v_;

            std::string feed_;
            std::string type_;
            std::string inert_;
            std::string odeSolver_;
            std::string daeSolver_;
            std::string kineticsPath_;

            bool constraints_;
            bool homo_;
            bool het_;

            const std::string& file_;

            std::vector<double>      inletValue_;

            std::vector<std::string> inputVector_;
            std::vector<std::string> inletName_;

            unsigned int reactorIndex_;
            unsigned int solverIndex_;
            unsigned int catalystIndex_;
            unsigned int fractionIndex_;
            unsigned int numericalIndex_;
            unsigned int kineticsIndex_;

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
            void velocity();
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
        L_                 = 0;
        av_                = 0;
        Dint_              = 0;
        Dext_              = 0;
        Dh_                = 0;
        alfa_              = 0;
        p_                 = 0;
        T_                 = 0;
        v_                 = 0;

        constraints_       = false;
        homo_              = false;
        het_               = false;

        reactorIndex_      = 0;
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
        reactor();
        catalyst();
        fractions();
        pressure();
        temperature();
        velocity();
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
        words[6] = "Volumetric flow rate || Velocity";
        words[7] = "Numerical solvers";
        words[8] = "Kinetics path";

        for (unsigned int i=0;i<inputVector_.size();i++)
        {
            if      (inputVector_[i]   == "Temperature")         {checkWord[0]  = true;}
            else if (inputVector_[i]   == "Pressure")            {checkWord[1]  = true;}
            else if (inputVector_[i]   == "Reactor")             {checkWord[2]  = true; reactorIndex_      = i;}
            else if (inputVector_[i]   == "Catalyst")            {checkWord[3]  = true; catalystIndex_     = i;}
            else if (inputVector_[i]   == "Mole" &&
                     inputVector_[i+1] == "fractions" )          {checkWord[4]  = true; fractionIndex_     = i;}
            else if (inputVector_[i]   == "Mass" &&
                     inputVector_[i+1] == "fractions" )          {checkWord[4]  = true; fractionIndex_     = i;}
            else if (inputVector_[i]   == "Solver" &&
                     inputVector_[i+1] == "options" )            {checkWord[5]  = true; solverIndex_       = i;}
            else if (inputVector_[i]   == "Volumetric" &&
                     inputVector_[i+1] == "flow" &&
                     inputVector_[i+2] == "rate" )               {checkWord[6]  = true;}
            else if (inputVector_[i]   == "Velocity")            {checkWord[6]  = true;}
            else if (inputVector_[i]   == "Numerical" &&
                     inputVector_[i+1] == "solvers" )            {checkWord[7]  = true; numericalIndex_    = i;}
            else if (inputVector_[i]   == "Kinetics" &&
                     inputVector_[i+1] == "path" )               {checkWord[8]  = true; kineticsIndex_     = i;}
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

            std::vector<bool>        checkWord(4);
            std::vector<std::string> words(4);

            double absIndex;
            double relIndex;
            double reactionIndex;
            double constraintIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "Absolute tollerance";
            words[1] = "Relative tollerance";
            words[2] = "Reactions";
            words[3] = "Constraints";

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if      (dummyVector[i] == "Absolute" &&
                         dummyVector[i+1] == "tollerance")         {checkWord[0] = true; absIndex         = i;}
                else if (dummyVector[i] == "Relative" &&
                         dummyVector[i+1] == "tollerance")         {checkWord[1] = true; relIndex         = i;}
                else if (dummyVector[i] == "Reactions")            {checkWord[2] = true; reactionIndex    = i;}
                else if (dummyVector[i] == "Constraints")          {checkWord[3] = true; constraintIndex  = i;}
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
                            else
                            {
                                error();
                                std::cout << "key word || specie || is MISSING in Solver options/Absolute tollerance sub-dictionary!\n" << std::endl;
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
                            else
                            {
                                error();
                                std::cout << "key word || specie || is MISSING in Solver options/Relative tollerance sub-dictionary!\n" << std::endl;
                                exit (EXIT_FAILURE);
                            }
                        }
                    }
                }
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
            
            std::vector<bool>        checkWord(4);
            std::vector<std::string> words(4);

            double typeIndex;
            double lengthIndex;
            double DintIndex;
            double DextIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;


            words[0] = "type";
            words[1] = "length";
            words[2] = "inner diameter";
            words[3] = "outer diameter";

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if      (dummyVector[i] == "type")              {checkWord[0] = true; typeIndex        = i;}
                else if (dummyVector[i] == "length")            {checkWord[1] = true; lengthIndex      = i;}
                else if (dummyVector[i] == "inner" &&
                         dummyVector[i+1] == "diameter")        {checkWord[2] = true; DintIndex        = i;}
                else if (dummyVector[i] == "outer" &&
                         dummyVector[i+1] == "diameter")        {checkWord[3] = true; DextIndex        = i;}
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
            
            {
                L_ = boost::lexical_cast<double>(dummyVector[lengthIndex+1]);
                std::string dim = dummyVector[lengthIndex+2];
                ConvertsToMeter(L_, dim);
            }
            
            {
                Dint_ = boost::lexical_cast<double>(dummyVector[DintIndex+2]);
                std::string dim = dummyVector[DintIndex+3];
                ConvertsToMeter(Dint_, dim);
            }

            {
                Dext_ = boost::lexical_cast<double>(dummyVector[DextIndex+2]);
                std::string dim = dummyVector[DextIndex+3];
                ConvertsToMeter(Dext_, dim);
            }

            type_ = dummyVector[typeIndex+1];
            if ( type_ != "annular")
            {
                error();
                std::cout << "key word ||" << " type " << "|| MUST be || annular || \n" << std::endl;
                exit (EXIT_FAILURE);
            }
            
            if ( Dint_ >= Dext_ )
            {
                error();
                std::cout << "key word || outer diameter || MUST be bigger than || inner dimeter || \n" << std::endl;
                exit (EXIT_FAILURE);
            }
            
            av_ = 4.*Dint_/(std::pow(Dext_,2.) - std::pow(Dint_,2.));
            Dh_ = (Dint_ + Dext_)*0.5;
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

                    double ReactorVolume = 3.14*0.25*(std::pow(Dext_,2.) - std::pow(Dint_,2.))*L_;

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

    void READinput::kinetics()
    {
        kineticsPath_ = inputVector_[kineticsIndex_+2];
    }

    void READinput::velocity()
    {
        for (unsigned int i=0;i<inputVector_.size();i++)
        {
            if (inputVector_[i] == "Velocity")
            {
                v_ = boost::lexical_cast<double>(inputVector_[i+1]);
                std::string dim = inputVector_[i+2];
                ConvertsToMeterPerSecond(v_,dim);
            }
            else if (inputVector_[i]   == "Volumetric" &&
                     inputVector_[i+1] == "flow" &&
                     inputVector_[i+2] == "rate" )
            {
                v_ = boost::lexical_cast<double>(inputVector_[i+3]);
                std::string dim = inputVector_[i+4];
                ConvertsToNm3perSecond(v_, dim);
                double Ain = 3.14*0.25*(std::pow(Dext_,2.) - std::pow(Dint_,2.));
                double T0 = 273.15;
                double P0 = 101325.;
                double nR = v_*P0/T0;
                double vnew = nR*T_/p_;
                v_ = vnew/Ain;
            }
        }
    }

    void READinput::recapOnScreen()
    {
        std::cout.precision(6);
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                          GENERAL INPUT                                         \n" << std::endl;
        std::cout << "Reactor type:                              " << type_ << std::endl;
        std::cout << "Reactor lenght                           = " << L_ << "\t[m]" << std::endl;
        std::cout << "Inner   diameter                         = " << Dint_ << "\t[m]" << std::endl;
        std::cout << "Outer   diameter                         = " << Dext_ << "\t[m]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                       REACTOR PROPERTIES                                       \n" << std::endl;
        std::cout << "Specific area (av)                       = " << av_ << "\t[1/m]" << std::endl;
        std::cout << "Specific catalytic area (alfa)           = " << alfa_ << "\t[1/m]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                      OPERATING CONDITIONS                                      \n" << std::endl;
        std::cout << "Feed velocity                            = " << v_ << "\t[m/s]" << std::endl;
        std::cout << "Feed pressure                            = " << p_ << "\t[Pa]" << std::endl;
        std::cout << "Feed temperature                         = " << T_ << "\t[K]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                         SOLVER OPTIONS                                         \n" << std::endl;
        std::cout << "Homogeneous  reactions:                    " << boolOnScreen(homo_) << std::endl;
        std::cout << "Heterogeneus reactions:                    " << boolOnScreen(het_) << std::endl;
        std::cout << "Feed in:                                   " << feed_ << " fractions" << std::endl;
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
