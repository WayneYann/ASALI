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

if ( input.getDAEsolver() == "BzzMath")
{
    #if ASALI_USE_BZZ == 1
    dae->resize(input.getModelName()[j]);

    BzzVector yFbzz;
    BzzVector x0bzz;
    ChangeDimensions( dae->NumberOfEquations(), &yFbzz, true);
    ChangeDimensions( dae->NumberOfEquations(), &x0bzz, true);
    ChangeDimensions( dae->NumberOfEquations(), &yF,    true);

    FromOSToBzz(x0,x0bzz);

    BzzVectorInt iAO;
    dae->AlgebraicEquations(iAO);

    BzzDaeObject o;
    o(x0bzz,0.,iAO,dae);

    {
        BzzVector yMax(dae->NumberOfEquations());
        if (input.getModelName()[j] == "heterogeneous")
        {
            unsigned int counter = 1;
            for(unsigned int s=1;s<=input.numberOfSpecies();s++)
            {
                yMax[counter++] = 1.;
            }

            for(unsigned int s=1;s<=input.numberOfSpecies();s++)
            {
                yMax[counter++] = 1.;
            }

            yMax[counter++] = 5000.;

            yMax[counter++] = 5000.;
        }
        else if (input.getModelName()[j] == "pseudoHomogeneous")
        {
            unsigned int counter = 1;
            for(unsigned int s=1;s<=input.numberOfSpecies();s++)
            {
                yMax[counter++] = 1.;
            }
            yMax[counter++] = 5000.;
        }

        BzzVector yMin(dae->NumberOfEquations());
        if (input.getModelName()[j] == "heterogeneous")
        {
            unsigned int counter = 1;
            for(unsigned int s=1;s<=input.numberOfSpecies();s++)
            {
                yMin[counter++] = 0.;
            }
            
            for(unsigned int s=1;s<=input.numberOfSpecies();s++)
            {
                yMin[counter++] = 0.;
            }
            
            yMin[counter++] = 0.;
            
            yMin[counter++] = 0.;
        }
        else if (input.getModelName()[j] == "pseudoHomogeneous")
        {
            unsigned int counter = 1;
            for(unsigned int s=1;s<=input.numberOfSpecies();s++)
            {
                yMin[counter++] = 0.;
            }
            yMin[counter++] = 0.;
        }
        o.SetMinimumConstraints(&yMin);
        o.SetMaximumConstraints(&yMax);
    }

    o.SetTolAbs(1.e-12);
    o.SetTolRel(1.e-07);

    dae->start();
    yFbzz = o(input.getReactorLength());

    FromBzzToOS(yFbzz,yF);

    dae->end();
    #endif
}
else if ( input.getDAEsolver() == "Sundials" )
{
    #if ASALI_USE_SUNDIALS == 1
    dae->resize(input.getModelName()[j]);

    OpenSMOKE::OpenSMOKEVectorDouble residuals(dae->NumberOfEquations());
    OpenSMOKE::OpenSMOKEVectorBool   algebraic(dae->NumberOfEquations());

    dae->AlgebraicEquations(algebraic);
    ChangeDimensions( dae->NumberOfEquations(), &yF, true);

    OpenSMOKE::DAESystem_IDA_Template *daeSystemObject;
    daeSystemObject = OpenSMOKE::DAESystem_IDA_Template::GetInstance(dae->NumberOfEquations(),algebraic.GetHandle());
    daeSystemObject->SetDaeSystem(dae);

    OpenSMOKE::OpenSMOKE_IDA_Sundials<OpenSMOKE::DAESystem_IDA_Template>   o(daeSystemObject);
    o.SetDimensions(dae->NumberOfEquations());
    o.SetAbsoluteTolerance(1e-12);
    o.SetRelativeTolerance(1e-07);
    o.SetMaximumNumberOfSteps(100000);
    o.SetAnalyticalJacobian(false);
    o.SetAlgebraicEquations(algebraic.GetHandle());
    o.SetInitialValues(0., x0.GetHandle(), residuals.GetHandle());


    dae->start();
    double tf = input.getReactorLength();
    o.Solve(tf);
    o.Solution(yF.GetHandle(), residuals.GetHandle());
    
    dae->end(]);
    #endif
}
