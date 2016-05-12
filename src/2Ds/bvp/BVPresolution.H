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

if ( input.getBVPsolver() == "BzzMath")
{
    #if ASALI_USE_BZZ == 1
    bvp->resize(input.getModelName()[j]);

    BzzVector yFbzz;
    BzzVector x0bzz;
    ChangeDimensions( bvp->NumberOfEquations(), &yFbzz, true);
    ChangeDimensions( bvp->NumberOfEquations(), &x0bzz, true);
    ChangeDimensions( bvp->NumberOfEquations(), &yF,    true);

    FromOSToBzz(x0,x0bzz);

    BzzVectorInt iAO;
    bvp->AlgebraicEquations(iAO);

    BzzDaeObject o;
    o(x0bzz,0.,iAO,bvp);

    {
        BzzVector yMax(ode->NumberOfEquations());
        if (input.getModelName()[j] == "honeyComb")
        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<input.axialPoints();a++)
            {
                for (unsigned int r=0;r<input.radialPoints();r++)
                {
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMax[counter++] = 1.;
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMax[counter++] = 1.;
                    yMax[counter++] = 5000.;
                    yMax[counter++] = 5000.;
                }
            }
        }
        else if (input.getModelName()[j] == "packedBed")
        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<input.axialPoints();a++)
            {
                for (unsigned int r=0;r<input.radialPoints();r++)
                {
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMax[counter++] = 1.;
                    yMax[counter++] = 5000.;
                }
            }
        }
        else if (input.getModelName()[j] == "microBed")
        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<input.axialPoints();a++)
            {
                for (unsigned int r=0;r<input.radialPoints();r++)
                {
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMax[counter++] = 1.;
                    yMax[counter++] = 5000.;
                    yMax[counter++] = 5000.;
                }
            }
        }
        BzzVector yMin(ode->NumberOfEquations());
        if (input.getModelName()[j] == "honeyComb")
        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<input.axialPoints();a++)
            {
                for (unsigned int r=0;r<input.radialPoints();r++)
                {
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMin[counter++] = 0.;
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMin[counter++] = 0.;
                    yMin[counter++] = 0.;
                    yMin[counter++] = 0.;
                }
            }
        }
        else if (input.getModelName()[j] == "packedBed")
        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<input.axialPoints();a++)
            {
                for (unsigned int r=0;r<input.radialPoints();r++)
                {
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMin[counter++] = 0.;
                    yMin[counter++] = 0.;
                }
            }
        }
        else if (input.getModelName()[j] == "microBed")
        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<input.axialPoints();a++)
            {
                for (unsigned int r=0;r<input.radialPoints();r++)
                {
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMin[counter++] = 0.;
                    yMin[counter++] = 0.;
                    yMin[counter++] = 0.;
                }
            }
        }
        o.SetMinimumConstraints(&yMin);
        o.SetMaximumConstraints(&yMax);
    }

    o.SetTolAbs(1.e-12);
    o.SetTolRel(1.e-09);

    bvp->start(input.getModelName()[j]);
    yFbzz = o(1e06);

    FromBzzToOS(yFbzz,yF);

    bvp->end(input.getModelName()[j]);
    #endif
}
else if ( input.getBVPsolver() == "Sundials" )
{
    #if ASALI_USE_SUNDIALS == 1
    bvp->resize(input.getModelName()[j]);

    OpenSMOKE::OpenSMOKEVectorDouble residuals(bvp->NumberOfEquations());
    OpenSMOKE::OpenSMOKEVectorBool   algebraic(bvp->NumberOfEquations());

    bvp->AlgebraicEquations(algebraic);
    ChangeDimensions( bvp->NumberOfEquations(), &yF, true);

    OpenSMOKE::DAESystem_IDA_Template *daeSystemObject;
    daeSystemObject = OpenSMOKE::DAESystem_IDA_Template::GetInstance(bvp->NumberOfEquations(),algebraic.GetHandle());
    daeSystemObject->SetDaeSystem(bvp);

    OpenSMOKE::OpenSMOKE_IDA_Sundials<OpenSMOKE::DAESystem_IDA_Template>   o(daeSystemObject);
    o.SetDimensions(bvp->NumberOfEquations());
    o.SetAbsoluteTolerance(1e-12);
    o.SetRelativeTolerance(1e-07);
    o.SetMaximumNumberOfSteps(100000);
    o.SetAnalyticalJacobian(false);
    o.SetAlgebraicEquations(algebraic.GetHandle());
    o.SetInitialValues(0., x0.GetHandle(), residuals.GetHandle());


    bvp->start(input.getModelName()[j]);
    double tf = 1e6;
    o.Solve(tf);
    o.Solution(yF.GetHandle(), residuals.GetHandle());
    
    bvp->end(input.getModelName()[j]);
    #endif
}
