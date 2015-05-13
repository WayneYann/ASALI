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
    class POSTprocessing
    {
        public:

            POSTprocessing(unsigned int& NC, unsigned int& SURF_NC);
            
            void resize(const unsigned int NP, const unsigned int NS);

            void setVariables  (const OpenSMOKE::OpenSMOKEVectorDouble y);
            void setGrid       (const OpenSMOKE::OpenSMOKEVectorDouble z);
            void setSampledGrid(const std::vector<double> sampled);

            void sampling();

            OpenSMOKE::OpenSMOKEVectorDouble getSampledVariables()   { return sampledVariables_;};

            ~POSTprocessing(){ };

        private:

            unsigned int NP_;
            unsigned int NS_;
            unsigned int NEs_;

            const unsigned int& NC_;
            const unsigned int& SURF_NC_;

            OpenSMOKE::OpenSMOKEVectorDouble* xBulk_;
            OpenSMOKE::OpenSMOKEVectorDouble* xWall_;
            OpenSMOKE::OpenSMOKEVectorDouble* teta_;
            OpenSMOKE::OpenSMOKEVectorDouble* sampledxBulk_;
            OpenSMOKE::OpenSMOKEVectorDouble* sampledxWall_;
            OpenSMOKE::OpenSMOKEVectorDouble* sampledteta_;

            OpenSMOKE::OpenSMOKEVectorDouble Tbulk_;
            OpenSMOKE::OpenSMOKEVectorDouble Twall_;
            OpenSMOKE::OpenSMOKEVectorDouble sampledTbulk_;
            OpenSMOKE::OpenSMOKEVectorDouble sampledTwall_;
            OpenSMOKE::OpenSMOKEVectorDouble z_;
            OpenSMOKE::OpenSMOKEVectorDouble sampledVariables_;
            OpenSMOKE::OpenSMOKEVectorDouble sampled_;

            void error() { std::cout << "\nASALI::READinput::ERROR\n" << std::endl;};
    };
    
    POSTprocessing::POSTprocessing(unsigned int& NC, unsigned int& SURF_NC):
    NC_(NC),
    SURF_NC_(SURF_NC)
    {
        NP_             = 0;
        NS_             = 0;
        NEs_            = 0;
    }

    void POSTprocessing::resize(const unsigned int NP, const unsigned int NS)
    {
        NP_   = NP;
        NS_   = NS;
        NEs_  = (NC_ + NC_ + SURF_NC_ + 1 + 1)*NS_;

        xBulk_         = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
        xWall_         = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
        teta_          = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
        sampledxBulk_  = new OpenSMOKE::OpenSMOKEVectorDouble[NS_];
        sampledxWall_  = new OpenSMOKE::OpenSMOKEVectorDouble[NS_];
        sampledteta_   = new OpenSMOKE::OpenSMOKEVectorDouble[NS_];

        for (unsigned int i=0;i<NP_;i++)
        {
            ChangeDimensions(NC_,      &xBulk_[i],        true);
            ChangeDimensions(NC_,      &xWall_[i],        true);
            ChangeDimensions(SURF_NC_, &teta_[i],         true);
        }
        
        for (unsigned int i=0;i<NS_;i++)
        {
            ChangeDimensions(NC_,      &sampledxBulk_[i], true);
            ChangeDimensions(NC_,      &sampledxWall_[i], true);
            ChangeDimensions(SURF_NC_, &sampledteta_[i],  true);
        }
        
        ChangeDimensions(NP_, &Tbulk_,        true);
        ChangeDimensions(NP_, &Twall_,        true);
        ChangeDimensions(NS_, &sampledTbulk_, true);
        ChangeDimensions(NS_, &sampledTwall_, true);

        ChangeDimensions(NEs_, &sampledVariables_, true);
    }

    void POSTprocessing::setVariables(const OpenSMOKE::OpenSMOKEVectorDouble  y)
    {
        unsigned int counter = 1;
        for (unsigned int i=0;i<NP_;i++)
        {
            for (unsigned int j=1;j<=NC_;j++)
                xBulk_[i][j] = y[counter++];
            for (unsigned int j=1;j<=NC_;j++)
                xWall_[i][j] = y[counter++];
            for (unsigned int j=1;j<=SURF_NC_;j++)
                teta_[i][j] = y[counter++];
            Tbulk_[i+1] = y[counter++];
            Twall_[i+1] = y[counter++];
        }
    }

    void POSTprocessing::setGrid(const OpenSMOKE::OpenSMOKEVectorDouble z)
    {
        ChangeDimensions(z.Size(), &z_, true);
        for (unsigned int k=1;k<=z_.Size();k++)
            z_[k] = z[k];
    }

    void POSTprocessing::setSampledGrid(std::vector<double> sampled)
    {
        ChangeDimensions(sampled.size(), &sampled_, true);
        for (unsigned int k=1;k<=sampled_.Size();k++)
            sampled_[k] = sampled[k-1];
    }

    void POSTprocessing::sampling()
    {
        std::vector<double> dzVector(NP_);
        double min;
        unsigned int counter = 1;
        for (unsigned int k=1;k<=NS_;k++)
        {
            for (unsigned int i=1;i<=NP_;i++)
            {
                dzVector[i-1] = std::fabs(sampled_[k] - z_[i]);
            }

            min = dzVector[0];
            for (unsigned int i=1;i<NP_;i++)
            {
                min = std::min(dzVector[i],min);
            }

            for (unsigned int i=1;i<=NP_;i++)
            {
                if ( std::fabs(sampled_[k] - z_[i]) == min)
                {
                    for (unsigned int j=1;j<=NC_;j++)
                        sampledxBulk_[k-1][j] = xBulk_[i-1][j];
                    for (unsigned int j=1;j<=NC_;j++)
                        sampledxWall_[k-1][j] = xWall_[i-1][j];
                    for (unsigned int j=1;j<=SURF_NC_;j++)
                        sampledteta_[k-1][j] = teta_[i-1][j];
                    sampledTbulk_[k] = Tbulk_[i];
                    sampledTwall_[k] = Twall_[i];

                    break;
                }
            }

            for (unsigned int j=1;j<=NC_;j++)
                sampledVariables_[counter++] = sampledxBulk_[k-1][j];
            for (unsigned int j=1;j<=NC_;j++)
                sampledVariables_[counter++] = sampledxWall_[k-1][j];
            for (unsigned int j=1;j<=SURF_NC_;j++)
                sampledVariables_[counter++] = sampledteta_[k-1][j];
            sampledVariables_[counter++] = sampledTbulk_[k];
            sampledVariables_[counter++] = sampledTwall_[k];
        }
    }
}
