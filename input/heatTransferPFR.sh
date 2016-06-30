touch heatTransfer.txt

echo '/*##############################################################################################' >> heatTransfer.txt
echo '#                                                                                              #' >> heatTransfer.txt
echo '#     #############       #############       #############       ####                ####     #' >> heatTransfer.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> heatTransfer.txt
echo '#    #    #####    #     #    #########      #    #####    #     #    #              #    #    #' >> heatTransfer.txt
echo '#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #' >> heatTransfer.txt
echo '#    #    #####    #     #    #              #    #####    #     #    #              #    #    #' >> heatTransfer.txt
echo '#    #             #     #    #########      #             #     #    #              #    #    #' >> heatTransfer.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> heatTransfer.txt
echo '#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #' >> heatTransfer.txt
echo '#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #' >> heatTransfer.txt
echo '#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #' >> heatTransfer.txt
echo '#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #' >> heatTransfer.txt
echo '#     ####     ####       #############       ####     ####       #############       ####     #' >> heatTransfer.txt
echo '#                                                                                              #' >> heatTransfer.txt
echo '#   Department of Energy                                                                       #' >> heatTransfer.txt
echo '#   Politecnico di Milano                                                                      #' >> heatTransfer.txt
echo '#   Author: Stefano Rebughini <stefano.rebughini@polimi.it>                                    #' >> heatTransfer.txt
echo '#                                                                                              #' >> heatTransfer.txt
echo '################################################################################################' >> heatTransfer.txt
echo ' ' >> heatTransfer.txt
echo 'Transport properties path           KineticsSchemes/CH4_UBI/kinetics' >> heatTransfer.txt
echo ' ' >> heatTransfer.txt
echo 'Reactor type            PackedBed' >> heatTransfer.txt
echo ' ' >> heatTransfer.txt 
echo 'Operating conditions' >> heatTransfer.txt
echo '{' >> heatTransfer.txt
echo '  specie name                 CO2' >> heatTransfer.txt
echo '  feed temperature            298 K' >> heatTransfer.txt
echo '  feed pressure               1.2 bar' >> heatTransfer.txt
echo '  solid temperature           773 K' >> heatTransfer.txt
echo '  velocity                    1.5 m/s' >> heatTransfer.txt
echo '}' >> heatTransfer.txt
echo ' ' >> heatTransfer.txt
echo 'Reactor properties' >> heatTransfer.txt
echo '{' >> heatTransfer.txt
echo '  reactor length        1.0  cm' >> heatTransfer.txt
echo '  shell diameter        3.0  mm' >> heatTransfer.txt
echo '  particle diameter     0.6  mm' >> heatTransfer.txt
echo '}' >> heatTransfer.txt
echo ' ' >> heatTransfer.txt
echo 'Solver options' >> heatTransfer.txt
echo '{' >> heatTransfer.txt
echo '  thermophysical properties    sutherland' >> heatTransfer.txt
echo '  pressure drops               Ergun' >> heatTransfer.txt
echo '  heat transfer                Yoshida' >> heatTransfer.txt
echo '}' >> heatTransfer.txt
echo ' ' >> heatTransfer.txt
echo 'Numerical solvers              OpenSMOKE' >> heatTransfer.txt
