touch massTransfer.txt

echo '/*##############################################################################################' >> massTransfer.txt
echo '#                                                                                              #' >> massTransfer.txt
echo '#     #############       #############       #############       ####                ####     #' >> massTransfer.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> massTransfer.txt
echo '#    #    #####    #     #    #########      #    #####    #     #    #              #    #    #' >> massTransfer.txt
echo '#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #' >> massTransfer.txt
echo '#    #    #####    #     #    #              #    #####    #     #    #              #    #    #' >> massTransfer.txt
echo '#    #             #     #    #########      #             #     #    #              #    #    #' >> massTransfer.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> massTransfer.txt
echo '#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #' >> massTransfer.txt
echo '#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #' >> massTransfer.txt
echo '#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #' >> massTransfer.txt
echo '#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #' >> massTransfer.txt
echo '#     ####     ####       #############       ####     ####       #############       ####     #' >> massTransfer.txt
echo '#                                                                                              #' >> massTransfer.txt
echo '#   Department of Energy                                                                       #' >> massTransfer.txt
echo '#   Politecnico di Milano                                                                      #' >> massTransfer.txt
echo '#   Author: Stefano Rebughini <stefano.rebughini@polimi.it>                                    #' >> massTransfer.txt
echo '#                                                                                              #' >> massTransfer.txt
echo '################################################################################################' >> massTransfer.txt
echo ' ' >> massTransfer.txt
echo 'Kinetics path           KineticsSchemes/CH4_UBI/kinetics' >> massTransfer.txt
echo ' ' >> massTransfer.txt
echo 'Reactor type            PackedBed' >> massTransfer.txt
echo ' ' >> massTransfer.txt 
echo 'Operating conditions' >> massTransfer.txt
echo '{' >> massTransfer.txt
echo '  feed temperature            298 K' >> massTransfer.txt
echo '  feed pressure               1.2 bar' >> massTransfer.txt
echo '  velocity                    1.5 m/s' >> massTransfer.txt
echo '}' >> massTransfer.txt
echo ' ' >> massTransfer.txt
echo 'Reactor properties' >> massTransfer.txt
echo '{' >> massTransfer.txt
echo '  reactor length        1.0  cm' >> massTransfer.txt
echo '  shell diameter        3.0  mm' >> massTransfer.txt
echo '  particle diameter     0.6  mm' >> massTransfer.txt
echo '}' >> massTransfer.txt
echo ' ' >> massTransfer.txt
echo 'Reaction properties' >> massTransfer.txt
echo '{' >> massTransfer.txt
echo '  mole fraction' >> massTransfer.txt
echo '  (' >> massTransfer.txt
echo '      CO     0.003' >> massTransfer.txt
echo '      O2     0.007' >> massTransfer.txt
echo '      N2     0.99' >> massTransfer.txt
echo '  )' >> massTransfer.txt
echo ' ' >> massTransfer.txt
echo '  inert specie  N2' >> massTransfer.txt
echo ' ' >> massTransfer.txt
echo '  results           mass' >> massTransfer.txt
echo '}' >> massTransfer.txt
echo ' ' >> massTransfer.txt
echo 'Solver options' >> massTransfer.txt
echo '{' >> massTransfer.txt
echo '  pressure drops    Eisfeld' >> massTransfer.txt
echo ' ' >> massTransfer.txt
echo '  mass transfer     Yoshida' >> massTransfer.txt
echo ' ' >> massTransfer.txt
echo '  numerical solver  OpenSMOKE' >> massTransfer.txt
echo '}' >> massTransfer.txt
