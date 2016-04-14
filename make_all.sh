source mybashrc
export ASALI_SUPPORT='-DASALI_USE_SUNDIALS=1 -DASALI_USE_BZZ=1'
mkdir exe
rm -f exe/heterogeneousPFR.sh
rm -f exe/idealPFR_heatTransfer.sh
rm -f exe/idealPFR_massTransfer.sh
rm -f exe/idealPFR.sh
rm -f exe/idealBATCH.sh
rm -f exe/catalyticSLAB.sh
rm -f exe/microPACKEDbed.sh

cd src
g++ -O3 -m64 -Wno-write-strings -fpermissive heterogeneousPFR/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IheterogeneousPFR -IheterogeneousPFR/utilities -IheterogeneousPFR/ode -IheterogeneousPFR/bvp -IheterogeneousPFR/post -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/heterogeneousPFR.sh

g++ -O3 -m64 -Wno-write-strings -fpermissive idealPFR_heatTransfer/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IidealPFR_heatTransfer -IidealPFR_heatTransfer/utilities -IidealPFR_heatTransfer/ode -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/idealPFR_heatTransfer.sh

g++ -O3 -m64 -Wno-write-strings -fpermissive idealPFR/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IidealPFR -IidealPFR/utilities -IidealPFR/ode -IidealPFR/dae -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/idealPFR.sh

g++ -O3 -m64 -Wno-write-strings -fpermissive idealPFR_massTransfer/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IidealPFR_massTransfer -IidealPFR_massTransfer/utilities -IidealPFR_massTransfer/ode -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/idealPFR_massTransfer.sh

g++ -O3 -m64 -Wno-write-strings -fpermissive idealBATCH/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IidealBATCH -IidealBATCH/utilities -IidealBATCH/ode -IidealBATCH/ic -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/idealBATCH.sh

g++ -O3 -m64 -Wno-write-strings -fpermissive catalyticSLAB/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IcatalyticSLAB -IcatalyticSLAB/utilities -IcatalyticSLAB/ode -IcatalyticSLAB/bvp -IcatalyticSLAB/post -I. I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/catalyticSLAB.sh

g++ -O3 -m64 -Wno-write-strings -fpermissive microPACKEDbed/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -ImicroPACKEDbed -ImicroPACKEDbed/utilities -ImicroPACKEDbed/ode -ImicroPACKEDbed/bvp -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/microPACKEDbed.sh

cd ..
