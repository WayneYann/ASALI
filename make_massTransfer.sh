source mybashrc
mkdir exe
rm -f exe/massTransfer.sh

cd src
g++ -O3 -m64 -Wno-write-strings -fpermissive massTransfer/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -ImassTransfer/. -ImassTransfer/utilities -ImassTransfer/dae -ImassTransfer/ode  -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/massTransfer.sh
cd ..


