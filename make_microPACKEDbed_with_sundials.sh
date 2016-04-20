source mybashrc
export ASALI_SUPPORT='-DASALI_USE_SUNDIALS=1 -DASALI_USE_BZZ=0'
mkdir exe
rm -f exe/microPACKEDbed_with_sundials.sh

cd src
g++ -O3 -m64 -Wno-write-strings -fpermissive microPACKEDbed/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -ImicroPACKEDbed -ImicroPACKEDbed/utilities -ImicroPACKEDbed/ode -ImicroPACKEDbed/bvp -I. -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$BOOST/lib/ -lgfortran -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/microPACKEDbed_with_sundialss.sh
cd ..
