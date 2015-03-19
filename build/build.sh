cmake -H.. -Blinux -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=fin
cd linux
make install
cd ..