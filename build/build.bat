cmake -H.. -Bwin32_msc -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=fin
cd win32_msc
nmake install
cd ..
