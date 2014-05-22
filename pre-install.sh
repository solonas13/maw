#! /bin/sh

sudo apt-get install g++
sudo apt-get install cmake
sudo apt-get install zlib1g-dev
wget https://libdivsufsort.googlecode.com/files/libdivsufsort-2.0.1.tar.bz2
tar -xvf libdivsufsort-2.0.1.tar.bz2
rm libdivsufsort-2.0.1.tar.bz2
cd libdivsufsort-2.0.1/
mkdir build
cd build
sudo cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX="/usr/local" ..
sudo make
cd ..
cd ..
sudo rm -r libdivsufsort-2.0.1/
wget https://github.com/simongog/sdsl/archive/master.zip
unzip master.zip
sudo rm master.zip
cd sdsl-master/
sudo ./install.sh $HOME/sdsl
cd ..
sudo rm -r sdsl-master/
