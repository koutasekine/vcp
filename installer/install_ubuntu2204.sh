#!/bin/bash
# sudo apt update
# yes | sudo apt upgrade
yes | sudo apt install build-essential
yes | sudo apt install libboost-all-dev
yes | sudo apt install libgmp-dev 
yes | sudo apt install libmpfr-dev
yes | sudo apt install curl
yes | sudo apt install wget

wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | sudo gpg --dearmor | sudo  tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt update
yes | sudo apt install intel-basekit
yes | sudo apt install intel-hpckit
echo "source /opt/intel/oneapi/setvars.sh" >> ~/.bashrc
source /opt/intel/oneapi/setvars.sh


echo 
echo "######################################################"
echo "######################################################"
echo "Install kv library:"
echo "Create a new folder in Home directory. "
echo -n "Please input a new folder name: "
read foldername
folderpath="${HOME}/${foldername}"

echo "Check: ${folderpath}"

mkdir "${folderpath}"
cd "${folderpath}"

kv=`curl http://verifiedby.me/kv/ | grep -o -E "download/(kv-[0-9]+\.[0-9]+\.[0-9]*\.tar\.gz)"`
kvver=${kv#*/}
kvver=${kvver%.tar.gz}

urlkv="http://verifiedby.me/kv/${kv}"
wget "${urlkv}"

if [ $? -ne 0 ]; then
    echo "[ERROR] Could not download..."
    exit 1
else
    tar -xvf "${kvver}.tar.gz"
    cp -r "${kvver}/kv/" ./
    cp -r "${kvver}/test/" ./
    cp -r "${kvver}/example/" ./
    rm "${kvver}.tar.gz"
    rm -r "${kvver}"
fi

echo "######################################################"
echo "Install VCP library:"
wget --no-check-certificate https://github.com/koutasekine/vcp/archive/master.tar.gz
if [ $? -ne 0 ]; then
    echo "[ERROR] Could not download..."
    exit 1
else
    tar -xvf master.tar.gz
    cp -r vcp-master/vcp/ ./
    cp -r vcp-master/test_matrix/ ./
    cp -r vcp-master/test_PDE/ ./
    cp -r vcp-master/vblas/ ./
    rm "master.tar.gz"
    rm -r vcp-master    
fi

echo 
echo "######################################################"
echo "Check for BLAS rounding mode changes: Please wait..."
source /opt/intel/mkl/bin/mklvars.sh intel64
cd "${folderpath}/test_matrix/"
g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 Check_pdblas_rounding.cpp -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lmpfr -fopenmp && ./a.out
g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 Check_pidblas_rounding.cpp -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lmpfr -fopenmp && ./a.out


echo
echo "######################################################"
echo "Check for Open MP rounding mode changes: Please wait..."
g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 Check_OpenMP.cpp -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lmpfr -fopenmp && ./a.out

echo "Finish!"