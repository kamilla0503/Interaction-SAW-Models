resultFolder=data/XY_MC/For_EE_3D11/  #SaveBeforeOptimizeOut40
rm -rf $resultFolder
mkdir $resultFolder

g++ *.cpp -O3 -I include -o testMC -std=c++17 -DREGIME_3D #-DSTARTHALF

for L in 4 5 6 7
do
  ./testMC $L > data/XY_MC/For_EE_3D/XY_Original_EE_$L.csv
done
wait
