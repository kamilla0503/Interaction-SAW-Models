resultFolder=SaveBeforeOptimize_3D   #data/XY_MC/For_EE_3D  #SaveBeforeOptimizeOut40
#rm -rf $resultFolder
#mkdir $resultFolder

#g++ *.cpp -O3 -I include -o testMC -std=c++17 -DREGIME_2D -DSTARTHALF
g++ *.cpp -O3 -I include -o testMC -std=c++17 -DREGIME_3D -DSTARTHALF

for L in 15 20 #4 5 6 7 #5 6 7 8 9 10 11
do
 # for J in 0.45 0.55 0.6 0.7 #0 0.25 0.5 0.75 1.0 1.25 1.5 #0 0.3 0.6 0.9 1.2#
  for J in 0.46 0.48 0.49 0.51 0.52 0.53 #0 0.25 0.5 0.75 1.0 1.25 1.5 #0 0.3 0.6 0.9 1.2#
  do
  ./testMC $L $J $resultFolder &
  done
#> data/XY_MC/For_EE/XY_Original_EE_$L.csv
done
wait
