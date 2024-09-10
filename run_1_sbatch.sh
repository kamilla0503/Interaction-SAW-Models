resultFolder=data/XY_MC/Run_Start/3D/
g++ *.cpp -O3 -I include -o testMC -std=c++17 -DREGIME_3D -DSTARTHALF -DSEED

for L in 100 200 300 #4 5 6 7 #5 6 7 8 9 10 11
do
 # for J in 0.45 0.55 0.6 0.7 #0 0.25 0.5 0.75 1.0 1.25 1.5 #0 0.3 0.6 0.9 1.2#
  for J in 0.25 0.23 0.21 0.19 #0.31 0.29 0.27 #0.43 0.4 0.37 0.34 #0.46 0.49 0.51 0.53 #0 0.25 0.5 0.75 1.0 1.25 1.5 #0 0.3 0.6 0.9 1.2#
  do
  sbatch --time=0-12:0 --wrap="./testMC $L $J $resultFolder"
  done
#> data/XY_MC/For_EE/XY_Original_EE_$L.csv
done

