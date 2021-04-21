#! /bin/sh
exe="ForceEquilibrium"
cp ~/../..//build/bin/$exe .
rm *_ChainExtensionDistribution.dat *_CrosslinkPosition.dat
./$exe -i LC_N32_B256_f3.bfm -d BondCreationBreaking.dat -u 0.70 -t 0.00001 -s 0.01
awk '{printf("%0.2f %.4f %d %d \n", $1,$2,$3,$4)> "tmp.dat" }' Modul.dat 
mv tmp.dat Modul.dat

