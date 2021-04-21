#! /bin/sh
###############################################################################
###############################################################################
###############################################################################
#simple function extracting substring in between an indicator and an delimiter
function getValue
{
local string=$1
local indicator=$2
local delimiter=_
local errorstatus
if [ "$#" -lt 2 ]; then
    errorstatus=0
elif [ "$#" -gt 2 ];then
    delimiter=$3
fi
if [ ! -z $errorstatus ]; then  
    echo $errorstatus
else 
    local tmpString=${string/*$indicator}
    local value=${tmpString/$delimiter*}
    echo $value
fi
}
###############################################################################
###############################################################################
###############################################################################
#calculate the modulus from the config and the bondcreation table
exe="ForceEquilibrium"
cp ~/../..//build/bin/$exe .
rm *_ChainExtensionDistribution.dat *_CrosslinkPosition.dat
./$exe -i LC_N32_B256_f3.bfm -d BondCreationBreaking.dat -u 0.70 -t 0.00001 -s 0.01
awk '{printf("%0.2f %.4f %d %d \n", $1,$2,$3,$4)> "tmp.dat" }' Modul.dat
mv tmp.dat Modul.dat
###############################################################################
###############################################################################
###############################################################################
#calculate the chan extension distribition
doStep=1
if [ "$doStep" == "1" ]; then 
    for file in *_ChainExtensionDistribution.dat ; do 
        awk -v size=2   '{b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; 
            bmin=b<bmin?b:bmin } 
            END { for(i=bmin;i<=bmax;++i) 
            print i*size,(i+1)*size,a[i] > "tmp.dat" }' $file
        mv tmp.dat ${file/.dat}"_histogram.dat"
    done 
fi 
doStep=1
if [ "$doStep" == "1" ]; then 
    gnu="ChainLength_Histogram.gp"
    out="ChainLength_Histogram.eps"
    echo "set terminal postscript enhanced color 'Helvetica, 24' eps  " > $gnu
    echo "set output '$out' " >> $gnu

    echo "b=2.688" >> $gnu
    echo "N=32" >> $gnu
    echo "Re=b*N**0.5" >> $gnu
    echo "set logscale y" >> $gnu
    echo "set ylabel 'absolute'" >> $gnu
    echo "set xlabel 'R'" >> $gnu
    echo "set label at 5,2E4 ' f=3, N=32'" >> $gnu
    echo "plot \\" >> $gnu
    for p in  0.79 0.81 0.83 0.85 0.87 0.89 0.91 0.93 0.95 0.97   ; do 
        for file in $(ls *ChainExtensionDistribution_histogram.dat | grep $p ); do 
            echo $file 
            c=$(getValue AA$file AAC _C )
            echo $c
            echo "'$file' u 1:3 w lp lw 2 title sprintf('p=%.3f',$c) ,\\" >> $gnu
        done 
    done 
    gnuplot $gnu
fi
###############################################################################
###############################################################################
###############################################################################
#plot the modulus
# data can be compared to Figure in 
# Figure 2 : 
# Lang, Michael, and T. Müller. "Analysis of the gel point of polymer model networks by computer simulations." Macromolecules 53.2 (2020): 498-512.
base=Modul
doStep=1
if [ "$doStep" == "1" ]; then 
    output=Modul.dat
    rm $output
    for file in *_ChainExtensionDistribution.dat ; do 
        c=$(getValue AA$file AAC _C )
        awk -v c=$c 'BEGIN{nEffective=0;nAll=0}
            {if($6>0){nEffective++};R2+=$6*$6; nAll++ }
            END{print c, R2, nEffective, nAll > "tmp.dat"}' $file 
        cat tmp.dat >> $output
    done
    sort -n -k1 $output > tmp.dat 
    mv tmp.dat $output 
fi
doStep=1
if [ "$doStep" == "1" ]; then 
    N=32
    gnu="$base.gp"
    out="$base.eps"
    echo "set terminal postscript enhanced color 'Helvetica, 24' eps  " > $gnu
    echo "set output '$out' " >> $gnu
    echo "b=2.68" >> $gnu
    echo "N=1" >> $gnu
    echo "set logscale yx" >> $gnu
    echo "set ylabel 'G/({/Symbol n}kT)={/Symbol G}=<R^2>/(Nb^2)'" >> $gnu
    echo "set xlabel '{/Symbol e}'" >> $gnu
    echo "pc=0.73" >> $gnu
    echo "set title 'modulus per chain vs rel. extent of reaction for N=32 f=3'" >> $gnu  
    echo "xsizeMin=0.002" >> $gnu
    echo "xsizeMax=0.5" >> $gnu
    echo "ysizeMin=0.0001" >> $gnu
    echo "ysizeMax=1." >> $gnu
    echo "set key bottom right " >> $gnu
    nChains=$(grep number_of_linear_chains LastConfig.bfm)
    nChains=${nChains/'#!number_of_linear_chains='}
    nMons=$(grep '!number_of_monomers=' LastConfig.bfm) 
    nMons=${nMons/'!number_of_monomers='}
   
    echo "plot [xsizeMin:xsizeMax][ysizeMin:ysizeMax]\\" >> $gnu
    # factor 2 for the density of phi=0.5
    echo "'$base.dat' u ((\$1-pc)/pc):(\$2/2./(b*b)/$N/\$3) w lp lw 2 lt 7 lc 7 notitle  " >> $gnu
    gnuplot $gnu
fi
###############################################################################
###############################################################################
###############################################################################
