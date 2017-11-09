#!/bin/bash


pred=(1.0 0.5 0.0)
food=(100 20)

for i in ${food[@]}
do

for j in ${pred[@]}
do

sed -i s/'i_s = [0-9].[0-9]'/'i_s = '$j/ params.h
sed -i s/'prey_pop_count = [0-9]*'/'prey_pop_count = '$i/ params.h
make

if [ $i == 100 ]
then
floc=highfood
elif [ $i == 20 ]
then
floc=lowfood
fi

if [ $j == '1.0' ]
then
ploc=high
elif [ $j == '0.5' ]
then
ploc=low
elif [ $j == '0.0' ]
then
ploc=none
fi

mv a.out ../RAND_ON/${floc}/${ploc}/${floc}-${ploc}.out

done
done


