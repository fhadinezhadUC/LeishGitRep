#!/bin/bash
folderpath="/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADBClustered/"
tsfmpath="/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfm-master/tsfm"
folders=$(ls -l $folderpath | grep "^d" | awk -F" " '{print $9}')


for name in $folders; do 
printf "goh"
if [ $name == "Logos" ]
then
printf "***************$name***************"

fi
done





