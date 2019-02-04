#!/bin/bash
# this script will run tsfm on all the genomes 
# for each genomes it will make a folder in each folder will make three folders: KLD,ID,Func, bubble_Table
# are we looking to find differences between TryTryp genes themself or between them and human?

folderpath="/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADBClustered/"
tsfmpath="/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfm-master/tsfm"
folders=$(ls -l $folderpath | grep "^d" | awk -F" " '{print $9}')


for name in $folders; do 
if [ $name != "Logos" ]
then
printf "***************$name***************"
	mkdir -p "$folderpath/Logos/$name/KLD"
	mkdir -p "$folderpath/Logos/$name/ID"
	mkdir -p "$folderpath/Logos/$name/Func"
	mkdir -p "$folderpath/Logos/$name/Bubble"
# run tsfm to make the logos 
python3 "$tsfmpath/tsfm.py" -c "$folderpath/struct_file.txt" --logo "$folderpath/$name/clustalW/$name"
mv -- *.eps "$folderpath/Logos/$name/Func"

python3 "$tsfmpath/tsfm.py" -c "$folderpath/struct_file.txt" --IDlogo "$folderpath/$name/clustalW/$name" "$folderpath/HOMO/clustalW/HOMO"
mv -- *.eps "$folderpath/Logos/$name/ID"

python3 "$tsfmpath/tsfm.py" -c "$folderpath/struct_file.txt" --KLDlogo "$folderpath/$name/clustalW/$name" "$folderpath/HOMO/clustalW/HOMO"
mv -- *.eps "$folderpath/Logos/$name/KLD"

python3 "$tsfmpath/tsfm.py" -c "$folderpath/struct_file.txt" --IDlogo --KLDlogo --bt "$folderpath/$name/clustalW/$name" "$folderpath/HOMO/clustalW/HOMO"
rm -- *.eps

mv *_Table.txt "$folderpath/Logos/$name/Bubble"
fi
done





