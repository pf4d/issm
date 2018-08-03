#!/bin/bash

#Synchronize Enums

#Get all lines of EnumDefinitions.h which hold Enum 
cat c.vim | sed "/ISSM's Enums begin/,/vim: ts=8/d" > temp
echo "\"ISSM's Enums begin" >> temp
cat ../../../../../src/c/shared/Enum/EnumDefinitions.h | grep -e "[0-9]Enum," -e "[a-z]Enum," -e "[A-Z]Enum," | grep -v StringToEnum | sed -e "s/,//g" | awk '{ printf "syn keyword cConstant " $1 "\n"}' >> temp
echo "\"ISSM's Enums end" >> temp
cat c.vim | sed "1,/ISSM's Enums end/d" >> temp
mv temp c.vim

#Synchronize objects
cat c.vim | sed "/ISSM's objects begin/,/vim: ts=8/d" > temp
echo "\"ISSM's objects begin" >> temp
find ../../../../../src/c/classes -name "*.cpp" -o -name "*.h" | sed -e "s/\// /g" -e "s/\.cpp//" -e "s/\.h//" | awk '{print  $(NF)}' | sort | uniq | awk '{ printf "syn keyword cType " $1 "\n"}'>> temp
find ../../../../../src/c/analyses -name "*Analysis.h" | sed -e "s/\// /g" -e "s/\.cpp//" -e "s/\.h//" | awk '{print  $(NF)}' | sort | uniq | awk '{ printf "syn keyword cType " $1 "\n"}'>> temp
echo "\"ISSM's objects end" >> temp
cat c.vim | sed "1,/ISSM's objects end/d" >> temp

mv temp c.vim
