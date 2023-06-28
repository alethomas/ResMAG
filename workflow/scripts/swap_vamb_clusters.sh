#!/bin/bash

input=$1 
output=$2
temp=$3

awk '{print $2 "\t" $1}' $input > $temp
sed 's/S[12]C//g' $temp > $output
rm $temp