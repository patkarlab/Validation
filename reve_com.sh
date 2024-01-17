#! /usr/bin/bash

input=$1

for i in `cat ${input}`
do

echo ${i} | tr 'ATCGatcg' 'TAGCtagc' | rev

done
