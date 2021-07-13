#!/bin/bash

src=$1
data=$2
out=$3
tissue=$4
sample=$5
celltype=$6

python $src/generateGenePrograms.py $data $out $tissue $sample $celltype
