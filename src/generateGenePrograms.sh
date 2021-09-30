#!/bin/bash

src=$1
programtype=$2
data=$3
out=$4
metadata=$5

python $src/generateGenePrograms.py $programtype $data $out $metadata
