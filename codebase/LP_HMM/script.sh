#!/bin/bash
datafile=$1;
lambda=$2;

echo $datafile;
echo $lambda;
java GenLP $datafile 15 $lambda 10000 1;
matlab < solveLP.m;
java NameMap sol varName;
