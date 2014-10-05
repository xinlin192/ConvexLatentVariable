#!/bin/bash
datafile=$1;
lambda=$2;

echo $datafile;
echo $lambda;

java GenLP $datafile $lambda;
matlab < solveLP.m;
java NameMap sol varName;
