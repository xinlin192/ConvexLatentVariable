#!/bin/bash
datafile=$1
dataset_assign=$2
theta=$3
lambda=$4

echo $datafile;
echo $dataset_assign;
echo $theta;
echo $lambda;

java GenLP $datafile $dataset_assign $theta $lambda;
matlab < solveLP.m;
java NameMap sol varName;
