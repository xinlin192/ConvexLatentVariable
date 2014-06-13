#!/bin/bash

fileName="$(./pathToFile $1)"

echo $fileName

matlab  -r "plotScatter('$fileName')" > log
mv $fileName.pdf ~/public_html/figures/
