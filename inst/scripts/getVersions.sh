#!/bin/bash

for pkg in oligoClasses Biobase affyio affxparser Biostrings BiocGenerics preprocessCore
do
    echo "==================================================="
    echo "Package: $pkg"
    echo `egrep -i "^Title" $pkg/DESCRIPTION`
    echo `egrep -i "^Version" $pkg/DESCRIPTION`
done
