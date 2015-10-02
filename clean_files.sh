#!/bin/bash

#if no time stamps originally
Date1=`date +"%m-%d-%Y"`

for i in *sub
do 
    echo "initial dir = $Date1" >> $i
done

