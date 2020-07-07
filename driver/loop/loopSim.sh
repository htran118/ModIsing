#!/bin/bash

stpNum=50000
calStp=1

while read lSze; do
  while read temp; do
    echo "size = $lSze, temp = $temp"
    ./isingModel.exe $lSze $temp $stpNum $calStp
    ./isingAnal.exe $lSze $temp $stpNum
  done < loop/temp.txt
done < loop/lSze.txt
exit 0
