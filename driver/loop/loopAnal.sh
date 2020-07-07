#!/bin/bash

stpNum=20000

rm -f stash.txt

while read lSze; do
  while read temp; do
    ./dataAnal.exe $lSze $temp $stpNum
  done < loop/tempLs.txt
done < loop/lSzeLs.txt
exit 0
