#!/bin/bash
file_autoc=w0_B64_autocorr_binN.dat
echo "1/N   w0  err " > $file_autoc
for i in $(seq 2 747)
do
	./w0  -p ../../data/  flow_B64.dat  -bin "$i" jack | grep ^w0 | awk -v n="$i" '{print n"  "$3"  "$4}' >> $file_autoc
done
./w0  -p ../../data/  flow_B64.dat  -bin 50 jack loop_B64.dat
