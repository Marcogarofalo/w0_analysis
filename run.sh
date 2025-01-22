#!/bin/bash
file_autoc=w0_B64_autocorr_bintoNb.dat
echo "1/N   w0  err " > $file_autoc
for i in $(seq 2 749)
do
	./w0  -p ../../data/  flow_B64.dat  -bin "$i" jack loop_B64.dat| grep "^w0 = " | awk -v n="$i" '{print 749/n"  "$3"  "$4}' >> $file_autoc
done
#./w0  -p ../../data/  flow_B64.dat  -bin 50 jack loop_B64.dat
