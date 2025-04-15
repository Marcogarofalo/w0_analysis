#!/bin/bash
file_autoc=w0_B64_autocorr_bintoNb.dat
file_autoc_l=w0+mul_B64_autocorr_bintoNb.dat
file_autoc_s=w0+mus_B64_autocorr_bintoNb.dat
file_autoc_c=w0+muc_B64_autocorr_bintoNb.dat
# echo "1/N   w0  err " > $file_autoc
# echo "1/N   w0  err " > $file_autoc_l
# echo "1/N   w0  err " > $file_autoc_s
# echo "1/N   w0  err " > $file_autoc_c
# is="5 6 7 8 9 10 12 15 17 20 30 40 50 60 70 80 90 100 200 300 500 749"
# # is="5 "
# for i in $is
# do
# 	./w0  -p ../../data/  flow_B64.dat  -bin "$i" jack loop_B64.dat > tmp
#     grep "^w0 = " tmp | awk -v n="$i" '{print 749/n"  "$3"  "$4}' >> $file_autoc
#     grep "^w0+mul_correction = " tmp | awk -v n="$i" '{print 749/n"  "$3"  "$4}' >> $file_autoc_l
#     grep "^w0+mus_correction = " tmp | awk -v n="$i" '{print 749/n"  "$3"  "$4}' >> $file_autoc_s
#     grep "^w0+muc_correction = " tmp | awk -v n="$i" '{print 749/n"  "$3"  "$4}' >> $file_autoc_c
# done
./w0  -p ../../data/  flow_B64.dat  -bin 20 jack loop_B64.dat
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_charm_OS_B64.dat  rewcOS
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_strange_OS_B64.dat  rewsOS
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_light_TM_B64.dat  rewlTM
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_charm_0.25_OS_B64.dat  rewcOS_0.25
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_charm_0.02_OS_B64.dat  rewcOS_0.02

./pion  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack loop_B64.dat
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_charm_OS_B64.dat rewcOS
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_strange_OS_B64.dat rewsOS
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_light_TM_B64.dat rewlTM

./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_sc_ndg_TM_B64.dat rew_scndgTM

./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_charm_0.25_OS_B64.dat rewcOS_0.25
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_charm_0.25_OS_50s_B64.dat rewcOS_0.25_50s
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_charm_0.25_OS_25s_B64.dat rewcOS_0.25_25s
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_charm_0.25_OS_1s_B64.dat rewcOS_0.25_1s


./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_light_TM_B64_1s.dat rewlTM_1s
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_light_TM_B64_25s.dat rewlTM_25s
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_light_TM_B64_50s.dat rewlTM_50s

./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_charm_0.02_OS_B64.dat rewcOS_0.02