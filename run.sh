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
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_pl0.3_OS_B64.dat  rewpl0.3OS
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_charm_OS_B64.dat  rewcOS
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_strange_OS_B64.dat  rewsOS
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_light_TM_B64.dat  rewlTM
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_charm_0.25_OS_B64.dat  rewcOS_0.25
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_charm_0.02_OS_B64.dat  rewcOS_0.02

## fit der w0 B64
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_charm_0.1_OS_B64.dat  rewcOS_0.1
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_step_0.1232229005_to_0.1234969805_OS_B64.dat  rew_0.1232229005_to_0.1234969805_OS
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_step_0.03125_to_0.03152408_OS_B64.dat  rew_0.03125_to_0.03152408_OS
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_strange_OS_B64.dat  rewsOS
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_step_0.01_to_0.01027408_OS_B64.dat  rew_0.01_to_0.01027408_OS
./w0_rew  -p ../../data/  flow_B64.dat  -bin 20 jack reweight_light_OS_B64.dat  rewlOS
###

./pion  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack loop_B64.dat
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_charm_OS_B64.dat rewcOS
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_strange_OS_B64.dat rewsOS
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_light_TM_B64.dat rewlTM
./pion_rew  -p ../../data/ onlinemeas_B64_LMA.dat  -bin 100 jack reweight_light_TM_B64_LMA.dat rewlTM
./pion_rew  -p ../../data/ onlinemeas_B64_LMA.dat  -bin 100 jack reweight_light_OS_B64_LMA.dat rewlOS

./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_charm_0.25_OS_B64.dat rewcOS_0.25
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_charm_0.25_OS_50s_B64.dat rewcOS_0.25_50s
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_charm_0.25_OS_25s_B64.dat rewcOS_0.25_25s
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_charm_0.25_OS_1s_B64.dat rewcOS_0.25_1s


./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_light_TM_B64_1s.dat rewlTM_1s
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_light_TM_B64_25s.dat rewlTM_25s
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_light_TM_B64_50s.dat rewlTM_50s

./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_charm_0.02_OS_B64.dat rewcOS_0.02
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_charm_0.1_OS_B64.dat rewcOS_0.1


./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_step_0.1232229005_to_0.1234969805_OS_B64.dat rew_0.1232229005_to_0.1234969805_OS
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_step_0.03125_to_0.03152408_OS_B64.dat rew_0.03125_to_0.03152408_OS
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_step_0.01_to_0.01027408_OS_B64.dat rew_0.01_to_0.01027408_OS

# ./fit_all_rew_der  jack ../../data/jackknife/  ../../data/fit_all
./fit_all_rew_der  jack ../../data/jackknife/  ../../data/fit_all_B64  ../fit_files_B64

## ndg
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_sc_ndg_TM_B64.dat rew_scndgTM
./pion_rew  -p ../../data/ onlinemeas_B64.dat  -bin 100 jack reweight_ndg_charm0.1_strange0_TM_B64.dat rew_charm0.1_strange0_ndgTM

### LMA
./pion_rew  -p ../../data/ onlinemeas_B64_LMA.dat  -bin 100 jack reweight_pl0.3_OS_B64_LMA.dat rewpl0.3OS
./pion_rew  -p ../../data/ onlinemeas_B64_LMA.dat  -bin 100 jack reweight_light_OS_B64_LMA.dat rewlOS
./pion_rew  -p ../../data/ onlinemeas_B64_LMA.dat  -bin 100 jack reweight_light_TM_B64_LMA.dat rewlTM

./pion_rew  -p ../../data/ onlinemeas_B64_LMA.dat  -bin 100 jack reweight_charm_0.1_OS_B64_LMA.dat rewcOS_0.1
./pion_rew  -p ../../data/ onlinemeas_B64_LMA.dat  -bin 100 jack reweight_strange_OS_B64_LMA.dat rewsOS
./pion_rew  -p ../../data/ onlinemeas_B64_LMA.dat  -bin 100 jack reweight_step_0.1232229005_to_0.1234969805_OS_B64_LMA.dat rew_0.1232229005_to_0.1234969805_OS
./pion_rew  -p ../../data/ onlinemeas_B64_LMA.dat  -bin 100 jack reweight_step_0.03125_to_0.03152408_OS_B64_LMA.dat rew_0.03125_to_0.03152408_OS
./pion_rew  -p ../../data/ onlinemeas_B64_LMA.dat  -bin 100 jack reweight_step_0.01_to_0.01027408_OS_B64_LMA.dat rew_0.01_to_0.01027408_OS

./fit_all_rew_der  jack ../../data/jackknife/  ../../data/fit_all_B64_LMA  ../fit_files_B64_LMA

./pion_rew  -p ../../data/ onlinemeas_B64_LMA.dat  -bin 100 jack reweight_charm_0.1_h10_OS_B64_LMA.dat rewcOS_0.1
./pion_rew  -p ../../data/ onlinemeas_B64_LMA.dat  -bin 100 jack reweight_charm_0.1_h5_OS_B64_LMA.dat rewcOS_0.1


## C80
# ./pion  -p ../../data/ onlinemeas_C80.dat  -bin 100 jack loop_C80.dat

./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_light_TM_C80_LMA.dat rewlTM

./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_light_OS_C80_LMA.dat rewlOS
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_charm_0.1_OS_C80_LMA.dat rewcOS_0.1
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_strange_OS_C80_LMA.dat rewsOS
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_strange_0.2_OS_C80_LMA.dat rewsOS_0.2
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_strange_pdmu_OS_C80_LMA.dat rewsOS_pdmu
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_strange_p2dmu_OS_C80_LMA.dat rewsOS_p2dmu
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_strange_mdmu_OS_C80_LMA.dat rewsOS_mdmu


./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_charm_0.1_h5_OS_C80_LMA.dat rewcOS_0.1_h5
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_strange_h5_OS_C80_LMA.dat rewsOS_h5
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_strange_pdmu_h5_OS_C80_LMA.dat rewsOS_pdmu_h5
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_strange_p2dmu_h5_OS_C80_LMA.dat rewsOS_p2dmu_h5
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_strange_mdmu_h5_OS_C80_LMA.dat rewsOS_mdmu_h5

./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_light_h10_OS_C80_LMA.dat rewlOS_h10
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_charm_0.1_h10_OS_C80_LMA.dat rewcOS_0.1_h10
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_strange_h10_OS_C80_LMA.dat rewsOS_h10
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_strange_pdmu_h10_OS_C80_LMA.dat rewsOS_pdmu_h10
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_strange_p2dmu_h10_OS_C80_LMA.dat rewsOS_p2dmu_h10
./pion_rew  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_strange_mdmu_h10_OS_C80_LMA.dat rewsOS_mdmu_h10

###### B32
./pion_rew  -p ../../data/ onlinemeas_B32.dat  -bin 100 jack reweight_light_OS_B32.dat rewlOS
./pion_rew  -p ../../data/ onlinemeas_B32.dat  -bin 100 jack reweight_strange_OS_B32.dat rewsOS
./pion_rew  -p ../../data/ onlinemeas_B32.dat  -bin 100 jack reweight_charm_OS_B32.dat rewcOS

./pion_rew  -p ../../data/ onlinemeas_mu1_B32.dat  -bin 100 jack reweight_light_OS_B32.dat rewlOS
./pion_rew  -p ../../data/ onlinemeas_mu1_B32.dat  -bin 100 jack reweight_strange_OS_B32.dat rewsOS
./pion_rew  -p ../../data/ onlinemeas_mu1_B32.dat  -bin 100 jack reweight_charm_OS_B32.dat rewcOS

./pion_rew  -p ../../data/ onlinemeas_B32.dat  -bin 100 jack reweight_charm_OS_more_B32.dat rewcOS
./pion_rew  -p ../../data/ onlinemeas_more_B32.dat  -bin 100 jack reweight_charm_OS_B32.dat rewcOS
./pion_rew  -p ../../data/ onlinemeas_more_B32.dat  -bin 100 jack reweight_charm_OS_more_B32.dat rewcOS

./rew_other_data  -p ../../data/ marcos_format/B32/B32_Marcos_P5P5_tm_mu2.5000e-03_mu2.5000e-03.bin  -bin 100 jack reweight_charm_OS_B32.dat rewcOS  onlinemeas_B32.dat


###### B24
./pion_rew  -p ../../data/ onlinemeas_B24.dat  -bin 100 jack reweight_light_OS_B24.dat rewlOS
./pion_rew  -p ../../data/ onlinemeas_B24.dat  -bin 100 jack reweight_strange_OS_B24.dat rewsOS
./pion_rew  -p ../../data/ onlinemeas_B24.dat  -bin 100 jack reweight_charm_OS_B24.dat rewcOS

./rew_other_data  -p ../../data/ Marcos_P5P5_tm_mu2.5000e-03_mu2.5000e-03.bin  -bin 100 jack reweight_charm_OS_B24.dat rewcOS  onlinemeas_B24.dat

###### C48
./pion_rew  -p ../../data/ onlinemeas_C48.dat  -bin 100 jack reweight_light_OS_C48.dat rewlOS
./pion_rew  -p ../../data/ onlinemeas_C48.dat  -bin 100 jack reweight_strange_OS_C48.dat rewsOS
./pion_rew  -p ../../data/ onlinemeas_C48.dat  -bin 100 jack reweight_charm_OS_C48.dat rewcOS
./pion_rew  -p ../../data/ onlinemeas_C48_confs434.dat  -bin 100 jack reweight_charm_OS_C48_confs434.dat rewcOS



./pion_rew  -p ../../data/ onlinemeas_C48.dat.not_thermalized  -bin 100 jack reweight_light_OS_C48.dat.not_thermalized rewlOS
./pion_rew  -p ../../data/ onlinemeas_C48.dat.not_thermalized  -bin 100 jack reweight_strange_OS_C48.dat.not_thermalized rewsOS
./pion_rew  -p ../../data/ onlinemeas_C48.dat.not_thermalized  -bin 100 jack reweight_charm_OS_C48.dat.not_thermalized rewcOS

./pion_rew  -p ../../data/ onlinemeas_C48.dat  -bin 100 jack reweight_charm_OS_C48_low_prec.dat rewcOS
##### E112
./pion_rew  -p ../../data/ onlinemeas_E112_LMA.dat  -bin 100 jack reweight_light_OS_E112_LMA.dat rewlOS
./pion_rew  -p ../../data/ onlinemeas_E112_LMA.dat  -bin 100 jack reweight_strange_OS_E112_LMA.dat rewsOS
# ./pion_rew  -p ../../data/ onlinemeas_E112_LMA.dat  -bin 100 jack reweight_charm_OS_E112_LMA.dat rewcOS

### A32
./pion_rew  -p ../../data/ onlinemeas_A32.dat  -bin 100 jack reweight_light_OS_A32.dat rewlOS
./pion_rew  -p ../../data/ onlinemeas_A32.dat  -bin 100 jack reweight_strange_OS_A32.dat rewsOS
./pion_rew  -p ../../data/ onlinemeas_A32.dat  -bin 100 jack reweight_charm_OS_A32.dat rewcOS

./rew_other_data  -p ../../data/  marcos_format/A32/A32_Marcos_P5P5_tm_mu3.0000e-03_mu3.0000e-03.bin  -bin 100 jack reweight_charm_OS_A32.dat rewcOS onlinemeas_A32.dat

./rew_other_data  -p ../../data/  marcos_format/A32/A32_Marcos_P5P5_tm_mu3.0000e-03_mu3.0000e-03.bin  -bin 100 jack reweight_strange_OS_A32.dat rewsOS onlinemeas_A32.dat
./rew_other_data  -p ../../data/  marcos_format/A32/A32_Marcos_P5P5_tm_mu5.0000e-03_mu5.0000e-03.bin  -bin 100 jack reweight_strange_OS_A32.dat rewsOS onlinemeas_A32.dat

## fake D96
# ./create_fake_jack  -p ../../data/ onlinemeas_D96_LMA.dat  -bin 100 jack reweight_charm_OS_D96_LMA.dat rewcOS
# ./create_fake_jack  -p ../../data/ onlinemeas_D96_LMA.dat  -bin 100 jack reweight_strange_p2dmu_OS_D96_LMA.dat rewsOS_p2dmu
# ./create_fake_jack  -p ../../data/ onlinemeas_D96_LMA.dat  -bin 100 jack reweight_strange_pdmu_OS_D96_LMA.dat rewsOS_pdmu
# ./pion_rew  -p ../../data/ onlinemeas_D96.dat  -bin 100 jack reweight_charm_dmu_OS_D96.dat rewcOS_dmu
./pion_rew  -p ../../data/ onlinemeas_D96.dat  -bin 100 jack reweight_charm_OS_D96.dat rewcOS
#./pion_rew  -p ../../data/ onlinemeas_D96.dat  -bin 100 jack reweight_strange_pdmu_OS_D96.dat rewsOS_pdmu
./pion_rew  -p ../../data/ onlinemeas_D96.dat  -bin 100 jack reweight_strange_OS_D96.dat rewsOS
# ./create_fake_jack  -p ../../data/ onlinemeas_D96_LMA.dat  -bin 100 jack reweight_strange_mdmu_OS_D96_LMA.dat rewsOS_mdmu
./pion_rew  -p ../../data/ onlinemeas_D96.dat  -bin 100 jack reweight_light_OS_D96.dat rewlOS

./fit_all_beta_rew_der  jack ../../data/jackknife/  ../../data/fit_all_beta  ../fit_files_all



############################


## w0 C80
./w0_rew  -p ../../data/ flow_C80_LMA.dat  -bin 20 jack reweight_charm_0.1_OS_C80_LMA.dat rewcOS_0.1
./w0_rew  -p ../../data/ flow_C80_LMA.dat  -bin 20 jack reweight_strange_OS_C80_LMA.dat rewsOS
./w0_rew  -p ../../data/ flow_C80_LMA.dat  -bin 20 jack reweight_strange_0.2_OS_C80_LMA.dat rewsOS_0.2
./w0_rew  -p ../../data/ flow_C80_LMA.dat  -bin 20 jack reweight_strange_pdmu_OS_C80_LMA.dat rewsOS_pdmu
./w0_rew  -p ../../data/ flow_C80_LMA.dat  -bin 20 jack reweight_strange_p2dmu_OS_C80_LMA.dat rewsOS_p2dmu
./w0_rew  -p ../../data/ flow_C80_LMA.dat  -bin 20 jack reweight_strange_mdmu_OS_C80_LMA.dat rewsOS_mdmu
./w0_rew  -p ../../data/ flow_C80_LMA.dat  -bin 20 jack reweight_light_OS_C80_LMA.dat rewlOS
##

# ./fit_all_rew_der_C80_LMA   jack ../../data/jackknife/  ../../data/fit_all_C80_LMA
# ./fit_all_rew_der  jack ../../data/jackknife/  ../../data/fit_all_C80_LMA  ../fit_files_C80_LMA


## w0 D96
# ./w0_rew  -p ../../data/ flow_D96.dat  -bin 20 jack reweight_charm_dmu_OS_D96.dat rewcOS_dmu
# ./w0_rew  -p ../../data/ flow_D96.dat  -bin 20 jack reweight_strange_pdmu_OS_D96.dat rewsOS_pdmu
./w0_rew  -p ../../data/ flow_D96.dat  -bin 20 jack reweight_charm_OS_D96.dat rewcOS
./w0_rew  -p ../../data/ flow_D96.dat  -bin 20 jack reweight_strange_OS_D96.dat rewsOS
./w0_rew  -p ../../data/ flow_D96.dat  -bin 20 jack reweight_light_OS_D96.dat rewlOS.1
##

## w0 B32
./w0_rew  -p ../../data/ flow_B32.dat  -bin 20 jack reweight_charm_OS_B32.dat rewcOS
./w0_rew  -p ../../data/ flow_B32.dat  -bin 20 jack reweight_strange_OS_B32.dat rewsOS
./w0_rew  -p ../../data/ flow_B32.dat  -bin 20 jack reweight_light_OS_B32.dat rewlOS
## w0 B24
./w0_rew  -p ../../data/ flow_B24.dat  -bin 20 jack reweight_charm_OS_B24.dat rewcOS
./w0_rew  -p ../../data/ flow_B24.dat  -bin 20 jack reweight_strange_OS_B24.dat rewsOS
./w0_rew  -p ../../data/ flow_B24.dat  -bin 20 jack reweight_light_OS_B24.dat rewlOS
##
./w0_rew  -p ../../data/ flow_C48.dat  -bin 20 jack reweight_charm_OS_C48.dat rewcOS
./w0_rew  -p ../../data_p/ flow_C48.dat  -bin 20 jack reweight_charm_OS_C48.dat rewcOS

./w0_rew  -p ../../data/ flow_C48.dat  -bin 20 jack reweight_strange_OS_C48.dat rewsOS
./w0_rew  -p ../../data/ flow_C48.dat  -bin 20 jack reweight_light_OS_C48.dat rewlOS

./w0_rew  -p ../../data/ flow_C48.dat.not_thermalized  -bin 20 jack reweight_charm_OS_C48.dat.not_thermalized rewcOS
./w0_rew  -p ../../data/ flow_C48.dat.not_thermalized  -bin 20 jack reweight_strange_OS_C48.dat.not_thermalized rewsOS
./w0_rew  -p ../../data/ flow_C48.dat.not_thermalized  -bin 20 jack reweight_light_OS_C48.dat.not_thermalized rewlOS


# ./w0_rew  -p ../../data/ flow_C48_therm.dat  -bin 20 jack reweight_charm_OS_C48_therm.dat rewcOS
# ./w0_rew  -p ../../data/ flow_C48_therm.dat  -bin 20 jack reweight_strange_OS_C48_therm.dat rewsOS
# ./w0_rew  -p ../../data/ flow_C48_therm.dat  -bin 20 jack reweight_light_OS_C48_therm.dat rewlOS
## E112
./w0_rew  -p ../../data/ flow_E112_LMA.dat  -bin 20 jack reweight_charm_OS_E112_LMA.dat rewcOS
./w0_rew  -p ../../data/ flow_E112_LMA.dat  -bin 20 jack reweight_strange_OS_E112_LMA.dat rewsOS
./w0_rew  -p ../../data/ flow_E112_LMA.dat  -bin 20 jack reweight_light_OS_E112_LMA.dat rewlOS

## A32
./w0_rew  -p ../../data/ flow_A32.dat  -bin 20 jack reweight_charm_OS_A32.dat rewcOS
./w0_rew  -p ../../data/ flow_A32.dat  -bin 20 jack reweight_strange_OS_A32.dat rewsOS
./w0_rew  -p ../../data/ flow_A32.dat  -bin 20 jack reweight_light_OS_A32.dat rewlOS


### rations
./pion_rew_ratio  -p ../../data/ onlinemeas_B64_LMA.dat  -bin 100 jack reweight_light_OS_B64_LMA.dat reweight_light_OS_B64_LMA.dat
./pion_rew_ratio  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_light_OS_C80_LMA.dat reweight_light_OS_C80_LMA.dat
./pion_rew_ratio  -p ../../data/ onlinemeas_C80_LMA.dat  -bin 100 jack reweight_light_OS_C80_LMA.dat reweight_strange_mdmu_OS_C80_LMA.dat

##### fit all beta w0
./fit_all_beta_w0rew_der  jack ../../data/jackknife/  ../../data/fit_all_beta  ../fit_files_all_w0


./set_w0_a  jack ../../data/jackknife/  ../../data/fit_all_beta  ../data_ens_E112   iso_param  E112