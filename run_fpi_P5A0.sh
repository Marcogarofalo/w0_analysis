# ./other_fpi  -p ../../data_fpi_A0   B64_TM_l1_l1_P5A0.dat  -bin 20  jack  B64_TM_l1_l1_P5P5.dat  B64_TM_l2_l2_P5A0.dat  B64_TM_deriv_P5P5.dat

res=()

run(){
    local -n tmp=$2
    local exe=$1
    $exe || tmp+=("$exe")
}
gammas=("A0P5"  "P5P5")
reg=("tm"  "OS")

#### B64
for ((ir=0; ir<${#reg[@]}; ir++)); do
for ((ig=0; ig<${#gammas[@]}; ig++)); do
run  "./derivative_corr  -p ../../data_fpi_A0/ B64_${reg[ir]}_l0_l0_${gammas[ig]}.dat -bin 1 jack B64_${reg[ir]}_l1_l1_${gammas[ig]}.dat  B64_${reg[ir]}_val_deriv_${gammas[ig]}.dat" res
done
done



run "./other_fpi  -p ../../data_fpi_A0   B64_tm_l0_l0_A0P5.dat  -bin 20  jack  B64_tm_l0_l0_P5P5.dat  B64_tm_val_deriv_A0P5.dat  B64_tm_val_deriv_P5P5.dat  onlinemeas_B64_LMA.dat" res
run "./other_fpi  -p ../../data_fpi_A0   B64_OS_l0_l0_A0P5.dat  -bin 20  jack  B64_OS_l0_l0_P5P5.dat  B64_OS_val_deriv_A0P5.dat  B64_OS_val_deriv_P5P5.dat  onlinemeas_B64_LMA.dat" res


run "./deriv_other_fpi  -p ../../data_fpi_A0   B64_tm_l0_l0_A0P5.dat  -bin 20  jack  B64_tm_l0_l0_P5P5.dat   reweight_light_OS_B64_LMA.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B64_OS_l0_l0_A0P5.dat  -bin 20  jack  B64_OS_l0_l0_P5P5.dat   reweight_light_OS_B64_LMA.dat  onlinemeas_B64_LMA.dat"  res

### B96
## the valence deivative is garbage, just to run the sim value
run "./other_fpi  -p ../../data_fpi_A0   B96_tm_l0_l0_A0P5.dat  -bin 20  jack  B96_tm_l0_l0_P5P5.dat  B96_tm_l0_l0_A0P5.dat  B96_tm_l0_l0_P5P5.dat  onlinemeas_B64_LMA.dat" res
run "./other_fpi  -p ../../data_fpi_A0   B96_OS_l0_l0_A0P5.dat  -bin 20  jack  B96_OS_l0_l0_P5P5.dat  B96_OS_l0_l0_A0P5.dat  B96_OS_l0_l0_P5P5.dat  onlinemeas_B64_LMA.dat" res

### C80

gammas=("A0P5"  "P5P5")
reg=("tm"  "OS")
for ((ir=0; ir<${#reg[@]}; ir++)); do
for ((ig=0; ig<${#gammas[@]}; ig++)); do
run  "./derivative_corr  -p ../../data_fpi_A0/ C80_${reg[ir]}_l0_l0_${gammas[ig]}.dat -bin 1 jack C80_${reg[ir]}_l1_l1_${gammas[ig]}.dat  C80_${reg[ir]}_val_deriv_${gammas[ig]}.dat" res
done
done

run "./other_fpi  -p ../../data_fpi_A0   C80_tm_l0_l0_A0P5.dat  -bin 20  jack  C80_tm_l0_l0_P5P5.dat  C80_tm_val_deriv_A0P5.dat  C80_tm_val_deriv_P5P5.dat  onlinemeas_C80_LMA.dat" res
run "./other_fpi  -p ../../data_fpi_A0   C80_OS_l0_l0_A0P5.dat  -bin 20  jack  C80_OS_l0_l0_P5P5.dat  C80_OS_val_deriv_A0P5.dat  C80_OS_val_deriv_P5P5.dat  onlinemeas_C80_LMA.dat" res


run "./deriv_other_fpi  -p ../../data_fpi_A0   C80_tm_l0_l0_A0P5.dat  -bin 20  jack  C80_tm_l0_l0_P5P5.dat   reweight_light_OS_C80_LMA.dat  onlinemeas_C80_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   C80_OS_l0_l0_A0P5.dat  -bin 20  jack  C80_OS_l0_l0_P5P5.dat   reweight_light_OS_C80_LMA.dat  onlinemeas_C80_LMA.dat"  res


## D96 

./derivative_corr  -p ../../data_fpi_A0/ D96_tm_smallstat_l0_l0_A0P5.dat -bin 1 jack D96_tm_smallstat_l1_l1_A0P5.dat  D96_tm_val_deriv_A0P5.dat
./derivative_corr  -p ../../data_fpi_A0/ D96_OS_smallstat_l0_l0_A0P5.dat -bin 1 jack D96_OS_smallstat_l1_l1_A0P5.dat  D96_OS_val_deriv_A0P5.dat

run "./other_fpi  -p ../../data_fpi_A0   D96_tm_l0_l0_A0P5.dat  -bin 20  jack  D96_tm_l0_l0_P5P5.dat  D96_tm_val_deriv_A0P5.dat  D96_tm_val_deriv_P5P5.dat  onlinemeas_D96.dat" res
run "./other_fpi  -p ../../data_fpi_A0   D96_OS_l0_l0_A0P5.dat  -bin 20  jack  D96_OS_l0_l0_P5P5.dat  D96_OS_val_deriv_A0P5.dat  D96_OS_val_deriv_P5P5.dat  onlinemeas_D96.dat" res

run "./deriv_other_fpi  -p ../../data_fpi_A0   D96_tm_l0_l0_A0P5.dat  -bin 20  jack  D96_tm_l0_l0_P5P5.dat   reweight_light_OS_D96.dat  onlinemeas_D96.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   D96_OS_l0_l0_A0P5.dat  -bin 20  jack  D96_OS_l0_l0_P5P5.dat   reweight_light_OS_D96.dat  onlinemeas_D96.dat"  res


## E112

./derivative_corr  -p ../../data_fpi_A0/ E112_tm_smallstat_l0_l0_A0P5.dat -bin 1 jack E112_tm_smallstat_l1_l1_A0P5.dat  E112_tm_val_deriv_A0P5.dat
./derivative_corr  -p ../../data_fpi_A0/ E112_OS_smallstat_l0_l0_A0P5.dat -bin 1 jack E112_OS_smallstat_l1_l1_A0P5.dat  E112_OS_val_deriv_A0P5.dat

# run "./other_fpi  -p ../../data_fpi_A0   E112_tm_l0_l0_A0P5.dat  -bin 20  jack  E112_tm_l0_l0_P5P5.dat  E112_tm_l1_l1_A0P5.dat  E112_tm_l1_l1_P5P5.dat  onlinemeas_E112_LMA.dat" res
# run "./other_fpi  -p ../../data_fpi_A0   E112_OS_l0_l0_A0P5.dat  -bin 20  jack  E112_OS_l0_l0_P5P5.dat  E112_OS_l1_l1_A0P5.dat  E112_OS_l1_l1_P5P5.dat  onlinemeas_E112_LMA.dat" res
run "./other_fpi  -p ../../data_fpi_A0   E112_tm_l0_l0_A0P5.dat  -bin 20  jack  E112_tm_l0_l0_P5P5.dat  E112_tm_val_deriv_A0P5.dat  E112_tm_val_deriv_P5P5.dat  onlinemeas_E112_LMA.dat" res
run "./other_fpi  -p ../../data_fpi_A0   E112_OS_l0_l0_A0P5.dat  -bin 20  jack  E112_OS_l0_l0_P5P5.dat  E112_OS_val_deriv_A0P5.dat  E112_OS_val_deriv_P5P5.dat  onlinemeas_E112_LMA.dat" res

run "./deriv_other_fpi  -p ../../data_fpi_A0   E112_tm_l0_l0_A0P5.dat  -bin 20  jack  E112_tm_l0_l0_P5P5.dat   reweight_light_OS_E112_LMA.dat  onlinemeas_E112_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   E112_OS_l0_l0_A0P5.dat  -bin 20  jack  E112_OS_l0_l0_P5P5.dat   reweight_light_OS_E112_LMA.dat  onlinemeas_E112_LMA.dat"  res


## small Volume

run "./deriv_other_fpi  -p ../../data_fpi_A0   B24_tm_l0_l0_A0P5.dat  -bin 20  jack  B24_tm_l0_l0_P5P5.dat   reweight_charm_OS_B24.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B24_tm_l1_l1_A0P5.dat  -bin 20  jack  B24_tm_l1_l1_P5P5.dat   reweight_charm_OS_B24.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B24_tm_l2_l2_A0P5.dat  -bin 20  jack  B24_tm_l2_l2_P5P5.dat   reweight_charm_OS_B24.dat  onlinemeas_B64_LMA.dat"  res
run "./fit_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../fit_smallV_B24_fpi_A0P5.txt  "  res

run "./deriv_other_fpi  -p ../../data_fpi_A0   B24_OS_l0_l0_A0P5.dat  -bin 20  jack  B24_OS_l0_l0_P5P5.dat   reweight_charm_OS_B24.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B24_OS_l1_l1_A0P5.dat  -bin 20  jack  B24_OS_l1_l1_P5P5.dat   reweight_charm_OS_B24.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B24_OS_l2_l2_A0P5.dat  -bin 20  jack  B24_OS_l2_l2_P5P5.dat   reweight_charm_OS_B24.dat  onlinemeas_B64_LMA.dat"  res
run "./fit_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../fit_smallV_B24_OS_fpi_A0P5.txt  "  res



run "./deriv_other_fpi  -p ../../data_fpi_A0   B32_tm_l0_l0_A0P5.dat  -bin 20  jack  B32_tm_l0_l0_P5P5.dat   reweight_charm_OS_B32.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B32_tm_l1_l1_A0P5.dat  -bin 20  jack  B32_tm_l1_l1_P5P5.dat   reweight_charm_OS_B32.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B32_tm_l2_l2_A0P5.dat  -bin 20  jack  B32_tm_l2_l2_P5P5.dat   reweight_charm_OS_B32.dat  onlinemeas_B64_LMA.dat"  res
run "./fit_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../fit_smallV_B32_fpi_A0P5.txt  "  res

run "./deriv_other_fpi  -p ../../data_fpi_A0   B32_OS_l0_l0_A0P5.dat  -bin 20  jack  B32_OS_l0_l0_P5P5.dat   reweight_charm_OS_B32.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B32_OS_l1_l1_A0P5.dat  -bin 20  jack  B32_OS_l1_l1_P5P5.dat   reweight_charm_OS_B32.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B32_OS_l2_l2_A0P5.dat  -bin 20  jack  B32_OS_l2_l2_P5P5.dat   reweight_charm_OS_B32.dat  onlinemeas_B64_LMA.dat"  res
run "./fit_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../fit_smallV_B32_OS_fpi_A0P5.txt  "  res


run "./deriv_other_fpi  -p ../../data_fpi_A0   C48_tm_l0_l0_A0P5.dat  -bin 20  jack  C48_tm_l0_l0_P5P5.dat   reweight_charm_OS_C48.dat  onlinemeas_C80_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   C48_tm_l1_l1_A0P5.dat  -bin 20  jack  C48_tm_l1_l1_P5P5.dat   reweight_charm_OS_C48.dat  onlinemeas_C80_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   C48_tm_l2_l2_A0P5.dat  -bin 20  jack  C48_tm_l2_l2_P5P5.dat   reweight_charm_OS_C48.dat  onlinemeas_C80_LMA.dat"  res
run "./fit_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../fit_smallV_C48_fpi_A0P5.txt  "  res

run "./deriv_other_fpi  -p ../../data_fpi_A0   C48_OS_l0_l0_A0P5.dat  -bin 20  jack  C48_OS_l0_l0_P5P5.dat   reweight_charm_OS_C48.dat  onlinemeas_C80_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   C48_OS_l1_l1_A0P5.dat  -bin 20  jack  C48_OS_l1_l1_P5P5.dat   reweight_charm_OS_C48.dat  onlinemeas_C80_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   C48_OS_l2_l2_A0P5.dat  -bin 20  jack  C48_OS_l2_l2_P5P5.dat   reweight_charm_OS_C48.dat  onlinemeas_C80_LMA.dat"  res
run "./fit_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../fit_smallV_C48_OS_fpi_A0P5.txt  "  res


run "./average_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../average_mul_over_fpi_df_dmuc_smallV_tm_fpi_A0P5.txt  "  res
run "./average_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../average_mul_over_fpi_df_dmuc_smallV_OS_fpi_A0P5.txt  "  res

##### strange

run "./deriv_other_fpi  -p ../../data_fpi_A0   B24_tm_l0_l0_A0P5.dat  -bin 20  jack  B24_tm_l0_l0_P5P5.dat   reweight_strange_OS_B24.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B24_tm_l1_l1_A0P5.dat  -bin 20  jack  B24_tm_l1_l1_P5P5.dat   reweight_strange_OS_B24.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B24_tm_l2_l2_A0P5.dat  -bin 20  jack  B24_tm_l2_l2_P5P5.dat   reweight_strange_OS_B24.dat  onlinemeas_B64_LMA.dat"  res
run "./fit_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../fit_smallV_B24_strange_fpi_A0P5.txt  "  res

run "./deriv_other_fpi  -p ../../data_fpi_A0   B24_OS_l0_l0_A0P5.dat  -bin 20  jack  B24_OS_l0_l0_P5P5.dat   reweight_strange_OS_B24.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B24_OS_l1_l1_A0P5.dat  -bin 20  jack  B24_OS_l1_l1_P5P5.dat   reweight_strange_OS_B24.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B24_OS_l2_l2_A0P5.dat  -bin 20  jack  B24_OS_l2_l2_P5P5.dat   reweight_strange_OS_B24.dat  onlinemeas_B64_LMA.dat"  res
run "./fit_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../fit_smallV_B24_OS_strange_fpi_A0P5.txt  "  res


run "./deriv_other_fpi  -p ../../data_fpi_A0   B32_tm_l0_l0_A0P5.dat  -bin 20  jack  B32_tm_l0_l0_P5P5.dat   reweight_strange_OS_B32.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B32_tm_l1_l1_A0P5.dat  -bin 20  jack  B32_tm_l1_l1_P5P5.dat   reweight_strange_OS_B32.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B32_tm_l2_l2_A0P5.dat  -bin 20  jack  B32_tm_l2_l2_P5P5.dat   reweight_strange_OS_B32.dat  onlinemeas_B64_LMA.dat"  res
run "./fit_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../fit_smallV_B32_strange_fpi_A0P5.txt  "  res

run "./deriv_other_fpi  -p ../../data_fpi_A0   B32_OS_l0_l0_A0P5.dat  -bin 20  jack  B32_OS_l0_l0_P5P5.dat   reweight_strange_OS_B32.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B32_OS_l1_l1_A0P5.dat  -bin 20  jack  B32_OS_l1_l1_P5P5.dat   reweight_strange_OS_B32.dat  onlinemeas_B64_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   B32_OS_l2_l2_A0P5.dat  -bin 20  jack  B32_OS_l2_l2_P5P5.dat   reweight_strange_OS_B32.dat  onlinemeas_B64_LMA.dat"  res
run "./fit_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../fit_smallV_B32_OS_strange_fpi_A0P5.txt  "  res


run "./deriv_other_fpi  -p ../../data_fpi_A0   C48_tm_l0_l0_A0P5.dat  -bin 20  jack  C48_tm_l0_l0_P5P5.dat   reweight_strange_OS_C48.dat  onlinemeas_C80_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   C48_tm_l1_l1_A0P5.dat  -bin 20  jack  C48_tm_l1_l1_P5P5.dat   reweight_strange_OS_C48.dat  onlinemeas_C80_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   C48_tm_l2_l2_A0P5.dat  -bin 20  jack  C48_tm_l2_l2_P5P5.dat   reweight_strange_OS_C48.dat  onlinemeas_C80_LMA.dat"  res
run "./fit_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../fit_smallV_C48_strange_fpi_A0P5.txt  "  res

run "./deriv_other_fpi  -p ../../data_fpi_A0   C48_OS_l0_l0_A0P5.dat  -bin 20  jack  C48_OS_l0_l0_P5P5.dat   reweight_strange_OS_C48.dat  onlinemeas_C80_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   C48_OS_l1_l1_A0P5.dat  -bin 20  jack  C48_OS_l1_l1_P5P5.dat   reweight_strange_OS_C48.dat  onlinemeas_C80_LMA.dat"  res
run "./deriv_other_fpi  -p ../../data_fpi_A0   C48_OS_l2_l2_A0P5.dat  -bin 20  jack  C48_OS_l2_l2_P5P5.dat   reweight_strange_OS_C48.dat  onlinemeas_C80_LMA.dat"  res
run "./fit_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../fit_smallV_C48_OS_strange_fpi_A0P5.txt  "  res


run "./average_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../average_mul_over_fpi_df_dmus_smallV_tm_fpi_A0P5.txt  "  res
run "./average_smallV_phys_point_fpi_A0P5  jack ../../data_fpi_A0/jackknife/  ../../data_fpi_A0/fit_smallV_phys_point  ../average_mul_over_fpi_df_dmus_smallV_OS_fpi_A0P5.txt  "  res


echo "===================="
echo "Results:"

for i in "${res[@]}"; do
    echo "Failed test: $i"
done
if [ ${#res[@]} -eq 0 ]; then
    echo "all good"
fi