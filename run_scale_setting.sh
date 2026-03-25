#!/bin/sh

res=()

run(){
    local -n tmp=$2
    local exe=$1
    $exe || tmp+=("$exe")
}


run "
./set_w0_a  jack ../../data/jackknife/  ../../data/fit_all_beta  ../data_ens_B64
" res
run "
./set_w0_a  jack ../../data/jackknife/  ../../data/fit_all_beta  ../data_ens_C80
" res
run "
./set_w0_a  jack ../../data/jackknife/  ../../data/fit_all_beta  ../data_ens_D96
" res
run "
./set_w0_a  jack ../../data/jackknife/  ../../data/fit_all_beta  ../data_ens_E112
" res


run "
./fit_cont_w0_fpi   jack ../../data/jackknife/  ../../data/fit_all_beta
" res
echo "===================="
echo "Results:"

for i in "${res[@]}"; do
    echo "Failed test: $i"
done
if [ ${#res[@]} -eq 0 ]; then
    echo "all good"
fi