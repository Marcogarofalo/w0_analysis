#!/bin/sh

res=()

run(){
    local -n tmp=$2
    local exe=$1
    $exe || tmp+=("$exe")
}

# ./w0_m0  -p ../../data/  flow_B64_m0.dat  -bin 20 jack  bubbles/B64/P5_mu7.2000e-04.bin || exit 1
# ./w0_m0  -p ../../data/  flow_B64_m0.dat  -bin 20 jack  bubbles/B64/P5_mu1.8250e-02.bin || exit 1
# ./w0_m0  -p ../../data/  flow_B64_m0.dat  -bin 20 jack  bubbles/B64/P5_mu2.3134e-01.bin || exit 1

# ./w0_m0  -p ../../data/  flow_C80_m0.dat  -bin 20 jack  bubbles/C80/P5_mu6.0000e-04.bin || exit 1
# ./w0_m0  -p ../../data/  flow_C80_m0.dat  -bin 20 jack  bubbles/C80/P5_mu1.9849e-01.bin || exit 1
# ./w0_m0  -p ../../data/  flow_C80_m0.dat  -bin 20 jack  bubbles/C80/P5_mu1.6044e-02.bin || exit 1



run "./w0_m0  -p ../../data/  flow_B64_LMA.dat -bin 20 jack  ../../gm2-mistuning/m0data/P5bub/B64/P5_mu7.2000e-04.bin"  res
run "./w0_m0  -p ../../data/  flow_B64_LMA.dat -bin 20 jack  ../../gm2-mistuning/m0data/P5bub/B64/P5_mu1.8250e-02.bin"  res
run "./w0_m0  -p ../../data/  flow_B64_LMA.dat -bin 20 jack  ../../gm2-mistuning/m0data/P5bub/B64/P5_mu2.3770e-01.bin"  res


run "./w0_m0  -p ../../data/  flow_C80_LMA.dat -bin 20 jack  ../../gm2-mistuning/m0data/P5bub/C80/P5_mu6.0000e-04.bin"  res
run "./w0_m0  -p ../../data/  flow_C80_LMA.dat -bin 20 jack  ../../gm2-mistuning/m0data/P5bub/C80/P5_mu1.6044e-02.bin"  res
run "./w0_m0  -p ../../data/  flow_C80_LMA.dat -bin 20 jack  ../../gm2-mistuning/m0data/P5bub/C80/P5_mu1.9849e-01.bin"  res


run "./w0_m0  -p ../../data/  flow_D96.dat -bin 20 jack  ../../gm2-mistuning/m0data/P5bub/D96/P5_mu5.4000e-04.bin"  res
run "./w0_m0  -p ../../data/  flow_D96.dat -bin 20 jack  ../../gm2-mistuning/m0data/P5bub/D96/P5_mu1.3576e-02.bin"  res
run "./w0_m0  -p ../../data/  flow_D96.dat -bin 20 jack  ../../gm2-mistuning/m0data/P5bub/D96/P5_mu1.6490e-01.bin"  res


run "./w0_m0  -p ../../data/  flow_E112_LMA.dat -bin 20 jack  ../../gm2-mistuning/m0data/P5bub/E112/P5_mu4.4000e-04.bin"  res
run "./w0_m0  -p ../../data/  flow_E112_LMA.dat -bin 20 jack  ../../gm2-mistuning/m0data/P5bub/E112/P5_mu1.1759e-02.bin"  res
run "./w0_m0  -p ../../data/  flow_E112_LMA.dat -bin 20 jack  ../../gm2-mistuning/m0data/P5bub/E112/P5_mu1.4150e-01.bin"  res

# garofalo@lattice03:~/analysis/flow/w0_analysis/build$ ls ../../../gm2-mistuning/m0data/P5bub/E112/P5_mu
# P5_mu1.1759e-02.bin  P5_mu1.4150e-01.bin  P5_mu4.4000e-04.bin  
# garofalo@lattice03:~/analysis/flow/w0_analysis/build$ ls ../../../gm2-mistuning/m0data/P5bub/B64/P5_mu
# P5_mu1.8250e-02.bin  P5_mu2.3770e-01.bin  P5_mu7.2000e-04.bin  
# garofalo@lattice03:~/analysis/flow/w0_analysis/build$ ls ../../../gm2-mistuning/m0data/P5bub/D96/P5_mu
# P5_mu1.3576e-02.bin  P5_mu1.6490e-01.bin  P5_mu5.4000e-04.bin  

./fit_all_w0_m0  jack ../../data/jackknife/  ../../data/fit_all_beta  ../fit_files_w0_m0.txt

echo "===================="
echo "Results:"

for i in "${res[@]}"; do
    echo "Failed test: $i"
done
if [ ${#res[@]} -eq 0 ]; then
    echo "all good"
fi