t_scale=1.0
alpha=0.2
JHund=10.0

t=$(echo "${t_scale}*(1.0 - ${alpha})" | bc -l)
lambda=$(echo "${alpha}*${t_scale}" | bc -l)



# along a1
gammax=1.0
gammay=0

Jxx=$(echo "(1.0/${JHund})*((${t}*${t}) - (${lambda}*${lambda}) + 2.0*(${lambda}*${lambda}*${gammay}*${gammay}))" | bc -l)
Jxy=$(echo "(1.0/${JHund})*(-2.0*(${lambda}*${lambda})*${gammax}*${gammay})" | bc -l)
Jxz=$(echo "(1.0/${JHund})*(2.0*${t}*${lambda}*${gammax})" | bc -l)


Jyx=$(echo "(1.0/${JHund})*(-2.0*(${lambda}*${lambda})*${gammax}*${gammay})" | bc -l)
Jyy=$(echo "(1.0/${JHund})*((${t}*${t}) - (${lambda}*${lambda}) + 2.0*(${lambda}*${lambda}*${gammax}*${gammax}))" | bc -l)
Jyz=$(echo "(1.0/${JHund})*(2.0*${t}*${lambda}*${gammay})" | bc -l)


Jzx=$(echo "(1.0/${JHund})*(-2.0*${t}*${lambda}*${gammax})" | bc -l)
Jzy=$(echo "(1.0/${JHund})*(-2.0*${t}*${lambda}*${gammay})" | bc -l)
Jzz=$(echo "(1.0/${JHund})*((${t}*${t}) - (${lambda}*${lambda}))" | bc -l)


echo "${Jxx} ${Jxy} ${Jxz}" > J1_a1.txt
echo "${Jyx} ${Jyy} ${Jyz}" >> J1_a1.txt
echo "${Jzx} ${Jzy} ${Jzz}" >> J1_a1.txt




# along a2
gammax=$(echo "0.5" | bc -l)
gammay=$(echo "sqrt(3.0)*0.5" | bc -l)

Jxx=$(echo "(1.0/${JHund})*((${t}*${t}) - (${lambda}*${lambda}) + 2.0*(${lambda}*${lambda}*${gammay}*${gammay}))" | bc -l)
Jxy=$(echo "(1.0/${JHund})*(-2.0*(${lambda}*${lambda})*${gammax}*${gammay})" | bc -l)
Jxz=$(echo "(1.0/${JHund})*(2.0*${t}*${lambda}*${gammax})" | bc -l)


Jyx=$(echo "(1.0/${JHund})*(-2.0*(${lambda}*${lambda})*${gammax}*${gammay})" | bc -l)
Jyy=$(echo "(1.0/${JHund})*((${t}*${t}) - (${lambda}*${lambda}) + 2.0*(${lambda}*${lambda}*${gammax}*${gammax}))" | bc -l)
Jyz=$(echo "(1.0/${JHund})*(2.0*${t}*${lambda}*${gammay})" | bc -l)


Jzx=$(echo "(1.0/${JHund})*(-2.0*${t}*${lambda}*${gammax})" | bc -l)
Jzy=$(echo "(1.0/${JHund})*(-2.0*${t}*${lambda}*${gammay})" | bc -l)
Jzz=$(echo "(1.0/${JHund})*((${t}*${t}) - (${lambda}*${lambda}))" | bc -l)


echo "${Jxx} ${Jxy} ${Jxz}" > J1_a2.txt
echo "${Jyx} ${Jyy} ${Jyz}" >> J1_a2.txt
echo "${Jzx} ${Jzy} ${Jzz}" >> J1_a2.txt






# along a3
gammax=$(echo "-0.5" | bc -l)
gammay=$(echo "sqrt(3.0)*0.5" | bc -l)

Jxx=$(echo "(1.0/${JHund})*((${t}*${t}) - (${lambda}*${lambda}) + 2.0*(${lambda}*${lambda}*${gammay}*${gammay}))" | bc -l)
Jxy=$(echo "(1.0/${JHund})*(-2.0*(${lambda}*${lambda})*${gammax}*${gammay})" | bc -l)
Jxz=$(echo "(1.0/${JHund})*(2.0*${t}*${lambda}*${gammax})" | bc -l)


Jyx=$(echo "(1.0/${JHund})*(-2.0*(${lambda}*${lambda})*${gammax}*${gammay})" | bc -l)
Jyy=$(echo "(1.0/${JHund})*((${t}*${t}) - (${lambda}*${lambda}) + 2.0*(${lambda}*${lambda}*${gammax}*${gammax}))" | bc -l)
Jyz=$(echo "(1.0/${JHund})*(2.0*${t}*${lambda}*${gammay})" | bc -l)


Jzx=$(echo "(1.0/${JHund})*(-2.0*${t}*${lambda}*${gammax})" | bc -l)
Jzy=$(echo "(1.0/${JHund})*(-2.0*${t}*${lambda}*${gammay})" | bc -l)
Jzz=$(echo "(1.0/${JHund})*((${t}*${t}) - (${lambda}*${lambda}))" | bc -l)


echo "${Jxx} ${Jxy} ${Jxz}" > J1_a3.txt
echo "${Jyx} ${Jyy} ${Jyz}" >> J1_a3.txt
echo "${Jzx} ${Jzy} ${Jzz}" >> J1_a3.txt


