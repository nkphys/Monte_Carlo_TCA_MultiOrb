for ms in {0..10}
do
seed=$(echo "1+${ms}"| bc -l)

cp input.inp input_run.inp
sed -i -e "s/VALUE_SEED/${seed}/g" input_run.inp

./sf input_run.inp > out_run.txt
mv ThetaPhi_Temp0.0001000000MicroState0.txt ThetaPhi_Temp0.0001000000MicroState${ms}_alongZ.txt

echo "${ms} done"
done

