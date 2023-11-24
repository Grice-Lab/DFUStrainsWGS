
echo "Missing" > gapped176.txt

for fname in /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/Patient176/GeneFastas/MAFFT/Mafft_*; do

	awk 'FNR==2{if(/-/) print FILENAME}' $fname >> gapped176.txt

done

