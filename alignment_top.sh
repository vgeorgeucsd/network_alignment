#~/bin/bash
OUTPUTDIR=/home/vikash/ssd2/SecondDrive/output_data
GAPS=0
WELL_NUMBER=2
POPULATION=1
ITERATIONS=0


well_time_path=/home/vikash/ssd2/SecondDrive/network_edge_lists/copy_test/
#for filename in ${well_time_path}well${WELL_NUMBER}*.txt; do
#    /home/vikash/ssd2/SecondDrive/network_alignment/Przulj_Code/a.out ${filename}
#done

if [ ! -d ${well_time_path}ODV/ ]
then
    mkdir ${well_time_path}ODV
fi
mv ${well_time_path}*.txt.* ${well_time_path}ODV

sed "s|CHANGE_INPUT_PATH|${well_time_path}ODV/|g;s|CHANGE_OUTPUT_PATH|${OUTPUTDIR}|g;s|CHANGE_GAP|${GAPS}|g;s|CHANGE_WELL_NUMBER|${WELL_NUMBER}|g" GDS_Vikash_old.py > GDS_Vikash_${WELL_NUMBER}.py

if [[ ${GAPS} -eq 0 ]]
then
    sed "s|WELL_NUMBER|${WELL_NUMBER}|g;s|POP_NUMBER|${POPULATION}|g;s|IT_NUMBER|${ITERATIONS}|g;s|GAP|${GAPS}|g;s|SIMILARITY_LIST_PATH|${OUTPUTDIR}/GDS_days/well${WELL_NUMBER}_lists|g;s|MAGNA_OUTPUT_PATH|${OUTPUTDIR}/magna/well${WELL_NUMBER}_outputs|g;s|NETWORK_PATH|${well_time_path}|g" process_files_Vikash.sh > process_files_Vikash_${WELL_NUMBER}.sh
else
    sed "s|WELL_NUMBER|${WELL_NUMBER}|g;s|POP_NUMBER|${POPULATION}|g;s|IT_NUMBER|${ITERATIONS}|g;s|GAP|${GAPS}|g;s|SIMILARITY_LIST_PATH|${OUTPUTDIR}/GDS_gap_days/GDS_gap_${GAPS}_days/well${WELL_NUMBER}_lists|g;s|MAGNA_OUTPUT_PATH|${OUTPUTDIR}/magna/well${WELL_NUMBER}_outputs|g;s|NETWORK_PATH|${well_time_path}|g" process_files_Vikash.sh > process_files_Vikash_${WELL_NUMBER}.sh
fi

python GDS_Vikash_${WELL_NUMBER}.py

chmod u+x process_files_Vikash_${WELL_NUMBER}.sh

./process_files_Vikash_${WELL_NUMBER}.sh

rm GDS_Vikash_${WELL_NUMBER}.py
rm process_files_Vikash_${WELL_NUMBER}.sh
