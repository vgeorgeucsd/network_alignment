#~/bin/bash

# magna++ settings
# the number of threads you wish to use
nThreads=4
# genetic algorithm settings
population=100
iterations=100

# balance between edge and node similarity, alpha as 0.5 is a balance between node and edge, alpha as 1 is only using node
alpha=1
# edge conservation options, they are: EC, ICS, S3
edgeConservation=EC

# The well you would like to analyze
wellname=well34
# your file paths
networkPath=/home/vivek/research/organoids/original_data/network_edge_lists
nodeSimilarityPath=/home/vivek/research/organoids/original_data/similarity_lists/${wellname}_lists
output_path=/home/vivek/research/organoids/output_data/magna/${wellname}_outputs

# your output prefixes
outFileName=${wellname}_p_${population}_i_${iterations}_a_${alpha}_temp_out

# graphBig=
# graphSmall=
# nodeSimilarityFile=

for filename in ${nodeSimilarityPath}/*.csv; do
  fname=`basename $filename`
  echo ${fname}
  d1=`echo $fname | cut -d'_' -f1 | cut -d'y' -f2`
  d2=`echo $fname | cut -d'_' -f2 | cut -d'.' -f1 | cut -d'y' -f2`
  pathN1=${networkPath}/${wellname}_day${d1}.txt
  pathN2=${networkPath}/${wellname}_day${d2}.txt
  v1=`awk '{for(i=1;i<=NF;i++) print $i}' ${pathN1} | sort -d | uniq | wc -w`
  v2=`awk '{for(i=1;i<=NF;i++) print $i}' ${pathN2} | sort -d | uniq | wc -w`
  echo "${d1}" >> temp.out
  echo "${d2}" >> temp.out

   if [ ${v1} -le ${v2} ]; then
     # G is the smaller network
     # H is the bigger network
     G=${pathN1}
     H=${pathN2}
     nSim=$filename
     # Removing un-needed info from files
     sed '1d' $filename | cut -d ',' -f 2- > tempNodeSimFile_${d1}_${d2}.csv

     ./magnapp_cli_linux64 -G ${G} -H ${H} -d tempNodeSimFile_${d1}_${d2}.csv -o ${outFileName}_${wellname}_${d1}_${d2} -m ${edgeConservation} -p ${population} -n ${iterations} -a ${alpha} -t ${nThreads} >> temp.out

     rm tempNodeSimFile_${d1}_${d2}.csv

   else
     G=${pathN2}
     H=${pathN1}
     nSim=$filename
     # Removing un-needed info from node similarity file and switching columns
     sed '1d' $filename | cut -d ',' -f 2- | awk -F, '{print $2,$1,$3}' OFS=, > tempNodeSimFile_${d2}_${d1}.csv

     ./magnapp_cli_linux64 -G ${G} -H ${H} -d tempNodeSimFile_${d2}_${d1}.csv -o ${outFileName}_${wellname}_${d2}_${d1} -m ${edgeConservation} -p ${population} -n ${iterations} -a ${alpha} -t ${nThreads} >> temp.out

     rm tempNodeSimFile_${d2}_${d1}.csv
   fi
done
mv temp* outputs/
# ./magnapp_cli_linux64 -G ${graphPath}/ex1.txt -H ${graphPath}/ex2.txt -d ${graphPath}/exgdvsim.csv -o ${outFileName} -m EC -p 10 -n 10 -a 0.5
# mv ${outFileName}* ${output_path}
