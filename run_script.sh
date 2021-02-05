#~/bin/bash
networkPath=/home/vivek/research/organoids/original_data/original_ucsd_well_networks
nodeSimilarityPath=/home/vivek/research/organoids/original_data/similarity_lists/well22_lists
output_path=/home/vivek/research/organoids/output_data/magna/well22_outputs

outFileName=test

graphBig=
graphSmall=
nodeSimilarityFile=

./magnapp_cli_linux64 -G ${graphPath}/ex1.txt -H ${graphPath}/ex2.txt -d ${graphPath}/exgdvsim.csv -o ${outFileName} -m EC -p 10 -n 10 -a 0.5
mv ${outFileName}* ${output_path}
