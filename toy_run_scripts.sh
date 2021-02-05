#~/bin/bash
graphPath=/home/vivek/research/organoids/sample_data/networks/toy/lists
output_path=/home/vivek/research/organoids/output_data/magna/toy_outputs

outFileName=no_double_quotes
./magnapp_cli_linux64 -G ${graphPath}/ex1.txt -H ${graphPath}/ex2.txt -d ${graphPath}/exgdvsim.csv -o ${outFileName} -m EC -p 10 -n 10 -a 0.5
mv ${outFileName}* ${output_path}
