import os 
edge_directory = '/home/vikash/ssd2/SecondDrive/network_edge_lists/'

well = 22

for i in range(250):
    edge_file = edge_directory + 'well%d_day%d.txt' % (well,i)
    
    if os.path.isfile(edge_file) == False:
        continue
    output_file = edge_directory + 'well%d_day%d_SANA.el' % (well,i)

    edge_list = list()

    with open(edge_file,'r') as file:
        for row in file:
            words=row.split()
            if [words[1],words[0]] not in edge_list:
                edge_list.append(words)
    f = open(output_file,'w')
    for i in range(len(edge_list)):
        f.write('%s %s\n' % (edge_list[i][0],edge_list[i][1]))
    f.close()
