#RECEP CAN ALTINBAĞ


f_gene = open('gene_data.txt', 'r')
readed_gene = f_gene.readlines()
f_gene.close()

my_dataset = []

for element in readed_gene:
    s = element.rstrip()
    my_dataset.append((s.split(' ')[0],s.split(' ')[1]))

print(my_dataset)

back_region = 0.0
gc_rich_region = 0.0
back_to_rich = 0.0
rich_to_back = 0.0

rich_region_nucleotides = []
back_region_nucleotides = []

if my_dataset[0][1] == 'N':
    back_region_nucleotides.append(my_dataset[0][0])
else:
    rich_region_nucleotides.append(my_dataset[0][0])

#COUNTS FROM DATA

for element in range(1,len(my_dataset)):
    if my_dataset[element][1] == 'N':
        if my_dataset[element][1] == my_dataset[element-1][1]:
            back_region = back_region + 1.0
        else:
            rich_to_back = rich_to_back + 1.0

        back_region_nucleotides.append(my_dataset[element][0])

    if my_dataset[element][1] == 'R':
        if my_dataset[element][1] == my_dataset[element-1][1]:
            gc_rich_region = gc_rich_region + 1.0
        else:
            back_to_rich = back_to_rich + 1.0

        rich_region_nucleotides.append(my_dataset[element][0])

back_region = back_region / (back_region + back_to_rich)
back_to_rich = 1.0 - back_region
gc_rich_region = gc_rich_region / (gc_rich_region + rich_to_back)
rich_to_back = 1.0 - gc_rich_region

print('back region prob:',back_region)
print('gc:',gc_rich_region)

print(rich_region_nucleotides)
print(back_region_nucleotides)

#GENERAL PROBS:
P_A_rich = rich_region_nucleotides.count('A')/len(rich_region_nucleotides)
P_T_rich = rich_region_nucleotides.count('T')/len(rich_region_nucleotides)
P_C_rich = rich_region_nucleotides.count('C')/len(rich_region_nucleotides)
P_G_rich = rich_region_nucleotides.count('G')/len(rich_region_nucleotides)

P_A_back = back_region_nucleotides.count('A')/len(back_region_nucleotides)
P_T_back = back_region_nucleotides.count('T')/len(back_region_nucleotides)
P_C_back = back_region_nucleotides.count('C')/len(back_region_nucleotides)
P_G_back = back_region_nucleotides.count('G')/len(back_region_nucleotides)

P_A_bayes_back = P_A_back / (P_A_back + P_A_rich)
P_A_bayes_rich = 1 - P_A_bayes_back
P_T_bayes_back = P_T_back / (P_T_back + P_T_rich)
P_T_bayes_rich = 1 - P_T_bayes_back
P_C_bayes_back = P_C_back / (P_C_back + P_C_rich)
P_C_bayes_rich = 1 - P_C_bayes_back
P_G_bayes_back = P_G_back / (P_G_back + P_G_rich)
P_G_bayes_rich = 1 - P_G_bayes_back

#BAYESIAN PROBS:
Bayesian_Dict = {
    "A" : P_A_bayes_back,
    "T" : P_T_bayes_back,
    "C" : P_C_bayes_back,
    "G" : P_G_bayes_back
}

#GENERAL PROBS:
being_back = rich_to_back / (1 - back_region + rich_to_back)
being_rich = 1 - being_back

print('Being Rich region: ',being_rich)
print('Being Back region',being_back)

#TESTING WITH TEST DATA
f_gene = open('test_data.txt', 'r')
readed_gene = f_gene.read()
f_gene.close()

print('TEST DATA: ',readed_gene)

#VITERBI ALGORITHM

#INIT VALUES
max_being_back = being_back
max_being_rich = being_rich

output_region = []
print(Bayesian_Dict)

temp_back_value = max_being_back * Bayesian_Dict[readed_gene[0]]
temp_rich_value = max_being_rich * (1 - Bayesian_Dict[readed_gene[0]])

if temp_rich_value > temp_back_value:
    output_region.append('R')
else:
    output_region.append('N')

#VITERBI ALGORITHM LOOP
for element in range(1, len(readed_gene)):
    b = []
    r = []

    b.append(temp_back_value * back_region * Bayesian_Dict[readed_gene[element]])
    b.append(temp_rich_value * rich_to_back * Bayesian_Dict[readed_gene[element]])
    r.append(temp_rich_value * gc_rich_region * (1 - Bayesian_Dict[readed_gene[element]]))
    r.append(temp_back_value * back_to_rich * (1 - Bayesian_Dict[readed_gene[element]]))


    if max(b) > max(r):
        output_region.append('N')
    else:
        output_region.append('R')

    temp_back_value = max(b)
    temp_rich_value = max(r)

output_gene = list(readed_gene)
print('\n Final Result:')
print(output_gene)
print(output_region)

#---------FILE WRITING
output_file1 =  '--- Final Results ---- ' + '\n' + ' R --> GC Rich region   N --> Not rich region'
output_file1 = str(output_file1) +'\n'  + str(output_gene) + '\n' + str(output_region)
f_output = open('output.txt', 'w+')
f_output.write(output_file1)
f_output.close()

#---------END OF FILE WRITING
