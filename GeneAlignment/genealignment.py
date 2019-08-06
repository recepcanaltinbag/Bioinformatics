#Gene Alignment Code
#
#Written by Recep Can Altınbağ
#
#1-2 August 2019
#PYTHON 3+

#README::::::::::::::::::::::::::::::::::::
#create input files gene1.txt and gene2.txt
#Then run the program
#It creates output.txt file

#Complexity: O(m+n)

import numpy as np

#---------FILE READING
f_gene1 = open('gene1.txt', 'r')
f_gene2 = open('gene2.txt', 'r')
readed_gene1 = f_gene1.read()
readed_gene2 = f_gene2.read()
f_gene1.close()
f_gene2.close()
print('Readed Genes: ')
print(readed_gene1,len(readed_gene1))
print(readed_gene2,len(readed_gene2))
#---------END OF FILE READING

#INITIALIZE
aligned_gene_1 = []
aligned_gene_2 = []
scores_array = np.zeros( (len(readed_gene1)+1, len(readed_gene2)+1 ))
udl_matrix = np.zeros( (len(readed_gene1)+1, len(readed_gene2)+1 ))
#Initialization: Max gene length to align(worst scenario)-> 10 000 000 000 (Nearly human genome)
scores_array[0][0]=0
left = -999999999
up = -999999999
diagonal = -999999999
#END OF THE INITILIZATION

#ALIGNMENT ALGORITHM SCORING MATRIX
# Scoring algorithm values:
# Mismatch = -1
# Match = 1
# Gaps = -2
for i in range(len(readed_gene1)+1):
    for j in range(len(readed_gene2)+1):
        dontlookleft = False #to safe in borders
        dontlookup = False

        if i==0 and j==0:
            continue
        if i-1 < 0:
            dontlookleft = True
        if j-1 < 0 :
            dontlookup = True

        upmatrixarray = []

        if dontlookleft == False:
            upmatrixarray.append((scores_array[i-1][j]-2.0,(i-1,j)))
            left = -2 + upmatrixarray[-1][0]
        else:
            l = -99999

        if dontlookup == False:
            upmatrixarray.append((scores_array[i][j-1]-2.0,(i,j-1)))
            up = -2 + upmatrixarray[-1][0]
        else:
            up = -99999

        if dontlookleft == False and dontlookup == False:
            if readed_gene1[i-1]==readed_gene2[j-1]:
                upmatrixarray.append((scores_array[i-1][j - 1] + 1.0,(i-1,j-1)))
                diagonal = 1 + upmatrixarray[-1][0]
            else:
                upmatrixarray.append((scores_array[i-1][j - 1] - 1.0,(i-1,j-1)))
                diagonal = -1 + upmatrixarray[-1][0]

        sortedlist=sorted(upmatrixarray,key=lambda l:l[0], reverse=True)
        udl_matrix[i][j] = np.argmax((up, diagonal, left))
        scores_array[i][j] = sortedlist[0][0]
#END OF THE SCORING PART


#ALIGNMENT ALGORITHM CREATING THE ALIGNED GENES WITH LOOKING UDL_MATRIX
# Meaning of the udl_matrix:
# 0 -> go up
# 1 -> go diagonal
# 2 -> go left
i = len(readed_gene1)
j = len(readed_gene2)
while i + j > 0 :

    if udl_matrix[i][j] == 0.0:
        aligned_gene_1.append('-')
        aligned_gene_2.append(readed_gene2[j-1])
        i = i
        j = j-1
    elif udl_matrix[i][j] == 1.0:
        if readed_gene1[i-1] == readed_gene2[j-1]:
            aligned_gene_1.append(readed_gene1[i-1])
            aligned_gene_2.append(readed_gene2[j-1])
        else:
            aligned_gene_1.append(readed_gene1[i-1]+'*')
            aligned_gene_2.append(readed_gene2[j-1]+'*')
        i = i-1
        j = j-1
    elif udl_matrix[i][j] == 2.0:
        aligned_gene_1.append(readed_gene1[i - 1])
        aligned_gene_2.append('-')
        i = i - 1
        j = j
#END OF THE ALGORITHM

#---------SCREEN WRITING
print('\n')
print('--- Final Results ---- ')
print(' * shows the mismatches')
aligned_gene_2.reverse()
aligned_gene_1.reverse()
print('\n')
print('Alignment of Gene 1 and Gene 2 :')
print(aligned_gene_1)
print(aligned_gene_2)
#---------END OF SCREEN WRITING


#---------FILE WRITING
output_file1 =  '--- Final Results ---- ' + '\n' + ' * shows the mismatches' + '\n' + 'Alignment Gene 1 and Gene 2 :'
output_file1 = str(output_file1) +'\n'  + str(aligned_gene_1) + '\n' + str(aligned_gene_2)
f_output = open('output.txt', 'w+')
f_output.write(output_file1)
f_output.close()

#---------END OF FILE WRITING
