# cook your dish here
import numpy as np
import matplotlib.pyplot as plt
import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.phylo as phylo
import biotite.sequence.graphics as graphics

# Obtain BLOSUM62
matrix = align.SubstitutionMatrix.std_protein_matrix()

type =input("enter type 1 or 2 \n 1-protein \n or \n 2-DNA ")
seq1 = input("enter the first sequence")
seq2 = input("enter the second sequence")
sz1 = len(seq1)  # virtical
sz2 = len(seq2)  # horizontal
algin_matrix = np.zeros((sz1 + 1, sz2 + 1))
for i in range(sz1 + 1):
    algin_matrix[i][0] = (-1 * i)

for i in range(sz2 + 1):
    algin_matrix[0][i] = (-1 * i)

# print(algin_matrix)
for i in range(1, sz1 + 1):
    for j in range(1, sz2 + 1):
        left = algin_matrix[i][j - 1] + (-1)
        up = algin_matrix[i - 1][j] + (-1)
        protein= -123456
        if type ==1:
            protein =matrix[seq1[i-1]][seq2[j-1]]
        diagonal = 0;
        if seq1[i - 1] == seq2[j - 1]:  # because we are one based
            diagonal = algin_matrix[i - 1][j - 1] + 1  # match
        else:
            diagonal = algin_matrix[i - 1][j - 1] - 2  # mismatch
        algin_matrix[i][j] = max(protein,left, up, diagonal)

#trce back

ans1 = ""
ans2 = ""
i = sz1
j = sz2
while(i>0 and j>0):
    if type ==1 and algin_matrix[i-1][j-1]+matrix[i-1][j-1]==algin_matrix[i-1][j-1]:
        ans1 += seq1[i - 1]
        ans2 += seq2[j - 1]
        i -= 1
        j -= 1
    elif seq1[i-1]==seq2[j-1]:#match
        ans1+=seq1[i-1]
        ans2+=seq2[j-1]
        i-=1
        j-=1
    elif algin_matrix[i-1][j-1]-2==algin_matrix[i-1][j-1]:#missmatch
        ans1+=seq1[i-1]
        ans2+=seq2[j-1]
        i-=1
        j-=1


    else:
        if algin_matrix[i-1][j]==algin_matrix[i][j]+1:#second gap
            ans1+=seq1[i-1]
            ans2+='_'
            i-=1
        else :
            ans1+='_'
            ans2+=seq2[j-1];
            j-=1


ans1 =ans1[::-1]
ans2=ans2[::-1]
#print()
#print(algin_matrix)
print('\n'+ans1+"\n"+ans2)