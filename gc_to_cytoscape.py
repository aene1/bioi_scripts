import os
import glob
import Bio
from Bio import SeqIO
import csv
import numpy as np

# summary_files = glob.glob("StaphAureusPhasters/*/summary.txt")
#

#">phage_i_j" where j is the jth gene in phage i

gc_folders = glob.glob("ecoli_clusters/*")

#ecoli

list1 = []
#matrix and edge list making code. First read the results in dictionary.
phage_dict = {}
for i in gc_folders:
    for j in SeqIO.parse(i,"fasta"):
        gene_cluste = i[15:]
        gene_cluster = gene_cluste[:(gene_cluste.find("."))]
        # phage_id = j.id[:str(j.id).find("_",8)]
        phage_id = j.id[:str(j.id).find("_",6)]
        list1.append(phage_id)
        genome_name = phage_id
        if genome_name not in phage_dict.keys():
            phage_dict[genome_name]=list()
        phage_dict[genome_name].append(gene_cluster)

# #pseduomona
# list_seq = list(SeqIO.parse('Paeruginosa_phage/pangenome_all.csv', "fasta"))
#
# phage_dict = {}
# for i in list_seq:
#     line = i.id.split(("|"))
#     gene_cluster = line[1][13:]
#     genome_name = line[2][12:]
#    # print(genome_name)
#     if genome_name not in phage_dict.keys():
#         phage_dict[genome_name]=list()
#     phage_dict[genome_name].append(gene_cluster)


genomes_list=list(phage_dict.keys())

values=list()
for i in range(len(genomes_list)):
    values.append(list())

outfile=open('edge_list_ecoli_final.csv','w')
for i in range(len(genomes_list)):
    genome_x=genomes_list[i]
    genes_for_x=phage_dict[genome_x]
    for j in range(len(genomes_list)):
        genome_y = genomes_list[j]
        genes_for_y = phage_dict[genome_y]
        intersection = [value for value in genes_for_x if value in genes_for_y]
        values[i].append(len(intersection))


#write an edge list
for i in range(len(genomes_list)):
    for j in range(len(genomes_list)):
        if (int(str(values[i][j])) != 0):
            if (str(genomes_list[i] != str(genomes_list[j]))):
        # print(genomes_list[i])
        # print(","+ genomes_list[j])
        # print(','+str(values[i][j])+ "\n")
                outfile.write(genomes_list[i])
                outfile.write(","+ genomes_list[j])
                outfile.write(','+str(values[i][j])+ "\n")


#write a matrix
# for i in range(len(genomes_list)):
#     outfile.write(genomes_list[i])
#     for j in range(len(genomes_list)):
#         outfile.write(','+str(values[i][j]))
#     outfile.write('\n')
#
# outfile.close()



