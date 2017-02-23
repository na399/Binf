import handin6

"""
file1_list = handin6.fasta_names_to_list("test1.fasta")
file2_list = handin6.fasta_names_to_list("test2.fasta")
"""
file1_list = handin6.fasta_names_to_list("uniprot_sprot_2015.fasta")
file2_list = handin6.fasta_names_to_list("uniprot_sprot_2007.fasta")


file2_list.sort()

#set the number of matches = 0
n = 0

#iterate the first list
for ID1 in file1_list:

    if handin6.binary_search(file2_list, ID1) is True:
        n += 1
        print n


print len(file1_list) - n

