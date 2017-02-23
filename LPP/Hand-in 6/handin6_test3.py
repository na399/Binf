import handin6

"""
file1_list = handin6.fasta_names_to_list("test1.fasta")
file2_dict = handin6.fasta_names_to_dict("test2.fasta")
"""
file1_list = handin6.fasta_names_to_list("uniprot_sprot_2015.fasta")
file2_dict = handin6.fasta_names_to_dict("uniprot_sprot_2007.fasta")


#set the number of matches = 0
n = 0

#iterate the first list
for ID1 in file1_list:
    if ID1 in file2_dict:
        #add 1 to the number of matches
        n += 1
        print n
            

#print the number of unique IDs in the first list (unmatched)
print len(file1_list) - n

