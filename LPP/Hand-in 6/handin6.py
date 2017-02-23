def fasta_label_to_id(label):
    """Extract meaningful ID from fasta label"""
    import re
    pattern = re.compile(r'>\w+[ |(]+(\w+)[ |)]')
    match = pattern.match(label)
    return match.group(1)

def fasta_names_to_list(path):
    """Take one argument, a string that is a path to a fasta file
    then create and return a list of protein ID"""

    file = open(path)
    fasta_list = file.readlines()
    ID_list = []

    for line in fasta_list:

        if line[0] == '>':
            ID_list.append(fasta_label_to_id(line))

    return ID_list



def fasta_names_to_dict(path):
    """Take one argument, a string that is a path to a fasta file
    then create and return a dictionary of protein ID"""

    file = open(path)
    fasta_list = file.readlines()
    ID_dict = dict()

    for line in fasta_list:

        if line[0] == '>':
            ID_dict[fasta_label_to_id(line)] = ''

    return ID_dict

def binary_search(sorted_list, element):
    """Search for element in list using binary search.
       Assumes sorted list"""
    # Current active list runs from index_start up to
    # but not including index_end
    index_start = 0
    index_end = len(sorted_list)
    while (index_end - index_start) > 0:
        index_current = (index_end-index_start)//2 + index_start
        if element == sorted_list[index_current]:
            return True
        elif element < sorted_list[index_current]:
            index_end = index_current
        elif element > sorted_list[index_current]:
            index_start = index_current+1
    return False

