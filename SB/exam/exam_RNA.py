def is_same_length(WT, MT):
    if len(WT) == len(MT):
        return True
    else:
        return False


def pair_list(seq):

    stack = []
    listPairs = []

    for i in range(0, len(seq)):

        if list(seq)[i] == '(':
            stack.append(i)
        if list(seq)[i] == ')':
            listPairs.append((stack[-1],i))
            del stack[-1]

    return listPairs


def hamming_distance(WT, MT):

    if is_same_length(WT,MT) == False:
        print("Lengths are different")
        return None
    else:
        distance = 0
        for i in range(0, len(WT)):
            if WT[i] != MT[i]:
                distance += 1

    return distance


def bp_distance(WT, MT):

    WT_pairs = pair_list(WT)
    MT_pairs = pair_list(MT)

    distance = 0

    for i in WT_pairs:
        if i not in MT_pairs:
            distance += 1

    for i in MT_pairs:
        if i not in WT_pairs:
            distance += 1

    return distance


def read_file(filename):

    # Open a  file
    the_file = open(filename)

    # Create a new dictionary
    the_list = []


    # Iterate over lines
    for line in the_file:
        # Replace whitespace, including /n, at the end of a line with a single space
        line = line.strip()

        # add line to the list
        the_list.append(line)

    # Close the file
    the_file.close()

    return the_list


def get_seq_list(eps_list):

    seq_list = []
    i = -1
    for item in eps_list:
        i += 1
        if item == '/sequence (\\':
            seq = eps_list[i+1]

    seq_list = list(seq[0:-1])

    return seq_list


def get_coor_list(eps_list):

    coor_list = []
    i = -1
    for item in eps_list:
        i += 1
        if item == '/coor [':
            while eps_list[i+1] != '] def':
                coor_list.append(eps_list[i+1])
                i += 1

    return coor_list


def get_prob_list(eps_list):

    prob_list = []
    i = -1
    for item in eps_list:
        i += 1
        if item == '/S [':
            while eps_list[i+1] != '] def':
                prob_list.append(float(eps_list[i+1]))
                i += 1

    return prob_list


def get_base_pair(eps_list):

    bp_list = []
    i = -1
    for item in eps_list:
        i += 1
        if item == '/pairs [':
            while eps_list[i + 1] != '] def':
                bp_list.append(eps_list[i + 1])
                i += 1

    base_pair_list = []
    for pair in bp_list:
        pair_mem_1 = int(pair.rsplit()[0][1:])
        pair_mem_2 = int(pair.rsplit()[1][:-1])
        base_pair_list.append((pair_mem_1, pair_mem_2))


    return base_pair_list


def get_processed_list(eps_list):

    seq_list = get_seq_list(eps_list)
    coor_list = get_coor_list(eps_list)
    prob_list = get_prob_list(eps_list)
    bp_list = get_base_pair(eps_list)

    processed_list = []


    for i in range(len(bp_list)):
        mem_1 = bp_list[i][0]
        mem_2 = bp_list[i][1]
        holder_dict = {mem_1: (seq_list[mem_1], coor_list[mem_1], prob_list[mem_1]),
                       mem_2: (seq_list[mem_2], coor_list[mem_2], prob_list[mem_2])}
        processed_list.append((0.5*(prob_list[mem_1]+prob_list[mem_2]), holder_dict))

    processed_list.sort(reverse=True)

    return processed_list


def dp_get_seq_list(eps_list):

    seq_list = []
    i = -1
    for item in eps_list:
        i += 1
        if item == '/sequence { (\\':
            seq = eps_list[i+1]

    seq_list = list(seq[0:-1])

    return seq_list





def dp_get_processed_list(eps_list):

    processed_list = []

    for line in eps_list:
        if line.endswith('ubox') and line[0] != '%':
            holder = line[:-5].rsplit()
            processed_list.append((float(holder[2]), holder[0], holder[1]))

    processed_list.sort(reverse=True)

    return processed_list


def mfe_question_1_2(filename):

    eps_list = read_file(filename)

    processed_list = get_processed_list(eps_list)

    print('=== start of the answer ===')

    print('\n')
    print('(average probability, {position in sequence: (letter, position in space, probability)})')

    for item in processed_list:
        print(item)

    print('\n')

    count = 0
    for item in processed_list:
        if item[0] > 0.8:
            count += 1

    print('The number of base pairs with an average probability of greater than 80%: ' + str(count))
    print('                                     from the total number of base pairs: ' + str(len(processed_list)))
    print('                                                                      or: ' + str(count/len(processed_list)*100) +' %')

    return print('=== end of the answer ===')

    print('\n')


def dp_question_1_2(filename):

    eps_list = read_file(filename)

    processed_list = dp_get_processed_list(eps_list)

    print('\n')
    print('=== start of the answer (dp) ===')
    print('(probability, base pair[0], base pair[1])')

    for item in processed_list:
        print(item)

    print('\n')

    count = 0
    for item in processed_list:
        if item[0] > 0.8:
            count += 1

    print('The number of base pairs with an average probability of greater than 80%: ' + str(count))
    print('                                     from the total number of base pairs: ' + str(len(processed_list)))
    print('                                                                      or: ' + str(count/len(processed_list)*100) +' %')

    return print('=== end of the answer ===')

    print('\n')




dp_question_1_2("RNA\sequence1_dp.eps")

dp_question_1_2("RNA\sequence2_dp.eps")


mfe_question_1_2("RNA\sequence1_pp.eps")


mfe_question_1_2("RNA\sequence2_pp.eps")



seq1    = 'AUCAGUUCUAGCAGGAGCUGUACUCAGAGACUCGGGAAAUUUUCCCGGAAUUUUACCCGGGUUUUUACGU'
seq1dot = '..(((((((....))))))).....(((((((((((.((..(((...)))..)).)))))))))))....'

seq2    = 'AUCGGUUCCAGCAGGAACUGUACUCGGGGGCUCGGGAAACCCUCCCGGGGUUUUACCCGGGUUUUUACGU'
seq2dot = '..(((((((....))))))).(((((((..(((((((.....)))))))......)))))))........'

print('\n')

print('Hamming distance    : ' + str(hamming_distance(seq1, seq2)))
print('Hamming distance(db): ' + str(hamming_distance(seq1dot, seq2dot)))
print('Base pair distance : ' + str(bp_distance(seq1dot, seq2dot)))
print('\n')


seqAli    = 'AUCAGUUCCAGCAGGAACUGUACUCAGAGACUCGGGAAACCCUCCCGGAAUUUUACCCGGGUUUUUACGU'

seqAliDot = '..(((((((....))))))).((..((((((((((((((((((...)))))))..)))))))))))..))'


print('Hamming distance    : ' + str(hamming_distance(seq1, seqAli)))
print('Hamming distance(db): ' + str(hamming_distance(seq1dot, seqAliDot)))
print('Base pair distance  : ' + str(bp_distance(seq1dot, seqAliDot)))
print('\n')

print('Hamming distance    : ' + str(hamming_distance(seq2, seqAli)))
print('Hamming distance(db): ' + str(hamming_distance(seq2dot, seqAliDot)))
print('Base pair distance  : ' + str(bp_distance(seq2dot, seqAliDot)))



