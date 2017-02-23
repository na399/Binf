#!/usr/bin/env python

import os
import Bio.PDB as PDB
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import itertools
import warnings


def save_dict(to_save_dict, filename):

    np.save(str(filename), to_save_dict)
    return None


def load_dict(filename):
    if filename[-4:] == ".npy":
        loaded_dict = np.load(filename).item()
    else:
        loaded_dict = np.load(filename+".npy").item()

    return loaded_dict


def load_list(filename):
    if filename[-4:] == ".npy":
        loaded_dict = np.load(filename)
    else:
        loaded_dict = np.load(filename+".npy")

    return loaded_dict

def find_k_mers(directory, length_aa=5, print_prot_name=0):
    """
    find k-mers

    :param directory: containing PDB files
    :param length_aa: length of sequence (k-mer) in amino acid residues DEFAULT = 5
    :param print_prot_name:
    :return CAs_coord_dict
    """

    # Make a list of the proteins from a collection (directory)
    protein_list = os.listdir(directory)

    # Make an empty dictionary to store the sequences of 5 amino acids (keys)
    CAs_coord_dict = {}


    # Iterate through all files in the directory
    for protein in protein_list:
        if print_prot_name == 1:
            print(protein)
        # Get a structure from file
        try:
            structure = parser.get_structure(protein, directory + "/" + protein)
        except:
            continue
        # Make an empty list to hold the residues for the structure
        residue_list = []

        # Iterate through all residues of a structure
        for residue in structure.get_residues():
            # Avoid non-standard residues with HET flags
            if residue.get_id()[0] == " ":
                # Add the standard residue to the list
                residue_list.append(residue)

        j = 1

        # Iterate through all sequences of length 5 amino acids
        for pos in range(1, len(residue_list) - (length_aa - 1)):
            # Get a sequence of length 5
            sequence_of_5 = residue_list[pos:(pos + length_aa)]
            sequence_of_5_name = ""

            # Iterate through all residues in the sequence
            for residue in sequence_of_5:
                sequence_of_5_name += residue.get_resname()

            # If the sequence is not already in the dictionary,
            # add it as a key with value of a dictionary with the protein name as a key
            if sequence_of_5_name not in CAs_coord_dict:
                CAs_coord_dict[sequence_of_5_name] = {protein: {"#1": {}}}
                j = 1
            # If the sequence is already existing in the dictionary, i.e. from other proteins,
            # add a new value to it as a dictionary with the protein name as a key
            else:
                if protein not in CAs_coord_dict[sequence_of_5_name]:
                    CAs_coord_dict[sequence_of_5_name][protein] = {}
                    CAs_coord_dict[sequence_of_5_name][protein]["#1"] = {}
                    j = 1
                else:
                    j = len(CAs_coord_dict[sequence_of_5_name][protein]) + 1
                    CAs_coord_dict[sequence_of_5_name][protein]["#" + str(j)] = {}

            i = 1

            for residue in sequence_of_5:
                for atom in residue.get_atom():
                    if atom.get_id() == "CA":
                        CAs_coord_dict[sequence_of_5_name][protein]["#" + str(j)]["CA" + str(i)] = atom.get_vector().get_array()

                i += 1

    return CAs_coord_dict




def get_k_mers(CAs_coord_dict, min_n=2):

    new_dict = {}

    for key in CAs_coord_dict:
        if len(CAs_coord_dict[key]) >= min_n:
            new_dict[key] = CAs_coord_dict[key]

    return new_dict


def print_dict(CAs_coord_dict):

    for key in CAs_coord_dict:
        print(key)
        print(len(CAs_coord_dict[key]))
        print(CAs_coord_dict[key])

    return None


def distance(CA1, CA2):

    dx = CA1[0] - CA2[0]
    dy = CA1[1] - CA2[1]
    dz = CA1[2] - CA2[2]

    distance = np.sqrt(dx*dx + dy*dy + dz*dz)

    return distance


def is_broken_chain(distance):
    """ Returns True if the two residues have
    more than 2 atoms that are less than 4 angstroms apart"""

    # Calculate distance using Atom class minus operator overload
    if distance > 4:
        return True

    # Finished comparing atoms and there were fewer than 4, so return False
    return False


def find_breaks_in_dict(dictionary, length_aa=5):

    to_remove_list = []

    count = 0
    count_break = 0

    for mer in dictionary:
        for protein in dictionary[mer]:
            for copy in dictionary[mer][protein]:
                if len(dictionary[mer][protein][copy]) != length_aa:
                    to_remove_list.append([mer, protein, copy])
                    break
                i = 1
                for CA in dictionary[mer][protein][copy]:
                    if i < len(dictionary[mer][protein][copy]):
                        if is_broken_chain(distance(dictionary[mer][protein][copy]["CA" + str(i)],
                                                    dictionary[mer][protein][copy]["CA" + str(i + 1)])):
                            count_break += 1
                            to_remove_list.append([mer, protein, copy])
                    i += 1
                    count += 1

    return to_remove_list


def remove_breaks(dictionary, to_remove_list):
    for item in to_remove_list:
        try:
            del dictionary[item[0]][item[1]][item[2]]
        except KeyError:
            continue
    return dictionary


def make_dict_without_breaks(dictionary, length_aa):

    to_remove_list = find_breaks_in_dict(dictionary, length_aa)
    new_dictionary = remove_breaks(dictionary, to_remove_list)

    return new_dictionary


def make_array_old(dictionary_copy):

    i = 1

    CAs_array = np.concatenate(([dictionary_copy["CA" + str(i)]],
                                [dictionary_copy["CA" + str(i + 1)]],
                                [dictionary_copy["CA" + str(i + 2)]],
                                [dictionary_copy["CA" + str(i + 3)]],
                                [dictionary_copy["CA" + str(i + 4)]],
                                ))

    return CAs_array


def make_tuple_CAs(dictionary_copy):

    list_CAs = []

    for i in range(len(dictionary_copy)):
        list_CAs.append([dictionary_copy["CA" + str(i + 1)]])

    tuple_CAs = tuple(list_CAs)


    return tuple_CAs


def make_array(tuple_CAs):
    try:
        CAs_array = np.concatenate(tuple_CAs)
    except ValueError:
        return None
    return CAs_array


def make_array_from_dict(dictionary_copy):

    array_CAs = make_array(make_tuple_CAs(dictionary_copy))

    return array_CAs


def calculate_rmsd(array_CAs_1, array_CAs_2, give_all_return=0):

    # register matrices
    X = array_CAs_1
    Y = array_CAs_2

    # Move X, Y to its center of mass
    X -= X.sum(0) / len(X)
    Y -= Y.sum(0) / len(Y)

    # transpose the matrices to get {3, 5} matrices
    X = np.transpose(X)
    Y = np.transpose(Y)

    # R=YXt
    R = np.dot(Y, np.transpose(X))

    ## Singular Value Decomposition
    # V, S, Wt=SVD(R)
    V, s, Wt = np.linalg.svd(R)  # S = np.diag(s) \ A = np.dot(np.dot(V, S), Wt)

    ## Check for reflection
    rotation = 0
    if (np.linalg.det(Wt) * np.linalg.det(np.transpose(V))) < 0.0:
        rotation = 1
        s[2] = -s[2]

    # Calculate RMSD (by formula)
    n = X.shape[1]
    E0 = sum(sum(np.square(X)) + sum(np.square(Y)))
    rmsd = np.sqrt(max([(1 / n) * (E0 - 2 * sum(s)), 0]))
    # do max() to make the inside of the log not negative
    # the goal is to maximise sum(s),
    # so it may produce the root of a negative value

    # Calculate U
    W = np.transpose(Wt)
    Vt = np.transpose(V)
    Z = np.identity(3)

    if rotation == 1:
        Z = np.diag([1,1,-1])

    U = np.dot(np.dot(W,Z),Vt)

    new_Y = np.dot(U,Y)

    if give_all_return == 1:
        return X, new_Y, rmsd, U

    return rmsd


def get_CAs_in_list_of_dicts(dictionary_mer, alt_return=0):
    """

    :param dictionary_mer:
    :return: list([{CA1:, CA2:, ..., CA5:}, {CA1:, CA2:, ..., CA5:}, ... ]) from the same mer
    """
    list_of_CAs = []

    for protein in dictionary_mer:
        for copy in dictionary_mer[protein]:
            list_of_CAs.append(dictionary_mer[protein][copy])

    return list_of_CAs


def make_mer_dict(directory, mer_length_aa, print_running_prot=1):

    output = str(directory) + '_' + str(mer_length_aa) + 'mers'

    CAs_coord_dict = find_k_mers(directory, mer_length_aa, print_running_prot)

    #save_dict(CAs_coord_dict, output)

    #CAs_coord_dict = load_dict(output)

    CAs_over2hits = get_k_mers(CAs_coord_dict)

    CAs_over2hits_no_breaks = make_dict_without_breaks(CAs_over2hits, mer_length_aa)

    save_dict(CAs_over2hits_no_breaks, output)

    return output


def make_rmsd_list(dict_npy_file, get_pdb=0):

    test_dict = load_dict(dict_npy_file)


    list_of_CAs = []
    rmsd_list = []
    rmsd_all_list = []
    count_mer = 0
    count_pair = 0
    count_random = 0

    if get_pdb != 0:
        if not os.path.exists('pdb_rmsd_'+ dict_npy_file):
            os.makedirs('pdb_rmsd_'+ dict_npy_file)

    for mer in test_dict:
        list_of_CAs = get_CAs_in_list_of_dicts(test_dict[mer])

        count_mer +=1

        for pair in itertools.combinations(range(len(list_of_CAs)), 2):

            count_pair += 1

            array_CAs_1 = make_array_from_dict(list_of_CAs[pair[0]])
            array_CAs_2 = make_array_from_dict(list_of_CAs[pair[1]])

            X, Y, rmsd, U = calculate_rmsd(array_CAs_1, array_CAs_2, 1)
            rmsd_list.append(rmsd)

            if get_pdb != 0 and count_pair < 10:

                structure = parser.get_structure("template", "template.pdb")

                io.set_structure(structure)

                i = 0
                for atom in structure.get_atoms():
                    new_pos = np.transpose(X[:, np.newaxis])[i][0]
                    atom.set_coord(new_pos)
                    i += 1

                io.save('pdb_rmsd_' + dict_npy_file + '\\'+ str(rmsd) + '_' + mer + '_X_' + str(test_dict[mer].keys())[9:] + '.pdb')

                structure = parser.get_structure("template", "template.pdb")

                io.set_structure(structure)

                i = 0
                for atom in structure.get_atoms():
                    new_pos = np.transpose(Y[:, np.newaxis])[i][0]
                    atom.set_coord(new_pos)
                    i += 1

                io.save('pdb_rmsd_' + dict_npy_file + '\\'+ str(rmsd) + '_' + mer + '_Y_' + str(test_dict[mer].keys())[9:] + '.pdb')

    if True:
        # Random pairs
        random_list_of_CAs = []
        list_holder = []

        for mer in test_dict:
            list_holder = get_CAs_in_list_of_dicts(test_dict[mer])
            for i in list_holder:
                random_list_of_CAs.append(i)

        to_run_pairs = np.random.choice(random_list_of_CAs, size=(len(rmsd_list),2))

        random_rmsd_list = []

        for pair in to_run_pairs:

            count_random +=1

            array_CAs_1 = make_array_from_dict(pair[0])
            array_CAs_2 = make_array_from_dict(pair[1])

            X, Y, rmsd, U = calculate_rmsd(array_CAs_1, array_CAs_2, 1)
            random_rmsd_list.append(rmsd)

            if get_pdb != 0 and count_random < 10:

                structure = parser.get_structure("template", "template.pdb")

                io.set_structure(structure)

                i = 0
                for atom in structure.get_atoms():
                    new_pos = np.transpose(X[:, np.newaxis])[i][0]
                    atom.set_coord(new_pos)
                    i += 1

                io.save('pdb_rmsd_' + dict_npy_file + '\\' + 'random_' + str(rmsd) + '_X' + '.pdb')

                structure = parser.get_structure("template", "template.pdb")

                io.set_structure(structure)

                i = 0
                for atom in structure.get_atoms():
                    new_pos = np.transpose(Y[:, np.newaxis])[i][0]
                    atom.set_coord(new_pos)
                    i += 1

                io.save('pdb_rmsd_' + dict_npy_file + '\\' + 'random_' + str(rmsd) + '_Y' + '.pdb')


    save_dict(rmsd_list, dict_npy_file + "_rmsd")
    save_dict(random_rmsd_list, dict_npy_file + "_rmsd_random")


    return dict_npy_file + "_rmsd"


def make_rmsd_list_old(dict_npy_file):

    test_dict = load_dict(dict_npy_file)
    print(len(test_dict))

    list_of_CAs = []
    rmsd_list = []
    count = 0

    for mer in test_dict:
        list_of_CAs = get_CAs_in_list_of_dicts(test_dict[mer])

        if False:
            for i in range(len(list_of_CAs)-1):
                array_CAs_1 = make_array_from_dict(list_of_CAs[i])
                array_CAs_2 = make_array_from_dict(list_of_CAs[i+1])

                rmsd = calculate_rmsd(array_CAs_1, array_CAs_2)
                rmsd_list.append(rmsd)

        if True:
            for pair in itertools.combinations(range(len(list_of_CAs)), 2):

                count += 1

                array_CAs_1 = make_array_from_dict(list_of_CAs[pair[0]])
                array_CAs_2 = make_array_from_dict(list_of_CAs[pair[1]])

                rmsd = calculate_rmsd(array_CAs_1, array_CAs_2)
                rmsd_list.append(rmsd)

    # Random pairs

    random_list_of_CAs = []
    list_holder = []

    for mer in test_dict:
        list_holder = get_CAs_in_list_of_dicts(test_dict[mer])
        for i in list_holder:
            random_list_of_CAs.append(i)

    to_run_pairs = np.random.choice(random_list_of_CAs, size=(len(rmsd_list),2))

    random_rmsd_list = []

    for pair in to_run_pairs:

        array_CAs_1 = make_array_from_dict(pair[0])
        array_CAs_2 = make_array_from_dict(pair[1])

        rmsd = calculate_rmsd(array_CAs_1, array_CAs_2)
        random_rmsd_list.append(rmsd)

    save_dict(rmsd_list, dict_npy_file + "_rmsd")
    save_dict(random_rmsd_list, dict_npy_file + "_rmsd_random")

    return dict_npy_file + "_rmsd"


def plot_rmsd_with_random_kde(rmsd_file):
    rmsd_list = load_list(rmsd_file)
    random_rmsd_list = load_list(rmsd_file + "_random")


    # Histograms
    # n, bins, patches = plt.hist(rmsd_list, normed=1)
    # n, bins, patches = plt.hist(random_rmsd_list, normed=1, alpha=0.5)

    fig, ax = plt.subplots(figsize=(8, 6))

    density = gaussian_kde(rmsd_list)
    xs = np.linspace(0, 6, 200)
    density.covariance_factor = lambda: .25
    density._compute_covariance()
    plt.plot(xs, density(xs), label="match (n="+ str(len(rmsd_list)) + ")")

    density = gaussian_kde(random_rmsd_list)
    xs = np.linspace(0, 6, 200)
    density.covariance_factor = lambda: .25
    density._compute_covariance()
    plt.plot(xs, density(xs), label="random (n="+ str(len(random_rmsd_list)) + ")")

    ax.grid(True)
    ax.legend(loc='right')
    ax.set_title(rmsd_file)
    ax.set_xlabel('RMSD')
    ax.set_ylabel('density')

    plt.show()

    return None

def plot_kde(rmsd_file):

    rmsd_list = load_list(rmsd_file)

    density = gaussian_kde(rmsd_list)
    xs = np.linspace(0, 6, 200)
    density.covariance_factor = lambda: .25
    density._compute_covariance()
    plt.plot(xs, density(xs), label= str(rmsd_file) + " (n=" + str(len(rmsd_list)) + ")")

    return None


warnings.filterwarnings("ignore")
plt.style.use('ggplot')

# Create parser
parser = PDB.PDBParser()

# Save pdb
io = PDB.PDBIO()


# Create plot from scratch, replace "top500H" with the directory containing PDB files
plot_rmsd_with_random_kde(make_rmsd_list(make_mer_dict("top500H", 5)))

# Take an existing dict of mers to calculate rmsd
# make_rmsd_list("500H_5mers", 0)

# .. as above with sample pdb files
# make_rmsd_list("500H_5mers", 1)

# Plot Kernel Density Estimation of pre-made lists of rmsd scores
# plot_rmsd_with_random_kde("500H_5mers_rmsd")

if False:

    fig, ax = plt.subplots(figsize=(8, 6))

    plot_kde("100H_3mers_rmsd")
    plot_kde("500H_4mers_rmsd")
    plot_kde("500H_5mers_rmsd")
    plot_kde("500H_6mers_rmsd")
    plot_kde("8000H_7mers_rmsd")
    plot_kde("8000H_8mers_rmsd")

    ax.grid(True)
    ax.legend(loc='right')
    ax.set_title('match')
    ax.set_xlabel('RMSD')
    ax.set_ylabel('density')

    plt.show()

if False:

    fig, ax = plt.subplots(figsize=(8, 6))

    plot_kde("500H_4mers_rmsd")
    plot_kde("500H_5mers_rmsd")
    plot_kde("500H_6mers_rmsd")

    ax.grid(True)
    ax.legend(loc='right')
    ax.set_title('match')
    ax.set_xlabel('RMSD')
    ax.set_ylabel('density')

    plt.show()



if False:

    fig, ax = plt.subplots(figsize=(8, 6))

    plot_kde("100H_3mers_rmsd_random")
    plot_kde("500H_4mers_rmsd_random")
    plot_kde("500H_5mers_rmsd_random")
    plot_kde("500H_6mers_rmsd_random")
    plot_kde("8000H_7mers_rmsd_random")
    plot_kde("8000H_8mers_rmsd_random")

    ax.grid(True)
    ax.legend(loc='right')
    ax.set_title('random')
    ax.set_xlabel('RMSD')
    ax.set_ylabel('density')

    plt.show()



if False:

    fig, ax = plt.subplots(figsize=(8, 6))

    plot_kde("500H_4mers_rmsd_random")
    plot_kde("500H_5mers_rmsd_random")
    plot_kde("500H_6mers_rmsd_random")

    ax.grid(True)
    ax.legend(loc='right')
    ax.set_title('random')
    ax.set_xlabel('RMSD')
    ax.set_ylabel('density')

    plt.show()
