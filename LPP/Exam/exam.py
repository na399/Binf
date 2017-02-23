# Natthawut Adulyanukosol's exam submission for Linux and Python Programming 2016

def read_speeches(filename):
    """
    Make a dictionary of speeches
    where keys are the title of speeches,
    and the values are the corresponding speeches

    Parameters
    ----------

    filename: str
        The name of the file (the path of the file)

    """

    # Open a speech file
    speech_file = open(filename)

    # Create a new dictionary
    speech_dict = {}

    # Iterate over lines
    for line in speech_file:
        # Replace whitespace, including /n, at the end of a line with a single space
        line = line.rstrip() + ' '

        # Given that a title begins with #
        if line.startswith('#'):
            # Remove '# ' at the beginning and ': ' at the end, to be used as a title
            title = line[2:-2]
            # Assign the tile as a key in the dictionary
            speech_dict[title] = ''
        # A speech line does not begins with #
        else:
            # Not begins with [ either
            if line.startswith('[') is False:
                # Append the speech line to the already existing string of the corresponding title
                # The tile variable is kept from the previous loop(s)
                speech_dict[title] += line

    # Close the file
    speech_file.close()

    return speech_dict


def merge_speeches(speeches_list):
    """Merge a list of speeches to a single string containing all speeches,
    separated by a single space

    Parameters
    ----------

    speeches_list: list
        A list with many string members

    """

    # Create a new string variable
    speeches_string = ''

    # Iterate over speeches in the given list
    for speech in speeches_list:
        # Append the speech and add a single space at the end
        speeches_string += speech + ' '

    return speeches_string


def count_words(text):
    """Count the number of words in the text

    Parameters
    ----------

    text: str
        A string of words

    """

    import re

    # Make a list of words (contiguous non-whitespace characters)
    word_list = re.findall(r'\S+', text)
    # Find the size of the list
    count = len(word_list)

    return count


def count_sentences(text):
    """Count the number of sentences in the text

    Parameters
    ----------

    text: str
        A string of sentences

    """

    import re

    # Make a list of sentences (separated by either '.', '!' or '?')
    sentence_list = re.split(r'[.!?]', text)
    # Find the size of the list
    count = len(sentence_list)

    return count


def count_syllables(text):
    """Count the number of syllables in the text

    Parameters
    ----------

    text: str
        A string of one or more English word(s)

    """

    import re

    # Make a list of vowel sounds presenting in the text (converted to lower-case letters)
    syllable_list = re.findall(r'[aiouy]+e*|e(?!d\b|ly)[aiouye]?|[td]ed|le\b', text.lower())
    # Find the size of the list
    count = len(syllable_list)

    return count


def calculate_flesch_score(text):
    """Calculate Flesch score from a given text

    Parameters
    ----------

    text: str
        A string in English

    """

    # Use the given formula and call other functions to get the values
    score = 206.835 - 1.015 * count_words(text) / count_sentences(text) - 84.6 * count_syllables(text) / count_words(
        text)

    return score


def plot_scores(scores,
                output_filename):
    """
    Plots score distributions

    Parameters
    ----------

    scores : dict object
        The values of this dictionary are lists of scores to be plotted. The keys of the dictionary
        will be used as labels for the plot.

    output_filename : str
        The name of the file used where the plot will be saved.

    """

    import matplotlib.pyplot as plt

    x_values = range(1, len(scores) + 1)
    y_values = scores.values()

    # Create a box plot
    plt.boxplot(y_values, showmeans=True)

    # Add labels on x-axis
    plt.setp(plt.gca(), xticks=x_values, xticklabels=scores.keys())

    # Add title and labels
    plt.title('Flesch scores of each candidate\'s speeches')
    plt.ylabel('Flesch scores')
    plt.xlabel('Candidate')

    # Add horizontal grid lines
    plt.grid(axis='y',  linestyle='--', which='major', color='grey', alpha=0.5)

    # Save figure
    plt.savefig(output_filename)

    return


def calculate_ngram_frequencies(text, n):
    """Calculate ngram frequencies and produce a dictionary of n-grams and their corresponding number of occurrences

    Parameters
    ----------

    text: str
        A string

    n: int
        The size of n-gram

    """

    import re

    # Create a new dictionary
    ngram_dict = {}

    # Find all sentences
    sentences_list = re.findall(r'[^\.\?!"]+', text)

    # Iterate over sentences in the list
    for sentence in sentences_list:
        # Split words by a whitespace character
        words_list = sentence.rsplit()

        # Iterate over ngrams in the sentence
        for i in range(len(words_list) - n + 1):

            # Join the words to size of n
            ngram = ' '.join(words_list[i:i + n])

            # Record the presence of a new ngram
            if not ngram in ngram_dict:
                ngram_dict[ngram] = 1

            # Add the number of occurrence of the ngram
            elif ngram in ngram_dict:
                ngram_dict[ngram] += 1

    return ngram_dict