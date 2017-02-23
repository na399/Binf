# Natthawut Adulyanukosol's exam submission for Linux and Python Programming 2016

import exam

# Q2.1
print '==================== Q2.1 ===================='

clinton_speeches_dict = exam.read_speeches('clinton_speeches.txt')
trump_speeches_dict = exam.read_speeches('trump_speeches.txt')

print 'The size of clinton_speeches_dict = ' + str(len(clinton_speeches_dict))
print 'The size of trump_speeches_dict = ' + str(len(trump_speeches_dict))

# Q2.2
print '==================== Q2.2 ===================='

clinton_speeches_list = clinton_speeches_dict.values()
trump_speeches_list = trump_speeches_dict.values()

clinton_speeches_all = exam.merge_speeches(clinton_speeches_list)
trump_speeches_all = exam.merge_speeches(trump_speeches_list)

print 'The length of the string of merged Clinton\'s speeches = ' + str(len(clinton_speeches_all))
print 'The length of the string of merged Trump\'s speeches = ' + str(len(trump_speeches_all))

# Q3.1
print '==================== Q3.1 ===================='

print 'The total number of words in Trump\'s speeches = ' + str(exam.count_words(trump_speeches_all))

# Q3.2
print '==================== Q3.2 ===================='

print 'The total number of sentences in Trump\'s speeches = ' + str(exam.count_sentences(trump_speeches_all))

# Q4.1
print '==================== Q4.1 ===================='

print 'The number of syllables in \"I eat apples\" = ' + str(exam.count_syllables("I eat apples"))

# Q4.2
print '==================== Q4.2 ===================='

clinton_scores = []
trump_scores = []

for title in clinton_speeches_dict:
    clinton_scores.append(exam.calculate_flesch_score(clinton_speeches_dict[title]))

for title in trump_speeches_dict:
    trump_scores.append(exam.calculate_flesch_score(trump_speeches_dict[title]))

print 'The Flesch scores of Clinton\'s speeches are'
print clinton_scores
print 'The Flesch scores of Trump\'s speeches are'
print trump_scores

# Q4.3

both_scores = {"Clinton": clinton_scores, "Trump": trump_scores}

exam.plot_scores(both_scores, "results.png")

# Q5
print '==================== Q5 ===================='

clinton_6gram_all = exam.calculate_ngram_frequencies(clinton_speeches_all, 6)
trump_6gram_all = exam.calculate_ngram_frequencies(trump_speeches_all, 6)


# define a function that returns the part of each item we want to sort by
def get_second(x):
    """returns the second item of x"""
    return x[1]


clinton_6gram_all_items = clinton_6gram_all.items()
clinton_6gram_all_items.sort(key=get_second)
print 'Clinton\'s 10 most frequent 6-grams are '
print clinton_6gram_all_items[::-1][0:10]

trump_6gram_all_items = trump_6gram_all.items()
trump_6gram_all_items.sort(key=get_second)
print 'Trump\'s 10 most frequent 6-grams are '
print trump_6gram_all_items[::-1][0:10]