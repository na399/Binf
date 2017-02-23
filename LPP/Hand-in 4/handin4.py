# coding=utf-8
# 1. Open this file, and read it in to a list of lists of floats (that is, convert the strings of each line to a list of two numeric values)

input = open('/home/lpp/experimental_results.txt', 'r')
inputlines = input.readlines()                                      # a list with each item from each line: ['a1 a2\n', 'b1 b2\n', ...]
input.close()

results = []
for line in inputlines:
    inputlines_split = line.replace("\n","").split(" ")             # create a list from a pair of str results from a row: [['a1','a2'],['b1','b2'], ...]
    holder = []                                                     # create a blank list to hold a pair of float results from the following for loop
    for item in inputlines_split:
        holder.append(float(item))                                  # convert each str result into float and put it in the holder in order: [a1, a2]
    results.append(holder)                                          # append the pair to the list of results, i.e., [[a1, a2], [b1 b2], ...]


# 2. Calculate and print the average value for each column by iterating over the rows.

sum_column = [0.0, 0.0]                                             # each member will be the sum of each respective column
average_column = [0.0, 0.0]                                         # each member will be the average of each respective column

for row in results:                                                 # for row by row
    for column in range(2):                                         # for column by column
        sum_column[column] = sum_column[column] + row[column]       # add the value to the sum

for column in range(2):                                             # for column by column
    average_column[column] = sum_column[column]/len(results)        # get the average by dividing sum with the number of values

print average_column

# 3. Write these two average values to a file called "results.txt" (remember to close the output file, since sometimes the values don't appear in the file until the file has been closed).  Include this file with your turn-in.

output = open('/home/lpp/results.txt','w')                          # open a new file with write permission
output.writelines("%s " % column for column in average_column)      # put average from each column in the file
output.close()



# 4. Now turn the list of lists around so that the outer list has two elements with each containing a list of all the values in the corresponding column. In other words, if I want to access the 26th value in the 2nd column, I would write: column_list[1][25].
# Hint: Note that these two nested lists are just two different ways of representing the same table of data. In the first case, we specify the row first, and each row is a list of two column values. In the second case, we specify the column first, and each column is a list of all the row values. Here is a sketch of the two variants:

column_list = [[],[]]                                               # create a blank list with 2 members representing the initial 2 columns

for row in range(len(results)):                                     # for row by row
    for column in range(2):                                         # for column by column
        column_list[column].append(results[row][column])            # append a value to its respective column in the column_list


# 5. Calculate the average value for each column again. Hint: you can now use the built-in sum function to do this easily.

average_column2 = [0.0, 0.0]                                        # each member will be the average of each respective column

for column in range(2):
    average_column2[column] = sum(column_list[column])/len(column_list[column])

print average_column2

# 6. Finally, turn this list of lists into a dictionary of lists so that you can access the element from before as: column_dict['column1'][25]

column_dict = {'column1':[],'column2':[]}                           # create a new dictionary with keys: 'column1' and 'column2', both associated with blank lists

for column in range(2):                                             # for row by row
    for row in range(len(column_list[column])):                     # for column by column
        column_dict['column'+str(column + 1)].append(column_list[column][row])  # append a new value to the appropriate list

print column_dict['column1'][25]