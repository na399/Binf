# Q1
input = open('/home/lpp/.bashrc','r')           # open the file for reading

# Q2
bashrc = input.readlines()                      # store the file content as a list 'bashrc'
line_6th = bashrc[5]                            # get the 6th line
char_2nd = line_6th[1]                          # get the 2nd character of the 6th line
print char_2nd

# Q3
print len(bashrc)

# Q4
line_5th = bashrc[4]                            # get the 5th line
line_5th_words = line_5th.split(' ')            # make a list by splitting the string with ' '
line_5th_words.remove('#')                      # remove '#'
print len(line_5th_words)

# Q5
output = open('/home/lpp/junk.txt','w')         # open the file for writing
output.writelines(bashrc[1:5])                  # put the 2th line the 5th line from .bashrc to the new file

# Q6
input.close()
output.close()
