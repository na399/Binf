# Natthawut Adulyanukosol's exam submission for Linux and Python Programming 2016

# Download the speech files
wget http://people.binf.ku.dk/wb/data/clinton_speeches.txt
wget http://people.binf.ku.dk/wb/data/trump_speeches.txt

# 1. Write a Unix command that counts the number of speeches in the Trump speech file.
echo ==================== Q1 ====================
echo The number of speeches in Trump file = 

# Search the Trump speech file for lines beginning with # (-e '^#') 
# Then print the number of such lines (-c)
# -c can also be replaced by piping the grep result to wc -l
grep -c -e '^#' trump_speeches.txt 


# 2.  Write a Unix command that counts the number of lines in the Trump speech file, 
# Ignoring empty lines and lines starting with # or [.
echo ==================== Q2 ====================
echo The number of lines in Trump file = 
 
# Search for lines NOT (-v) beginning with either # (-e '^#') or [ (-e '^\[') 
# Then print the number of such lines (-c) 
grep -v -c -e '^#' -e '^\[' trump_speeches.txt 


# 3. Write a Unix command that uses grep to show only the first line in each speech in the Trump file 
# (i.e. the line after the title).
echo ==================== Q3 ====================
echo The first line in each speech in Trump file:

# Search for lines beginning with # (-e '^#')
# Print them along with each line's the line after it (-A 1)
# Discard the lines beginning with either # or group separator (--) (-v -e '^#' -e '^--$')
grep -A 1 -e '^#' trump_speeches.txt | grep -v -e '^#' -e '^--$' 

# 4. Write Unix commands that counts how many times Trump mentions Clinton, 
# And how many times Clinton mentions Trump, respectively  (both in upper and lower case). 
# Remember that there might be multiple occurrences in a line.
echo ==================== Q4 ====================
echo The number of times Trump mentions Clinton =

# Edit a line with Clinton (both cases) by having the rest of line after Clinton as a new line (\n)
# Count how many lines that have Clinton (ignoring case)
sed 's/Clinton/Clinton\n/g;s/clinton/clinton\n/g' trump_speeches.txt | grep -c -i -e 'Clinton'


echo The number of times Clinton mentions Trump =

# As the above commnand
sed 's/Trump/Trump\n/g;s/trump/trump\n/g' clinton_speeches.txt | grep -c -i 'Trump'



# 5. Write a Unix command that counts the number of sentences in the Trump speech file, 
# Ignoring lines starting with # or [. 
# You can define a sentence as Any number of characters that is not a period, followed by a period.
echo ==================== Q5 ====================
echo The number of sentences in Trump file = 

# Find a period preceded by alphanumeric character(s) followed by a space ([[:alnum:]]\+\. ) 
# and then add a new line after them (&\n)
grep -v -e '^#' -e '^\[' trump_speeches.txt | sed 's/[[:alnum:]]\+\. /&\n/g' | wc -l


# Assumption: a period must be followed by a space and an uppercase letter ([[:alnum:]]\+\. \)\([[:upper:]]\)
echo with improved commands \(see codes and comments\)
grep -v -e '^#' -e '^\[' trump_speeches.txt | sed 's/\([[:alnum:]]\+\. \)\([[:upper:]]\)/\1\n\2/g' | wc -l


# Assumption: not only a period can end a sentence, also a question mark (?) and an exclamation mark (!)
grep -v -e '^#' -e '^\[' trump_speeches.txt | sed 's/\([[:alnum:]]\+[.?!] \)\([[:upper:]]\)/\1\n\2/g' | wc -l


# The regular expression can be further improved to recognise true sentences. 

