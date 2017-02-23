#!/usr/bin/python


import os
import sys


original_file = sys.argv[1]

# add .gz to the file name
compressed_file = original_file + '.gz'


# check whether the file name is correct
def checkfilename():
    global compressed_file, correct_file
    print 'The compressed file is ' + compressed_file + '?'
    correct_file = raw_input('y/n: ')

    if correct_file != 'y':
        compressed_file = raw_input('Type your file name \nThe compressed file is ')
        checkfilename()

    if correct_file == 'y':
        pass


checkfilename()


# calculation
original_size = float(os.path.getsize(original_file))
compressed_size = float(os.path.getsize(compressed_file))

compression = ((original_size - compressed_size) / original_size) * 100

print 'The original size is ' + str(original_size) + ' bytes.'
print 'The compressed size is ' + str(compressed_size) + ' bytes.'
print 'The compression percentage is ' + str(compression) + '%.'
