#! /usr/bin/python

import sys

current_bin = None
current_freq = 0
bin = None

# Input comes from STDIN
for line in sys.stdin:
    #remove trailing '\n'
    line = line.strip()
    # extract key value from the output of Mapper
    bin, freq = line.split('\t')

    # Convert freq (currently string)into int
    try:
        freq = int(freq)
    except ValueError:
        # count was not a number, so silently
        # ignore/discard this line
        continue

    # Assuming that Hadoop has sorted the output from Mapper before passing it to reducer
    if current_bin == bin:
        current_freq += freq
    else:
        if current_bin:
            x_lo, y_lo = bin.split(',')
            x_hi = str(float(x_lo) + 0.1)
            y_hi = str(float(y_lo) + 0.1)
            print '%s,%s,%s,%s,%s' % (x_lo,x_hi,y_lo,y_hi,current_freq)
        current_freq = freq
        current_bin = bin

# Printing the last one
if current_bin == bin:
    x_lo, y_lo = bin.split(',')
    x_hi = str(float(x_lo) + 0.1)
    y_hi = str(float(y_lo) + 0.1)
    print '%s,%s,%s,%s,%s' % (x_lo,x_hi,y_lo,y_hi,current_freq)

