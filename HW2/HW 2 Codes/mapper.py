#! /usr/bin/python

import sys
import math

# input comes from STDIN (standard input)
for line in sys.stdin:
    # remove leading and trailing whitespace
    line = line.strip()
    # split the line into words
    numbers = line.split('\t')

    n = float((math.floor(10*float(numbers[0])))/10), float((math.floor(10*float(numbers[1])))/10)
    num_st1 = str(n[0])
    num_st2 = str(n[1])

    print '%s,%s\t%s' % (num_st1, num_st2, 1)