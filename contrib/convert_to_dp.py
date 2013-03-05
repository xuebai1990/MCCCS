#!/usr/bin/env python

import re

def convert_into_dp(line):
    floatRe=re.compile(r"([0-9]+\.?[0-9]*|\.[0-9]+)[eEdD]([-+]?[0-9]+)")
    line=floatRe.sub(r"\1E\2_dp",line)
    return line

if __name__ == '__main__':
    import sys,os

    if len(sys.argv)<2:
        print "Usage: ",sys.argv[0]," file1 [file2 ...]"
        quit()

    for fileName in sys.argv[1:]:
        try:
           tmpName=fileName+'.tmp'
           outfile=open(tmpName,'w')
           for line in open(fileName,'r'):
               outfile.write(convert_into_dp(line))
           outfile.close()
           os.rename(tmpName, fileName)
        except:
           print "Error for ", fileName

