import os
import sys

thelist = os.listdir('Z:\Compounds')
thefile = open('test.txt', 'w') 
for item in thelist:
  print>>thefile, item
#for item in thelist:
#  thefile.write("%s\n" % item)