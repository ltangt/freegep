import sys
import os
import random
def feature(maxv):
    rand = random.Random()
    return [rand.randint(0,maxv), rand.randint(0,maxv), rand.randint(0,maxv)]

def measure(fs):
    return fs[0]*2-fs[1]-fs[2]

if __name__ == '__main__' :
    fh = file('dataset.txt','w')
    for i in range(100) :
        row = feature(10)
        row.append(measure(row))
        for elem in row:
            fh.write(str(elem)+' ')
        fh.write('\n')
    fh.close()
