import os
import sys

def main(n):
    with open ('file' + str(n) + '.txt','w') as filehd:
        for i in range(10000000):
            filehd.write(str(i))

if __name__ == '__main__':
    if len(sys.argv) < 1:
        print 'no parameter given'
    main(sys.argv[1]) 
