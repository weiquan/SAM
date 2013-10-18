import sys
def diffSam_main(opt, arg): 
    Distance = opt.max_diff
    
    filename1, filename2 = arg[0], arg[1]

    fp1 = open(filename1, 'r')
    fp2 = open(filename2, 'r')

    leftPair, rightPair = ['', ''], ['', '']

    #read header from file1
    line1 = fp1.readline()
    while line1[0] == '@':
        line1 = fp1.readline()
    fp1.seek(-len(line1),1)
    #read header from file2
    line2 = fp2.readline()
    while line2[0] == '@':
        line2 = fp2.readline()
    fp2.seek(-len(line2),1)
    #core loop
    fileEndFlag1, fileEndFlag2 = False, False
    while True:
        leftPair[0], leftPair[1], rightPair[0], rightPair[1] = '', '', '', ''
        
        #read 2 lines from file1 and file2
        for i in range(2):
            
            line1 = fp1.readline()
            if len(line1) == 0:
                fileEndFlag1 = True
                break   #break for
            leftPair[i] = line1[:-1]
            
            line2 = fp2.readline()
            if len(line2) == 0:
                fileEndFlag2 = True
                break   #break for
            rightPair[i] = line2[:-1]
        
        if fileEndFlag1 == True or fileEndFlag2 == True:
            break #break while 
        
        l0, l1, r0, r1 = leftPair[0].split('\t'), leftPair[1].split('\t'), rightPair[0].split('\t'), rightPair[1].split('\t')
        
        if l0[0] == l1[0] == r0[0] == r1[0]:
            pass
        else:
            print >>sys.stderr, '>nameError' 
            print >>sys.stderr, leftPair[0]
            print >>sys.stderr, leftPair[1]
            print >>sys.stderr, rightPair[0]
            print >>sys.stderr, rightPair[1]		
            continue

        leftChr0, leftPos0 = l0[2], int(l0[3]) 
        leftChr1, leftPos1 = l1[2], int(l1[3]) 
        rightChr0, rightPos0 = r0[2], int(r0[3]) 
        rightChr1, rightPos1 = r1[2], int(r1[3]) 
        if (leftChr0, leftChr1 == rightChr0, rightChr1 and abs(leftPos0 - rightPos0) < Distance and abs(leftPos1 - rightPos1) < Distance):
            pass
        elif(leftChr0, leftChr1 == rightChr1, rightChr0 and abs(leftPos0 - rightPos1) < Distance and abs(leftPos1 - rightPos0) < Distance):
            pass
        else:
            print >>sys.stderr, '>posDiff'
            print >>sys.stderr, leftPair[0]
            print >>sys.stderr, leftPair[1]
            print >>sys.stderr, rightPair[0]
            print >>sys.stderr, rightPair[1]		

    fp1.close()
    fp2.close()

import optparse

if __name__ == '__main__':
    #options init
    usage = "usage: %prog [Options] <file1> <file2>"
    parser = optparse.OptionParser(usage)
    parser.add_option('-d', '--diff', action = 'store', type = 'int',  dest='max_diff', help = 'max diff between different SAM', default= 4)

    #get opt
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        exit(1)
    diffSam_main(options, args)







