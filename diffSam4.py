import sys
def diffSam_main(opt, arg): 
    Distance = opt.max_diff
    filename1 = arg[0]
    filename2 = arg[1]
    fp1 = open(filename1, 'r')
    fp2 = open(filename2, 'r')
    #l1 = len(fp1.readlines())
    #l2 = len(fp2.readlines())
    #if l1 != l2 or (l1-4)%2 != 0:
    #	print >>sys.stderr, 'file1 line num != file line num'
    #	sys.exit(1)
    leftPair = [0,0]
    rightPair = [0, 0]
    #fp1.seek(0)
    #fp2.seek(0)

    line1 = fp1.readline()
    while line1[0] == '@':
        line1 = fp1.readline()
    fp1.seek(-len(line1),1)
    line2 = fp2.readline()
    while line2[0] == '@':
        line2 = fp2.readline()
    fp2.seek(-len(line2),1)
    flag = True
    while flag:
        leftPair[0] = ''
        leftPair[1] = ''
        rightPair[0] = ''
        rightPair[1] = ''
        for i in range(2):
            line1 = fp1.readline()
            if len(line1) == 0:
                flag = False
                break
            leftPair[i] = line1[:-1]
            
            line2 = fp2.readline()
            if len(line2) == 0:
                flag = False
                break
            rightPair[i] = line2[:-1]

        leftPos0 = int(leftPair[0].split('\t')[3])
        leftPos1 = int(leftPair[1].split('\t')[3])
        rightPos0 = int(rightPair[0].split('\t')[3])
        rightPos1 = int(rightPair[1].split('\t')[3])
        if (abs(leftPos0 - rightPos0)<Distance and abs(leftPos1 - rightPos1) <Distance) or (abs(leftPos0 - rightPos1)<Distance and abs(leftPos1 - rightPos0)<Distance):
            pass
        else:
            print >>sys.stderr, '>>diff!'
            print >>sys.stderr, leftPair[0]
            print >>sys.stderr, leftPair[1]
            print >>sys.stderr, rightPair[0]
            print >>sys.stderr, rightPair[1]		
    fp1.close()
    fp2.close()

import optparse

if __name__ == '__main__':
    #options Init
    usage = "usage: %prog [Options] <file1> <file2>"
    parser = optparse.OptionParser(usage)
    parser.add_option('-d', '--diff', action = 'store', type = 'int',  dest='max_diff', help = 'max diff between different SAM', default= 4)

    #get opt
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        exit(1)
    diffSam_main(options, args)







