import sys
class samTrack:
    def __init__(self, line):
        self.line = line.strip()
        split_line = line.split('\t')
        self.qname = split_line[0]
        self.flag = int(split_line[1]) 
        self.rname = split_line[2] 
        self.pos = int(split_line[3])
        self.mapq = int(split_line[4])
        self.cigar = split_line[5]
        self.mrname = split_line[6]
        self.mpos = int(split_line[7])
        self.isize = int(split_line[8])
        self.seq = split_line[9]
        self.qual = split_line[10]
        if len(split_line) >11:
            self.opt= split_line[11]
        else:
            self.opt = None
def skipSamHeader(fp):
    line = fp.readline()
    while line[0] == '@':
        line = fp.readline()
    fp.seek(-len(line),1)

def diffSam_main(opt, arg): 
    Distance = int(opt.max_diff)
    
    filename1, filename2 = arg[0], arg[1]

    fp1 = open(filename1, 'r')
    fp2 = open(filename2, 'r')

   # leftPair, rightPair = ['', ''], ['', '']

    #read header from file1
    #read header from file2
    skipSamHeader(fp1)
    skipSamHeader(fp2)
    #core loop
    n = 0
    fileEndFlag1, fileEndFlag2 = False, False
    while True:
    
        samTrack0 , samTrack1= [None, None], [None, None]
        #read 2 lines from file1 and file2
        for i in range(2):
            
            line1 = fp1.readline()
            if len(line1) == 0:
                fileEndFlag1 = True
                break   #break for
            samTrack0[i] = samTrack(line1[:-1])
            
            line2 = fp2.readline()
            if len(line2) == 0:
                fileEndFlag2 = True
                break   #break for
            samTrack1[i] = samTrack(line2[:-1])
        
        if fileEndFlag1 == True or fileEndFlag2 == True:
            break #break while 
        
        if samTrack0[0].qname == samTrack0[1].qname == samTrack1[0].qname == samTrack1[1].qname:
            #print >>sys.stderr, '>nameSame' 
            pass
        else:
            print '@nameError' 
            print samTrack0[0].qname
            print samTrack0[1].qname
            print samTrack1[0].qname
            print samTrack1[1].qname
            if opt.print_SamTrack:
                print samTrack0[0].line
                print samTrack0[1].line
                print samTrack1[0].line
                print samTrack1[1].line
            continue

        
        #print samTrack0[0].pos
        #print samTrack0[1].pos
        #print samTrack1[0].pos
        #print samTrack1[1].pos

        if (samTrack0[0].rname, samTrack0[1].rname) == (samTrack1[0].rname, samTrack1[1].rname) and \
            abs(samTrack0[0].pos - samTrack1[0].pos) < Distance and \
            abs(samTrack0[1].pos - samTrack1[1].pos) < Distance:
            #print samTrack0[0].pos, samTrack1[0].pos, samTrack0[1].pos, samTrack1[1].pos
            pass
        elif (samTrack0[0].rname, samTrack0[1].rname) == (samTrack1[1].rname, samTrack1[0].rname) and \
             abs(samTrack0[0].pos - samTrack1[1].pos) < Distance and \
             abs(samTrack0[1].pos - samTrack1[0].pos) < Distance:
            #print samTrack0[0].pos, samTrack1[1].pos, samTrack0[0].pos, samTrack1[1].pos
            pass
        else:
            print '@posDiff'
            print samTrack0[0].rname, samTrack0[0].pos
            print samTrack0[1].rname, samTrack0[1].pos
            print samTrack1[0].rname, samTrack1[0].pos
            print samTrack1[1].rname, samTrack1[1].pos
            if opt.print_SamTrack:
                print samTrack0[0].line
                print samTrack0[1].line
                print samTrack1[0].line
                print samTrack1[1].line
        n += 2
        if n %100000 == 0:
            print >>sys.stderr, '%d lines finished!'%(n)
    fp1.close()
    fp2.close()

import optparse

if __name__ == '__main__':
    #options init
    usage = "usage: %prog [Options] <file1> <file2>"
    parser = optparse.OptionParser(usage)
    parser.add_option('-d', '--diff', action = 'store', type = 'int',  dest='max_diff', help = 'max diff between different SAM', default= 4)
    parser.add_option('-p', '--print', action = 'store_true', dest='print_SamTrack', help = 'print sam track', default= False)

    #get opt
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        exit(1)
    diffSam_main(options, args)







