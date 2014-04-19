import sys
class samTrack:
    def __init__(self, line):
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
def parseWgsimAnswer(seqName):
    import re
    pattern = re.compile(r'^(\S+)_(\d+)_(\d+)')
    chr, pos0, pos1 = pattern.search(seqName).groups() 
    pos0 = int(pos0)
    pos1 = int(pos1)
    return chr, pos0, pos1
    
def diffSam_wgsim_main(opt, arg):
    Distance = opt.max_diff
    Length = opt.read_length
    filename1, filename2 = arg[0], arg[1]
    
    fp1 = open(filename1, 'r')
    fp2 = open(filename2, 'r')
    if options.ans_file != '':
        fp_ans = open(options.ans_file)
    leftPair, rightPair = ['', ''], ['', '']
    #read header from file1
    #read header from file2
    skipSamHeader(fp1)
    skipSamHeader(fp2)
    #core loop
    fileEndFlag1, fileEndFlag2 = False, False
    while True:
        answerChr, answerpos0, answerpos1 = None, None, None
        samTrack0 , samTrack1= [None, None], [None, None]
        #read paired-end reads track from file1 and file2
        for i in range(2):
            line1 = fp1.readline()
            if len(line1) == 0:
                fileEndFlag1 = True
                break   #break for
            samTrack0[i] = samTrack(line1[:-1])
            leftPair[i] = line1[:-1]
            line2 = fp2.readline()
            if len(line2) == 0:
                fileEndFlag2 = True
                break   #break for
            samTrack1[i] = samTrack(line2[:-1])
            rightPair[i] = line2[:-1]
        if fileEndFlag1 == True or fileEndFlag2 == True:
            break #break while
        #parse answer
        if samTrack0[0].qname == samTrack0[1].qname == samTrack1[0].qname == samTrack1[1].qname:
            #print >>sys.stderr, '>nameSame'
            if options.ans_file == '':#parse answer from read name
                answerChr, answerpos0, answerpos1 = parseWgsimAnswer(samTrack0[0].qname) 
                answerpos1 -= Length-1
            else:#parse answer from answer file
                answerChr, answerpos0, answerpos1 = fp_ans.readline.strip().split('\t')
                answerpos1 -= Length-1
        else:
            print >>sys.stderr, '>nameError'
            print >>sys.stderr, samTrack0[0].qname
            print >>sys.stderr, samTrack0[1].qname
            print >>sys.stderr, samTrack1[0].qname
            print >>sys.stderr, samTrack1[1].qname
            continue
        #Check whether pos is correct
        if (samTrack0[0].rname, samTrack0[1].rname) == (samTrack1[0].rname, samTrack1[1].rname) and \
            abs(samTrack0[0].pos - samTrack1[0].pos) < Distance and \
            abs(samTrack0[1].pos - samTrack1[1].pos) < Distance:
            pass
        elif (samTrack0[0].rname, samTrack0[1].rname) == (samTrack1[1].rname, samTrack1[0].rname) and \
            abs(samTrack0[0].pos - samTrack1[1].pos) < Distance and \
            abs(samTrack0[1].pos - samTrack1[0].pos) < Distance:
            pass
        else:
            if samTrack0[0].rname == samTrack0[1].rname == answerChr and \
                ( abs(min(samTrack0[0].pos, samTrack0[1].pos) - answerpos0) < Distance and
                  abs(max(samTrack0[0].pos, samTrack0[1].pos) - answerpos1) < Distance):
                print >>sys.stderr, '>posDiff   1T_2F'
            elif samTrack1[0].rname == samTrack1[1].rname == answerChr and \
                ( abs(min(samTrack1[0].pos, samTrack1[1].pos) - answerpos0) < Distance and
                  abs(max(samTrack1[0].pos, samTrack1[1].pos) - answerpos1) < Distance):
                print >>sys.stderr, '>posDiff   1F_2T'
            else:
                print >>sys.stderr, '>posDiff   1F_2F'
            print >>sys.stderr, leftPair[0]
            print >>sys.stderr, leftPair[1]
            print >>sys.stderr, rightPair[0]
            print >>sys.stderr, rightPair[1]
    fp1.close()
    fp2.close()
    if options.ans_file != '':
        fp_ans.close()
import optparse

if __name__ == '__main__':
    #options init
    usage = "usage: %prog [Options] <file1> <file2>"
    parser = optparse.OptionParser(usage)
    parser.add_option('-d', '--diff', action = 'store', type = 'int',  dest='max_diff', help = 'max diff between different SAM', default= 4)
    parser.add_option('-l', '--length', action = 'store', type = 'int',  dest='read_length', help = 'read length', default= 100)
    parser.add_option('-a', '--answer', action = 'store', type = 'string',  dest='ans_file', help = 'answer file name', default= '')

    #get options
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        exit(1)
    diffSam_wgsim_main(options, args)







