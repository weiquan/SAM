import sys

def readfq(fp):
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break
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
        if len(split_line) > 11:
            self.opt = split_line[11:]
        else:
            self.opt = []
    def parseXA(self):
        aln_list = []
        for opt in self.opt:
            if opt[:2] == 'XA':
                aln_list = opt.split(':')[2].split(';')[:-1]
        return aln_list
def skipSamHeader(fp):
    line = fp.readline()
    while line[0] == '@':
        line = fp.readline()
    fp.seek(-len(line), 1)

def parseCigar(cigar):
    _encode = 'MIDNSHP'

    result = []
    n = ''
    for c in cigar:

        if c.isdigit():
            n += c
        elif c in _encode:
            if n == '':
                print cigar
                raise ValueError("end of CIGAR string reached, but an operator was expected")
            result.append((c, int(n)))
            n = ''  
    return result
def parseWgsimAnswer(seqName):
    import re
    pattern = re.compile(r'^(\S+)_(\d+)_(\d+)_')
    chr, pos0, pos1 = pattern.search(seqName).groups() 
    pos0 = int(pos0)
    pos1 = int(pos1)
    return chr, pos0, pos1
def lev1(seq1, seq2):
    oneago = None
    thisrow = range(1, len(seq2) + 1) + [0]
    for x in xrange(len(seq1)):
        twoago, oneago, thisrow = oneago, thisrow, [0] * len(seq2) + [x + 1]
        for y in xrange(len(seq2)):
            delcost = oneago[y] + 1
            addcost = thisrow[y - 1] + 1
            subcost = oneago[y - 1] + (seq1[x] != seq2[y])
            thisrow[y] = min(delcost, addcost, subcost)
    return thisrow[len(seq2) - 1]

def diffSam_wgsim_main(opt, arg):
    dictRef = {}
    if opt.ref != '':
        fp_ref = open(opt.ref)
        for name, seq, qual in readfq(fp_ref):
            if name in dictRef:
                print >>sys.stderr, 'more than one seq in %s named with %s'%(opt.ref, name)
                sys.exit(1)
            dictRef[name] = seq
        fp_ref.close()


    

    Distance = opt.max_diff
    Length = opt.read_length
    filename1, filename2 = arg[0], arg[1]

    fp1 = open(filename1, 'r')
    fp2 = open(filename2, 'r')
    if opt.ans_file != '':
        fp_ans = open(opt.ans_file)
    leftPair, rightPair = ['', ''], ['', '']
    #read header from file1
    #read header from file2
    skipSamHeader(fp1)
    skipSamHeader(fp2)
    #core loop
    fileEndFlag1, fileEndFlag2 = False, False
    while True:
        answerChr, answerpos0, answerpos1 = None, None, None
        samTrack0, samTrack1 = [None, None], [None, None]
        #read paired-end reads track from file1 and file2
        for i in range(2):
            line1 = fp1.readline()
            if len(line1) == 0:
                fileEndFlag1 = True
                break  # break for
            samTrack0[i] = samTrack(line1[:-1])
            leftPair[i] = line1[:-1]
            line2 = fp2.readline()
            if len(line2) == 0:
                fileEndFlag2 = True
                break  # break for
            samTrack1[i] = samTrack(line2[:-1])
            rightPair[i] = line2[:-1]
        if fileEndFlag1 is True or fileEndFlag2 is True:
            break   # break while
        #parse answer
        if samTrack0[0].qname == samTrack0[1].qname ==\
           samTrack1[0].qname == samTrack1[1].qname:
            #print >>sys.stderr, '>nameSame'
            if opt.ans_file == '':  # parse answer from read name
                answerChr, answerpos0, answerpos1 =\
                    parseWgsimAnswer(samTrack0[0].qname)
                answerpos1 -= Length-1
            else:  # parse answer from answer file
                answerChr, answerpos0, answerpos1 =\
                    fp_ans.readline().strip().split('\t')
                answerpos1 -= Length-1
        else:
            print >>sys.stderr, '>nameError'
            print >>sys.stderr, samTrack0[0].qname
            print >>sys.stderr, samTrack0[1].qname
            print >>sys.stderr, samTrack1[0].qname
            print >>sys.stderr, samTrack1[1].qname
            continue
        #Check whether pos is correct
        if (samTrack0[0].rname, samTrack0[1].rname) ==\
           (samTrack1[0].rname, samTrack1[1].rname) and\
           abs(samTrack0[0].pos - samTrack1[0].pos) < Distance and\
           abs(samTrack0[1].pos - samTrack1[1].pos) < Distance:
            pass
        elif(samTrack0[0].rname, samTrack0[1].rname) ==\
            (samTrack1[1].rname, samTrack1[0].rname) and\
            abs(samTrack0[0].pos - samTrack1[1].pos) < Distance and\
                abs(samTrack0[1].pos - samTrack1[0].pos) < Distance:
            pass
        else:
            flag_map_mate0 = 0
            flag_map_mate1 = 0

            if samTrack0[0].rname == samTrack0[1].rname == answerChr and \
                (abs(min(samTrack0[0].pos, samTrack0[1].pos) - answerpos0) < Distance and
                    abs(max(samTrack0[0].pos, samTrack0[1].pos) - answerpos1) < Distance):
                aln_list0 = samTrack1[0].parseXA()
                for string in aln_list0:
                    chrom, pos, cigar, NM = string.split(',')
                    pos = int(pos)
                    if(chrom == answerChr and (abs(pos - answerpos0) < Distance or abs(pos- answerpos1) < Distance)):
                        flag_map_mate0 = 1
                        break
                aln_list1 = samTrack1[0].parseXA()
                for string in aln_list1:
                    chrom, pos, cigar, NM = string.split(',')
                    pos = int(pos)
                    if(chrom == answerChr and (abs(pos - answerpos0) < Distance or abs(pos- answerpos1) < Distance)):
                        flag_map_mate1 = 1
                        break
                print >>sys.stderr, '>posDiff   1T_2F %s_%s'%('01'[flag_map_mate0], '01'[flag_map_mate1])

            elif samTrack1[0].rname == samTrack1[1].rname == answerChr and \
                (abs(min(samTrack1[0].pos, samTrack1[1].pos) - answerpos0) < Distance and
                  abs(max(samTrack1[0].pos, samTrack1[1].pos) - answerpos1) < Distance):
                aln_list0 = samTrack0[0].parseXA()
                for string in aln_list0:
                    chrom, pos, cigar, NM = string.split(',')
                    pos = int(pos)
                    if(chrom == answerChr and (abs(pos - answerpos0) < Distance or abs(pos- answerpos1) < Distance)):
                        flag_map_mate0 = 1
                        break
                aln_list1 = samTrack0[1].parseXA()
                for string in aln_list1:
                    chrom, pos, cigar, NM = string.split(',')
                    pos = int(pos)
                    if(chrom == answerChr and (abs(pos - answerpos0) < Distance or abs(pos- answerpos1) < Distance)):
                        flag_map_mate1 = 1
                        break
                print >>sys.stderr, '>posDiff   1F_2T %s_%s'%('01'[flag_map_mate0], '01'[flag_map_mate1])
            else:
                flag_map_mate0 = 0
                flag_map_mate1 = 0
                aln_list0 = samTrack0[0].parseXA()
                for string in aln_list0:
                    chrom, pos, cigar, NM = string.split(',')
                    pos = int(pos)
                    if(chrom == answerChr and (abs(pos - answerpos0) < Distance or abs(pos- answerpos1) < Distance)):
                        flag_map_mate0 = 1
                        break
                aln_list1 = samTrack0[1].parseXA()
                for string in aln_list1:
                    chrom, pos, cigar, NM = string.split(',')
                    pos = int(pos)
                    if(chrom == answerChr and (abs(pos - answerpos0) < Distance or abs(pos- answerpos1) < Distance)):
                        flag_map_mate1 = 1
                        break
                str0 = '01'[flag_map_mate0]+'_'+'01'[flag_map_mate1]
                flag_map_mate0 = 0
                flag_map_mate1 = 0
                aln_list0 = samTrack1[0].parseXA()
                for string in aln_list0:
                    chrom, pos, cigar, NM = string.split(',')
                    pos = int(pos)
                    if(chrom == answerChr and (abs(pos - answerpos0) < Distance or abs(pos- answerpos1) < Distance)):
                        flag_map_mate0 = 1
                        break
                aln_list1 = samTrack1[1].parseXA()
                for string in aln_list1:
                    chrom, pos, cigar, NM = string.split(',')
                    pos = int(pos)
                    if(chrom == answerChr and (abs(pos - answerpos0) < Distance or abs(pos- answerpos1) < Distance)):
                        flag_map_mate1 = 1
                        break
                str1 = '01'[flag_map_mate0]+'_'+'01'[flag_map_mate1]

                print >>sys.stderr, '>posDiff   1F_2F %s %s'%(str0, str1)
            
            if opt.ref != '':
                d = [0,0,0,0]
                for i, track in enumerate((samTrack0[0], samTrack0[1], samTrack1[0], samTrack1[1])):
                   
                    if track.flag & 4 != 0:
                        d[i] = -1
                        continue
                    cigar_list = parseCigar(track.cigar)
                    l_ref, read_start, read_end = 0, 0, len(track.seq) 
                    for c, j in cigar_list:
                        if c == 'M' or c == 'D': l_ref += j
                    if cigar_list[0][0] == 'S': read_start = cigar_list[0][1]
                    if cigar_list[-1][0] == 'S': read_end -= cigar_list[-1][1] 
                    ref = dictRef[track.rname][track.pos-1:track.pos-1+l_ref].upper()
                    read = track.seq[read_start:read_end].upper()
                    if opt.print_seq:
                        print>>sys.stderr, ref
                        print>>sys.stderr, read
                    d[i] = lev1(ref, read)
                print>>sys.stderr, '>edit distance %d %d %d %d'%(d[0], d[1], d[2], d[3])
                print>>sys.stderr, '>diff distance %d %d'%(d[0]-d[2], d[1]- d[3])

            print >>sys.stderr, leftPair[0]
            print >>sys.stderr, leftPair[1]
            print >>sys.stderr, rightPair[0]
            print >>sys.stderr, rightPair[1]
    fp1.close()
    fp2.close()
    if opt.ans_file != '':
        fp_ans.close()
import optparse

if __name__ == '__main__':
    #options init
    usage = "usage: %prog [Options] <file1> <file2>"
    parser = optparse.OptionParser(usage)
    parser.add_option('-d', '--diff', action = 'store', type = 'int',  dest='max_diff', help = 'max diff between different SAM', default= 4)
    parser.add_option('-l', '--length', action = 'store', type = 'int',  dest='read_length', help = 'read length', default= 100)
    parser.add_option('-a', '--answer', action = 'store', type = 'string',  dest='ans_file', help = 'answer file name', default= '')
    parser.add_option('-r', '--reference', action = 'store',  type = 'string', dest='ref', help = 'if reference file name is given, the edit distance will be computed', default= '')
    parser.add_option('-p', '--print', action = 'store_true',  dest='print_seq', help = 'print reference and read', default= False)

    #get options
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        exit(1)
    diffSam_wgsim_main(options, args)







