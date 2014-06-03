import sys



class samTrack:
    FLAG_PAIRED = 0x0001
    FLAG_UNMAPPED = 0x0004
    FLAG_READ0 = 0x0040
    FLAG_READ1 = 0x0080
    FLAG_REVERSE = 0x0010
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
    def isMapped(self):
        return self.flag & samTrack.FLAG_UNMAPPED ==0
    def isReverse(self):
        return self.flag & samTrack.FLAG_REVERSE != 0
    def n_hits(self):
        if not self.isMapped():
            return 0
        return len(self.parseXA())+1
    def mapPrimary(self, chrom, left, right, d):
        if self.rname != chrom:
            return False
        if self.isReverse():
            pos = pos2rightcoord(self.pos, self.cigar)
            return abs(pos-right) < d 
        else:
            pos = pos2leftcoord(self.pos, self.cigar)
            return abs(pos-left) <d
    def mapAlternative(self, chrom, left, right, d):
        alns = self.parseXA()
        if len(alns) == 0:
            return False
        for aln in alns:
            aln_chrom, aln_pos, aln_cigar, aln_nm = aln.split(',')
            aln_pos = int(aln_pos)
            if aln_chrom != chrom:
                continue
            if aln_pos <0:#reverse
                pos = pos2rightcoord(abs(aln_pos), '*', len(self.seq))
                return abs(pos-right) < d 
            else:
                pos = pos2leftcoord(abs(aln_pos), '*')
                return abs(pos-left) <d
        return False


def pos2rightcoord(pos, cigar =None, readLen = None):
    if (cigar is None or cigar == '*') and readLen == None:
        return None
    if (cigar == '*' or cigar is None )and readLen is not None:
        return pos + readLen-1
    if cigar is not None:#use cigar to convert right coord to left coord
        cigarList =  __parseCigar(cigar)
        if cigarList[-1][1] in 'SH':
            pos += int(cigarList[-1][0])
            cigarList = cigarList[:-1]
        for n, c in cigarList:
            if c in 'MDN':
                pos += int(n)
        pos -= 1
        return pos
    
    else:
        print >>sys.stderr, '[ERROR]: prog should not be here!'
        return None
def pos2leftcoord(pos, cigar ='*'):
    if cigar == '*':
        return pos
    cigarList =  __parseCigar(cigar)
    if cigarList[0][1] in 'SH':
        pos -= int(cigarList[0][0])
    return pos
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
def __parseCigar(cigar):
    import re
    p = re.compile(r'(\d+)([MIDNSHPX=])')
    return p.findall(cigar)
def __parseCigarIter(cigar):
    import re
    p = re.compile(r'(\d+)([MIDNSHPX=])')
    return p.finditer(cigar)    
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


def alneval_wgsim_main(opt, arg):
    
    d = opt.max_shift
    n_reads = 0
    n_mapped = 0
    n_unmapped = 0
    n_primary = 0
    n_alternative = 0
    n_hits = 0


    fp = open(arg[0], 'r')
    for line in fp:
        line = line.strip()
        #skip header
        if line[0] == '@':
            continue
        #main 
        track = samTrack(line)
        chrom, left, right = parseWgsimAnswer(track.qname)
        if not track.isMapped():#unmapped
            n_reads += 1
            n_unmapped += 1
            if opt.print_seq:
                if track.qual != '*':
                    print '@'+track.qname
                    print track.seq
                    print '+'
                    print track.qual
                else:
                    print '>'+track.qname
                    print track.seq
            continue
        n_reads +=1
        n_mapped +=1
        n_hits += track.n_hits()
        if track.mapPrimary(chrom, int(left), int(right), d):
            n_primary += 1
            continue
        if track.mapAlternative(chrom, left, right, d):
            n_alternative +=1
            continue

    fp.close()    
    print >>sys.stderr, 'Reads Number:{0}'.format(n_reads)
    print >>sys.stderr, 'unmapped Number:{0}'.format(n_unmapped)
    print >>sys.stderr, 'Mapped Number:{0}'.format(n_mapped)
    print >>sys.stderr, 'Primary mapped Number:{0}'.format(n_primary)
    print >>sys.stderr, 'Alternative mapped Number:{0}'.format(n_alternative)
    print >>sys.stderr, 'Primary + alternative :{0}'.format(n_alternative+n_primary)
    print >>sys.stderr, 'Hits Number:{0}'.format(n_hits)




import optparse

if __name__ == '__main__':
    #options init
    usage = "usage: %prog [Options] <file>"
    parser = optparse.OptionParser(usage)
    parser.add_option('-s', '--shift', action = 'store', type = 'int',  dest='max_shift', help = 'max shift between answer and alignment pos', default= 4)
    parser.add_option('-p', '--print', action = 'store_true',  dest='print_seq', help = 'print unmapped reads', default= False)

    #get options
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        exit(1)
    alneval_wgsim_main(options, args)







