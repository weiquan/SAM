#!/usr/bin/python2.7
import sys
from multiprocessing import Process, Array


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
                if abs(pos-right) < d: return True
            else:
                pos = pos2leftcoord(abs(aln_pos), '*')
                if abs(pos-left) <d: return True
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
n_reads = 0
n_unmapped =1
n_mapped =2
n_hits =3 
n_primary =4
n_alternative =5

def alneval_wgsim_worker(n_processers, process_id, opt, buf, data):
    d = opt.max_shift
    print process_id, data
    for i, line in enumerate(buf):
        if i % n_processers !=process_id:
            continue
        line = line.strip()
        #skip header
        if line[0] == '@':
            continue
         
        #main 
        track = samTrack(line)
        if opt.debug:
           print track.qname, track.rname, track.pos    
        chrom, left, right = parseWgsimAnswer(track.qname)
        if not track.isMapped():#unmapped
            data[n_reads] += 1
            data[n_unmapped] += 1
            if opt.debug: print 'Unmapped'
            continue
        data[n_reads] +=1
        data[n_mapped] +=1
        data[n_hits] += track.n_hits()
        if track.mapPrimary(chrom, int(left), int(right), d):
            data[n_primary] += 1
            if opt.debug: print 'Primary'
            continue
        if track.mapAlternative(chrom, left, right, d):
            data[n_alternative] +=1
            if opt.debug: print 'Alternative'
            continue
        if opt.debug :print 'bad map'
    print process_id, data
def readBuf(fp, n_lines):
    buf = []
    for i in xrange(n_lines):
        line = fp.readline()
        if len(line) == 0:
            break
        buf.append(line)
    return buf
def alneval_wgsim_main(opt, arg):
    
    n_lines = opt.n_lines
    #INIT data
    data = []
    for i in xrange(opt.n_processer):
        data.append(Array('i', [0 for i in xrange(6)]))
    fp = open(arg[0] ,'r')
    buf = readBuf(fp, n_lines)    
    while len(buf) > 0: 
        processes = []

        for i in xrange(opt.n_processer):
            t = Process(target=alneval_wgsim_worker, args =(opt.n_processer, i, opt, buf, data[i]))
            processes.append(t)
        for i in xrange(opt.n_processer):
            processes[i].start()
        for i in xrange(opt.n_processer):
            processes[i].join()   
        buf = readBuf(fp, n_lines)    
    print data[0]
    print data[1]
    for i in xrange(1, opt.n_processer):
        for key in xrange(6):
            data[0][key] += data[i][key]
    fp.close()

    print >>sys.stderr, 'Reads Number:{0}'.format(data[0][n_reads])
    print >>sys.stderr, 'unmapped Number:{0}'.format(data[0][n_unmapped])
    print >>sys.stderr, 'Mapped Number:{0}'.format(data[0][n_mapped])
    print >>sys.stderr, 'Primary mapped Number:{0}'.format(data[0][n_primary])
    print >>sys.stderr, 'Alternative mapped Number:{0}'.format(data[0][n_alternative])
    print >>sys.stderr, 'Primary + alternative :{0}'.format(data[0][n_primary]+data[0][n_alternative])
    print >>sys.stderr, 'Hits Number:{0}'.format(data[0][n_hits])




import optparse

if __name__ == '__main__':
    #options init
    usage = "usage: %prog [Options] <file>"
    parser = optparse.OptionParser(usage)
    parser.add_option('-s', '--shift', action = 'store', type = 'int',  dest='max_shift', help = 'max shift between answer and alignment pos', default= 4)
    parser.add_option('--debug', action = 'store_true',  dest='debug', help = 'print seq with alignment tag(not support multi processs)', default= False)
    parser.add_option('--processers', action = 'store',  type = 'int', dest='n_processer', help ='processers num', default= 1)
    parser.add_option('-n', action = 'store',  type = 'int', dest='n_lines', help ='buf lines', default= 1000000)

    #get options
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        exit(1)
    alneval_wgsim_main(options, args)







