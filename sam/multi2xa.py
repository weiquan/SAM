import sys
import itertools
#convert multi to xa
FLAG_PAIRED = 0x0001
FLAG_READ0 = 0x0040
FLAG_READ1 = 0x0080
class samTrack:
    @staticmethod
    def isPaired(flag):
        return flag & FLAG_PAIRED !=0
    @staticmethod
    def isRead0(flag):
        return flag &FLAG_READ0 != 0
    @staticmethod
    def isRead1(flag):
        return flag &FLAG_READ1 != 0
    class hit:
        def __init__(self, chrom ='', pos=0, cigar='', nm=-1):
            self.chrom = chrom
            self.pos = pos
            self.cigar = cigar
            self.nm = nm
    def __init__(self, line= None):
        if line == None:
            print >>sys.stderr, '[ERROR]: line == None'
            sys.exit(1)
            return
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
        self.xa = []
        if len(split_line) > 11:
            self.opt = split_line[11:]
        else:
            self.opt = []
    def isSameRead(self, samTrack):

        if self.seq == samTrack.seq:
            return True
        nt2bt = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3}
        seq0 = self.seq
        seq1 = samTrack.seq[::-1]
        for i, j in itertools.izip(seq0, seq1):
            x = nt2bt[i]
            y = nt2bt[j]
            if x+y != 3:
                return False
        return True
    def mergeHit(self, samtrack):
        pos = ''
        if samtrack.flag & 0x10 != 0:
            pos = '-'+str(samtrack.pos)
        else:
            pos = '+' +str(samtrack.pos)
        hit = self.hit(samtrack.rname, pos, samtrack.cigar)
        self.xa.append(hit)
    def parseXA(self):
        aln_list = []
        for opt in self.opt:
            if opt[:2] == 'XA':
                aln_list = opt.split(':')[2].split(';')[:-1]
        return aln_list
    def printSAM(self):
        if len(self.xa) != 0:
            xa_string = 'XA:Z:'
            for hit in self.xa:
                xa_string += hit.chrom+','+str(hit.pos)+','+hit.cigar+','+str(hit.nm)+';'
            print '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}'.format(
                self.qname, self.flag, self.rname, self.pos, self.mapq, self.cigar, 
                self.mrname, self.mpos, self.isize, self.seq, self.qual, xa_string)
            return
        print '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}'.format(
            self.qname, self.flag, self.rname, self.pos, self.mapq, self.cigar, self.mrname, self.mpos, self.isize, self.seq, self.qual)
import optparse
if __name__ == '__main__':
    usage = "usage: %prog [Options] <file>"

    parser = optparse.OptionParser(usage)
    parser.add_option('-p', '--pair', action = 'store_true',   dest='paired', help = 'SE/PE', default= False)

    opts, arg = parser.parse_args()

    if len(arg) != 1:
        parser.print_help()
        exit(1)
    fn = arg[0]
    fp = open(fn, 'r')
    buffer = []
    last_sam = None
    sam0, sam1 = None, None
    for line in fp:
        line = line.strip()
        if line[0] == '@':#skip header
            print line
            continue
        sam = samTrack(line)
        if last_sam == None:
            if samTrack.isPaired(sam.flag):
                if samTrack.isRead0(sam.flag):
                    sam0 = sam
                elif samTrack.isRead1(sam.flag):
                    sam1 = sam
                else:
                    print >>sys.stderr, 'Read {0} is not paired'.format(sam.qname)
            else:
                sam0 = sam
        elif sam.qname == last_sam.qname:
            if sam.isPaired(sam.flag):
                if sam.isRead0(sam.flag):
                    if sam0 == None:
                        sam0 = sam
                    else: sam0.mergeHit(sam)
                elif sam.isRead1(sam.flag):
                    if sam1 == None:
                        sam1 = sam
                    else:sam1.mergeHit(sam)
                else:
                    print >>sys.stderr, 'Read {0} is not paired'.format(sam.qname)
            else:
                sam0.mergeHit(sam)
        else:
            if sam0 != None: sam0.printSAM()
            if sam1 != None: sam1.printSAM()
            if samTrack.isPaired(sam.flag):
                if samTrack.isRead0(sam.flag):
                    sam0 = sam
                    sam1 = None
                elif samTrack.isRead1(sam.flag):
                    sam0 = None
                    sam1 = sam
                else:
                    print >>sys.stderr, 'Read {0} is not paired'.format(sam.qname)
            else: 
                sam0 = sam
                sam1 = None
  
        last_sam = sam
    if sam0 != None: sam0.printSAM()
    if sam1 != None: sam1.printSAM()
    
    fp.close()
