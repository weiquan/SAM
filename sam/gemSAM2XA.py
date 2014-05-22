import sys
import itertools

class samTrack:
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
        if len(buffer) == 0:
            buffer.append(sam)
            last_sam = sam
            continue
        if sam.qname == last_sam.qname:
            buffer.append(sam)
            last_sam = sam
        elif len(buffer) >0:
            sam0, sam1 = None, None
            i = 0
            while i < len(buffer):
                if sam0 == None:
                    sam0 = buffer[i]
                    i += 1

                elif buffer[i].seq == '*':
                    if i < len(buffer) and buffer[i+1].seq == '*':
                        sam0.mergeHit(buffer[i])
                        if sam1 == None:
                            print >>sys.stderr, '[ERROR]: sam1 should not be None!'
                            exit(1)
                        sam1.mergeHit(buffer[i+1])
                    else:
                        print >>sys.stderr, '[ERROR]: sam0 = * sam1 != *'
                        exit(1)
                    i += 2

                elif sam0.isSameRead(buffer[i]):
                    sam0.mergeHit(buffer[i])
                    i += 1
                elif sam1 == None:
                    sam1 = buffer[i]
                    i += 1
                elif sam1.isSameRead(buffer[i]):
                    sam1.mergeHit(buffer[i])
                    i += 1
                else:
                    print >>sys.stderr, '[ERROR]: Neither SE nor PE!'
                    sys.exit(1)
            
            if sam0 != None:
                sam0.printSAM()
            if sam1 != None:
                sam1.printSAM()

  

            buffer = [sam]
            last_sam = sam
    if len(buffer) != 0:
        
        sam0, sam1 = None, None
        i = 0
        while i < len(buffer):
            if sam0 == None:
                sam0 = buffer[i]
                i += 1

            elif buffer[i].seq == '*':
                if i < len(buffer) and buffer[i+1].seq == '*':
                    sam0.mergeHit(buffer[i])
                    if sam1 == None:
                        print >>sys.stderr, '[ERROR]: sam1 should not be None!'
                        exit(1)
                    sam1.mergeHit(buffer[i+1])
                else:
                    print >>sys.stderr, '[ERROR]: sam0 = * sam1 != *'
                    exit(1)
                i += 2

            elif sam0.isSameRead(buffer[i]):
                sam0.mergeHit(buffer[i])
                i += 1
            elif sam1 == None:
                sam1 = buffer[i]
                i += 1
            elif sam1.isSameRead(buffer[i]):
                sam1.mergeHit(buffer[i])
                i += 1
            else:
                print >>sys.stderr, '[ERROR]: Neither SE nor PE!'
                sys.exit(1)

        if sam0 != None:
            sam0.printSAM()
        if sam1 != None:
            sam1.printSAM()
    fp.close()