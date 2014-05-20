import optparse
import re
import itertools

FLAG_PAIRED = 0x0001
FLAG_PROPERPAIRED = 0x0002
FLAG_UNMAPPED0 = 0x0004
FLAG_UNMAPPED1= 0x0008
FLAG_STRANDBCKWARD0=0x0010
FLAG_STRANDBCKWARD1=0x0020
FLAG_READ0=0x0040
FLAG_READ1=0x0080
FLAG_NOTPRIMARY=0x0010


class gem ():
    """docstring for gem """
   
    class hit():
        def __init__(self, string, mate=None):
            self.chrom, self.strand, self.pos, self.cigar = string.split(':')
            self.pos = int(self.pos)
            self.mate = mate
    def __init__(self, is_pe, track):
        words = track.split('\t')
        self.qname = words[0]
        self.is_pe = is_pe
        self.properPaired = True #PE
        if re.match(r'(\S+)[/|][12]',self.qname) != None:#SE
            self.properPaired = False
            self.qname = self.qname[:-2] #trim '/1'
        if is_pe and self.properPaired:
            self.seq0, self.seq1 = words[1].split()
            self.qual0, self.qual1 = words[2].split()
            self.map = words[3]
            self.hit0 = []
            self.hit1 = []
            for string in words[4].split(','):
                if string.find(':::') > 0:
                    string, ann = string.split(':::')
                aln0, aln1 = string.split('::')
                hit0 = self.hit(aln0)
                hit1 = self.hit(aln1)
                hit0.mate = hit1
                hit1.mate = hit0
                self.hit0.append(hit0)
                self.hit1.append(hit1)
        else:
            self.seq0 = words[1]
            self.qual0 = words[2]
            self.seq1 = None
            self.qual1 = None
            self.map = words[3]
            self.hit0 = []
            self.hit1 = None
            for string in words[4].split(','):
                if string.find(':::') >0:
                    string, ann = string.split(':::')
                aln0 = string
                hit0 = self.hit(aln0)
                self.hit0.append(hit0)
    def printHitSAM(self):
        flag = 0
        if self.is_pe:
            flag &=FLAG_PAIRED
            
        if self.is_pe and self.properPaired:
            flag &=FLAG_PROPERPAIRED #to be edited
            flag0 = flag1 = flag
            flag0 &= FLAG_READ0
            flag1 &= FLAG_READ1
            isize0 = isize1 = abs(self.hit0[0].pos - self.hit1[0].pos)+len(self.seq1)

            if self.hit0[0].strand == '+' and self.hit1[0].strand == '-':
                isize1 = -isize1
                flag1 &= FLAG_STRANDBCKWARD1
            elif self.hit0[0].strand == '-' and self.hit1[0].strand == '+':
                isize0 = - isize0
                flag0 &= FLAG_STRANDBCKWARD0
            xa0 = xa1 = 'XA:Z:'
            for hit0, hit1 in itertools.izip(self.hit0[1:], self.hit1[1:]):
                xa0 += hit0.chrom+','+str(hit0.pos)+','+hit0.cigar+','+'-1;'
                xa1 += hit1.chrom+','+str(hit1.pos)+','+hit1.cigar+','+'-1;'

            print '{0}  {1} {2} {3} 255 {4} {5} {6} {7} {8} {9}  {10}'.format(
                    self.qname, flag0, self.hit0[0].chrom, self.hit0[0].pos, self.hit0[0].cigar, 
                    self.hit0[0].mate.chrom, self.hit0[0].mate.pos, isize0, self.seq0, self.qual0, xa0)
            print '{0}  {1} {2} {3} 255 {4} {5} {6} {7} {8} {9}  {10}'.format(
                    self.qname, flag1, self.hit1[0].chrom, self.hit1[0].pos, self.hit1[0].cigar, 
                    self.hit1[0].mate.chrom, self.hit1[0].mate.pos, isize1, self.seq1, self.qual1, xa1)
        else:
            if self.hit0[0].strand == '-':
                flag &= FLAG_STRANDBCKWARD0
            if self.is_pe and re.match(r'(\S+)[/|][12]', self.qname) != None:
                flag &= FLAG_READ0
            if self.is_pe and re.match(r'(\S+)[/|][12]', self.qname) != None:
                flag &= FLAG_READ1
            opt = 'XA:Z:'#to be edited
            for hit in self.hit0[1:]:
                opt += hit.chrom+','+str(hit.pos)+','+hit.cigar+','+'-1;'
            print '{0}  {1} {2} {3} 255 {4} * 0 0 {5} {6}  {7}'.format(
                    self.qname, flag, self.hit0[0].chrom, self.hit0[0].pos, self.hit0[0].cigar, self.seq0, self.qual0, opt)
if __name__ == '__main__':
    usage = "usage: %prog [Options] <file>"

    parser = optparse.OptionParser(usage)
    parser.add_option('-p', '--pair', action = 'store_true',   dest='paired', help = 'SE/PE', default= False)

    opts, arg = parser.parse_args()
    if len(arg) != 1:
        parser.print_help()
        exit(1)
    fn = arg[0]
    fp = open(fn ,'r')
    for track in fp:
        track = track.strip()
        g = gem(opts.paired, track)
        g.printHitSAM()
    fp.close()

   