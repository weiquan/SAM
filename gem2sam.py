import optparse
import re

class gem ():
    """docstring for gem """
    FLAG_PAIRED = 0x0001
    FLAG_PROPERPAIRED = 0x0002
    FLAG_UNMAPPED0 = 0x0004
    FLAG_UNMAPPED1= 0x0008
    FLAG_STRANDBCKWARD0=0x0010
    FLAG_STRANDBCKWARD1=0x0020
    FLAG_READ0=0x0040
    FLAG_READ1=0x0080
    FLAG_NOTPRIMARY=0x0010
    class hit():
        def __init__(self, string, mate=None):
            self.chrom, self.strand, self.pos, self.cigar = string.split(':')
            self.pos = int(self.pos)
            self.mate = mate
    def __init__(self, is_pe, track):
        words = track.split('\t')
        self.qname = words[0]
        self.pe = is_pe
        self.properPaired = True #PE
        if re.find(r'[/|][12]',self.qname) > 0:#SE
            self.paired = False
            self.qname = self.qname[:-2] #trim '/1'
        if is_pe and self.properPaired:
            self.seq0, self.seq1 = words[1].split()
            self.qual0, self.qual1 = words[2].split()
            self.map = words[3]
            self.hit0 = []
            self.hit1 = []
            for string in words[4].split(','):
                if string.find(':::'):
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
                if string.find(':::'):
                    string, ann = string.split(':::')
                aln0 = string
                hit0 = self.hit(aln0)
                self.hit0.append(hit0)
    def printHitSAM(self):
        flag = 0
        if self.is_pe:
            flag &=gem.FLAG_PAIRED
            
        if self.properPaired:
            flag &=gem.FLAG_PROPERPAIRED #to be edited
        else:
            if self.strand == '-':
                flag &= gem.FLAG_STRANDBCKWARD0
            if self.is_pe and re.find(r'[/|][1]', self.qname) >0:
                flag &= gem.FLAG_READ0
            if self.is_pe and re.find(r'[/|][2]', self.qname) >0:
                flag &= gem.FLAG_READ1
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
    fn = arg[0]
    fp = open(fn ,'r')
    for track in fp:
        track = track.strip()
        
    fp.close()

   