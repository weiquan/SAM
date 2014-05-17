import optparse
import re
class gem (object):
    """docstring for gem """
    def __init__(self, track):
        words = track.split('\t')
        self.qname = words[0]
        self.paired = True #PE
        if re.find(r'[/|][12]',self.qname) > 0:#SE
            self.paired = False
            self.qname = self.qname[:-2] #trim '/1'
        if self.paired:
            self.seq0 = words[1]
            self.seq1 = words[2]
            self.qual0 = words[3]
            self.qual1 = words[4]
        else:
            self.seq0 = words[1]
            self.qual0 = words[2]
            self.seq1 = None
            self.qual1 = None

        
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

   