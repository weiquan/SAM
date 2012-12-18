import sys;
import pysam;
usage = """ statSam <prefix.bam> 
            
            Note: SAM file should be converted to Bam file and BAI file by SAMtools first!
            
            samtools view -bS <prefix.sam> > <prefix.bam>
            samtools sort <prefix.bam> <prefix.sorted>
            samtools index <prefix.bam>
            """

def main(argv):
    if len(sys.argv) < 2:
            print usage;
    fn_sam = argv[1];
    fp_sam = pysam.Samfile(fn_sam);
    if fp_sam == False:
        print fn_sam + " file open fail!"
        exit(1)
    statMappedNum = fp_sam.mapped;
    statUnMappedNum = fp_sam.unmapped; 
    statPerfectHit = 0;
    statGoodHit = 0
    for read in fp_sam.fetch():
        realPos = [0, 0]
        realPos[0] = int(read.qname.split('_')[2]);#parse pos answer
        realPos[1] = int(read.qname.split('_')[3]);
        
        realPos[1] = realPos[1] - read.qlen
        pos = read.pos +1#0-based to 1-based
        if  pos in realPos:
            statPerfectHit = statPerfectHit +1;
        else:
            if pos in read.positions:
                statGoodHit = statGoodHit +1;
    print "Num of Mapped Reads: %d"%(statMappedNum);
    print "Num of Unmapped Reads: %d"%(statUnMappedNum);
    print "Num of Perfect Hits: %d"%(statPerfectHit);
    print "Num of Good Hits: %d"%(statGoodHit);
        
    fp_sam.close();

if __name__ == '__main__':
        main(sys.argv)
