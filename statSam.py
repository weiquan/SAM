#*
# * =====================================================================================
# *
# *       Filename:  statSam.py
# *
# *    Description:  
# *
# *        Version:  1.0
# *        Created:  12/22/2012 03:41:52 PM
# *       Revision:  none
# *   
# *
# *         Author:  Wei Quan (mn), wquanhit@gmail.com
# *        Company:  BIC, HIT
# *
# * =====================================================================================
# */
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
    totReads = fp_sam.mapped + fp_sam.unmapped;
    statMappedNum = fp_sam.mapped;
    statUnMappedNum = fp_sam.unmapped; 
    statPerfectHit = 0;
    statGoodHit = 0
    for read in fp_sam.fetch():
        #parse pos answer
        realPos = [0, 0]
        realPos[0] = int(read.qname.split('_')[2]);
        realPos[1] = int(read.qname.split('_')[3]);
        #end pos to start pos
        if realPos[1] > realPos[0]:
            realPos[1] = realPos[1] - read.qlen
        else:
            realPos[0] = realPos[0] - read.qlen
        
        pos = read.pos +1#0-based leftmost pos to 1-based leftmost pos
        if  pos in realPos:
            statPerfectHit = statPerfectHit +1;
        else:
            if realPos[0] in read.positions or realPos[1] in read.positions:
                statGoodHit = statGoodHit +1;

    print "Num of total Reads: %d"%(totReads);
    print "Num of Mapped Reads: %d, percent %f %%"%(statMappedNum, float(statMappedNum)/float(totReads) * 100);
    print "Num of Unmapped Reads: %d, percent %f %%"%(statUnMappedNum, float(statUnMappedNum)/float(totReads) * 100);
    print "Num of Perfect Hits: %d, percent %f %%"%(statPerfectHit, float(statPerfectHit)/float(totReads) * 100);
    print "Num of Good Hits: %d, percent %f %%"%(statGoodHit, float(statGoodHit)/float(totReads) * 100);
        
    fp_sam.close();

if __name__ == '__main__':
        main(sys.argv)
