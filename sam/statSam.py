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
import os;
import sys;
import pysam;
usage = """ statSam [opt] <prefix.sam> 
            
            Note: SAM file should be converted to Bam file and BAI file by SAMtools first!
            
            samtools view -bS <prefix.sam> > <prefix.bam>
            samtools sort <prefix.bam> <prefix.sorted>
            samtools index <prefix.sorted.bam>

            eg: samtools view -bS 123.sam > 123.bam
                samtools sort 123.bam 123.sorted
                samtools index 123.sorted.bam
                python statSam [opt] 123.sorted.bam
        """
def samtools_pipeline(fn_sam):
    cmd1 = 'samtools view -bS ' + fn_sam +' >tmp.bam'
    print cmd1
    os.system(cmd1)
    cmd2 = 'samtools sort tmp.bam tmp.sorted'
    print cmd2
    os.system(cmd2)
    cmd3 = 'samtools index tmp.sorted.bam'
    print cmd3
    os.system(cmd3)
def dump_tmp():
    os.system('rm -f tmp.bam')
    os.system('rm -f tmp.sorted.bam')
    os.system('rm -f tmp.sorted.bam.bai')

def parse_answer(str, realPos):
    realPos[0] = int(str.split('_')[1]);
    realPos[1] = int(str.split('_')[2]);
    return realPos
def main(argv):
    if len(sys.argv) < 2:
            print usage;
            exit(1)
    
    fn_sam = argv[1];
    samtools_pipeline(fn_sam)        
    print 60
    fp_sam = pysam.Samfile('tmp.sorted.bam');
    

    
    if fp_sam == False:
        print fn_sam + " file open fail!"
        exit(1)
    totReads = fp_sam.mapped + fp_sam.unmapped;
    statMappedNum = fp_sam.mapped;
    statUnMappedNum = fp_sam.unmapped; 
    statPerfectHit = 0;
    statBadHits = 0
    for read in fp_sam.fetch():
        #parse pos answer
        realPos = [0, 0]
        parse_answer(read.qname, realPos)
        #end pos to start pos
        if realPos[1] > realPos[0]:
            realPos[1] = realPos[1] - read.qlen+1
        else:
            realPos[0] = realPos[0] - read.qlen+1
        
        pos = read.pos +1#0-based leftmost pos to 1-based leftmost pos
        #print realPos[0], realPos[1], pos
        if  pos in realPos:
            statPerfectHit = statPerfectHit +1;
        else:
            statBadHits = statBadHits +1;

    print "Num of total Reads: %d"%(totReads);
    print "Num of Mapped Reads: %d, percent %f %%"%(statMappedNum, float(statMappedNum)/float(totReads) * 100);
    print "Num of Unmapped Reads: %d, percent %f %%"%(statUnMappedNum, float(statUnMappedNum)/float(totReads) * 100);
    print "Num of Perfect Hits: %d, percent %f %%"%(statPerfectHit, float(statPerfectHit)/float(totReads) * 100);
    print "Num of Bad Hits: %d, percent %f %%"%(statBadHits, float(statBadHits)/float(totReads) * 100);
        
    fp_sam.close();
    dump_tmp()
if __name__ == '__main__':
        main(sys.argv)
