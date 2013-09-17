#*
# * =====================================================================================
# *
# *       Filename:  diffSam.py
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
#!/usr/bin/python2.7
import sys
import os
usage = """ diffSam2 <file1.sam> <file2.sam>
        print alned reads which are alned better in file1 than in file2 """
#from samtools
#static inline int strnum_cmp(const char *a, const char *b)
#{
#   char *pa, *pb;
#   pa = (char*)a; pb = (char*)b;
#   while (*pa && *pb) {
#       if (isdigit(*pa) && isdigit(*pb)) {
#           long ai, bi;
#           ai = strtol(pa, &pa, 10);
#           bi = strtol(pb, &pb, 10);
#           if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
#       } else {
#           if (*pa != *pb) break;
#           ++pa; ++pb;
#       }
#   }
#   if (*pa == *pb)
#       return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;
#   return *pa<*pb? -1 : *pa>*pb? 1 : 0;
#}


def comp_seq_name(seqName1, seqName2):
    start1 = 0
    start2 = 0
    end1 = 0
    end2 = 0
    while start1 < len(seqName1) and start2 < len(seqName2):
        if str.isdigit(seqName1[start1]) and str.isdigit(seqName2[start2]):
            end1 = start1 +1
            end2 = start2 +1
            while end1 < len(seqName1) and str.isdigit(seqName1[end1]):
                end1 = end1 +1
            
            while end2 < len(seqName2) and str.isdigit(seqName2[end2]):
                end2 = end2 +1
            
            #if(start1 == end1):
            #   print >>sys.stderr, '%s\t%d\t%d'%(seqName1, start1, start2)
            #if(start1 >= len(seqName1)):
            #print >>sys.stderr, '>%s\t%d\t%d'%(seqName1, start1, end1)
            num1 = long(seqName1[start1:end1])
            start1 = end1
            num2 = long(seqName2[start2:end2])
            start2 = end2
            
            if num1 < num2:
                return -1
            elif num1 > num2:
                return 1
        else:
            if seqName1[start1] < seqName2[start2]:
                return -1
            elif seqName1[start1] > seqName2[start2]:
                return 1
            else:
                start1 = start1+1
                start2 = start2+1
    if start1 == len(seqName1) or start1 == len(seqName2):
        if len(seqName1) < len(seqName2):
            return -1
        elif len(seqName1) > len(seqName2):
            return 1
        else:
            return 0


            

SAM_FLAG_PAIRED  = 0x0001
SAM_FLAG_PROPER_PAIRED  = 0x0002
SAM_FLAG_UNMAPED  = 0x0004
SAM_FLAG_MATE_UNMAPPED = 0x0008
SAM_FLAG_STRAND =0x0010
SAM_FLAG_MATE_STRAND =0x0020
SAM_FLAG_MATE_READ1 =0x0040
SAM_FLAG_MATE_READ2 =0x0080
SAM_FLAG_SECANDARY_ALN =0x0100
SAM_FLAG_QC =0x0200
SAM_FLAG_PCR_DUP =0x0400

        
def diffSam_main(opt, args):
    #print opt
    prefix1 = args[0][:-4]
    prefix2 = args[1][:-4]
    Bam1 = prefix1 +".bam"
    Bam2 = prefix2 + '.bam'
    sortedPrefix1 = prefix1 + '.sorted'
    sortedPrefix2 = prefix2 + '.sorted'
    sortedBam1 = sortedPrefix1 + '.bam'
    sortedBam2 = sortedPrefix2 + '.bam'
    #convert sam to bam
    cmd1 = "samtools view -bS "+sys.argv[1]+'>'+Bam1
    cmd2 = "samtools view -bS "+sys.argv[2]+'>'+Bam2
    print>>sys.stderr, cmd1
    os.system(cmd1)
    print>>sys.stderr, cmd2
    os.system(cmd2)
    #sort bam by name
    cmd1 = "samtools sort -n "+Bam1+' '+sortedPrefix1
    cmd2 = "samtools sort -n "+Bam2+' '+sortedPrefix2
    print>>sys.stderr, cmd1
    os.system(cmd1)
    print>>sys.stderr, cmd2
    os.system(cmd2)
    #convert sorted bam to sam
    cmd1 = "samtools view "+sortedBam1+'>'+prefix1+'.sorted.sam'
    cmd2 = "samtools view "+sortedBam2+'>'+prefix2+'.sorted.sam'
    print>>sys.stderr, cmd1
    os.system(cmd1)
    print>>sys.stderr, cmd2
    os.system(cmd2)
    fn1 = prefix1+'.sorted.sam'
    fn2 = prefix2+'.sorted.sam'
    fp1 = open(fn1, 'r')
    fp2 = open(fn2, 'r')

    line1 = fp1.readline()
    line2 = fp2.readline()

    while len(line1) and len(line2):
        split_line1 = line1.split('\t')
        split_line2 = line2.split('\t')

        if comp_seq_name(split_line1[0], split_line2[0]) < 0:#seq in line1 not in line2
            #print>>sys.stderr, line1
            #print>>sys.stderr, '++++++++'
            #print>>sys.stderr, line2
            if opt.more or opt.map or opt.correct:
                print '>track in %s not in %s'%(fn1, fn2)
            print line1,
            line1 = fp1.readline()
        elif comp_seq_name(split_line1[0], split_line2[0]) > 0:#seq in line2 not in line1
            line2 = fp2.readline()
        elif comp_seq_name(split_line1[0], split_line2[0]) == 0:#seq both in line1 and line2
            if opt.map and\
                    int(split_line1[1])&SAM_FLAG_PROPER_PAIRED != 0 and\
                    int(split_line2[1])&SAM_FLAG_PROPER_PAIRED == 0:#alned in line1 not alned in line2
                print '>proper paired in %s not in %s'%(fn1, fn2)
                print line1,
            if opt.correct:
                pass
            line1 = fp1.readline()
            line2 = fp2.readline()
        else:#some unknown problem
            print>>sys.stderr, 'erro!!!'
    while len(line1):
        print line1,
        line1 = fp1.readline()

    fp1.close()
    fp2.close() 

import optparse

if __name__ == '__main__':
    #options Init
    usage = "usage: %prog [Options] <file1> <file2>"
    parser = optparse.OptionParser(usage)
    #parser.add_option('-1', dest='filename1', help = '*SAM file 1')
    #parser.add_option('-2', dest='filename2', help = '*SAM file 2')
    parser.add_option('--more', action = 'store_true', dest='more', help = 'print track in SAM file1 not in SAM file2', default= False)
    parser.add_option('--map', action = 'store_true', dest='map', help = 'print track mapped in SAM file1 not mapped in SAM file2', default= False)
    parser.add_option('--correct', action = 'store_true', dest='correct', help = 'print track with a right pos in SAM file1 not with it in SAM file2[not supported]', default= False) 
    #get opt
    options, args = parser.parse_args()
    #print options['more']
    #print options['map']
    #print options['correct']

    if len(args) != 2:
        parser.print_help()
        exit(1)
    diffSam_main(options, args)
