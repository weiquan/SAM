#*
# * =====================================================================================
# *
# *       Filename:  RepeatSam.py
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
import sys
import os
usage = '''*****Test repeat alignment in samfile****	
		Usage:	RepeatSam.py <file1.sam> '''

#from samtools
#static inline int strnum_cmp(const char *a, const char *b)
#{
#	char *pa, *pb;
#	pa = (char*)a; pb = (char*)b;
#	while (*pa && *pb) {
#		if (isdigit(*pa) && isdigit(*pb)) {
#			long ai, bi;
#			ai = strtol(pa, &pa, 10);
#			bi = strtol(pb, &pb, 10);
#			if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
#		} else {
#			if (*pa != *pb) break;
#			++pa; ++pb;
#		}
#	}
#	if (*pa == *pb)
#		return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;
#	return *pa<*pb? -1 : *pa>*pb? 1 : 0;
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
			#	print >>sys.stderr, '%s\t%d\t%d'%(seqName1, start1, start2)
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


			

	
		
def testRepeat_main(argv):
	if len(argv) < 2:
		print usage
		exit(1)
	
	prefix1 = sys.argv[1][:-4]
	Bam1 = prefix1 +".bam"
	sortedPrefix1 = prefix1 + '.sorted'
	sortedBam1 = sortedPrefix1 + '.bam'
	#convert sam to bam
	cmd1 = "samtools view -bS "+sys.argv[1]+'>'+Bam1
	print>>sys.stderr, cmd1
	os.system(cmd1)
	#sort bam by name
	cmd1 = "samtools sort -n "+Bam1+' '+sortedPrefix1
	print>>sys.stderr, cmd1
	os.system(cmd1)
	#convert sorted bam to sam
	cmd1 = "samtools view "+sortedBam1+'>'+prefix1+'.sorted.sam'
	print>>sys.stderr, cmd1
	os.system(cmd1)
	print>>sys.stderr, "test repeat in file "+ sys.argv[1]
	fp = open(prefix1+'.sorted.sam', 'r')
	last_line = ''
	line = fp.readline()
	while len(line):
		if comp_seq_name(line.split('\t')[0], last_line.split('\t')[0]) != 0:
			last_line = line
		else:
			print line
		line = fp.readline()
	fp.close()



if __name__ == '__main__':
	testRepeat_main(sys.argv)
