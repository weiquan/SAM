import sys
usage = 'diffSam4.py <file1> <file2>'
filename1 = sys.argv[1]
filename2 = sys.argv[2]
fp1 = open(filename1, 'r')
fp2 = open(filename2, 'r')
l1 = len(fp1.readlines())
l2 = len(fp2.readlines())
if l1 != l2 or l1%2 != 0:
	print >>sys.stderr, 'file1 line num != file line num'
	sys.exit(1)
leftPair = [0,0]
rightPair = [0, 0]
while i in range(l1/2):
	leftPair[0] = fp1.readline()[:-1]
	leftPair[1] = fp1.readline()[:-1]
	rightPair[0] = fp2.readline()[:-1]
	rightPair[1] = fp2.readline()[:-1]
	leftPos0 = leftPair[0].split('\t')[3]
	leftPos1 = leftPair[1].split('\t')[3]
	rightPos0 = rightPair[0].split('\t')[3]
	rightPos1 = rightPair[1].split('\t')[3]
	if (leftPos0 == rightPos0 and leftPos1 == rightPos1) or (leftPos0 == rightPos1 and leftPos1 == rightPos0):
		pass
	else:
		print >>sys.stderr, '>>diff!'
		print >>sys.stderr, leftPair0		
		print >>sys.stderr, leftPair1		
		print >>sys.stderr, rightPair0		
		print >>sys.stderr, rightPair1		
fp1.close()
fp2.close() 