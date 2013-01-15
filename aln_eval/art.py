import sys
import getopt
import const

def usage():
    usage_wgsim_alnEval = '''*************************************************************************
Usage:          aln_eval.py     art     [opt]   <samfile.sam>   <answer.sam>

opt:            -h/--help       help 
                -u  str         print unmapped reads to file with name str
                -m  str         print mapped reads to file with name str
                -a              if both unmapped in alignment and answer, the alignment is considered as perfectAln
                                eg: alignment answer    default -a mode
                                    unmap     map       unmap    unmap
                                    unmap     unmap     unmap    map
*************************************************************************'''
    print >>sys.stderr, usage_wgsim_alnEval
def check_sorted_by_name(fileName):
    fp = open(fileName, 'r')
    last_str = ''
    cur_str = ''
    line = fp.readline()
    while len(line):
        cur_str = line.split('\t')[0]
        if(comp_seq_name(cur_str, last_str) < 0):
            print >>sys.stderr, "name %s > name %s"%(cur_str, last_str)
            print >>sys.stderr, "file %s is not sorted by seq name"%(fileName)
            return False
        last_str = cur_str
        line = fp.readline()
    fp.close()
    return True


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

def art_alnEval(arg):

    optDic = {'OutPutMapped':False, 'OutPutUnMapped':False, 'SameMatch':False}
    statDic = {'unmappedReadsNum':0, 'mappedReadsNum':0, 'perfectAlnNum':0}
    #parse opt
    try:
        opts, args = getopt.getopt(arg[1:], "huma",['help'])
    except getopt.GetoptError as err:
        print >>sys.stderr, str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in('-u'):
            optDic['OutPutUnMapped'] = True #output unmapped
            fpUnmapped = open(a, 'w')
        elif o in('-m'):
            optDic['OutPutMapped'] = True #output mapped
            fpMapped = open(a, 'w')
        elif o in('-a'):
            optDic['SameMatch'] = True
        else:
            assert False, "unhandled option"
    # ...
    if len(args) != 2:
        usage()
        exit(1)
    samFileName, alnFileName = args[0:2]
    if not ( check_sorted_by_name(samFileName) and check_sorted_by_name(alnFileName) ):
        exit(1)

    totReadsNum = 0
    print >>sys.stderr, 'Begin evaluation...'
    fp = open(samFileName, 'r')
    fpAnswer = open(alnFileName, 'r')
    #skip head
    line0 = fp.readline()
    while len(line0) and line0[0] == '@':
        line0 = fp.readline()
    line1 = fpAnswer.readline()
    while len(line1) and line1[0] == '@':
        line1 = fpAnswer.readline()
    #core loop
    while len(line0) and len(line1):
        seqName0, flag0, chr0, pos0 = line0.split('\t')[0:4]
        seqName1, flag1, chr1, pos1 = line1.split('\t')[0:4]
        flag0, flag1, pos0, pos1 = (int(flag0), int(flag1), long(pos0), long(pos1))
        
        if flag0 & const.FLAG_UNMAP and comp_seq_name(seqName0, seqName1) != 0: #un map
            optDic['unmappedReadsNum'] +=1
            totReadsNum +=1
            if optDic['OutPutUnMapped']:
                print >>fpUnmapped, line0
            line0 = fp.readline()
            continue
        elif flag0 & const.FLAG_UNMAP and comp_seq_name(seqName0, seqName1) == 0:
            if optDic['SameMatch'] and flag1 & const.FLAG_UNMAP:
                statDic['mappedReadsNum'] += 1
                statDic['perfectAlnNum'] += 1
                if optDic['OutPutMapped']:
                    print >>fpMapped, line0
            else:    
                statDic['unmappedReadsNum'] += 1
                if optDic['OutPutUnMapped']:
                    print >>fpUnmapped, line0
            totReadsNum += 1
           
            line0 = fp.readline()
            line1 = fpAnswer.readline()
            continue
        else:# map
            statDic['mappedReadsNum'] += 1
            totReadsNum += 1 

        if comp_seq_name(seqName0, seqName1) < 0:
            print >>sys.stderr, "seq %s is not in answer file"%(seqName0)
            exit(1)
        elif comp_seq_name(seqName0, seqName1) > 0:
            statDic['unmappedReadsNum'] += 1
            if optDic['OutPutUnMapped']:
                print >>fpUnmapped, line1
            line1 = fpAnswer.readline()

        else:    
            if (flag0 & const.FLAG_REVERSE) == (flag1 & const.FLAG_REVERSE) and chr0 == chr1 and pos0 == pos1:
                statDic['perfectAlnNum'] += 1
                if optDic['OutPutMapped']:
                    print >>fpMapped, line0
            line0 = fp.readline()#may have bug when a seq have multi alignment track
            line1 = fpAnswer.readline()
    
    if len(line0) > 0:
        print >>sys.stderr, "seq %s is not in answer file"%(seqName0)
        exit(1)
    while len(line1):
        statDic['unmappedReadsNum'] += 1
        totReadsNum += 1
        if optDic['OutPutUnMapped']:
            print >>fpUnmapped, line1
        line1 = fpAnswer.readline()
    
    if optDic['OutPutUnMapped']:
        fpUnmapped.close()
    if optDic['OutPutMapped']:
        fpMapped.close()
    fp.close()
    fpAnswer.close()
    
    
    print '**********************************************'
    print 'Evaluation of the alignment result!'
    print 'Total reads Num : %u'%(totReadsNum)
    print "Unmapped reads Num : %u, percent %f"%(statDic['unmappedReadsNum'], (statDic['unmappedReadsNum'] *1.0)/totReadsNum)
    print "Mapped reads Num : %u, percent %f"%(statDic['mappedReadsNum'], (statDic['mappedReadsNum'] * 1.0)/totReadsNum)
    print "Prefect alignment Num : %u, percent %f"%(statDic['perfectAlnNum'], (statDic['perfectAlnNum'] * 1.0)/totReadsNum)
    print '**********************************************'