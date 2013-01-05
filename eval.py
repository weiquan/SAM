# *
# * =====================================================================================
# *
# *       Filename:  eval.py
# *
# *    Description:  eval aln result  
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
import re
import sys
import getopt
usage = '''
Usage:           eval.py <cmd> [opt]
Cmd:                      aln        evaluation alignment result in sam file'''
usage_alnEval = '''
Usage:       eval.py aln [opt] <samfile>
opt:                                     '''
def parseWgsimAnswer(seqName):
    pattern = re.compile(r'_\d+_\d+_')
    match = pattern.search(seqName)
    if match == None:
        print >>sys.err, 'seq ' + seqName + "is not generated by wgsim!"
        exit(1);
    chrome = seqName[:match.start()]
    posLeft, posRight = seqName[match.start()+1:match.end()-1].split('_')
    posLeft = long(posLeft)
    posRight = long(posRight)
    return chrome, posLeft, posRight
def alnEval(arg):
    FLAG_UNMAP = 0x4
    FLAG_REVERSE = 0x10

    unmappedReadsNum = 0
    mappedReadsNum = 0
    totReadsNum = 0
    perfectAlnNum = 0

    OutPutUnMappedReads = False
    fp_unmapped = 0
    #parse opt
    try:
        opts, args = getopt.getopt(arg[1:], "hu:",['help'])
    except getopt.GetoptError as err:
        # print help information and exit:
        print >>sys.stderr, str(err) # will print something like "option -a not recognized"
        print >>sys.stderr,usage_alnEval
        sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            print >>sys.stderr,usage_alnEval
            sys.exit()
        elif o in('-u'):
            OutPutUnMappedReads = True
            fp_unmapped = open(a, 'w')
        else:
            assert False, "unhandled option"
    # ...
    if len(args) != 1:
        print >>sys.stderr,usage_alnEval
        exit(1)

    samFileName = args[0]

    print >>sys.stderr, 'Begin evaluation...'
    fp = open(samFileName, 'r')
    line = fp.readline()
    while len(line):
        #skip the head 
        if line[0] == '@':
            line = fp.readline()
            continue
        #
        samFormat = line.split('\t')
        #print samFormat

        seqName, flag, ref, pos, Mapq, cigar, mateRef, matePos, templateLen, seq, qual = samFormat[:11]
        flag = int(flag)
        pos = long(pos)
        #templateLen = int(templateLen)
        templateLen = 100

        answerChr, answerLeftPos, answerRightPos = parseWgsimAnswer(seqName)
        if flag & FLAG_UNMAP != 0:
            unmappedReadsNum = unmappedReadsNum +1
            if(OutPutUnMappedReads):
                print >>fp_unmapped, line
        else:
            mappedReadsNum = mappedReadsNum +1
            if ref == answerChr:
                if flag & FLAG_REVERSE == 0:
                    if pos == answerLeftPos:
                        perfectAlnNum = perfectAlnNum +1
                else:
                    if pos == answerRightPos - templateLen+1:
                        perfectAlnNum = perfectAlnNum +1
        totReadsNum = totReadsNum +1
        if totReadsNum % 100000 == 0:
            print >>sys.stderr, 'eval %u reads...'%(totReadsNum)

        line = fp.readline()
    if(OutPutUnMappedReads):
        fp_unmapped.close()
    fp.close()
    print '**********************************************'
    print 'Evaluation of the alignment result!'
    print 'Unmapped reads Num : %u'%(unmappedReadsNum)
    print 'Mapped reads Num : %u'%(mappedReadsNum)
    print 'Total reads Num : %u'%(totReadsNum)
    print 'Prefect alignment Num : %u'%(perfectAlnNum)
    print '**********************************************'

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print >>sys.stderr, usage;
        exit(1)
    cmd = sys.argv[1]
    if(cmd == 'aln'):
        alnEval(sys.argv[1:])
    else:
        print >>sys.stderr, 'no cmd called '+cmd
        print >>sys.stderr, usage;



