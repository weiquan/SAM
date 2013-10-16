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
import sys
import wgsim
import art
usage ='usage: %prog [Opt] <file1> <file2>' 
UHIT = 0 #un hit
PHIT = 1 #perfect hit
GHIT = 2 #good hit
BHIT = 3 #bad hit
FLAG_UNMAP = 0x4
FLAG_REVERSE = 0x10
import optparse
if __name__ == '__main__':
    parser = optparse.OptionParser(usage)
    parser.add_option('--wgsim', action = 'store_true', dest='wgsim', help='reads are generated by wgsim', default =False)
    parser.add_option('--art', action = 'store_true', dest='art', help='reads are generated by art', default=False)
    opt, args = parser.parse_args()
    if(opt.wgsim):
        wgsim.wgsim_alnEval(sys.argv[1:])
    elif(opt.art):
        art.art_alnEval(sys.argv[1:])
    else:
        print >>sys.stderr, '[CMD ERROR]:%s'%(sys.argv)
        parser.print_help()
        exit(1)



