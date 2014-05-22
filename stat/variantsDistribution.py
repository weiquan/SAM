def readfq(fp):
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break
import optparse
def stat_wgsim(opt, fn):
    MAX_DIFF = 10
    fp = open(fn)
    stat_table=[0 for i in xrange(MAX_DIFF)]
    stat_err=[0 for i in xrange(MAX_DIFF)]
    stat_sub =[0 for i in xrange(MAX_DIFF)]
    stat_indel = [0 for i in xrange(MAX_DIFF)]
    for name, seq, qual in readfq(fp):
        mutations_str0 = name.split('_')[4]
        n_err, n_sub, n_indel = mutations_str0.split(':')
        n_err = min(int(n_err), MAX_DIFF-1)
        stat_err[n_err] += 1
        n_sub = min(int(n_sub), MAX_DIFF-1)
        stat_sub[n_sub] += 1
        n_indel = min(int(n_indel), MAX_DIFF-1)
        stat_indel[n_indel] += 1
        stat_table[min(n_err+n_sub+n_indel, 9)] += 1
        mutations_str1 = name.split('_')[5]
        n_err, n_sub, n_indel = mutations_str1.split(':')
        n_err = min(int(n_err), MAX_DIFF-1)
        stat_err[n_err] += 1
        n_sub = min(int(n_sub), MAX_DIFF-1)
        stat_sub[n_sub] += 1
        n_indel = min(int(n_indel), MAX_DIFF-1)
        stat_indel[n_indel] += 1
        stat_table[min(n_err+n_sub+n_indel, 9)] += 1
    print '===mutations stat==='
    print stat_table

    print '===sequencing erros stat==='
    print stat_err
    
    print '===SNP stat==='
    print stat_sub

    print '===Indel stat==='
    print stat_indel
    fp.close()

def __parseCigarIter(cigar):
    import re
    p = re.compile(r'(\d+)([MIDNSHPX=])')
    return p.finditer(cigar)    



def stat_art(opt, fn):
    MAX_DIFF = 10
    #fn must be sam format
    stat_table=[0 for i in xrange(MAX_DIFF)]
    stat_sub =[0 for i in xrange(MAX_DIFF)]
    stat_indel = [0 for i in xrange(MAX_DIFF)]

    fp = open(fn, 'r')
    for line in fp:
        line = line.strip()
        if line[0] == '@':#skip header
            continue
        words = line.split('\t')
        cigar = words[5]
        n_diff = 0
        n_sub = 0
        n_indel = 0
        for i, c in __parseCigarIter(cigar):
            if c != '=':
                n_diff += i

            elif c == 'X':
                n_sub += i
            elif c == 'I' or c == 'D':
                n_indel += i
        n_diff = min(int(n_diff), MAX_DIFF-1)
        n_sub = min(int(n_sub), MAX_DIFF-1)
        n_indel = min(int(n_indel), MAX_DIFF-1)
        stat_table[n_diff] += 1
        stat_sub[n_sub] += 1
        stat_indel[n_indel] += 1

    print '===mutations stat==='
    print stat_table

    print '===SNP stat==='
    print stat_sub

    print '===Indel stat==='
    print stat_indel
    fp.close()
import sys
if __name__ == '__main__':
    usage = "usage: %prog [Options] <Reads.fq> "
    parser = optparse.OptionParser(usage)
    parser.add_option('-s', '--simulator', action = 'store', type = 'str', dest='simulator', help = 'reads simulaotr generator', default= 'wgsim')
    #parser.add_option('-p', '--print', action = 'store_true', dest='print', help = 'print wrong alignment place', default= True)
    options, args = parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        sys.exit(1)
    if options.simulator == 'wgsim':
        stat_wgsim(options, args[0])
    elif options.simulator == 'art':
        stat_art(options, args[0])