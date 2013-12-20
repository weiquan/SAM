INDEL_penalty = 1
Mismatch_penalty = 1
Match_penalty = 0
def Match(a, b):
	if a == b:
		return Match_penalty
	else:
		return Mismatch_penalty
#http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance
def lev0(a, b):
    if not a: return len(b)
    if not b: return len(a)
    return min(lev0(a[1:], b[1:])+(Match(a[0], b[0])), lev0(a[1:], b)+INDEL_penalty, lev0(a, b[1:])+INDEL_penalty)
def lev1(seq1, seq2):
    oneago = None
    thisrow = range(1, len(seq2) + 1) + [0]
    for x in xrange(len(seq1)):
        twoago, oneago, thisrow = oneago, thisrow, [0] * len(seq2) + [x + 1]
        for y in xrange(len(seq2)):
            delcost = oneago[y] + 1
            addcost = thisrow[y - 1] + 1
            subcost = oneago[y - 1] + (seq1[x] != seq2[y])
            thisrow[y] = min(delcost, addcost, subcost)
    return thisrow[len(seq2) - 1]

complNT = {'A':'T', 'T':'A', 'C':'G', 'G':'c'}
import optparse
import sys

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
                
def parseCigar(cigar):
    _encode = 'MIDNSHP'

    result = []
    for c in cigar:
        n = ''
        if c.isdigit():
            n += c
        elif c in _encode:
            if n == '':
                raise ValueError("end of CIGAR string reached, but an operator was expected")
            result.append((c, int(n)))
    return result

if __name__ == '__main__':
    usage = "usage: %prog [Options] <testFile.SAM> <Ref>"
    parser = optparse.OptionParser(usage)
    parser.add_option('-d', '--diff', action = 'store', type = 'int',  dest='max_diff', help = 'max diff between different SAM', default= 4)
    parser.add_option('-p', '--print', action = 'store_true', dest='print', help = 'print wrong alignment place', default= True)

    #get opt
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        exit(1)
    fn_sam, fn_ref = args[0], args[1] 
    fp_ref = open(fn_ref, 'r')
    index = {}
    for name, seq, qual in readfq(fp_ref):
        if name in index:
            print >>sys.stderr, 'more than one seq in %s named with %s'%(fn_ref, name)
            sys.exit(1)
        index[name] = seq
    fp_ref.close()

    fp_sam = open(fn_sam, 'r')
    for line in fp_sam:
        #skip the header
        if line[0] == '@':
            continue
        words = line.strip().split('\t')
        ans_chrom = words.split('_')[0]
        ans_pos0 = int(words.split('_')[1])
        ans_pos1 = int(words.split('_')[2])


        name = words[0]
        flag = int(words[1])
        chrom = words[2]
        pos = int(words[3])
        seq = words[9]
        rev_seq = [complNT[c] for c in seq]
        rev_seq = rev_seq[::-1]

        cigar = words[5]
        cigar_list = parseCigar(cigar)
        read_start = 0
        read_end = len(seq)
        if cigar_list[0][0] == 'S':
            read_start = cigar_list[0][1] 
        if cigar_list[-1][0] == 'S':
            read_end -= cigar_list[-1][1]
        read = seq[read_start:read_end]
        ref = index[chrom][pos:pos+read_end-read_start]
        d = lev1(ref, read)

        l = len(seq)
        ref0 = index[ans_chrom][ans_pos0:ans_pos0+l]
        read0 = seq
        if flag & 0x0010 != 0:
            read0 = rev_seq
        d1 = lev1(ref0, read0)

        ref1 = index[ans_chrom][ans_pos1-l+1:ans_pos1]
        read1 = rev_seq
        if flag & 0x0010 == 0:
            read0 = seq
        d2 = lev1(ref1, read1)
        if d < min(d1, d2):
            print '[0]: >d_ans < d_mapped'
            print line.strip()
        elif d == min(d1, d2):
            print '[1]: >ans and mapped pos have the same distance'
            print line.strip()
        else:
            print '[0]: >d_ans > d_mapped'
            print line.strip()

    fp_sam.close()
    



