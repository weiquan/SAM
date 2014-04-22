#stat diff from reads files


def readfq(fp):
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break
import re


def wgsim_statdiff(fn, stat):
    stat['indel'][0] = 0
    stat['indel'][1] = 0
    stat['indel'][2] = 0
    stat['diff'][0] = 0
    stat['diff'][1] = 0
    stat['diff'][2] = 0
    stat['diff'][3] = 0
    stat['diff'][4] = 0
    stat['diff'][5] = 0
    stat['mismatch'][0] = 0
    stat['mismatch'][1] = 0
    stat['mismatch'][2] = 0
    stat['mismatch'][3] = 0
    stat['mismatch'][4] = 0
    stat['mismatch'][5] = 0

    fp = open(fn)

    for name, seq, qual in readfq(fp):
        r = re.findall('\d+_\d+_\d+', name)
        for i in xrange(2):
            n_err, n_sub, n_indel = r[i].split('_')
            n_err = int(n_err)
            n_sub = int(n_sub)
            n_indel = int(n_indel)
            n_mismatch = n_err + n_sub
            n_diff = n_mismatch + n_indel
            if n_indel == 0:
                stat['index'][0] += 1
                if n_mismatch >= 5:
                    stat['mismatch'][5] += 1
                else:
                    stat['mismatch'][n_mismatch] += 1
            else:
                stat['indel'][n_indel] += 1
            if n_diff >= 5:
                stat['diff'][5] += 1
            else:
                stat['diff'][n_diff] += 1
    fp.close()
    return 1


def art_statdiff(fn, stat):
    print 'wgsim'
    return 1


def stat_wrap(simulater, fn, stat):
    if simulater == 'wgsim':
        return wgsim_statdiff(fn, stat)
    elif simulater == 'art':
        return art_statdiff(fn, stat)
if __name__ == '__main__':
    fn = 'in.txt'
    stat = {'indel': [0, 0, 0],
            'diff': [0, 0, 0, 0, 0, 0],
            'mismatch': [0, 0, 0, 0, 0, 0]}

    stat_wrap('wgsim', fn, stat)
