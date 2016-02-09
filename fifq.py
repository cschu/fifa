#!/usr/bin/env python
import sys

import argparse
from itertools import izip_longest

def anabl_getSeqsFromFastX(fn, X=2):
    it = open(fn)
    for head in it:
        try:
            yield (head.strip(), tuple(map(lambda x:x.strip(),
                                       ([it.next() for i in xrange(X - 1)]))))
        except:
            break
        pass
    it.close()
    pass

def getFastqIdentifier(string):
    if string.endswith('/1') or string.endswith('/2'):
        return string[1:-2]
    else:
        return string.split()[0][1:]

def extractSequences(inR1, outR1, keepSequences, inR2=None, outR2=None, fmt='fq', overwrite='ab'):
    fastx = 4 if fmt == 'fq' else 2
    fwdOut = open(outR1, overwrite)
    fwdGen = anabl_getSeqsFromFastX(inR1, X=fastx)

    revOut, revGen = None, None
    revSid, revSeq = None, None
    if inR2 is not None and outR2 is not None:
        revOut = open(outR2, overwrite)
        revGen = anabl_getSeqsFromFastX(inR2, X=fastx)

    fqid1, fqid2 = None, None

    #print fwdOut is None, fwdGen is None, revOut is None, revGen is None
    #sys.exit(1)
    seqCounter = 0
    while True:
        try:
            fwdSid, fwdSeq = fwdGen.next()
            fqid1 = getFastqIdentifier(fwdSid)
            if revGen is not None:
                revSid, revSeq = revGen.next()
                fqid2 = getFastqIdentifier(revSid)
        except:
            sys.stderr.write("Broke out of main loop after %i sequences.\n%i sequences still not found." % (seqCounter, len(keepSequences)))
            break
        seqCounter += 1
        if fqid1 != fqid2 and fqid2 is not None:
            sys.stderr.write('Error: fqid-mismatch %s %s.\n' % (fqid1, fqid2))
            sys.exit(1)
        if fqid1 in keepSequences:
            keepSequences.difference_update(set([fqid1]))
            fwdOut.write(('%s\n' * fastx) % ((fwdSid,) + fwdSeq))
            if revOut is not None:
                revOut.write(('%s\n' * fastx) % ((revSid,) + revSeq))
        else:
            pass
    fwdOut.close()
    if revOut is not None:
        revOut.close()

    #open(inR1 + '.missing', 'wb').write('\n'.join([id_ for id_ in keepSequences]))

    pass



def main(argv):

    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--in1', help='The r1-file (single-end reads or left paired-end reads).')
    parser.add_argument('--in2', help='The r2-file (right paired-end reads)')
    parser.add_argument('--input-format', help='Input sequences stored in fa or fq file(s).', default='fq')
    parser.add_argument('--out1', type=str, help='The r1-output file.')
    parser.add_argument('--out2', type=str, help='The r2-output file.')
    parser.add_argument('--keep', type=str, help='List (file) of sequence identifies to keep.')
    args = parser.parse_args()

    #print args.keep
    #print keepSequences
    #return None

    inputs1, inputs2 = args.in1.split(','), [None]
    if 'in2' in args:
        inputs2 = args.in2.split(',')
    # print inputs1, inputs2
    # return None

    keepSequences = set([sid.strip() for sid in open(args.keep)])
    for in1, in2 in izip_longest(inputs1, inputs2):
        print in1, in2
        extractSequences(in1, args.out1, keepSequences, inR2=in2, outR2=args.out2)

    open(args.keep + '.missing', 'wb').write('\n'.join([id_ for id_ in keepSequences]))



"""
wanted = set([sid.strip() for sid in open(sys.argv[2])])
total = len(wanted)
extracted = set()

for sid, seq, sep, qual in anabl_getSeqsFromFastX(sys.argv[1]):
	found


for sid, seq in anabl_getContigsFromFASTA(sys.argv[1]):
	found = False
	for wid in wanted:
		if sid.startswith(wid):
			found = True
			break
	if found:
		wanted.difference_update(set([wid]))
		extracted.add(wid)
		sys.stdout.write('>%s\n%s\n' % (sid, seq))


if wanted:
	sys.stderr.write('%i sequence(s) were not found in %s.\n' % (len(wanted), sys.argv[1]))
	for sid in wanted:
		sys.stderr.write('%s\n' % sid)
"""
main(sys.argv[1:])
