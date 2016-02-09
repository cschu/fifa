#!/usr/bin/env python
import sys


def anabl_getContigsFromFASTA(fn):
    """
    Returns generator object to access sequences from a multi-FASTA file.
    Originates from 'anabl' - BLAST analysing tool, hence the prefix.
    """
    head, seq = None, ''
    for line in open(fn):
        if line[0] == '>':
            if head is not None:
                yield (head, seq) 
            head, seq = line.strip().strip('>'), ''
        else:
            seq += line.strip()
    yield (head, seq) 


wanted = set([sid.strip() for sid in open(sys.argv[2])])
total = len(wanted)
extracted = set()

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
