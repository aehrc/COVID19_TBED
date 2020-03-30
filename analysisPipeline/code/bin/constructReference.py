import sys
import numpy as np
import collections
import psutil

data = {}
with open(sys.argv[1]) as f:
	for line in f:
		line = line.rstrip()
		if line[0] == ">":
			name = line
			if name not in data:
				data[name] = ""
		else:
			data[name] = data[name] + str(line)

list_data = [list(data[i]) for i in data]
array_data = np.array(list_data)

refseq = []
for column in array_data.T:
	seq = collections.Counter(column).most_common(1)[0][0]
	refseq.append(seq)

refstring = ''.join(refseq)

with open(sys.argv[2], 'w') as outfile:
	outfile.write('>Reference Sequence\n')
	outfile.write("%s" % refstring)
