#!/usr/bin/python3
import sys
import math
ifile = sys.argv[1]
total_length = 0
alphabet = {}
i = 0
for line in open(ifile, 'r'):
	if i%2 != 0:
		total_length += len(line.strip())
		for s in line.strip():
			if s in alphabet:
				alphabet[s] += 1
			else:
				alphabet[s] = 1
	i += 1
print("total length ", total_length)
print("alphabet", alphabet)
print("LG ", math.log(total_length, len(alphabet)))

print(alphabet)

k = 0
size = sum(alphabet.values())
for count in alphabet.values():
	k += -math.log(count / size, len(alphabet)) * (count / size)

print("uk = ", math.log(total_length, len(alphabet)))
print("fk = ", math.log(total_length, len(alphabet)) / k)
print("k = ", math.floor(math.log(total_length, len(alphabet)) / k))
