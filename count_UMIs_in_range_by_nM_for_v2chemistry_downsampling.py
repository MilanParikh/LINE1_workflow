import pysam
import sys
import random

bamfilename = sys.argv[1]
range = sys.argv[2]
max_edit_distance = int(sys.argv[3])
strand = sys.argv[4]
clip_max = int(sys.argv[5])
downsample_to = float(sys.argv[6])

bam = pysam.AlignmentFile(bamfilename,'rb')

UMI_dict = dict()

for read in bam.fetch(region=range):
	if random.random() < downsample_to and read.has_tag("nM") and read.get_tag("nM") <= max_edit_distance and read.has_tag("CB") and read.has_tag("UB") and read.is_reverse == (strand == '+') and 'N' not in read.cigarstring:
		if (read.cigartuples[0][0] != 4 or (read.cigartuples[0][1] <= clip_max)) and (read.cigartuples[-1][0] != 4 or (read.cigartuples[-1][1] <= clip_max)):
			CB = read.get_tag("CB")
			UB = read.get_tag("UB")
			if CB not in UMI_dict:
				UMI_dict[CB] = set()
			UMI_dict[CB].add(UB)

for CB in UMI_dict:
	print(CB+'\t'+str(len(UMI_dict[CB])))
