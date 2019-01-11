import gzip

from Bio import SeqIO

selected_reads = []

sample_name = snakemake.wildcards.sample
counter = 1
with open(snakemake.input[0], 'rt') as handle:
	for record in SeqIO.parse(handle, 'fasta'):
		if len(record.seq) >= int(snakemake.params.min_len):
			record.id = f"{sample_name}_{counter:07d}_length_{len(record.seq)}"
			record.description = ''
			counter += 1
			selected_reads.append(record)

with gzip.open(snakemake.output[0], 'wt') as handle:
	SeqIO.write(selected_reads, handle, 'fasta')
