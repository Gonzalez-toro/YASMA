#!/usr/bin/env python3 

import sys
import os
from pprint import pprint

from subprocess import PIPE, Popen, call
from pathlib import Path

from os.path import isfile

from time import time
from collections import Counter, deque
from itertools import count, chain

from statistics import median, mean
import argparse



parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam_file', 
	nargs="?",
	required=True,
	help='bamfile of aligned small RNAs (tested with shortstack)')

parser.add_argument('-o', '--output_directory', 
	nargs="?",
	default=f'SmoothLoc_{round(time())}',
	help='working folder for locus analysis')


parser.add_argument('-d', '--dicercall', 
	nargs='+', 
	help='list of sRNA sizes to be analyzed separately', 
	default = [20,21,22,23,24])


parser.add_argument('-r', '--readgroups', 
	nargs='+', 
	help='list of read groups (libraries) to be considered for the annotation. All present read groups will be considered for counting.', 
	default = ['all'])


parser.add_argument('-f', '--force', 
	action='store_true',
	default=False,
	help='force remake of supporting files')


parser.add_argument('--partial_wigs', 
	action='store_true',
	default=False,
	help='only make wiggle files associated with essential functions (ignoring size and strand specific coverages)')


parser.add_argument('--window', 
	nargs="?",
	default=100,
	type=int,
	help='kernel desity bandwith window')


parser.add_argument('--merge_dist', 
	nargs="?",
	default=150,
	type=int,
	help='maximum gap size between valid regions to merge to a single locus')

parser.add_argument('--pad', 
	nargs="?",
	default=10,
	type=int,
	help='number of bases added to either end of a defined locus (arbitrary)')

parser.add_argument('--rpm', 
	nargs="?",
	default=1.0,
	type=float,
	help='rpm depth threshold for nucleating a locus')

parser.add_argument('--dicer_ratio', 
	nargs="?",
	default=3,
	type=float,
	help='ratio of dicer to non-dicer reads to be considered for a locus region')


parser.add_argument('--extension_ratio', 
	nargs="?",
	default=0.5,
	type=float,
	help="fraction of rpm threshold to be considered for extending a locus' boundaries")


def check_reqs():
	tool_responses = {
	'samtools version' : 'Samtools compilation details:',
	'gt --version' : 'gt (GenomeTools)',
	'bgzip --version' : 'bgzip (htslib)',
	'tabix --version' : 'tabix (htslib)',
	'wigToBigWig' : 'wigToBigWig v 2.8',
	}


	fails = []

	for tool, response in tool_responses.items():
		p = Popen(tool.split(), stdout=PIPE, stderr=PIPE, encoding='utf-8')

		out, err = p.communicate()

		merged = out + err

		tool = tool.split()[0]


		# print(out)
		# print(err)
		if response in merged:
			pass_str = "[x]"
		else:
			pass_str = "[ ]"
			fails.append(tool)




		print(" ", pass_str, tool)
		# sys.exit()

	do_not_prepare_gff = False
	do_not_make_bigwig = False

	if 'samtools' in fails:
		sys.exit("Error: samtools not found in PATH (required)")

	for tool in ['gt','bgzip','tabix']:
		if tool in fails:
			do_not_prepare_gff = True
			break

	if 'wigToBigWig' in fails:
		do_not_make_bigwig = True

	if do_not_prepare_gff:
		print("Warning: will not prepare indexed gff for jbrowse due to missing reqs")
	if do_not_make_bigwig:
		print("Warning: will not prepare bigwig files due to missing reqs")

	return(do_not_prepare_gff, do_not_make_bigwig)


class Logger(object):
	def __init__(self, file_name):
		self.terminal = sys.stdout
		self.file_name = file_name
		self.log = open(file_name, "w")
		# with open(file_name, "w") as outf:
		# 	outf.write("")

	def clear_ansi(self, message):
		return(message.replace("\033[1m", "").replace("\033[0m",""))

	def write(self, message):
		self.terminal.write(message)
		# with open(self.file_name, 'a') as outf:
		# 	outf.write(message)  
		self.log.write(self.clear_ansi(message))

	def flush(self):
		self.terminal.flush()
		self.log.flush()

def get_library_depth(out_dir, file):
	depth_file = f"./{out_dir}/library_depth.txt"
	if isfile(depth_file) and not force:
		with open(depth_file, 'r') as f:
			line = f.readline()
			line = line.split("\t")

			if line[0] == file:
				# print(f'\nread depth from {depth_file}...')
				return(int(line[1]))



	print('reading annotation RG depth...')

	call = ['samtools', 'view', '-F', '4']

	for rg in annotation_readgroups:
		call += ['-r', rg]

	call += [file]

	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

	depth = 0
	for line in p.stdout:
		depth += 1

	p.wait()


	with open(depth_file, 'w') as outf:
		print(file, depth, sep='\t', file=outf)

	return(depth)

def samtools_view(bam):

	if not isfile(f"{bam}.bai"):
		# call = f"samtools index {bam}"
		call= ['samtools','index',bam]
		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
		out,err=p.communicate()
		print(out)
		print(err)


	# call = f"samtools view -@ 4 -F 4 {bam}"
	call = ['samtools', 'view', '-F', '4', bam]
	# print(call)
	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

	for i,line in enumerate(p.stdout):

		line = line.strip().split()

		# read_id, flag, sam_chrom, sam_pos, _, length, _, _, _, _,_,_,_,_,_,_,_,_,rg= line

		read_id = line[0]
		flag = line[1]
		seq = line[9]

		if flag == "16":
			strand = '-'
		elif flag == '0':
			strand = "+"
		else:
			strand = False

		# print(line)
		length = int(line[5].rstrip("M"))
		# sam_pos = int(sam_pos)

		# length = len(line[9])

		sam_pos = int(line[3])
		sam_chrom = line[2]

		rg = line[18].lstrip("RG:Z:")


		if length in dcr_range:
			size = 'dcr'
		elif length in non_range:
			size = 'non'
		else:
			size = False



		yield(strand, length, size, sam_pos, sam_chrom, rg, seq)

	p.wait()

def get_chromosomes(file):
	chromosomes = []
	rgs = []
	# call = f"samtools view -@ 4 -H {file}"
	call = ['samtools','view','-H', file]
	# print(call)

	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
	out, err = p.communicate()

	for o in out.strip().split("\n"):
		o = o.strip().split('\t')
		if o[0] == "@SQ":
			name = o[1].split(":")[-1]
			length = int(o[2].split(":")[-1])
			chromosomes.append((name,length))
		if o[0] == "@RG":
			rgs.append(o[1].split(":")[-1])

	return(chromosomes, rgs)


def prepare_gff(gff_input):

	sorted_input = gff_input.replace(".gff3", ".sorted.gff3")
	zipped_input = sorted_input.replace(".gff3", ".gff3.gz")


	# print(gff_input)
	# print(sorted_input)
	# print(zipped_input)



	print("  sorting...")
	c2 = f"gt gff3 -retainids -sortlines -tidy {gff_input}"
	with open(sorted_input, 'w') as f:
		c2 = Popen(c2.split(), stdout=f)

	print("  zipping...")
	c3 = f"bgzip -f {sorted_input}"
	call(c3.split())

	print("  indexing...")
	c4 = f"tabix -f -p gff {zipped_input}"
	call(c4.split())



def inf_counter():
	i = 1
	while True:
		yield(i)
		i += 1

def test1(rpm):
	if rpm == None:
		return('-')

	if rpm >= rpm_cutoff:
		return('n')
	elif rpm >= ext_cutoff:
		return('e')
	else:
		return('-')

def test2(dcr, non):
	if dcr == None or non == None:
		return('-')

	if dcr >= non * dicer_ratio and dcr > 0:
		return('x')
	else:
		return('-')



class coverageClass():
	def __init__(self, bandwidth=0):
		self.ls = deque()
		self.bandwidth = bandwidth

	def get(self):
		try:
			d = self.ls.popleft()
		except IndexError:
			d = 0

		return(d)

	def add(self, length):
		for r in range(length + self.bandwidth + 1):
			try:
				self.ls[r] += 1
			except IndexError:
				self.ls.append(1)


def flatten_lol(ls):
	out = set()
	for l in ls:
		for e in l:
			out.add(e)

	return(out)





class locusClass():
	def __init__(self):
		self.reads = deque([])

		self.in_locus = False
		self.nucleated = False
		self.last_hit_pos = 0
		self.last_end = 0
		self.start = False
		self.stop  = False

	def hit(self, pos, hit_type):
		self.last_hit_pos = pos

		if not self.in_locus:
			self.start = pos

		self.in_locus = True


		if hit_type == 'nuc':
			self.nucleated = True


	def add(self, read):
		self.reads[-1].append(read)

	def check(self, pos):
		# print(self.reads)

		if self.in_locus:

			if pos - self.last_hit_pos > merge_dist:
				self.in_locus = False

				self.stop = self.last_hit_pos

				if self.nucleated:
					self.nucleated = False
					return(True)

		else:

			while True:

				if len(self.reads) == 0:
					break


				if len(self.reads[0]) == 0:
					self.reads.popleft()
				elif self.reads[0][0][0] + 35 < pos - merge_dist:
					self.reads.popleft()

				else:
					break
					
		self.reads.append([])






	def summarize(self):

		# start = self.reads[0][1]
		start = self.start - pad
		stop  = self.stop + pad


		# pprint(self.reads)
		# print(self.reads.values())

		reads = chain.from_iterable(self.reads)
		# reads = [r for r in flatten_lol(self.reads) if r[1] >= start and r[1] + r[2] <= stop]


		reads = [r for r in reads if r[0] + r[1] >= start and r[0] <= stop]

		# reads = flatten_lol([])

		len_c    = Counter()
		size_c   = Counter()
		strand_c = Counter()
		rg_c     = Counter()
		read_c   = Counter()

		read_starts = []
		read_stops  = []


		for r in reads:
			sam_pos, sam_length, sam_size, sam_strand, sam_rg, sam_read = r

			read_starts.append(sam_pos)
			read_stops.append(sam_pos + sam_length)

			len_c.update([sam_length])
			size_c.update([sam_size])
			strand_c.update([sam_strand])
			rg_c.update([sam_rg])
			read_c.update([sam_read])


		# print(start, stop)
		start = min(read_starts) - pad
		stop  = max(read_stops) + pad


		name = f"Cl_{next(cluster_counter)}"


		if len(reads) == 0:
			print(f"WARNING: {name} detected no reads. Likely an error. Skipping.")
			return(0,0)

		# pprint(len_c)
		# pprint(read_c)


		dist_to_last = start - self.last_end
		self.last_end = stop

		coords = f"{chrom}:{start}..{stop}"
		length = stop - start
		n_reads = len(reads)
		frac_top = round(strand_c["+"] / n_reads,3)
		frac_dicercall = round(size_c['dcr'] / n_reads, 3)
		rpm = round(n_reads / library_depth * 1000000, 2)


		cum_count = 0
		top_reads = read_c.most_common(100)
		with open(reads_file, 'a') as outf:
			for rank, read in enumerate(top_reads):

				seq, dep = read
				rpm = round(dep * read_equivalent, 4)

				cum_count += dep

				loc_prop = round(cum_count / n_reads, 4)

				print(name, seq, rank, dep, rpm, loc_prop, file=outf, sep='\t')

				if loc_prop >= 0.3:
					break




		predominant_length, predominant_length_depth = len_c.most_common(1)[0]
		predominant_length_depth = round(predominant_length_depth/n_reads,3)
		# print(predominant_length, predominant_length_depth)

		if frac_top >= 0.8:
			strand = '+'
		elif frac_top <= 0.2:
			strand = "-"
		else:
			strand = '.'

		depth_by_length = round(n_reads / length, 3)


		to_print = [name, coords, length, dist_to_last, n_reads, rpm, depth_by_length, frac_top, strand]
		to_print += [frac_dicercall, size_c['dcr'],  size_c['non']]
		to_print += [predominant_length, predominant_length_depth]
		to_print += [len_c[d] for d in dcr_range]
		to_print = "\t".join(map(str, to_print))


		with open(results_file, 'a') as outf:
			print(to_print, sep='\t', file=outf)

		# print(name, coords, length, dist_to_last, n_reads, rpm, depth_by_length, frac_top, strand, 
		# 		frac_dicercall, size_c['dcr'],  size_c['non'], 
		# 		predominant_length, round(predominant_length_depth/n_reads,3), "\t".join([str(len_c[d]) for d in dcr_range]), sep='\t')


		# sys.exit()

		with open(gff_file, 'a') as outf:
			print(f"{chrom}\tsmoothLocus\tnc_RNA\t{start}\t{stop}\t.\t.\t.\tID={name};dicercall={predominant_length};frac_dicercall={predominant_length_depth}", file=outf)


		to_print = [name]

		to_print.append(sum([rg_c[rg] for rg in annotation_readgroups]))
		to_print.append(sum([rg_c[rg] for rg in bam_rgs]))
		to_print += [rg_c[rg] for rg in bam_rgs]
		to_print = "\t".join(map(str, to_print))

		with open(count_file, 'a') as outf:
			print(to_print, file=outf)


		# self.reads = {}

		return(length, dist_to_last)



class wiggleClass():
	def __init__(self, file):
		self.file = f"./{out_dir}/coverages/{file}.wig"
		self.outf = open(self.file, 'w')
		self.reset()


	def reset(self):
		self.val = 0
		self.start_pos = 1


	def add(self, val, pos):


		if val != self.val:
			span = pos - self.start_pos

			if span > 0:

				print(f"variableStep chrom={chrom} span={span}", file=self.outf)
				print(f"{self.start_pos} {self.val}", file=self.outf)

				self.val = val
				self.start_pos = pos

	def convert(self, cleanup=False):

		self.outf.close()

		wig = self.file

		bigwig = wig.replace(".wig", ".bigwig")
		print(f"  {wig} -> {bigwig}", flush=True)

		call = f"wigToBigWig {wig} ./{out_dir}/chrom.sizes.txt {bigwig}"

		p = Popen(call.split(), stdout=PIPE, stderr=PIPE, encoding='utf-8')

		out, err= p.communicate()

		if out.strip() + err.strip() != "":

			print(out)
			print(err)

		if cleanup:
			os.remove(wig)


def main():
	global dcr_range, non_range, half_window, window, pad, rpm_cutoff, ext_cutoff, merge_dist, chrom, library_depth, cluster_counter, annotation_readgroups, bam_rgs, force, out_dir, dicer_ratio, results_file, gff_file, count_file, reads_file, read_equivalent



	args = parser.parse_args()
	input_bam  = args.bam_file

	out_dir    = args.output_directory
	force = args.force
	# rgs        = args.readgroups
	dicercall      = args.dicercall
	dicer_ratio    = args.dicer_ratio
	window = args.window
	annotation_readgroups = args.readgroups

	pad = args.pad
	merge_dist = args.merge_dist

	rpm_cutoff = args.rpm
	partial_wigs = args.partial_wigs


	extension_multiplier = args.extension_ratio

	assert isfile(input_bam), f"{input_bam} does not exist"

	Path(out_dir).mkdir(parents=True, exist_ok=True)
	Path(f'./{out_dir}/coverages').mkdir(parents=True, exist_ok=True)


	log_file = f"{out_dir}/Log.txt"

	if isfile(log_file) and not force:
		sys.exit(f"Error: Log file is already exists ({log_file}). The annotator will not over-write by default (use --force to override). Be warned: this will trigger the overwrite of all files in this folder!")

	sys.stdout = Logger(log_file)


	dicercall = [int(d) for d in dicercall]
	dcr_range = set([r for r in range(min(dicercall), max(dicercall) + 1)])
	non_range = set([r for r in range(15,30) if r not in dcr_range])


	assert window % 2 == 0, "Window must be an even number!"
	half_window = int(window / 2)




	print()
	print()
	print("\033[1m-- annotator v0.3x --\033[0m")

	print()
	print()
	print(f"\033[1m[Prerequisites]\033[0m")

	do_not_prepare_gff, do_not_make_bigwig = check_reqs()

	# rpm_cutoff = round(rpm_cutoff / window, 6)


	chromosomes, bam_rgs = get_chromosomes(input_bam)



	with open(f"./{out_dir}/chrom.sizes.txt", 'w') as outf:
		for chrom, size in chromosomes:
			print(chrom, size, sep='\t', file=outf)


	if annotation_readgroups[0].lower() == 'all':
		annotation_readgroups = bam_rgs
	else:
		for rg in annotation_readgroups:
			if rg not in bam_rgs:
				sys.exit(f"Error: submitted readgroup '{rg}' not found within bamfile header:\n{bam_rgs}")

	annotation_readgroups = set(annotation_readgroups)
		# def get(self):





	## initiating output files
	gff_file   = f"{out_dir}/Annotation.gff3"

	with open(gff_file, 'w') as outf:
		print("##gff-version 3", file=outf)

		for c, l in chromosomes:
			print(f"##sequence-region   {c} 1 {l}", file=outf)


	count_file = f'{out_dir}/Counts.txt'
	with open(count_file, 'w') as outf:
		print("cluster", 'ann_depth', 'tot_depth', "\t".join(bam_rgs), sep='\t', file=outf)



	results_file = f"{out_dir}/Results.txt"
	with open(results_file, 'w') as outf:
		print("#name\tlocus\tlength\tgap\tdepth\trpm\tdepth:length\tfrac_top\tstrand\tfrac_dicer\tdcr_reads\tnon_reads\tdicercall\tfrac_dicercall\t" + "\t".join(map(str, dcr_range)), file=outf)


	reads_file = f"{out_dir}/TopReads.txt"
	with open(reads_file, 'w') as outf:
		print("cluster\tseq\trank\tdepth\trpm\tlocus_prop", file=outf)




	# wig_densities = {'dcr' : wiggleClass(f"./{out_dir}/dcr.dens.wig"), 'non' : wiggleClass(f"./{out_dir}/non.dens.wig")}
	# wig_rpm_pass = wiggleClass(f"./{out_dir}/rpm_passing.wig")
	# wig_ratio_pass = wiggleClass(f"./{out_dir}/ratio_passing.wig")
	# wig_pass = wiggleClass(f"./{out_dir}/passing_all.wig")


	print()
	print(f"\033[1m[General settings]\033[0m")
	print(f"             input_bam: {input_bam}")
	print(f"      output_directory: {out_dir}")
	print(f" annotation_readgroups: {list(annotation_readgroups)}")
	print(f"           dicer_sizes: {list(dcr_range)}")
	print(f"                 force: {force}")
	print(f"          partial_wigs: {partial_wigs}")
	print(f"              log_file: {log_file}")
	print()


	print(f"\033[1m[Annotation settings]\033[0m")
	print(f"     window: {window}")
	print(f" merge_dist: {merge_dist}")
	print(f"        pad: {pad}")



	library_depth = get_library_depth(out_dir, input_bam)


	read_equivalent = 1 / library_depth * 1000000
	depth_cutoff = library_depth / rpm_cutoff / 1000000
	ext_cutoff = rpm_cutoff * extension_multiplier


	print()
	print('\033[1m[Depth settings]\033[0m')
	print(f'    ann_rg_depth: {library_depth:,} reads')
	print(f'          1 read: {round(read_equivalent,5)} rpm')
	print(f"      rpm_cutoff: {rpm_cutoff} rpm -> {round(depth_cutoff,2)} reads")
	print(f"      ext_cutoff: {ext_cutoff} rpm > {round(depth_cutoff*extension_multiplier,2)} reads")
	print(f" extension_ratio: {extension_multiplier}")
	print(f"     dicer_ratio: {dicer_ratio}")



	coverage_names = ['dcr','non']

	if not partial_wigs:
		for s in ["+","-"]:
			for l in list(dcr_range) + ['non']:
				coverage_names.append(f"{l}{s}")


	wig_d = {c : wiggleClass(c) for c in coverage_names + ['rpm_passing', 'ratio_passing', 'passing_all']}


	cluster_counter = inf_counter()

	sam_iter = samtools_view(input_bam)
	read = next(sam_iter)
	sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read = read


	total_locus_count = 0
	chrom_count = 0

	for chrom, chrom_length in chromosomes:
		chrom_count += 1


		print()
		print()
		print(f"{chrom_count} / {len(chromosomes)}")
		print(f"chrom: {chrom}")
		print(f"       {chrom_length:,} bp")
		pos = 0

		print("  ", end='')

		locus_lengths = []
		locus_gaps = []



		coverages = {c : coverageClass() for c in coverage_names}
		coverage_buffer = {c : deque([0]*half_window) for c in coverage_names}

		window_coverages = {'dcr' : coverageClass(window), 'non' : coverageClass(window)}
		locus = locusClass()



		while sam_chrom == chrom:

			# if pos == 393471:
			# 	sys.exit()

			# if pos == 1000000:
			# 	sys.exit("timeup!")


			corrected_pos = pos - half_window


			if pos == sam_pos:

				if sam_size:

					if sam_rg in annotation_readgroups:

						coverages[sam_size].add(sam_length)
						window_coverages[sam_size].add(sam_length)

						if not partial_wigs:
							if sam_size == 'non':
								coverages[f'{sam_size}{sam_strand}'].add(sam_length)
							else:
								coverages[f'{sam_length}{sam_strand}'].add(sam_length)



					locus.add((sam_pos, sam_length, sam_size, sam_strand, sam_rg, sam_read))

				try:
					read = next(sam_iter)
					sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read = read
					# print(pos, sam_id, sam_size, sep='\t')
				except StopIteration:
					break



			elif pos < sam_pos:

				read_count = {}
				# dens_rpm = {}
				# dens = {}
				win_cov = {}
				rpms = {}



				for size in coverage_names:

					cov = coverages[size].get()
					coverage_buffer[size].append(cov)
					coverage_buffer[size].popleft()
					wig_d[size].add(round(coverage_buffer[size][0] * read_equivalent,4), corrected_pos)

					if size in ['dcr','non']:
						win_cov[size] = window_coverages[size].get()
						rpms[size] = round(win_cov[size] * read_equivalent, 4)


				t1 = test1(round(coverage_buffer['dcr'][0] * read_equivalent, 4))
				t2 = test2(win_cov['dcr'], win_cov['non'])

				tests = f"{t1}{t2}"


				# print(tests)

				if t1 == 'n':
					wig_d['rpm_passing'].add(1, corrected_pos)

				elif t1 == 'e':
					wig_d['rpm_passing'].add(0.5, corrected_pos)

				else:
					wig_d['rpm_passing'].add(0, corrected_pos)


				if not win_cov['dcr'] or win_cov['dcr'] == 0:
					ratio = 0
				else:
					try:
						ratio = round(win_cov['dcr'] / win_cov['non'], 2)
					except ZeroDivisionError:
						ratio = dicer_ratio


				wig_d['ratio_passing'].add(ratio, corrected_pos)


				if tests == "nx":
					locus.hit(corrected_pos, 'nuc')
					wig_d['passing_all'].add(1, corrected_pos)

				elif tests == "n-" or t1 == 'e':
					locus.hit(corrected_pos, 'ext')
					wig_d['passing_all'].add(0.5, corrected_pos)

				else:
					wig_d['passing_all'].add(0, corrected_pos)



				# if tests and "-" not in tests:
				# if rds['dcr'] and rds['dcr'] > 0 and pos-half_window > 0:
				# if pos >
				# if "n" in tests or 'e' in tests:
				# 	print(chrom, pos-half_window, 
				# 		"||", coverage_buffer['dcr'][0], round(dens['dcr'], 4), dens_rpm['dcr'], rds['dcr'], 
				# 		"||", coverage_buffer['non'][0], round(dens['non'], 4), dens_rpm['non'], rds['non'],
				# 		'||', tests,
				# 		sep='\t')




				if locus.check(corrected_pos):
					length, gap = locus.summarize()

					locus_lengths.append(length)
					locus_gaps.append(gap)
					total_locus_count += 1


				pos += 1


				if pos % 100000 == 0:
					print(".", end='', flush=True)
				if pos % 1000000 == 0:
					print(" ", end='', flush=True)




		for key in wig_d.keys():
			wig_d[key].add(0,pos)
			wig_d[key].reset()
		# for size in ['dcr','non']:
		# 	wig_densities[size].add(0, pos)
		# wig_rpm_pass.add(0, pos)
		# wig_pass.add(0, pos)

		locus_count = len(locus_gaps)
		med_length = median(locus_lengths)
		med_gap = median(locus_gaps)
		print()
		print(f"  {locus_count:,} loci found")
		print(f"  {med_length:,} median length")
		print(f"  {med_gap:,} median gap")

		# break


	print()
	print(f"{total_locus_count:,} loci found in total")

	print()

	if not do_not_make_bigwig:
		print("converting wigs to bigwigs...")
		for key in wig_d.keys():
			wig_d[key].convert()
	else:
		print("Not making bigwig files due to missing req...")


	print()
	if not do_not_prepare_gff:
		print("indexing GFF for jbrowse...")
		prepare_gff(gff_file)
	else:
		print("Not preparing indexed gff due to missing reqs...")


if __name__ == "__main__":
	main()








