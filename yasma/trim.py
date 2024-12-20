

from .generics import *

from shutil import rmtree










@cli.command(group='Processing', help_priority=2)



@optgroup.group('\n  Basic options',
				help='')


@optgroup.option("-ul", "--untrimmed_libraries", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_glob_path,
	multiple=True,
	help='Path to untrimmed libraries. Accepts wildcards (*).')



@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")


@optgroup.option("-a", "--adapter", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	type=str,
	help="Adapter sequence which is meant to be trimmed.")




@optgroup.group('\n  Cutadapt options',
				help='')

@optgroup.option("--min_length",
	default = 15,
	help= 'Minimum allowed size for a trimmed read. (default 10)')


@optgroup.option("--max_length",
	default = 50,
	help= 'Maxiumum allowed size for a trimmed read. (default 50)')


@optgroup.option("--cores",
	default = 1,
	type=int,
	help= 'Number of CPU cores to use. 0 has cutadapt "autodetect" the number of cores (default 1)')


@optgroup.option('--cleanup', is_flag=True, default=False, help='Removes download and untrimmed data to save space')


def trim(**params):
	'''Wrapper for trimming using cutadapt.'''


	rc = requirementClass()
	rc.add_cutadapt()
	rc.check()

	ic = inputClass(params)
	ic.check(['untrimmed_libraries', 'adapter'])

	ic.inputs['trimmed_libraries'] = []


	output_directory        = ic.output_directory
	untrimmed_libraries     = ic.inputs['untrimmed_libraries']
	adapter                 = ic.inputs['adapter']

	max_length              = params['max_length']
	min_length              = params['min_length']


	Path(output_directory, "trim").mkdir(parents=True, exist_ok=True)

	log_file = Path(output_directory,"trim/log.txt")
	sys.stdout = Logger(log_file)



	if adapter == "PRE-TRIMMED":
		ic.inputs['trimmed_libraries'] = untrimmed_libraries.copy()
		ic.write()
		print('trimming complete')
		sys.exit("Libraries PRE-TRIMMED... Skipping trim.")
	elif adapter == 'None':
		sys.exit("Adapter identification failed. Trimming not possible. Are you sure this is sRNA-seq?")


	trimmed_libraries = []

	for file in untrimmed_libraries:
		print(f"trimming: {file}", flush=True)


		path = Path(file)



		suffixes = file.suffixes
		for i,s in enumerate(suffixes):
			if s.endswith(('.fastq', '.fq')):
				suffixes[i] = '.t.fq'
			elif s.endswith(('.fasta', '.fa')):
				suffixes[i] = '.t.fa'

		out_file = Path(output_directory, 'trim', file.stem.split('.')[0] + ".t.fq.gz")

		# cutadapt -a [adapter] --minimum-length 15--maximum-length 50-O 4 --max-n 0 --trimmed-only -o [out_file] [file]

		call = ["cutadapt", "-a", adapter, "--minimum-length", str(min_length), "--maximum-length", str(max_length), '-j', str(params['cores']), "-O", "4", "--max-n", "0", "--trimmed-only", "-o", out_file, file]

		p = Popen(call, stdout=PIPE, encoding=ENCODING)

		for line in p.stdout:
			print(line.strip())


		p.wait()

		trimmed_libraries.append(Path(out_file))



	if params['cleanup']:

		print("Warning: flag 'cleanup' activated. SRR downloads and untrimmed files will be deleted in 10 seconds...")
		time.sleep(10)

		for srr in ic.inputs['srrs']:

			try:
				rmtree(Path(output_directory, 'download', srr))
			except FileNotFoundError:
				pass

		for file in ic.inputs['untrimmed_libraries']:

			try:
				Path(file).unlink()
			except FileNotFoundError:
				pass

		print("  -> deletion successful. Scrubbing inputs.json")

		ic.inputs['untrimmed_libraries'] = []
		ic.write()



	ic.inputs['trimmed_libraries'] = trimmed_libraries
	ic.write()


	print("trimming complete!")

