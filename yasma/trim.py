
import sys
import os

import click
from click_option_group import optgroup

from pathlib import Path, PurePath
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint

# from random import sample

# import numpy as np
# from statistics import quantiles
# import math
import shutil
# import re

from .generics import *
from .cli import cli








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




def trim(**params):
	'''Wrapper for trimming using cutadapt.'''


	rc = requirementClass()
	rc.add_cutadapt()
	rc.check()

	ic = inputClass(params)
	ic.check(['untrimmed_libraries', 'adapter'])


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
		sys.exit("Libraries PRE-TRIMMED... Skipping trim.")


	trimmed_libraries = []

	for file in untrimmed_libraries:


		path = PurePath(file)



		suffixes = file.suffixes
		for i,s in enumerate(suffixes):
			if s.endswith(('.fastq', '.fq')):
				suffixes[i] = '.t.fq'
			elif s.endswith(('.fasta', '.fa')):
				suffixes[i] = '.t.fa'

		out_file = Path(output_directory, 'trim', file.stem.split('.')[0] + "".join(suffixes))

		call = ["cutadapt", "-a", adapter, "--minimum-length", str(min_length), "--maximum-length", str(max_length), "-O", "4", "--max-n", "0", "--trimmed-only", "-o", out_file, file]

		p = Popen(call, stdout=PIPE, encoding=ENCODING)

		for line in p.stdout:
			print(line.strip())


		p.wait()

		trimmed_libraries.append(Path(out_file))



	ic.inputs['trimmed_libraries'] = trimmed_libraries
	ic.write()




