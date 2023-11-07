

import sys
import os

import click
from click_option_group import optgroup

from pathlib import Path
import shutil
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint
from random import sample, seed

import numpy as np
# from statistics import quantiles
import math

from .generics import *
from .cli import cli

from statistics import mean, median, StatisticsError

from time import time

import re







@cli.command(group='Utilities', help_priority=4)



@optgroup.group('\n  Required options',
                help='')
@optgroup.option('-d', '--depth',
	required=True,
	help="Allows the user to subsample alignments for the annotation to a defined depth. Accepts an integer number of reads, which can be modified with a 10^3 prefix (ex. 10M).")

@optgroup.group('\n  Input options',
                help='')

@optgroup.option("-a", "--alignment_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
	help='Alignment file input (bam or cram).')

@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")

@optgroup.option('-r', '--annotation_readgroups', 
	required=False,
	multiple=True,
	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")




@optgroup.group('\n  Other options',
                help='')
@optgroup.option('--seed',
	default=0,
	type=int,
	help="Seed used for sampling.")

@optgroup.option('--compression',
	default='bam',
	type=click.Choice(['bam', 'cram']),
	help="Compression used for output sampled alignments (bam or cram).")
@optgroup.option('--override', is_flag=True, default=False, help='Overrides config file changes without prompting.')


@optgroup.option('--force', is_flag=True, default=False, help='force resubsample')

@optgroup.option('--init_dir', is_flag=True, default=False, 
	help='Initiate a new directory for alignments. This overrides the normal behavior to perform subsets in the source folder for an alignment.')


# def replace_parent(path, old_pattern, new_pattern):

# 	for p in path.parents:
# 		if str(p) == old_pattern:




def subsample(**params):
	'''Utility to subsample libraries to a specific read-count.'''

	rc = requirementClass()
	rc.add_samtools()
	rc.check()

	ic = inputClass(params)
	ic.check(['alignment_file', 'annotation_readgroups'])

	output_directory        = str(ic.output_directory)
	alignment_file          = ic.inputs["alignment_file"]
	annotation_readgroups   = ic.inputs['annotation_readgroups']
	project_name            = ic.inputs['project_name']

	target_depth            = params['depth']
	compression             = params['compression']
	force                   = params['force']
	override                = params['override']
	init_dir                = params['init_dir']



	chromosomes, bam_rgs = get_chromosomes(alignment_file)


	# print(alignment_file)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)

	# chromosomes = [c for c in chromosomes if c[0] == 'NC_037320.1']

	chrom_depth_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg','chrom'])

	keys = list(chrom_depth_c.keys())
	for key in keys:
		if key[0] in annotation_readgroups:
			chrom_depth_c[key[1]] += chrom_depth_c[key]

		del chrom_depth_c[key]



	subsample = parse_subsample(target_depth, alignment_file, compression, sum(chrom_depth_c.values()), seed=params['seed'])


	if init_dir:
		Path(output_directory, 'align').mkdir(parents=True, exist_ok=True)
		for i, file in enumerate(subsample.files):
			subsample.files[i] = Path(output_directory, 'align', file.name)


	perform_subsample(subsample, force=force)

