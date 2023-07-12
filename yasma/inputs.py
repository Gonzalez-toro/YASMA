


import sys

import click
from click_option_group import optgroup

from pathlib import Path
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint


from .generics import *
from .cli import cli

from shutil import copyfile






@cli.command(group='Utilities', help_priority=3)



@optgroup.group('\n  Required',
				help='')


@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")




@optgroup.group('\n  Inputs which may be logged',
				help='')

@optgroup.option("-a", "--alignment_file", 
	required=False, 
	default=None,
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')

@optgroup.option('-r', '--annotation_readgroups', 
	required=False,
	default=None,
	multiple=True,
	# type=list,
	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")



@optgroup.option("-g", "--genome_file", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	default=None,
	type=click.Path(exists=True),
	help='Genome or assembly which was used for the original alignment.')


@optgroup.option("-j", "--jbrowse_directory", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	default=None,
	type=click.Path(exists=True),
	help='A path to a working directory for a jbrowse2 instance.')



@optgroup.option("--gene_annotation_file", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	default=None,
	type=click.Path(exists=True),
	help='Annotation file for genes (gff3) matching the included genome.')


@optgroup.group('\n  Options',
				help='')

@optgroup.option("--overwrite", 
	is_flag=True,
	default=False,
	help='Allowing the overwrite of inputs already logged.')


def inputs(**params):
	'''A tool to log inputs, which will be referenced by later tools.'''


	rc = requirementClass()
	rc.add_samtools()
	rc.check()

	ic = inputClass(params)
	ic.check_chromosomes()

	print()
	print(color.BOLD + "Inputs:" + color.END)
	pprint(ic.inputs)





























