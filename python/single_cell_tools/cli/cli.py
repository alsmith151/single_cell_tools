import click
import json

import single_cell_tools


@click.command()
@click.option("-i",'--bam', required=True, help='BAM file')
@click.option("-b", '--barcodes', required=True, help='JSON file with barcodes')
@click.option("-o", '--output', required=True, help='Output directory')
@click.option("-p", '--threads', default=1, help='Number of threads')
def split_bam(bam: str, barcodes: str, output: str, threads: int):

    with open(barcodes, 'r') as f:
        barcodes = json.load(f)
    
    single_cell_tools.libsct.split_bam_file(bam, barcodes, threads, output)


split_bam()