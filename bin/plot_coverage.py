import subprocess
import matplotlib.pyplot as plt
import numpy as np
from os import path
import click
import pysam


@click.command()
@click.argument('bam_fn')
@click.option('--logscale',is_flag=True)
def plot_depth(bam_fn,logscale):
    """ Plot coverage """
    # Check that the bam file exists.
    if not (path.exists(bam_fn) and path.isfile(bam_fn)) :
        print(f"Error! {bam_fn} is an invalid path or not a file!")
        return

    # output = subprocess.check_output(["samtools","depth",bam_fn])

    pos_depth = []
    pos_depth.append(list())
    pos_depth.append(list())
    pos_depth.append(list())

    samfile = pysam.AlignmentFile(bam_fn,"rb")
    for pileupcolumn in samfile.pileup("MN908947.3"):
        pos_depth[0].append(pileupcolumn.pos)
        pos_depth[1].append(pileupcolumn.n)

    fig, ax = plt.subplots()
    fig.set_figwidth(20)
    ax.bar(pos_depth[0],pos_depth[1],width=1.0)
    if logscale:
        ax.set_yscale('log')
    plt.xlabel('Position')
    plt.ylabel('Depth')
    plt.savefig("coverage.png")


if __name__ == "__main__":
    plot_depth()
