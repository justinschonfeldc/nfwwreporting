import os
import click
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace

@click.command()
@click.argument('tree_fn')
@click.argument('output_fn')
@click.option('--samples',default="",help="Specify a file which contains the names of the samples, one per line.")
def plot_trees(tree_fn,output_fn,samples):
    """ Plot a tree """

    # Read in the samples
    sampleList = []
    inf = open(samples,"r")
    for line in inf:
        sampleList.append(line.strip())
    inf.close()

    print("Sample List:")
    print(sampleList)

    os.environ['QT_QPA_PLATFORM']='offscreen'
    style = TreeStyle()
    style.mode = "c"
    style.show_leaf_name = True

    # Set Node Styles
    rstyle = NodeStyle()
    rstyle["shape"] = "sphere"
    rstyle["size"] = 15
    rstyle["fgcolor"]="blue"
    nstyle = NodeStyle()
    nstyle["shape"] = "sphere"
    nstyle["size"] = 0
    sstyle = NodeStyle()
    sstyle["shape"] = "sphere"
    sstyle["size"] = 15
    sstyle["fgcolor"]="red"

    #Create the tree object
    ct = Tree(tree_fn)

    for n in ct.traverse():
        if n.name == "Wuhan" or n.name == "Wuhan|402124":
            n.set_style(rstyle)
        elif n.name in sampleList:
            n.set_style(sstyle)
        elif not n.is_leaf():
            n.set_style(nstyle)
    
    #Draw the tree and write to a file
    ct.render(output_fn, tree_style=style)

if __name__ == "__main__":
    plot_trees()

