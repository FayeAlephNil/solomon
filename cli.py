import click
from solomon.examples import alternating_genus_1
from solomon.modular import mod_rep_orb
from solomon.utils import primes
import networkx as nx

@click.command()
@click.option("-n", default=1, help="Number of punctures")
@click.option("-N", "N", default=None, help="Up to This many Punctures")
@click.option("-p", default=2, help="Prime p")
@click.option("-P", "P", default=None,help="Up to This number of primes")
@click.option("--output", default="out.graph", help="File to Output Graph to")
@click.option("--dir", "dirr", default=None,help="Directory for Graphs")
@click.option("--projective", "proj", default=False, help="Projectivize")
def gen_graph(*,n,N,p,P,output,dirr,proj):
    if (N != None or P != None) and (dirr == None):
        print("Must specify a directory with --dir if using -N or -P")
    elif (N != None or P != None):
        our_range = range(2,int(N)+1) if N != None else [n]
        our_primes = primes(int(P)) if P != None else [p]
        headers = ["n = " + str(x) for x in our_range]
        summary = [[]] * len(our_range)
        for x,y in [(x,y) for x in our_range for y in our_primes]:
            outt = dirr + "/" + str(x) + "-" + str(y) + "-" + output
            num = gen_one_graph(x,y,outt,proj=proj)
            print(str((x,y)) + " has size " + str(num))
    gen_one_graph(n,p,output,proj=proj)

def gen_one_graph(n,p,output,proj=False):
    (S,pure_rep,pure_orb) = alternating_genus_1(n,proj=proj)
    (rep_p, orb_p) = mod_rep_orb(pure_rep,pure_orb, p)
    orb_p.advance_until()
    G = nx.convert_node_labels_to_integers(orb_p.digraph)
    nx.drawing.nx_pydot.write_dot(G, output)
    return len(orb_p.digraph.nodes)

def print_table(data, headers=None, padding=2):
    """Prints a list of lists as a table.

    Args:
        data: The list of lists to print.
        headers: An optional list of headers for the table.
        padding: The amount of padding to add to each cell.
    """

    if not data:
        return

    # Calculate column widths
    column_widths = [max(len(str(item)) for item in col) for col in zip(*(data + ([headers] if headers else [])))]

    # Print header
    if headers:
        print("".join(str(item).ljust(width + padding) for item, width in zip(headers, column_widths)))
        print("-" * (sum(column_widths) + padding * len(column_widths)))

    # Print data
    for row in data:
        print("".join(str(item).ljust(width + padding) for item, width in zip(row, column_widths)))

if __name__ == '__main__':
    gen_graph()
