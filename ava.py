import csv
import argparse
from Bio import SeqIO

def get_lengths(filepath):
    """Reads the file located at the --seq
    argument, and returns a dictionary containing
    their lengths.
    """

    lengths = {}

    with open(filepath, 'r') as f:
        for s in SeqIO.parse(f, 'fasta'):
            lengths[s.id] = float(len(s.seq))

    return lengths

def file_len(resultfile):

    with open(resultfile, 'r') as f:
        i = f.readlines()
        return len(i)

def find_homologues(resultfile, lengths, pid_thresh, length_thresh):
    """Finds homologous genes within the input --seq file.

    If both percent identity and percent length are both over the
    threshold, the longer variant is kept.
    """

    homologues = {}

    ignore = set([])
    
    total_lines = file_len(resultfile)

    processed_lines = 0 

    with open(resultfile, 'r') as f:
        rdr = csv.reader(f, delimiter = ",")
        reported = False
        for line in rdr:
            
            if round(100.0 * processed_lines / total_lines, 0) in range(0,100,5):
                if reported is False:
                    print("{}% complete".format(round(100 * processed_lines / total_lines), 0)      )
                    reported = True
            else:
                reported = False
            processed_lines += 1

            name1 = line[0]
            name2 = line[1]
            pid   = line[2]
            aln   = float(line[3])

            if float(pid) > pid_thresh \
                and aln / lengths[name2] > (length_thresh / 100) \
                and name1 not in ignore \
                and name2 not in ignore:

                if lengths[name1] >= lengths[name2]:
                    
                    try:
                        homologues[name1].append(name2)
                    except KeyError:
                        homologues[name1] = [name2]

                    ignore.add(name2)
                
                else:

                    if name1 in homologues:
                        try:
                            homologues[name2].extend(homologues[name1])
                            homologues[name2].append(name1)

                        except KeyError:

                            homologues[name2] = homologues[name1] + [name1]
                    
                        ignore.add(name1)
                        del homologues[name1]
                    
                    else:
                        homologues[name2] = [name1]
    return homologues

def extract_non_redundant(infilepath, homologues, outfilepath):
    """ Writes out the longest variant of each homologue found
    by find_homologues().

    Results are written to a fasta file.
    """
    print("{} non-redundant genes found".format(len(homologues)))

    print("Extracting non-redundant")

    nr = []
    with open(infilepath, 'r') as f:

        for s in SeqIO.parse(f,'fasta'):
            if s.id in homologues.keys():
                nr.append(s)

    with open(outfilepath, 'w') as o:
        SeqIO.write(nr, o, 'fasta')

def arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument('--seq', required = True, help = "Path to source FFN  file.")
    parser.add_argument('--result', required = True, help = "Path to BLAST output file.")
    parser.add_argument('--out', required = True, help = "Output path.")
    parser.add_argument('--identity', default = 90.0, type = float, help = "Percent ID threshold.")
    parser.add_argument('--length', default = 50.0, type = float, help = "Percent length threshold.")
    return parser.parse_args()

def main():

    args = arguments()
    
    lengths = get_lengths(args.seq)

    homologues = find_homologues(args.result, lengths, args.identity, args.length)

    extract_non_redundant(args.seq, homologues, args.out)

if __name__ == '__main__':
    main()
