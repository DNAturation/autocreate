import argparse
import csv
import os
import json


def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-j', '--jsons', help = 'Path to JSON directory')
    parser.add_argument('-o', '--out', help = 'CSV output file path.')
    parser.add_argument('-t', '--test', help = 'Name of the MIST test name')
    
    return parser.parse_args()

def load_json(filepath):

    with open(filepath, 'r') as f:
        data = json.load(f)
    return data

def allele_calls(genes):
    '''Read a MIST JSON and gathers the allele calls from it.
    
    Will state if the gene was found to be missing (0), truncated (-1), or
    a previously unobserved allele (?). The latter case is an indication
    that you should run update_definitions.py
    '''


    calls = {}

    for gene in genes:
        
        if genes[gene]['BlastResults'] is None:
            calls[gene] = '0'
        
        elif genes[gene]['IsContigTruncation']:
            calls[gene] = '-1'
        
        elif genes[gene]['CorrectMarkerMatch'] is False:
            calls[gene] = '?'
            print("WARNING, '?'s are present, did you run update_definitions.py?")
        else:
            calls[gene] = genes[gene]['MarkerCall']

    return calls

def write_csv(results, outpath):
    
    '''Writes results as a CSV file'''

    # enforce alphanumeric ordering
    genome_order = sorted(results.keys())
    gene_order = sorted(results[genome_order[0]].keys())

    header = ['genomes'] + gene_order

    with open(outpath, 'w') as f:
        out = csv.writer(f)
        
        out.writerow(header)

        for genome in genome_order:

            line = [genome]

            for gene in gene_order:
                
                line.append(results[genome][gene])
            
            out.writerow(line)

def main():
    
    args = arguments()

    jsons = [os.path.join(args.jsons, x) for x in os.listdir(args.jsons)]

    results = {}

    for j in jsons:

        data = load_json(j)
        
        strain = data['Results'][0]['Strain']

        genes = data['Results'][0]['TestResults'][args.test]
       
        results[strain] = allele_calls(genes)

    write_csv(results, args.out)

if __name__ == '__main__':
    main()

