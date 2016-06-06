import subprocess
import os
import argparse
import re
import glob
import multiprocessing
import string
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord



### TODO Add filter to get genes for a cgMLST, or accessory scheme
### TODO Change names/args to reflect that the scheme initially created
        ### is NOT a cgmlst scheme nor a pangenome scheme
        ### nor is it a pangenome scheme
### TODO Fix divide_schemes.R so that threshold for absence is no longer hardcoded



# SRC_DIR=os.path.pardir() #goes up one level
SCRIPT_DIR=os.getcwd()
T="$(date +%s)"

def arguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('--workdir', required=True, help='Working directory for this script.')
    parser.add_argument('--reference', required=True, help='Fasta filename (not path) within --genomes path to be used as reference.')
    parser.add_argument('--genomes', required=True, help='Path to directory containing genomes as FASTAS.')
    parser.add_argument('--prokkaout', default='prokka_out/')
    return parser.parse_args()




def fasta_rename(file, dir):
    with open(file, 'r') as f:
        FASTANAME = f.readlines()[0]
        FASTANAME = FASTANAME.split()[0].replace('>', '')
        FASTANAME = FASTANAME+'.fasta'
        newfile = re.sub(r'(^>.*)', r'>1\1', f.read())
    with open(os.path.join(dir, FASTANAME), 'w') as g:
        g.write(newfile)



### Set up directories ###
def mkdir(WORKDIR):
    direc = ['alleles/', 'blast_out/', 'jsons/', 'msa/', 'temp/']
    if not os.access(WORKDIR, os.F_OK):
        os.mkdir(WORKDIR)
        for item in direc:
            os.mkdir(os.path.join(WORKDIR, item))
    else:
        for item in direc:
            if not os.access(os.path.join(WORKDIR, item), os.F_OK):
                os.mkdir(os.path.join(WORKDIR, item))

### Get non-redundant gene set ###
def prefixget(REFERENCE):
    PROKKA_PREFIX=os.path.splitext(REFERENCE)
    return PROKKA_PREFIX[0]

def run_prokka(prokkaout, prefix, genomes, reference, workdir):

    print ("\nRunning Prokka\n")
    prokargs= ('prokka',
               '--outdir', os.path.join(workdir, prokkaout),
               '--prefix', prefix,
               '--locustag', prefix,
               '--cpus', str(0),
               os.path.join(genomes, reference))
    subprocess.call(prokargs)

### all-vs-all BLAST search to filter homologues
def run_blastdb(prokkaout, prefix, workdir):
    print("\nStarting all-vs-all BLAST search of CDS\n")
    dbblastargs = ('makeblastdb',
                 '-in', os.path.join(workdir, prokkaout, prefix + '.ffn'),
                 '-dbtype', 'nucl',
                 '-out', '{}temp/{}_db'.format(workdir, prefix))
    subprocess.call(dbblastargs)

def run_blastn(prokkaout, prefix, workdir):
    if not os.access('{}blast_out/'.format(workdir), os.F_OK):
        os.mkdir('{}blast_out/'.format(workdir))
    blastargs = ('blastn',
                 '-query', os.path.join(workdir, prokkaout, prefix)+'.ffn',
                 '-db', '{}temp/{}_db'.format(workdir, prefix),
                 '-num_threads', str(multiprocessing.cpu_count()),
                 '-outfmt', str(10),
                 '-out', '{}blast_out/all_vs_all.csv'.format(workdir))
    subprocess.call(blastargs)


# AVA filters the all-vs-all search for homologues
    # Thresholds are defaulted as 90% PID and 50% length
    # Only the longest variant is kept
def run_ava(scriptdir, prokkaout, prefix, workdir):
    print("\nFiltering out homologues\n")
    if not os.access('{}blast_out/'.format(workdir), os.F_OK):
        os.mkdir('{}blast_out/'.format(workdir))
    avargs=('python3', '{}/ava.py'.format(scriptdir),
            '--seq', os.path.join(workdir, prokkaout, prefix)+'.ffn',
            '--result', '{}blast_out/all_vs_all.csv'.format(workdir),
            '--out', '{}blast_out/non_redundant.fasta'.format(workdir))
    subprocess.call(avargs)



### Create .markers file for MIST ###
def markers(prokkaout, prefix, workdir):
    print("\nSplitting to discrete fastas.\n")
    duplist = ['']+list(string.ascii_lowercase)
    if not os.listdir('{}alleles/'.format(workdir)):
        with open(workdir+prokkaout+prefix+'.ffn', 'r') as f:
            for rec in SeqIO.parse(f, 'fasta'):
                header = rec.id
                sequence = rec.seq
                for dup in duplist:
                    if os.path.isfile('{}alleles/'.format(workdir)+header+'z.fasta'):
                        print('Error, too many duplicate genes')
                        subprocess.call(exit(1))
                    elif os.path.isfile('{}alleles/'.format(workdir)+header+dup+'.fasta'):
                        continue

                    else:
                        with open('{}alleles/'.format(workdir)+header+dup+'.fasta', 'w+') as g:
                            g.seek(0)
                            dat = g.readlines()
                            g.seek(0)
                            num = len(dat)/2
                            line = '>{}\n'.format(int(num)+1)
                            g.writelines(line)
                            g.write(str(sequence).strip())
                        break
        # for item in d:
        #     with open('{}alleles/'.format(workdir)+item+'.fasta', 'w') as f:
        #         record = SeqRecord(d[item], id=1)
        #         SeqIO.write(record, f, 'fasta')


        # markargs=('csplit', '--quiet',
        #           '--prefix', '{}alleles/'.format(workdir),
        #           '-z', os.path.join(workdir, prokkaout, prefix)+'.ffn',
        #           '/>/', '{*}')
        # subprocess.call(markargs)
        # for afile in os.listdir('{}alleles/'.format(workdir)):
        #     with open ('{}alleles/'.format(workdir) + afile, 'r') as f:
        #         lines = f.readlines()
        #         fullname = lines[0]
        #         name = fullname.split()[0][1:]
        #
        #
        #         for line in lines:
        #             with open ('{}alleles/'.format(workdir)+name+'.fasta', 'a+') as g:
        #                 g.seek(0)
        #                 dat = g.readlines()
        #                 g.seek(0)
        #                 if '>' in line:
        #                     num = len(dat)/2
        #                     if num == 0:
        #                         line = '>{}\n'.format(int(num)+1)
        #                         g.writelines(line)
        #                     else:
        #                         line = '\n>{}\n'.format(int(num)+1)
        #                         g.writelines(line)
        #                 else:
        #                     g.writelines(line.strip())
        #
        #     os.remove('{}alleles/'.format(workdir) + afile)


def renamer():
    ls = glob.glob('alleles/*')
    for i in ls:
        fasta_rename(i, 'alleles/')


def build(scriptdir, workdir):
    print("\nBuilding reference genome .markers file\n")
    bargs = ('python3', os.path.join(workdir, scriptdir, 'marker_maker.py'),
             '--fastas', '{}alleles/'.format(workdir),
             '--out', '{}wgmlst.markers'.format(workdir),
             '--test', 'wgmlst')
    subprocess.call(bargs)

### run MIST ###
def run_mist(genomes, workdir):
    os.chdir(workdir)
    pool = multiprocessing.Pool(int(multiprocessing.cpu_count()/2))
    print("\nRunning MIST in parallel.\n")
    files = glob.glob(genomes+'*.fasta')
    for file in files:
# A hack that will find the shortest path to MIST
    # The notion is that it won't accidentally find
    # debugging binaries buried in the project directory
        mistargs = ('/usr/local/bin/MIST/MIST.exe',
                 '-t', 'wgmlst.markers',
                 '-T', 'temp/',
                 '-a', 'alleles/',
                 '-b', '-j', 'jsons/'+ str(os.path.splitext(os.path.basename(file))[0]) +'.json',
                 file)
        pool.apply_async(subprocess.call, args=(mistargs,))
    pool.close()
    pool.join()

### Update allele definitions ###
def update(scriptdir, workdir):
    print("\nUpdating allele definitions.\n")
    # os.chdir(workdir)
    upargs = ('python3', os.path.join(scriptdir, 'update_definitions.py'),
              '--alleles', 'alleles/',
              '--jsons', 'jsons/',
              '--test', 'wgmlst')
    subprocess.call(upargs)

### Align genes with clustalo ###
def align(workdir):
    print("\nAligning genes.\n")
    # os.chdir(workdir)
    pool = multiprocessing.Pool(int(multiprocessing.cpu_count()/2))
    pathlist = glob.glob('alleles/*.fasta')
    for path in pathlist:
        alargs = ('clustalo',
                  '-i', path,
                  '-o', 'msa/{}'.format(os.path.basename(path)))
        pool.apply_async(subprocess.call, args=(alargs,))
    pool.close()
    pool.join()

### Divide Reference-based calls into core, genome, accessory schemes ###
def divvy(scriptdir, workdir, prefix):
    print("\nParsing JSONs.\n")
    # os.chdir(workdir)
    aargs=('python3', os.path.join(scriptdir, 'json2csv.py'),
             '--jsons', 'jsons/',
             '--test', 'wgmlst',
             '--out', prefix+'_calls.csv')
    subprocess.call(aargs)
    print("\nDividing markers into core and accessory schemes.\n")
    bargs=('Rscript', os.path.join(scriptdir, 'divide_schemes.R'),
           prefix+'_calls.csv', 'wgmlst.markers')
    subprocess.call(bargs)
    print("\nScript complete at `date`\n")


def main():
    args = arguments()
    if not os.access(os.path.join(args.genomes, args.reference), os.F_OK):
        print (os.path.join(args.genomes, args.reference)+'does not exist.')
        subprocess.call('exit')
    mkdir(args.workdir)
    prefix = prefixget(args.reference)
    run_prokka(args.prokkaout, prefix, args.genomes, args.reference, args.workdir)
    run_blastdb(args.prokkaout, prefix, args.workdir)
    run_blastn(args.prokkaout, prefix, args.workdir)
    run_ava(SCRIPT_DIR, args.prokkaout, prefix, args.workdir)
    markers(args.prokkaout, prefix, args.workdir)
    renamer()
    build(SCRIPT_DIR, args.workdir)
    run_mist(args.genomes, args.workdir)
    update(SCRIPT_DIR, args.workdir)
    align(args.workdir)
    divvy(SCRIPT_DIR, args.workdir, prefix)

if __name__ == '__main__':
    main()


T="$(($(date +%s)-T))"
print("Total run time: %02d:%02d:%02d:%02d\n" "$((T/86400))" "$((T/3600%24))" "$((T/60%60))" "$((T%60))")

