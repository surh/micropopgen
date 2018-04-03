#!/usr/bin/env python
# Copyright (C) 2018 Sur Herrera Paredes
# This file is based on the MarkerScanner.pl script of the AMPHORA
# pipeline by Martin Wu.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import fyrd
import argparse
import os
from Bio import SearchIO, SeqIO
import time


def fasta_seq_lenghts(fasta_file, split=False):
    """Read sequences in fasta file and obtain sequence lengths"""

    fasta = SeqIO.parse(fasta_file, 'fasta')
    Sequences = dict()
    for s in fasta:
        if split:
            key = s.description.split()[1]
        else:
            key = s.id
        Sequences[key] = [str(s.seq), len(s.seq)]

    return Sequences


def get_hmm_hits(hmmfile, query_fasta, dbfile, name, outdir='./'):
    """Read HMMER files and get hits"""

    # Read query fasta
    queries = fasta_seq_lenghts(query_fasta)
    # db = fasta_seq_lenghts(db_fasta, split=True)
    db = read_marker_list(dbfile)

    # Find hits and save tophit for every query
    hmmsearch = SearchIO.parse(hmmfile, 'hmmer3-text')
    print("==Read==")
    hmm_hits = {k: [] for k in db}
    for query in hmmsearch:
        for hit in query:
            hit_span, query_span = hit_and_query_span(hit)
            query_cov = query_span / queries[query.id][1]
            hit_cov = hit_span / db[hit.id]
            # print("\t{}\t{}\t{}".format(query.id, query_cov, hit_cov))
            if query_cov > 0.7 and hit_cov > 0.7:
                hmm_hits[hit.id].append(query.id)
                break

    # print(hmm_hits)
    # Write file per marker
    res = dict()
    for marker in hmm_hits:
        # print(marker)
        # marker_file = strip_right(hmmfile, '.hmms')
        # marker_file = marker_file + '.' + marker + '.faa'
        marker_file = outdir + '/' + name + '.faa'
        with open(marker_file, mode='w') as out:
            i = 0
            for hit in hmm_hits[marker]:
                out.write(">" + name + "_" + str(i) + "\n")
                out.write(queries[hit][0] + "\n")
                i = i + 1
            res[marker] = i

    return res


def hit_and_query_span(hit):
    """Calculate total hit and query span"""

    hit_span = 0
    query_span = 0
    for hsp in hit.fragments:
        query_span = query_span + hsp.query_span
        hit_span = hit_span + hsp.hit_span

    return(hit_span, query_span)


def hmmscan_file(filename, db, job_name=None, outpath='./logs/',
                 scriptpath='./scripts/', partition='',
                 time='00:30:00', hmmscan='hmmscan',
                 indir='', outdir='', fasta_suffix='.faa',
                 out_suffix='.hmms', memory='300mb',
                 maxjobs=1000):
    """Use fyrd to run hmmscan on a given file"""

    # Get basename
    basename = strip_right(filename, fasta_suffix)

    # Build hmmscan filenames
    infile = '/'.join([indir, filename])
    outfile = '/'.join([outdir, basename])
    outfile = ''.join([outfile, out_suffix])

    # Build hmmscan command
    command = ' '.join([hmmscan,
                        "-Z", str(5000),
                        "-E", str(1e-3),
                        db,
                        infile,
                        ">", outfile])

    # Build fyrd filenames
    if job_name is None:
        job_name = '.'.join(['hmmscan', basename])
    print(job_name)

    print("\tCreating fyrd.Job")
    fyrd_job = fyrd.Job(command,
                        runpath=os.getcwd(), outpath=outpath,
                        scriptpath=scriptpath,
                        clean_files=False, clean_outputs=False,
                        mem=memory, name=job_name,
                        outfile=job_name + ".log",
                        errfile=job_name + ".err",
                        partition=args.queue,
                        nodes=1, cores=1, time=time)

    # Submit joobs
    print("\tSubmitting job")
    fyrd_job.submit(max_jobs=maxjobs)

    return job_name, outfile, fyrd_job


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Uses HMMER hmmscan search all protein fasta "
                          "sequences in all fasta files in a directory, "
                          "against a database of profiles o marker genes")

    # Define required arguments
    required.add_argument("--indir", help=("Directory where fasta files are"),
                          required=True, type=str)
    required.add_argument("--outdir", help=("Directory where to write output "
                                            "files"),
                          required=True, type=str)
    required.add_argument("--db", help=("Database of hmm profiles for maker "
                                        "sequences"),
                          required=True, type=str)

    # Define other arguments
    # parser.add_argument("--markers_pep", help=("Location of fasta file of "
    #                                            "marker genes. If left empty "
    #                                            "it will default to "
    #                                            "'markers.fas' in the same "
    #                                            "directory as --db."),
    #                     type=str, default='')
    parser.add_argument("--fasta_suffix", help=("Suffix of fasta files in "
                                                "indir"),
                        type=str, default='.faa')
    parser.add_argument("--out_suffix", help=("Suffix of output files"),
                        type=str, default='.hmms')
    parser.add_argument("--hmmscan", help=("Binary executable of hmmscan "
                                           "tool of the HMMER package"),
                        type=str, default='hmmscan')
    parser.add_argument("--logs", help=("Directory where to write logfiles "
                                        "from fyrd."),
                        type=str, default='logs/')
    parser.add_argument("--scripts", help=("Directory where to write scripts "
                                           "from fyrd."),
                        type=str, default='scripts/')
    parser.add_argument("--maxjobs", help=("Maximum number of fyrd jobs "
                                           "running simultaneously"),
                        type=int, default='500')
    parser.add_argument("--queue", help=("Queue (partition) to use for "
                                         "submitting jobs"),
                        type=str)
    parser.add_argument("--memory", help=("Amunt of memory to reserve per "
                                          "job"),
                        type=str, default="300mb")
    parser.add_argument("--time", help=("Amunt of time to reserve per "
                                        "job"),
                        type=str, default="1:00:00")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    # Check if exectuable exists
    if which(args.hmmscan) is None:
        raise FileNotFoundError("Executable for hmmscan not found")
    else:
        args.hmmscan = which(args.hmmscan)

    # Check if markers.pep is passed and exists, if not
    # set default and check if exists
    # if args.markers_pep == '':
    #     args.markers_pep = "/".join([os.path.dirname(args.db), 'markers.fas'])
    # if not os.path.isfile(args.markers_pep):
    #     raise FileExistsError("Markers fasta does not exist")

    return args


def read_marker_list(infile):
    """Read hmm profiles file and return length of profiles"""

    hmm_lenghts = dict()
    with open(infile) as fh:
        for l in fh:
            if l.startswith('NAME'):
                name = l.split()[1]
            elif l.startswith('LENG'):
                hmm_lenghts[name] = int(l.split()[1])

    return hmm_lenghts


def strip_right(text, suffix):
    # tip from http://stackoverflow.com/questions/1038824
    # MIT License
    if not text.endswith(suffix):
        return text
    # else
    return text[:len(text)-len(suffix)]


def which(program):
    """Check if executable exists. Returns path of executable."""

    # From:
    # https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    # Under MIT license

    def is_exe(fpath):
        """Check if path is executable"""

        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def submit_hmmscan_file(f, args, name=None):
    """Create a job for rinnung hmmscan and submit it"""

    # Creat subdirectory for hmmscan files
    hmmscandir = args.outdir + '/' + 'hmms/'
    if not os.path.isdir(hmmscandir):
        os.mkdir(hmmscandir)

    job_name, hmmfile, job = hmmscan_file(filename=f,
                                          db=args.db,
                                          outpath=args.logs,
                                          scriptpath=args.scripts,
                                          partition=args.queue,
                                          time=args.time,
                                          hmmscan=args.hmmscan,
                                          indir=args.indir,
                                          outdir=hmmscandir,
                                          fasta_suffix=args.fasta_suffix,
                                          out_suffix=args.out_suffix,
                                          memory=args.memory,
                                          maxjobs=args.maxjobs,
                                          job_name=name)

    return hmmfile, job


def submit_get_hmm_hits(hmmfile, job, fasta_file, args):
    """Create a job for rinnung hmmscan and submit it"""

    # Creat subdirectory for fasta files
    markersdir = args.outdir + '/' + 'markers/'
    if not os.path.isdir(markersdir):
        os.mkdir(markersdir)

    # Get strain name
    strain_name = os.path.basename(hmmfile)
    strain_name = strip_right(strain_name, args.out_suffix)

    # Fasta File
    fasta_file = args.indir + '/' + fasta_file
    if not os.path.isfile(fasta_file):
        raise FileNotFoundError("Fasta file not found")

    res = get_hmm_hits(hmmfile, query_fasta=fasta_file,
                       dbfile=args.db, name=strain_name,
                       outdir=markersdir)
    Res = {strain_name: res}
    # job = fyrd.Job(get_hmm_hits, f, {'query_fasta': o[1]},
    #                depends=o[0], runpath=os.getcwd(),
    #                outpath=args.logs,
    #                scriptpath=args.scripts)
    # job.submit(max_jobs=args.maxjobs)

    return Res


if __name__ == "__main__":
    args = process_arguments()

    # Get list of fasta files from indir
    fasta_files = os.listdir(args.indir)
    fasta_files = list(filter(lambda f: f.endswith(args.fasta_suffix),
                              fasta_files))
    # print(fasta_files)

    # Make directories
    if os.path.isdir(args.logs):
        raise FileExistsError("Directory for fyrd logs ({}) "
                              "already exists".format([args.logs]))
    else:
        os.mkdir(args.logs)
    if os.path.isdir(args.scripts):
        raise FileExistsError("Directory for fyrd scripts ({}) "
                              "already exists".format([args.scripts]))
    else:
        os.mkdir(args.scripts)
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    # Submit hmmscan jobs
    print("===hmmscan===")
    hmm_files = dict()
    for f in fasta_files:
        hmmfile, job = submit_hmmscan_file(f=f, name=None, args=args)
        hmm_files[hmmfile] = [job, f]
        print(f)

    # Submit hits_job
    print("===hits===")
    time.sleep(10)
    for f, o in hmm_files.items():
        print(f)
        submit_get_hmm_hits(hmmfile=f, job=o[0], fasta_file=o[1], args=args)
