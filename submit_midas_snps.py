#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes
import sutilspy
import os

if __name__ == "__main__":
    import argparse
    
    # Parse arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Required arguments
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--samples", help = "File with samples (one per line) to be processed",
                          type = str, required = True)
    required.add_argument("--indir", help ="Directory where sample files are located",
                          type = str, required = True)
    required.add_argument("--outdir", help = "Directory where to save output",
                          type = str, required = True)
    
    # Optional arguments
    parser.add_argument("--sample_col",help = "Column where sample id is located in --samples",
                        default = 1, type = int)
    parser.add_argument("--method", help = "Method to use for submissions", type = str,
                        default = 'qsub', choices = ['qsub','slurm'])
    parser.add_argument("--logdir", help = "If method is cluster-based, where to store the logfiles",
                         type = str, default = "logs")
    parser.add_argument("--submissions_dir", help = "Directory where to store submission dirs",
                        type = str, default = "submissions")
    parser.add_argument("--queue", help = "If method is  'slurm, the partition to use", default = "hbfraser",
                        choices = ['hbfraser','owners'], type = str)
    parser.add_argument("--memory", help = "Amount of memory to request",
                        default = "10G", type = str)
    parser.add_argument("--time", help = "If method is slurm, amount of time to reserve",
                        type = str, default = "4:00:00")
    parser.add_argument("--species_cov", help = "Include species with greated coverage than identified",
                        type = float, default = 3.0)
    parser.add_argument("--mapid", help = "Minimum mapping identity", type = float,
                        default = 94.0)
    parser.add_argument("--mapq", help = "Minumum mappping quality",
                        type = int, default = 20)
    parser.add_argument("--baseq", help = "Discard bases with quality under the specified value",
                        type = int, default = 30)
    parser.add_argument("--readq", help = "Minimum read quality",
                        type = int, default = 30)
    parser.add_argument("--trim", help = "Trim N base-pairs from 3' end of the read",
                        type = int, default = 0)
    parser.add_argument("--discard", help = "Flag to discard discordant read pairs",
                        action = "store_true")
    parser.add_argument("--baq", help = "Flag to enable per-base alignment quality (BAQ)",
                        action = "store_true")
    parser.add_argument("--adjust_m1", help = "Adjust MAPQ",
                        action = "store_true")
    args = parser.parse_args()
    #args.sample_col -= 1
    
    
    # Read samples
    samples = sutilspy.io.return_column(infile = args.samples,
                                        col = args.sample_col,
                                        separator = '\t',
                                        header = False)
    # Prepare directories
    if args.method in ['qsub','slurm']:
        if not os.path.isdir(args.logdir):
            print("Creating directory {}".format(args.logdir))
            os.mkdir(args.logdir)
        if not os.path.isdir(args.submissions_dir):
            print("Creating directory {}".format(args.submissions_dir))
            os.mkdir(args.submissions_dir)
    
    
    #print(samples)
    
    pre_commands = []
    
    # Add module dependencies
    pre_commands.append("module load MIDAS/1.2.1")
    pre_commands.append("echo MIDAS database is $MIDAS_DB")
    bin = "run_midas.py"
    
    for sample in samples:
        print("== Processing sample {}".format(sample))
        sample_file_base = args.indir + "/" + sample
        read1 = sample_file_base + "_read1.fastq.bz2"
        read2 = sample_file_base + "_read2.fastq.bz2"
        
        # We need to check if the read files exist
        if not os.path.isfile(read1):
            raise FileNotFoundError("File {} not found".format(read1))
        if not os.path.isfile(read2):
            raise FileNotFoundError("File {} not found".format(read2))
        
        # We also need to check if the species coverage file exists
        species_file = args.outdir + "/" + sample + "species/species_profile.txt"
        if not os.path.isfile(species_file):
            raise FileNotFoundError("File {} not found".format(species_file))
        
        # Build MIDAS comand
        midas_command = [bin,"snps",args.outdir + "/" + sample,
                         "-1", read1, "-2", read2,"-t","8",
                         "--remove_temp",
                         "--species_cov",args.species_cov,
                         "--mapid", args.mapid,
                         "--mapq", args.mapq,
                         "--baseq", args.baseq,
                         "--readq", args.readq]
        if args.trim > 0:
            midas_command.extend(["--trim",args.trim])
        if args.discard:
            midas_command.append("--discard")
        if args.baq:
            midas_command.append("--baq")
        if args.adjust_mq:
            midas_command.append("--adjust_mq")
        
        midas_command = " ".join(midas_command)
        #print("#######")
        #print(midas_command)
        #print("#######")
        
        commands = pre_commands[:]
        commands.append(midas_command)
        
        job_name = sample + ".midas"     
        logfile = args.logdir + "/midas.species." + sample + ".log"
        errorfile = args.logdir + "/midas.species." + sample + ".err"
   
        submission_file = args.submissions_dir + "/midas.species." + sample + ".bash"
        
        with open(submission_file,'w') as fh:
            if args.method == 'qsub':
                #memory = "16000mb"
                nodes = "nodes=1:ppn=8"
                sutilspy.io.write_qsub_submission(fh = fh, commands = commands,
                                                  name = job_name,
                                                  memory = args.memory,
                                                  logfile = logfile,
                                                  errorfile = errorfile,
                                                  nodes = nodes)
            elif args.method == 'slurm':
                #memory = "16G"
                nodes = "1"
                cpus = "8"
                sutilspy.io.write_slurm_submission(fh = fh,
                                                   commands = commands,
                                                   name = job_name,
                                                   memory = args.memory,
                                                   logfile = logfile,
                                                   errorfile = errorfile,
                                                   queue = args.queue,
                                                   nodes = '1',
                                                   cpus = '8',
                                                   time = args.time)
            else:
                raise ValueError("Invalid method {}".format(args.method))
        fh.close()
        os.chmod(submission_file, 0o744)
        
        # Submit submission file
        if args.method == 'qsub':
            print(submission_file)
            #sutilspy.io.qsub_submissions([submission_file],args.logdir)
        elif args.method == 'slurm':
            print(submission_file)
            #sutilspy.io.sbatch_submissions([submission_file], args.logdir)
        elif args.method == 'bash':
            sutilspy.io.run_command(submission_file)
        else:
            raise ValueError("Incorrect method supplie ({})".format(args.method))
    