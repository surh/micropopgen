#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes

def check_set_of_runs(runs, dir):
    print("\t Checking runs..,")
    for run in runs:
        run_sra = dir + "/" + run + ".sra"
        check = 0
        if os.path.exists(run_sra):
            command = 'vdb-validate ' + run_sra
            check = download_runs.run_command(command)
            #check = subprocess.run('vdb-validate ' + run_sra + " &", shell = True)
            if check.returncode != 0:
                raise IntegrityError("\rRun {} did not pass the validation".format(run))
        else:
            raise MissingFileError("\tRun {} file does not exist in {}".format(run,outdir))
        
    return(check)

def fastq_dump_runs(runs,indir,outdir,keep):
    if not os.path.isdir(indir):
        raise FileNotFoundError("Input directory {} does not exist".format(indir))
    if not os.path.isdir(outdir):
        print("\tCreating output directory {}".format(outdir))
        os.mkdir(outdir)
        
    FILES = [[], []]
    for run in runs:
        run_sra = indir + "/" + run + ".sra"
        
        if os.path.exists(run_sra):
            command = 'fastq-dump -I -O ' + outdir + ' --split-files --bzip2 ' + run_sra
            check = download_runs.run_command(command)
            #check = subprocess.run('fastq-dump -O ' + outdir + ' --split-files ' + run_sra + " &", shell = True)
            if check.returncode != 0:
                raise ProcessError("\tRun {} could not be processed by fastq-dump".format(run))
            else:
                read1 = outdir + "/" + run + "_1.fastq.bz2"                
                FILES[0].append(read1)
                read2 = outdir + "/" + run + "_2.fastq.bz2"                
                FILES[1].append(read2)
        else:
            raise MissingFileError("\tRun {} file does not exist in {}".format(run,outdir))
        
    return(FILES)

def concatenate_run(file_sets,outdir,name_prefix, extension = ".fastq"):
    if not os.path.isdir(outdir):
        print("\tCreating output directory {}".format(outdir))
        os.mkdir(outdir)
    
    i = 1
    FILES = []
    for files in file_sets:
        newfile = outdir + "/" + name_prefix + "_read" + str(i) + extension
        try:
            concatenate_files(files, newfile)
            i += 1
            FILES.append(newfile)
        except (ProcessError):
            raise ProcessError("Could not concatenate files from read {}".format(i))
    
    return(FILES)

def process_ebi_metadata(infile,accession,accession_col = 4):
    accession_col -= 1
    with open(infile,'r') as meta_file:
        reader = csv.reader(meta_file,delimiter = '\t')
        header = next(reader)
        header.append('subject_id')
        print("\tSeaching for accession {} in column {}".format(accession,header[accession_col]))
        #print(header[accession_col])
        nruns = 0
        res = []
        for row in reader:
            if row[accession_col] == accession:
                nruns += 1
                #print(row[22]) 
                
                # Search for subject id
                title = row[22]
                m = re.search('containing sample (\d+) from participant (\d+)', title)
                if m is not None:
                    #print(m.group(1))
                    row.append(m.group(1))
                else:
                    row.append('NA')
                
                res.append(row)
        meta_file.close()
    return(header,res,nruns)

def process_sample(sample,runs,indir,fastqdir,outdir,keep = False):
    
    # Validate files
    try:
        check_set_of_runs(runs,indir)
    except IntegrityError as error:
        print("\tWARNING: Run(s) in sample {} did not pass integrity check. SKIPPING".format(sample))
        raise ProcessError("\tSample didn't pass check")
    except MissingFileError as error:
        print("\tWARNING: Missing file(s) for run(s) in sample {}. SKIPPING".format(sample))
        raise ProcessError("\tSample didn't pass check")    
    
    # Proceed to fastq-dump
    try:
        run_fastq = fastq_dump_runs(runs,indir,fastqdir,keep)
    except FileNotFoundError as error:
        print("\tERROR: Input directory {} does not exist. TERMINATING".format(indir))
        raise FileNotFoundError("Input directory {} does not exist".format(indir))
    except ProcessError as error:
        print("\tWARNING: Run(s) in sample {} could not be processed by fastq-dump. SKIPPING".format(sample))
        raise ProcessError("\tSample could not be processed  with fastq-dump")
    except MissingFileError as error:
        print("\tWARNING:Run(s) file(s) for sample {} missing".format(sample))
        raise ProcessError("Run(s) file(s) for sample {} missing".format(sample))
    
    # Proceed to concatenate
    try:
        concatenated_files = concatenate_run(run_fastq, outdir, sample, ".fastq.bz2")
    except ProcessError as error:
        print("\tWARNING. Could not concatenate files from sample {}. SKIPPING")
        raise ProcessError("Could not concatenate files from sample {}".format(sample))
    
    return(concatenated_files)

def create_single_submission(name, group,submission_dir,outdir,logdir):
    submission_file = submission_dir + "/" + name
    with open(submission_file,'w') as fh:
        fh.write("#!/bin/bash\n")
        fh.write("#PBS -N download." + name + "\n")
        fh.write("#PBS -d " + outdir + "\n")
        fh.write("#PBS -o " + logdir + "/download." + name + ".log\n")
        fh.write("#PBS -e " + logdir + "/download." + name + ".err\n")
        fh.write("#PBS -l mem=1000mb\n")
        
        # Add lines for every run in sample
        ascp_command = 'ascp -i /godot/hmp/aspera/asperaweb_id_dsa.openssh -k 1 -T -l200m'
        sra_prefix = 'anonftp\@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/'
        for run in group:
            run_location = run[0:6] + "/" + run + "/" + run + ".sra"
            command = " ".join([ascp_command, sra_prefix + "/" + run_location, outdir + "\n"])
            fh.write(command)
    fh.close()
    os.chmod(submission_file, 0o744)
    
    return(submission_file)

def aspera_download(groups,outdir):
    
    ascp_command = 'ascp -i /godot/hmp/aspera/asperaweb_id_dsa.openssh -k 1 -T -l200m'
    sra_prefix = 'anonftp\@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/'
    
    FAILED = []
    for name, runs in groups.items():
        for run in runs:
            run_location = run[0:6] + "/" + run + "/" + run + ".sra"
            command = " ".join([ascp_command, sra_prefix + "/" + run_location, outdir + "\n"])
            print(command)
            try:
                check = run_command(command)
                if check.returncode != 0:
                    raise CalledProcessError("Aspera download failed for run {}".format(run))
            except (CalledProcessError):
                print("\tWARNING: Failed downloading run {}".format(run))
                FAILED.append([run, 'dummy'])
    #FAILED = [[1,'a'],[2,'a']]
    return(FAILED)

def create_submission_sets(runs_per_sample, split_by, ngroups):
    print("\n=============================================")
    GROUPS = dict()
    if split_by == 'sample':
        print("== Entering splity by sample")
        GROUPS = runs_per_sample
    elif split_by == 'groups':
        samples = runs_per_sample.keys()
        total_samples = len(samples)
        samples_per_submission = ceil(total_samples / ngroups)
        
        print("== Splitting {} samples into {} submissions".format(total_samples,ngroups))
        i = 0
        group_i = 0
        for sample, runs in runs_per_sample.items():
            if (i % samples_per_submission) == 0:
                #GROUPS.append([])
                id = 'group' + str(group_i)
                GROUPS[id] = []
                group_i += 1
            GROUPS[id].extend(runs)
            i += 1
    else:
        raise ValueError("Unrecognized split_by value")
    
    print("\tSplitted runs")
    print("=============================================")
    return(GROUPS)

def create_submission_files(groups, outdir, logdir):
    print("\n=============================================")
    # Create output directory
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        print("== Outdir ({}) already exists. Using it.".format(outdir))
    
    submission_dir = tempfile.mkdtemp(suffix = None, prefix = 'submissions',
                                      dir = outdir)
    
    SUBMISSIONS = []
    #i = 1
    for name, runs in groups:
        #name = str(i)
        #name = "group" + name + ".bash"
        #print(name)
        newfile = create_single_submission(name,runs,
                                           submission_dir,
                                           outdir,logdir)
        SUBMISSIONS.append(newfile)
        #print(newfile)
        #print(group)
        #i += 1
    
    print("\tCreated submision files")
    print("=============================================")
    return(SUBMISSIONS)
