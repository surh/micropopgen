def call_cmd(cmd):
    """Call aspera to download runs from SRA"""
    
    print("\n=============================================")
    command = cmd
    print(command)
    try:
        check = run_command(command)
        if check.returncode != 0:
            print("Raising error")
            raise CalledProcessError("Aspera download failed for run {}".format(cmd))
    except CalledProcessError:
        print("\tWARNING: Failed downloading run {}".format(run))
        FAILED.append([1, 'dummy'])
    print("=============================================\n\n")
    return(FAILED)

def run_command(command):
    status = 0;
    print("Executing:\n>{}".format(command))
    status = subprocess.run(command, shell = True)
    print("Status={}\n\n".format(status));

    return(status)

if __name__ == "__main__":
    cmd = "ls -la"
    fail = call_cmd(cmd)
    print(fail)
    
