import subprocess
from subprocess import CalledProcessError

def call_cmd(cmd):
    """Call aspera to download runs from SRA"""
    
    print("\n=============================================")
    command = cmd
    print(command)
    FAILED = []
    try:
        check = run_command(command)
        if check.returncode != 0:
            print("Raising error")
            raise CalledProcessError(check.returncode,command,check.stdout,check.stderr)
    except CalledProcessError:
        print("\tWARNING: Failed command {}".format(command))
        FAILED.append([1, 'dummy'])
    print("=============================================\n\n")
    return(FAILED)

def run_command(command):
    status = 0;
    print("Executing:\n>{}".format(command))
    try:
        status = subprocess.run(command, shell = True)
        print("Status={}\n\n".format(status));
    except:
        print("\tWARNING")
        
    print("run_command done")
    return(status)

if __name__ == "__main__":
    cmd = "ls -la"
    cmd = 'a'
    fail = call_cmd(cmd)
    print(fail)
    
