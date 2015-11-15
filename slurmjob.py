import sys
from subprocess import Popen, PIPE

class SLURMJob:

    """
    Job for working with SLURM at Princeton
    """
    def __init__(self, name='', queue='default', nodes=1, ppn=1,
                 walltime = '01:00:00', mail='a', addr=None,
                 cwd = True, command = None, array=None, depends=None, workdir=None):
        self.name = name
        self.queue = queue
        self.nodes = nodes
        self.ppn = ppn
        self.walltime = walltime
        self.mail = mail
        self.array = array
        self.depends = depends
        self.addr = addr
        self.cwd = cwd
        if command is None:
            sys.stderr.write('A command is REQUIRED')
            return None
        self.command = command
        if workdir is None:
            sys.stderr.write('A command is REQUIRED')
            return None
        self.workdir = workdir

    def set_command(self, command):
        self.command = command

    def set_name_command(self, name, command):
        self.name = name
        self.command = command

    #depends should be a list of jobs, uses afterok rule by default
    def set_depends(self, depends):
        self.depends = depends

    def submit(self, filename):
        self.write(filename)
        stdout = Popen('sbatch ' + filename, shell=True, stdout=PIPE).stdout.read()
        jobid = stdout.strip().split()[-1]
        return jobid

    def write(self, filename):
        ofile = open(filename, 'w')
        ofile.write('#!/bin/bash\n')
        ofile.write('#SBATCH --ntasks=1\n')
        ofile.write('#SBATCH --ntasks-per-node=1\n')
        ofile.write('#SBATCH --nodes=' + str(self.nodes) + '\n')
        ofile.write('#SBATCH --cpus-per-task=' + str(self.ppn) + '\n')
        ofile.write('#SBATCH --time=' + str(self.walltime) + '\n')
        #ofile.write('#SBATCH -M ' + str(self.addr) + '\n')
        #ofile.write('#SBATCH -m=' + str(self.mail) + '\n')
        ofile.write('#SBATCH --output=' + self.workdir + '/' + self.name + '.out\n')
        ofile.write('#SBATCH --error=' + self.workdir + '/' + self.name + '.err\n')
        ofile.write('#SBATCH --mem-per-cpu=8192\n')
        if self.depends is not None:
            ofile.write('#SBATCH -d afterok:' + ':'.join(self.depends) + '\n')
        if self.cwd:
            ofile.write('cd $GRID_HOME\n')
        ofile.write(str(self.command) + '\n')
        ofile.write('exit 0\n')
        ofile.close()

if __name__ == "__main__":
    job = SLURMJob(name='test', command='ls')
    job.submit('test.sh')
