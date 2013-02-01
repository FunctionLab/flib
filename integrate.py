from optparse import OptionParser
from subprocess import Popen, PIPE
import re
import os
import sys

binaries = { 'Counter':'Counter', 'NetworkCombiner':'NetworkCombiner', 'DChecker':'DChecker', 'Dat2Dab':'Dat2Dab' }

job_queue_id = re.compile("Your job (?P<job_id>\d+) .*")

usage = "usage: %prog [options]"
parser = OptionParser(usage, version="%prog dev-0.0.1")
parser.add_option("-S", "--split-gold-standard",
                        dest="split",
                        action="store_true",
                        default=False,
                        help="split gold standard? if not passed, assumed gold standard is already split by previous run or k=1")
parser.add_option("-k", "--k-intervals",
                        dest="intervals",
                        type="int",
                        default=1,
                        help="Number of cross validation intervals? 1 = no cross validation.")
parser.add_option("-a", "--answers-file",
                        dest="answers",
                        help="Gold standard relationships file",
                        metavar="FILE")
parser.add_option("-b", "--sleipnir-binaries-dir",
                        dest="sleipnir",
                        help="Directory containing sleipnir binaries. If not passed, binaries used must be available in the path.",
                        metavar="FILE")
parser.add_option("-m", "--memory-per-integration",
                        dest="gb_per",
                        type="int",
                        default=3,
                        help="Integer memory requirement per context integration (in GB).")
"""
parser.add_option("-l", "--holdout-file",
                        dest="holdout",
                        help="Holdout file (used for dchecking)",
                        metavar="FILE")
parser.add_option("-z", "--zeros-file",
                        dest="zeros",
                        help="Zeros file",
                        metavar="FILE")
parser.add_option("-d", "--data-directory",
                        dest="datadir",
                        help="Directory where data .qdabs are located.",
                        metavar="string")
parser.add_option("-e", "--genes-file",
                        dest="gfile",
                        help="Tab-delimited text file containing two " \
                        "columns, the first a one-based integer index and " \
                        "the second the unique identifier of each gene " \
                        "in the genome.",
                        metavar="string")
parser.add_option("-c", "--contexts-directory",
                        dest="contdir",
                        help="Directory where contexts files are located",
                        metavar="string")
parser.add_option("-t", "--threads-per",
                        dest="job_threads",
                        help="Number of threads per job",
                        type="int",
                        default=1,
                        metavar="int")
parser.add_option("-w", "--working-directory",
                        dest="workdir",
                        help="Directory for integration",
                        type="string",
                        default=".",
                        metavar="DIRECTORY")
parser.add_option("-r", "--alphas-file",
                        dest="alphas",
                        help="Alphas file for bayesian regularization",
                        metavar="FILE")
parser.add_option("-p", "--pseudo-count",
                        dest="pseudo",
                        help="Pseudo count to be used with alphas file " \
                        "(default 1)",
                        type="int",
                        default=1,
                        metavar="int")
parser.add_option("-C", "--combiner-flag", dest = "combiner", help = "Combine " \
            "contexts into one large network.", 
            action="store_true", default=False)
parser.add_option("-K", "--dcheck", dest = "dcheck", help = "DCheck " \
            "contexts after integration.", 
            action="store_true", default=False)
parser.add_option("-D", "--delete-contexts-flag", dest = "del_ctxt", 
            help = "Delete context integrations after combinining",
            action="store_true", default=False)

"""
(options, args) = parser.parse_args()

#You always need a gold standard, and CV has to be 1 or greater
if options.answers is None:
    sys.stderr.write("--answers-file file is required.\n")
    sys.exit()
if options.intervals < 0:
    sys.stderr.write("--k-intervals must be 1 or greater.\n")
    sys.exit()

#Add paths if required.
if options.sleipnir is not None:
    for tool in binaries:
        binaries[tool] = options.sleipnir + '/' + binaries[tool]


"""
if options.datadir is None:
    sys.stderr.write("--data-directory is required.\n")
    sys.exit()
if options.gfile is None:
    sys.stderr.write("--genes-file is required.\n")
    sys.exit()
"""

"""
#Make global network.
cmdline = binaries['Counter'] + ' -w ' + options.answers + \
                 ' -d ' + options.datadir + \
                 ' -o ' + options.workdir + '/' + \
                 ' -Z ' + options.zeros
stdout = Popen('qsub -wd ' + options.workdir + \
               ' -N GlobalCounts -j y ' \
               '-o ' + options.workdir + '/ '\
               '-l gb=' + str(gb_per) + \
               ' "' + cmdline + '"',\
               shell=True, stdout=PIPE).stdout.read()
glob_count_job = job_queue_id.search(stdout).group('job_id')


ctxt_jobs = []
contexts = None
#Make context networks
#  helps if counts exists even w/o contexts for networks.bin
try:
    os.mkdir(options.workdir + '/counts')
except OSError:
    pass
if options.contdir is not None:
    contexts = os.listdir(options.contdir)
    cmdline = binaries['Counter'] + ' -w ' + options.answers + \
              ' -d ' + options.datadir + \
              ' -o ' + options.workdir + '/counts/' + \
              ' -Z ' + options.zeros
    contexts_submitted = 0
    while contexts_submitted < len(contexts):
        job_ctxts = contexts[contexts_submitted:(contexts_submitted + options.job_threads)]
        job_thds = len(job_ctxts)
        job_cmd = cmdline + ' -t ' + str(job_thds) + ' ' + ' '.join([options.contdir + '/' + context for context in job_ctxts])
        stdout = Popen('qsub -wd ' + options.workdir + ' -N CtxtCounts -j y -o ' + options.workdir + '/ -l gb=' + str(gb_per * job_thds) + ' -l job_capacity=' + str(job_thds) + ' "' + job_cmd + '"', shell=True, stdout=PIPE).stdout.read()
        ctxt_jobs.append(job_queue_id.search(stdout).group('job_id'))
        contexts_submitted += job_thds

#Make networks.bin file
os.system("ls " + options.datadir + "/*dab | perl -pe 's/.*\/(.*)\.q?dab/$.\t$+/' > " + options.workdir + "/datasets.txt")
cmdline = binaries['Counter'] + " -k " + options.workdir + '/counts/ -o ' + options.workdir + '/networks.bin -s ' + options.workdir + '/datasets.txt -b ' + options.workdir + '/global.txt'
#Regularization
if options.alphas is not None:
    cmdline = cmdline + " -a " + options.alphas + " -p " + str(options.pseudo)
depends = glob_count_job
if options.contdir is not None:
    depends += ',' + ','.join(ctxt_jobs)
    os.system("ls " + options.contdir + "/* | perl -pe 's/.*\/(.*)/$.\t$+\t$+/' > " + options.workdir + "/contexts.txt")
    cmdline = cmdline + " -X " + options.workdir + "/contexts.txt"
stdout = Popen('qsub -wd ' + options.workdir + ' -N NetworksBin -j y -o ' + options.workdir + '/ -hold_jid ' + depends + ' -l 1hr "' + cmdline + '"', shell=True, stdout=PIPE).stdout.read()
networks_job = job_queue_id.search(stdout).group('job_id')

#Make global predictions
try:
    os.mkdir(options.workdir + '/predictions')
except OSError:
    pass
cmdline = binaries['Counter'] + " -n " + options.workdir + '/networks.bin -s ' + options.workdir + '/datasets.txt -d ' + options.datadir + ' -e ' + options.gfile + ' -o ' + options.workdir + '/predictions' + ' -Z ' + options.zeros
stdout = Popen('qsub -wd ' + options.workdir + ' -N GlobalPred -j y -o ' + options.workdir + '/ -hold_jid ' + networks_job + ' -l gb=' + str(gb_per) + ' "' + cmdline + '"', shell=True, stdout=PIPE).stdout.read()
glob_pred_job = job_queue_id.search(stdout).group('job_id')
if options.contdir is not None:
    contexts_submitted = 0
    ctxt_jobs = []
    while contexts_submitted < len(contexts):
        job_ctxts = contexts[contexts_submitted:(contexts_submitted + options.job_threads)]
        job_thds = len(job_ctxts)
        job_cmd = cmdline + ' -t ' + str(job_thds) + ' ' + ' '.join([options.contdir + '/' + context for context in job_ctxts])
        stdout = Popen('qsub -wd ' + options.workdir + ' -N CtxtPred -j y -o ' + options.workdir + '/ -hold_jid ' + networks_job + ' -l gb=' + str(gb_per * job_thds) + ' -l job_capacity=' + str(job_thds) + ' "' + job_cmd + '"', shell=True, stdout=PIPE).stdout.read()
        ctxt_jobs.append(job_queue_id.search(stdout).group('job_id'))
        contexts_submitted += job_thds

dcheck_jobs = []
if options.holdout is not None and options.dcheck:
    try:
        os.mkdir(options.workdir + '/dcheck')
    except OSError:
        pass
    depends = glob_pred_job + ',' + ','.join(ctxt_jobs)
    for context in contexts:
        job_cmd = binaries['DChecker'] + ' -i ' + options.workdir + '/predictions/' + context + '.dab -w ' + options.holdout + ' -l ' + options.contdir + '/' + context + ' > ' + options.workdir + '/dcheck/' + context + '.auc'
        stdout = Popen('qsub -wd ' + options.workdir + ' -N DChecker -j y -o ' + options.workdir + '/ -hold_jid ' + depends + ' -l 1hr "' + job_cmd + '"', shell=True, stdout=PIPE).stdout.read()
        dcheck_jobs.append(job_queue_id.search(stdout).group('job_id'))

if options.contdir is not None and options.combiner:
    depends = None
    if options.holdout is not None and options.dcheck:
        depends = ','.join(dcheck_jobs)
    else:
        depends = glob_pred_job + ',' + ','.join(ctxt_jobs)
    job_cmd = 'mv ' + options.workdir + '/predictions/global.dab ' + options.workdir + ' ; ' + binaries['NetworkCombiner'] + ' -o ' + options.workdir + '/global_average.dab -d ' + options.workdir + '/predictions/' 
    stdout = Popen('qsub -wd ' + options.workdir + ' -N CtxtAvg -j y -o ' + options.workdir + '/ -hold_jid ' + depends + ' -l gb=' + str(gb_per) + ' "' + job_cmd + '"', shell=True, stdout=PIPE).stdout.read()
    combine_job = job_queue_id.search(stdout).group('job_id')

    if options.del_ctxt:
        job_cmd = 'rm -fr ' + options.workdir + '/predictions/'
        stdout = Popen('qsub -wd ' + options.workdir + ' -N CtxtRm -j y -o ' + options.workdir + '/ -hold_jid ' + combine_job + ' "' + job_cmd + '"', shell=True, stdout=PIPE).stdout.read()
"""
