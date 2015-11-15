from optparse import OptionParser
from subprocess import Popen, PIPE
import re
import os
import sys

gb_per = 4

job_queue_id = re.compile("Your job (?P<job_id>\d+) .*")

usage = "usage: %prog [options]"
parser = OptionParser(usage, version="%prog dev-0.0.1")
parser.add_option("-a", "--answers-file",
                        dest="answers",
                        help="Answers file",
                        metavar="FILE")
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
parser.add_option("-q", "--context-pos",
			dest="ctxtpos",
			help="Use positive edges between context genes (default=on)",
			action="store_false", default=True)
parser.add_option("-Q", "--context-neg",
			dest="ctxtneg",
			help="Use negative edges between context genes (default=on)",
			action="store_false", default=True)
parser.add_option("-j", "--bridge-pos",
			dest="bridgepos",
			help="Use bridging positives between context and non-context genes (default=off)",
			action="store_true", default=False)
parser.add_option("-J", "--bridge-neg",
			dest="bridgeneg",
			help="Use bridging negatives between context and non-context genes (default=on)",
			action="store_false", default=True)
parser.add_option("-W", "--weights-file",
                        dest="weights",
                        help="Weighted context file (for weighted Counter)",
                        metavar="FILE")
parser.add_option("-F", "--flipneg-flag",
                        dest="flipneg",
                        help="Flip weights (one minus original) for negative standards (default=on)",
                        action="store_false", default=True)
parser.add_option("-N", "--noweightneg-flag",
                        dest="noweightneg",
                        help="Use weight one for all negative standards (when using weight file, default=off)",
                        action="store_true", default=False)
parser.add_option("-m", "--multiplier",
                        dest="multiplier",
                        help="multiplier used for weighting (default = 1)",
                        type="int", default=1)
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
parser.add_option("-b", "--sleipnir-binaries-dir",
                        dest="sleipnir",
                        help="Directory containing sleipnir binaries",
                        metavar="FILE")
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
parser.add_option("-L", "--log-dir", dest="log_dir",
                        help="directory to put qsub files",
                        type="string",metavar="DIRECTORY", default=".")


(options, args) = parser.parse_args()

if options.answers is None:
    sys.stderr.write("--answers-file file is required.\n")
    sys.exit()
if options.datadir is None:
    sys.stderr.write("--data-directory is required.\n")
    sys.exit()
if options.gfile is None:
    sys.stderr.write("--genes-file is required.\n")
    sys.exit()

binaries = {'Counter':'Counter', 'NetworkCombiner':'NetworkCombiner', 'DChecker':'DChecker'}
if options.sleipnir is not None:
    for tool in binaries:
        binaries[tool] = options.sleipnir + '/' + binaries[tool]

# Establish output file directory
out_dir = options.workdir
if options.log_dir is not None:
    out_dir = options.log_dir

#Make global network.
cmdline = binaries['Counter'] + ' -w ' + options.answers + \
                 ' -d ' + options.datadir + \
                 ' -o ' + options.workdir + '/' + \
                 ' -Z ' + options.zeros
if options.weights is not None:
   cmdline = cmdline + ' -W ' + options.weights
   if options.flipneg is not True:
       cmdline = cmdline + ' -F '
   if options.multiplier is not None:
       cmdline = cmdline + ' -f ' + str(options.multiplier)
   if options.noweightneg is True:
       cmdline = cmdline + ' -N '
stdout = Popen('qsub -wd ' + options.workdir + \
               ' -N GlobalCounts -j y ' \
               '-o ' + out_dir + '/ '\
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
    if not options.ctxtpos:
        cmdline = cmdline + ' -q '
    if not options.ctxtneg:
        cmdline = cmdline + ' -Q '
    if options.bridgepos:
        cmdline = cmdline + ' -j '
    if not options.bridgeneg:
        cmdline = cmdline + ' -J ' 
    contexts_submitted = 0
    while contexts_submitted < len(contexts):
        job_ctxts = contexts[contexts_submitted:(contexts_submitted + options.job_threads)]
        job_thds = len(job_ctxts)
        job_cmd = cmdline + ' -t ' + str(job_thds) + ' ' + ' '.join([options.contdir + '/' + context for context in job_ctxts])
        stdout = Popen('qsub -wd ' + options.workdir + ' -N CtxtCounts -j y -o ' + out_dir + '/ -l gb=' + str(gb_per * job_thds) + ' -l job_capacity=' + str(job_thds) + ' "' + job_cmd + '"', shell=True, stdout=PIPE).stdout.read()
        ctxt_jobs.append(job_queue_id.search(stdout).group('job_id'))
        contexts_submitted += job_thds

#Make networks.bin file
os.system("ls " + options.datadir + "/*dab | perl -pe 's/.*\/(.*)\.q?dab/$.\t$+/' > " + options.workdir + "/datasets.txt")
cmdline = binaries['Counter'] + " -k " + options.workdir + '/counts/ -o ' + options.workdir + '/networks.bin -s ' + options.workdir + '/datasets.txt' 
if options.weights:
    weights_file = os.path.basename(options.weights)
    cmdline = cmdline + ' -b ' + options.workdir + '/' + weights_file[0:len(weights_file)-4] + '.txt'
else:
    cmdline = cmdline + ' -b ' + options.workdir + '/global.txt'
#Regularization
if options.alphas is not None:
    cmdline = cmdline + " -a " + options.alphas + " -p " + str(options.pseudo)
depends = glob_count_job
if options.contdir is not None:
    depends += ',' + ','.join(ctxt_jobs)
    os.system("ls " + options.contdir + "/* | perl -pe 's/.*\/(.*)/$.\t$+\t$+/' > " + options.workdir + "/contexts.txt")
    cmdline = cmdline + " -X " + options.workdir + "/contexts.txt"
stdout = Popen('qsub -wd ' + options.workdir + ' -N NetworksBin -j y -o ' + out_dir + '/ -hold_jid ' + depends + ' -l 1hr "' + cmdline + '"', shell=True, stdout=PIPE).stdout.read()
networks_job = job_queue_id.search(stdout).group('job_id')

#Make global predictions
try:
    os.mkdir(options.workdir + '/predictions')
except OSError:
    pass
cmdline = binaries['Counter'] + " -n " + options.workdir + '/networks.bin -s ' + options.workdir + '/datasets.txt -d ' + options.datadir + ' -e ' + options.gfile + ' -o ' + options.workdir + '/predictions' + ' -Z ' + options.zeros
stdout = Popen('qsub -wd ' + options.workdir + ' -N GlobalPred -j y -o ' + out_dir + '/ -hold_jid ' + networks_job + ' -l gb=' + str(gb_per) + ' "' + cmdline + '"', shell=True, stdout=PIPE).stdout.read()
glob_pred_job = job_queue_id.search(stdout).group('job_id')
if options.contdir is not None:
    contexts_submitted = 0
    ctxt_jobs = []
    while contexts_submitted < len(contexts):
        job_ctxts = contexts[contexts_submitted:(contexts_submitted + options.job_threads)]
        job_thds = len(job_ctxts)
        job_cmd = cmdline + ' -t ' + str(job_thds) + ' ' + ' '.join([options.contdir + '/' + context for context in job_ctxts])
        stdout = Popen('qsub -wd ' + options.workdir + ' -N CtxtPred -j y -o ' + out_dir + '/ -hold_jid ' + networks_job + ' -l gb=' + str(gb_per * job_thds) + ' -l job_capacity=' + str(job_thds) + ' "' + job_cmd + '"', shell=True, stdout=PIPE).stdout.read()
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
        job_cmd = binaries['DChecker'] + ' -i ' + options.workdir + '/predictions/' + context + '.dab -w ' + options.holdout + ' -c ' + options.contdir + '/' + context
        if not options.ctxtpos:
            job_cmd = job_cmd + ' -q '
        if not options.ctxtneg:
            job_cmd = job_cmd + ' -Q '
        if options.bridgepos:
            job_cmd = job_cmd + ' -j '
        if not options.bridgeneg:
            job_cmd = job_cmd + ' -J ' 
        job_cmd = job_cmd + ' > ' + options.workdir + '/dcheck/' + context + '.auc'
        stdout = Popen('qsub -wd ' + options.workdir + ' -N DChecker -j y -o ' + out_dir + '/ -hold_jid ' + depends + ' -l 1hr "' + job_cmd + '"', shell=True, stdout=PIPE).stdout.read()
        dcheck_jobs.append(job_queue_id.search(stdout).group('job_id'))

if options.contdir is not None and options.combiner:
    depends = None
    if options.holdout is not None and options.dcheck:
        depends = ','.join(dcheck_jobs)
    else:
        depends = glob_pred_job + ',' + ','.join(ctxt_jobs)
    job_cmd = 'mv ' + options.workdir + '/predictions/global.dab ' + options.workdir + ' ; ' + binaries['NetworkCombiner'] + ' -v0 -o ' + options.workdir + '/global_average.dab -d ' + options.workdir + '/predictions/' 
    stdout = Popen('qsub -wd ' + options.workdir + ' -N CtxtAvg -j y -o ' + options.workdir + '/ -hold_jid ' + depends + ' -l gb=' + str(gb_per) + ' "' + job_cmd + '"', shell=True, stdout=PIPE).stdout.read()
    combine_job = job_queue_id.search(stdout).group('job_id')

    if options.del_ctxt:
        job_cmd = 'rm -fr ' + options.workdir + '/predictions/'
        stdout = Popen('qsub -wd ' + options.workdir + ' -N CtxtRm -j y -o ' + out_dir + '/ -hold_jid ' + combine_job + ' "' + job_cmd + '"', shell=True, stdout=PIPE).stdout.read()
