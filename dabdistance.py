from optparse import OptionParser
import random
import os
import sys

#CONSTANTS
binaries = { 'Counter':'Counter', 'NetworkCombiner':'NetworkCombiner', 'DChecker':'DChecker', 'Dat2Dab':'Dat2Dab', 'Distancer': 'Distancer' }

#FUNCTIONS
#functions that call Distancer
#Counter Learn
def distance(job, job_name, distancer, measure, inputfile, outputfile, threads, pbs_folder):
    #Make global network.
    job.ppn = threads
    cmdline = distancer +   ' -i ' + inputfile + \
                            ' -o ' + outputfile + \
                            ' -d ' + measure + \
                            ' -t ' + str(threads) + \
                            ' -z '
    job.set_name_command(job_name + '-Distance', cmdline)
    job.submit(pbs_folder + '/' + job_name + 'Distance.pbs')

#CMD LINE PROCESSING
usage = "usage: %prog [options]"
parser = OptionParser(usage, version="%prog dev-0.0.1")
#CONTROL PARAMETERS
#Cross validation
parser.add_option("-i", "--input-directory",
                        dest="input",
                        type="string",
                        metavar="DIR",
                        help="Directory of input DABs.")
parser.add_option("-o", "--output-directory",
                        dest="output",
                        type="string",
                        metavar="DIR",
                        help="Directory for output DABs.")
#Queue details
parser.add_option("-Q", "--queue",
                        dest="queue",
                        type="string",
                        default="discovery",
                        help="What cluster parameters should be used (e.g. discovery, cetus)?")
parser.add_option("-E", "--email",
                        dest="email",
                        type="string",
                        help="An e-mail address is required for submission to some clusters (e.g. discovery).")
#Where do programs exist
parser.add_option("-B", "--sleipnir-binaries-dir",
                        dest="sleipnir",
                        help="DIR containing sleipnir binaries. If not passed, binaries used must be available in the path.",
                        metavar="DIR")
#How should this be run?
parser.add_option("-T", "--threads-per",
                        dest="threads",
                        help="Number of threads per job",
                        type="int",
                        default=1,
                        metavar="int")
parser.add_option("-m", "--measure",
                        dest="measure",
                        help="Which distancer measure should be used?.",
                        type="string",
                        metavar="string")
#Minor Details
parser.add_option("-w", "--working-directory",
                        dest="workdir",
                        help="Perform integration in DIR.",
                        type="string",
                        metavar="DIR")
(options, args) = parser.parse_args()

#Add paths if required.
if options.sleipnir is not None:
    for tool in binaries:
        binaries[tool] = options.sleipnir + '/' + binaries[tool]

job = None
if options.queue == 'discovery':
    from pbsjob import PBSJob
    job = PBSJob(addr=options.email, command="echo test", walltime="47:59:00", queue="largeq")
else:
    sys.stderr.write("Unknown queue '" + options.queue + "' -- Is there an implemented job interface for this queue?")
    sys.exit()

#find networks
networks = os.listdir(options.input)

for network in networks:
    distance(job, network.split('.')[0], binaries['Distancer'], options.measure, os.path.join(options.input, network), os.path.join(options.output, network), options.threads, options.workdir)


