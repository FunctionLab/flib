from optparse import OptionParser
import re
import sys
import os
import commands
from idmap import idmap
from ortho import ortho
from go import go
from math import sqrt
from math import exp,log
from itertools import count, izip
from operator import itemgetter, attrgetter


"""
input decision_values, real_labels{1,-1}, #positive_instances, #negative_instances
output [A,B] that minimize sigmoid likilihood
"""
def SigmoidTrain(deci, label, prior1=None, prior0=None):
	#Count prior0 and prior1 if needed
	if prior1==None or prior0==None:
		prior1, prior0 = 0, 0
		for i in range(len(label)):
			if label[i] > 0:
				prior1+=1
			else:
				prior0+=1
	
	#Parameter Setting
	maxiter=100	#Maximum number of iterations
	minstep=1e-10	#Minimum step taken in line search
	sigma=1e-12	#For numerically strict PD of Hessian
	eps=1e-5
	
	#Construct Target Support
	hiTarget=(prior1+1.0)/(prior1+2.0)
	loTarget=1/(prior0+2.0)
	length=prior1+prior0
	t=[]
        
	for i in range(length):
		if label[i] > 0:
			t.append(hiTarget)
		else:
			t.append(loTarget)

	#Initial Point and Initial Fun Value
	A,B=0.0, log((prior0+1.0)/(prior1+1.0))
	fval = 0.0

	for i in range(length):
		fApB = deci[i]*A+B
		if fApB >= 0:
			fval += t[i]*fApB + log(1+exp(-fApB))
		else:
			fval += (t[i] - 1)*fApB +log(1+exp(fApB))

	for it in range(maxiter):
		h11=h22=sigma
		h21=g1=g2=0.0
		for i in range(length):
			fApB = deci[i]*A+B
			if (fApB >= 0):
				p=exp(-fApB)/(1.0+exp(-fApB))
				q=1.0/(1.0+exp(-fApB))
			else:
				p=1.0/(1.0+exp(fApB))
				q=exp(fApB)/(1.0+exp(fApB))
			d2=p*q
			h11+=deci[i]*deci[i]*d2
			h22+=d2
			h21+=deci[i]*d2
			d1=t[i]-p
			g1+=deci[i]*d1
			g2+=d1

		#Stopping Criteria
		if abs(g1)<eps and abs(g2)<eps:
			break

		det=h11*h22-h21*h21
		dA=-(h22*g1 - h21 * g2) / det
		dB=-(-h21*g1+ h11 * g2) / det
		gd=g1*dA+g2*dB

		#Line Search
		stepsize = 1
		while stepsize >= minstep:
			newA = A + stepsize * dA
			newB = B + stepsize * dB

			#New function value
			newf = 0.0
			for i in range(length):
				fApB = deci[i]*newA+newB
				if fApB >= 0:
					newf += t[i]*fApB + log(1+exp(-fApB))
				else:
					newf += (t[i] - 1)*fApB +log(1+exp(fApB))

			#Check sufficient decrease
			if newf < fval + 0.0001 * stepsize * gd:
				A, B, fval = newA, newB, newf
				break
			else:
				stepsize = stepsize / 2.0

		if stepsize < minstep:
			print "line search fails",A,B,g1,g2,dA,dB,gd
			return [A,B]

	if it>=maxiter-1:
		print "reaching maximal iterations",g1,g2
	return [A,B]

# score & parameter array [A, B]
def sigmoid2prob(score, go_param):
    fApB = score * go_param[0] + go_param[1]
    if (fApB >= 0):
        return exp(-fApB)/(1.0+exp(-fApB))
    else:
        return 1.0/(1+exp(fApB)) 
    

usage = "usage: %prog [options]<Prediction files>+"
parser = OptionParser(usage, version = "%prog dev-unreleased")

parser.add_option("-o", "--output-dir", dest="output", help="output directory", metavar="STRING")

(options, args) = parser.parse_args()

if options.output is None:
    sys.stderr.write("--output-dir is required.\n")
    sys.exit()
    
for fname in args:
    print >> sys.stderr, fname
    
    predictf = open(fname)
    deci = []
    ans = []
    pairs = []
    for line in predictf:       
        fields = line.split()
	if len(fields) > 0 and line[0] == '#':
		continue
        
        gene = fields[0]
        label = float(fields[1])
        score = float(fields[2])
        
        if label != 0:
            deci.append(score)
            ans.append(label)
        
        pairs.append([gene, label, score])
        
    predictf.close()
    params = SigmoidTrain(deci, ans)
    
    outf = open(options.output + '/'  + os.path.basename(fname) + '.prob', 'w')
    for pair in pairs:
        prob = sigmoid2prob(pair[2], params)
        
        outf.write('\t'.join([ pair[0], str(pair[1]), str(prob)])+'\n')
    outf.close()
    