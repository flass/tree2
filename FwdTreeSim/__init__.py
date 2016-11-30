#!/usr/bin/python
# -*- coding: utf-8 -*-

__all__ = ["models", "simulators", "genome"]

__author__ = "Florent Lassalle <f.lassalle@imperial.ac.uk>"
__date__ = "22 August 2016"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""

import os
import cPickle as pickle
#~ import pickle

nodelabelprefix = dict(livetip='S', deadtip='E', node='N')

def loadpickle(fileorpath, autoclosefile=True):
	if type(fileorpath)==str and os.path.exists(os.path.dirname(fileorpath)):
		fpickle = open(fileorpath, 'r')
	elif isinstance(fileorpath, file):
		if not os.access(os.R_OK):
			raise ValueError, "argument file object is not readable"
		else:
			fpickle = fileorpath
	else:
		raise ValueError, "argument should be a file path or file-like object; '%s' is not"%(repr(fileorpath))
	obj = pickle.load(fpickle)
	if autoclosefile: fpickle.close()
	return obj

def deleteGenneratorAttr(obj, genattrname='eventidgen'):
	"""delete object attributes which are generators in order to pickle the object"""
	if genattrname in obj.__dict__:
		del obj.__dict__[genattrname]
	
		
def dumppickle(obj, fileorpath, autoclosefile=True, genattrname='eventidgen'):
	if type(fileorpath)==str and os.path.exists(os.path.dirname(fileorpath)):
		fpickle = open(fileorpath, 'w')
		fpathstr = " in file '%s'"%fileorpath
	elif isinstance(fileorpath, file):
		if not os.access(os.W_OK):
			raise ValueError, "argument file object is not writeable"
		else:
			fpickle = fileorpath
		fpathstr = ""
	else:
		raise ValueError, "argument should be a file path or file-like object; '%s' is not"%(repr(fileorpath))
	# check that objects don't have generators and delete them if there are
	if isinstance(obj, simulators.BaseTreeSimulator):
		deleteGenneratorAttr(obj, genattrname=genattrname)
	elif isinstance(obj, list):
		for subobj in obj:
			if isinstance(subobj, simulators.BaseTreeSimulator):
				deleteGenneratorAttr(subobj, genattrname=genattrname)
	elif isinstance(obj, dict):
		for subobj in obj.values():
			if isinstance(subobj, simulators.BaseTreeSimulator):
				deleteGenneratorAttr(subobj, genattrname=genattrname)
			
	pickle.dump(obj, file=fpickle, protocol=2)
	print "saved %s in binary format"%repr(obj) + fpathstr
	if autoclosefile: fpickle.close()
