#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Simulate phylogenetic tree(s) forward in time, notably under models recording events of duplication, transfer and loss for gene trees relative to a species tree."""

__author__ = "Florent Lassalle <florent.lassalle@imperial.ac.uk>"
__date__ = "27 July 2016"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""

import os
import copy
import random
from tree2.FwdTreeSim import models, nodelabelprefix, loadpickle, dumppickle
import tree2

######################################
# Tree Simulators
# wrappers for the simulation models; perfom the simulation from t0 to the end
######################################

class BaseTreeSimulator(object):
	"""
	Base wrapper class for forward tree simulations.
	
	Main attributes:
	- self.model           : a model instance from *Model classes from 'models' module.
	- self.tree            : a (list of) tree2.Node (or derived class) phylogenetic tree object; 
	or self.trees             uses the tree provided through 'starttree' argument or (default) a single-node tree (set) as a seed for the simulation.
	- self.t               : the number of the last simulation iteration / time slice.
	
	Attributes than can be provided at instance intitiation through kyword arguments
	- self.times           : list of evolutionary time spent within each time slice.
	- self.extincts        : list of leaf nodes representing lineages considered extinct.
	- self.eventsrecord    : record events by time slice where they occurred;
	                          structured as a list (which indexes correspond to the time slice number) of lists of event objects (models.BaseEvent children classes).
	- self.eventsmap       : record events by tree branch on which they occurred (in case of transfer events, donor node is used);
	                          structured as a dictionary (which keys are reference tree node labels) of lists of event objects (models.BaseEvent children classes).
	- self.contempbranches : record by time slice of the extant branches present at this time that were subject to the evolutionary processes;
	                          structured as a list of list of of node objects.
	
	Falcultative keyword arguments:
	- ngen                 : the end number of generation to simulate (count includes potential previous iterations).
	- noTrigger            : when initiating an instance, automatic checkdata() and evolve(ngen) triggers from parent classes
	                          are deactivated when noTrigger=True.
	"""
	
	def __init__(self, model, **kwargs):
		print 'invoke _BaseTreeSimulator__init__'
		print 'kwargs:', kwargs
		self.model = model
		self.t = 0
		self.times = kwargs.get('times', [])
		self.eventsrecord = kwargs.get('eventsrecord', [])
		self.eventsmap = kwargs.get('eventsmap', {})
		self.extincts = kwargs.get('extincts', [])
		self.contempbranches =  kwargs.get('contempbranches', [])
		self.nodeClass = kwargs.get('nodeClass', tree2.Node)
		self.ngen = kwargs.get('ngen')
		self.eventidgen = models.eventIdGen()	# generator object providing serial unique identifiers for event objects
		
	def __getattr__(self, name):
		if name=='__dict__': return self.__dict__
		try:
			# search normal attributes
			return self.__dict__[name]
		except KeyError:
			try:
				# search private cached attirbutes
				return self.__dict__['_'+name]
			except KeyError:
				try:
					# look up an attribute-generating method
					return type(self).__dict__['get_'+name](self)
				except KeyError:
					raise AttributeError, "'%s' is not an attribute of %s instance of %s"%(name, repr(self), type(self))
		
	def checkdata(self):
		"""placeholder; method should be called from descendant classes"""
		for attr in ['times', 'eventsrecord', 'extincts', 'contempbranches']:
			assert len(getattr(self, attr)) == self.t
		return 0
		
	def get_timeslices(self):
		ts = []
		for t in range(self.t):
			low = 0 if t==0 else ts[t-1][1] # t low bound is t-1 up bound
			up = low + self.times[t]
			ts.append((low, up))
		return ts
		
	@staticmethod
	def prune_dead_lineages(tree, extincts):
		livetree = copy.deepcopy(tree)
		for leaf in extincts:
			livetree.pop(leaf.label())
		# set root branch lenght to zero
		livetree.set_lg(0)	
		return livetree
		
	############################
	## miscealaneous methods
		
	def dumppickle(self, fileorpath, autoclosefile=True):
		# cannot pickle generator objects, have to delete it
		if 'eventidgen' in self.__dict___:
			prompt =  "cannot pickle generator objects, have to delete event id generator 'self.eventidgen' first\n"
			prompt += "(will discontinue numeration of events for further simulation).\n"
			prompt += "Delete generator and pickle? (y/n) "
			doit = raw_input(prompt)
			while not (doit in ['y', 'n']):
				print "answer 'y' (for yes) or 'n' (for no)"
				doit = raw_input(prompt)
			if doit=='y':
				del self.eventidgen
				dumppickle(self, fileorpath, autoclosefile=autoclosefile)
			else:
				print "DID NOT save %s object"%repr(self)
		else:
			dumppickle(self, fileorpath, autoclosefile=autoclosefile)
			
	# loaddpickle directly inherited from FwdTreeSim through import; can be called either as FwdTreeSim.loaddpickle() or [FwdTreeSim.]simulators.loaddpickle()


################################################
# Classes describing how wether a single or many
# trees are simulated within the population
################################################
		
class SingleTreeSimulator(BaseTreeSimulator):
	def __init__(self, model, starttree=None, **kwargs):
		super(SingleTreeSimulator, self).__init__(model=model, **kwargs)
		self.tree = starttree or self.nodeClass(l=float(0), lab="Root")
		
		if not kwargs.get('noTrigger'):
			self.checkdata()
			# if ngen specified, launch simulation for ngen iterations
			if self.ngen: self.evolve(self.ngen)
		
	def checkdata(self):
		"""assert (self-)consistency of data, i.e. support of tree object class and that nodes in eventsrecord are included in the tree"""
		assert isinstance(self.tree, tree2.Node) # any descendant class of tree2.Node
		lerrnodes = []
		for ex in self.extincts:
			if ex not in self.tree:
				#~ raise TreeReferenceError(ex, self.tree)
				lerrnodes.append(ex)
		if lerrnodes: raise tree2.AggregateTreeReferenceError(lerrnodes, self.tree)
		for t in self.eventsrecord:
			evt = self.eventsrecord[t]
			for u in evt:
				evtu = evt[u]
				for v in evtu:
					if evtu[v] not in self.tree:
						#~ raise TreeReferenceError(evtu[v], self.tree)
						lerrnodes.append(evtu[v])
		if lerrnodes: raise tree2.AggregateTreeReferenceError(lerrnodes, self.tree)
		return 0

	def evolve(self, ngen, verbose=False, nodeathspan=[], stopcondition=(lambda x: (None, None)), **kwargs):
		"""simulation engine, iterates step of the simulation model"""
		evtidgen = kwargs.get('evtidgen', self.eventidgen)
		while self.t < ngen:
			self.t += 1
			levents, devents, contempbranches, brlen = self.model.stepforward(self, allowdeath=(self.t not in nodeathspan), evtidgen=evtidgen)
			# record what happened
			self.eventsrecord.append(levents)
			self.eventsmap.update(devents)
			self.contempbranches.append(contempbranches)
			self.times.append(float(brlen))
			if verbose: self.verbevolve(devents)
			# check if tree got full extinct
			s, v = stopcondition(self)
			if s:
				print "Stop at time %d because: %s"%(self.t, v)
				rval = 0
				break
			elif set(self.tree.get_leaves())==set(self.extincts):
				print "Extinction of all lineages at time", self.t
				rval = 1
				break
		else:
			# set root branch lenght to zero
			self.tree.set_lg(0)
			# destroys cache outdated attributes
			if hasattr(self, '_extanttree'): self._extanttree = None
			rval = 0
		return rval
	
	def verbevolve(self, devents):
		print "\ntime:", self.t,
		nextanttm1 = len(self.contempbranches[-1])
		#~ if self.t > 1:
			#~ nextanttm1 = len(self.contempbranches[self.t-2])	# aims at timeslice t-1, but shifted due to Python 0-based numbering ;
			#~ nextanttm1 = len(self.contempbranches[-1])	# aims at timeslice t-1, but shifted due to Python 0-based numbering ;
			#~ # equivalent to self.contempbranches[-1], but only if we keep calling this function above the update of contempbranches in evolve()
		#~ else:
			#~ nextanttm1 = 1
		print "\tstarted with %d live branches"%(nextanttm1)
		if (isinstance(self.model, models.UniformDiscreetBirthDeathModel) or isinstance(self.model, models.GenericDiscreetBirthDeathModel)):
			for evtype in devents:
				netgrowth = sum(node.nb_children() for node in devents[evtype]) - len(devents[evtype])
				print "  %s:\t%d\t(%.2g per lineage)"%(evtype, netgrowth, float(netgrowth)/nextanttm1)
		elif isinstance(self.model, PartialMoranProcess):
			for evtype in devents:
				print "\t%s:\t%d"%(evtype, sum(int(not (type(e) is int)) for e in devents[evtype])),
	
	################################	
	## post-simulation tree-cleaning
	
	def labeltreenodes(self, dictprefix=nodelabelprefix, onlyExtants=True):
		"""puts distinctive labes at internal nodes and extinct and extant leaves (default prefixes are N, E, S,respectively). Labelling follows increasing order from root"""
		self.tree.complete_internal_labels(prefix=dictprefix['livetip'], onlyLeaves=True, exclude=self.extincts)
		if not onlyExtants:
			self.tree.complete_internal_labels(prefix=dictprefix['deadtip'], onlyLeaves=True)
			self.tree.complete_internal_labels(prefix=dictprefix['node'], excludeLeaves=True)
		
	def scaletree(self, relToExtantMrca=True):
		if relToExtantMrca: t = self.get_extanttree(compute=True)
		else: t = self.tree
		root2tipdist = self.get_extants()[0].distance_root()
		self.tree /= root2tipdist
		self.extanttree /= root2tipdist
		self.timeslices = [float(ts)/root2tipdist for ts in self.timeslices]
		
		
	#############################
	## result description methods
			
	def get_extanttree(self, compute=True):
		"""returns phylogenetic tree 'cleaned' of its dead branches ; NB: original tree will have all its nodes labelled afterward"""
		if not compute: return self._extanttree
		# only works with all nodes being labelled, to use labels as references rather than the node objects (which refer to the original tree)
		self.labeltreenodes()
		livetree = self.prune_dead_lineages(self.tree, self.extincts)
		self._extanttree = livetree	# create or update cache attribute
		return livetree
				
	def nb_extant(self):
		return self.tree.nb_leaves() - len(self.extincts)
		
	def get_extants(self, depthsorted=False):
		e = list(set(self.tree.get_leaves()) - set(self.extincts))
		if depthsorted: e.sort(key=lambda x: x.depth())
		return e
	
class MultipleTreeSimulator(BaseTreeSimulator):
	def __init__(self, model, **kwargs):
		print 'invoke _MultipleTreeSimulator__init__'
		print 'kwargs:', kwargs
		super(MultipleTreeSimulator, self).__init__(model=model, **kwargs)
		self.popsize = kwargs.get('popsize', 100)
		self.trees = [self.nodeClass(l=float(0), lab="Root_%d"%i) for i in range(self.popsize)]
		
		if not kwargs.get('noTrigger'):
			self.checkdata()
			# if ngen specified, launch simulation for ngen iterations
			if self.ngen: self.evolve(self.ngen)
		
	def checkdata(self):
		"""assert (self-)consistency of data, i.e. support of tree object class and that nodes in eventsrecord are included in the trees"""
		for i in range(len(self.trees)):
			if not isinstance(self.trees[i], tree2.Node): raise TypeError, "element %d of 'trees' attribute is not a tree2.Node instance"%(i)
		lerrnodes = []
		allnodes = sum((t.get_all_children() for t in self.trees), [])
		for ex in self.extincts:
			if ex not in allnodes:
				lerrnodes.append(ex)
		if lerrnodes: raise tree2.AggregateTreeReferenceError(lerrnodes, self.trees)
		for t in self.eventsrecord:
			evt = self.eventsrecord[t]
			for u in evt:
				evtu = evt[u]
				for v in evtu:
					if evtu[v] not in allnodes:
						lerrnodes.append(evtu[v])
		if lerrnodes: raise tree2.AggregateTreeReferenceError(lerrnodes, self.trees)
		return 0
		
		
	def get_extants(self, depthsorted=False):
		"""generate a single list of all extant leaves across all trees in  the self.trees set"""
		sleave = set(sum((tree.get_leaves() for tree in self.trees), []))
		extants = list(sleave - set(self.extincts))
		if depthsorted: extants.sort(key=lambda x: x.depth())
		if isinstance(self.model, models.MoranProcess): assert len(extants) == self.popsize
		return extants
		
	def evolve(self, ngen, verbose=False, nodeathspan=[], stopcondition=(lambda x: (None, None)), **kwargs):
		"""simulation engine, iterates step of the simulation model"""
		evtidgen = kwargs.get('evtidgen', self.eventidgen)
		while self.t < ngen: 
			self.t += 1
			levents, devents, contempbranches, brlen = self.model.stepforward(self, allowdeath=(self.t not in nodeathspan), evtidgen=evtidgen)
			# record what happened
			self.eventsrecord.append(levents)
			self.eventsmap.update(devents)
			self.contempbranches.append(contempbranches)
			self.times.append(float(brlen))
			if verbose: self.verbevolve(devents)
			# check if simulation should stop, e.g. because tree got full extinct
			s, v = stopcondition(self)
			if s:
				print "Stop at time %d because: %s"%(self.t, v)
				rval = 0
				break
		else:
			# destroys cache outdated attributes
			if hasattr(self, '_extanttrees'): self._extanttrees = None
			rval = 0
		# give labels to extant tips
		for i, extant in enumerate(self.get_extants()): extant.edit_label("%s%d"%(nodelabelprefix['livetip'], i))
		return rval
		
	def get_most_extant_tree(self):
		"""choose the tree with most extant leaves for the focal output of simulation"""
		pass
	
	def verbevolve(self, devents):
		print "\ntime:", self.t
		
	def get_extants(self, depthsorted=False):
		# generate a single list of all extant leaves across all trees in  the currtees set
		sleave = set(sum((tree.get_leaves() for tree in self.trees), []))
		extants = list(sleave - set(self.extincts))
		if depthsorted: extants.sort(key=lambda x: x.depth())
		return extants
		
	#############################
	## result description methods
			
	def get_extanttrees(self, compute=True, use_copy=True):
		"""returns phylogenetic tree 'cleaned' of its dead branches ; NB: original tree will have all its nodes labelled afterward"""
		if not compute: return self._extanttrees
		livetrees = []
		if use_copy: trees = [t.deepcopybelow(keep_lg=True, add_ref_attr=True) for t in self.trees]
		else: trees = self.trees
		for tree in trees:
			# only works with all nodes being labelled, to use labels as references rather than the node objects (which refer to the original tree) in Node.pop() 
			self.labeltreenodes()
			livetrees.append(self.prune_dead_lineages(tree, self.extincts))
		self._extanttrees = livetrees	# create or update cache attribute
		return livetrees

class DTLtreeSimulator(MultipleTreeSimulator):
	"""Simulates gene trees from (species) reference trees under a DTL model (see Szollosi et al. 2013, Lateral Gene Transfer from the Dead, Systematic Biology 62(3):386–397)
	
	Attribute 'genetrees' (these gene trees could as well be a locus/plasmid tree as well as a gene tree) is generated orignially by deep-copying the reference tree set, 
	sampling a reference tree at proba 'rootfreq' (to model the heterogeneous ocurrence of accessory genes in genomes of a prokaryotic population).
	A reference to the the original node is kept for each node the gene tree copies under the additional attribute 'ref'.
	During the simulation proceess, modifications are made on gene trees, which can be: copying a subtree from itself and grafting it next to the original (duplication), 
	removing a subtree (loss), or removing a subtree and replacing it with a copy of a subtree __from the reference tree set__ (transfer).
	
	As for parents classes, events are recorded in the 'eventsrecord' and 'eventsmap' attributes, storing events object indexed by time slice and reference tree node origin, respectively. 
	For transfers, an additional dictionary attribute 'transferrec' record the event objects indexed by reference tree node destination.
	"""
	def __init__(self, model, rootfreq=0.5, noTrigger=False, **kwargs):
		print 'invoke _DTLtreeSimulator__init__'
		print 'kwargs:', kwargs
		super(MultipleTreeSimulator, self).__init__(model=model, noTrigger=True, **kwargs) 
		# automatic checkdata() and evolve(ngen) triggers from parent classses are deactivated by noTrigger=True
		self.rootfreq = rootfreq # frequency a which the gene family is found at the root of each tree in the multiple reference tree set (species lineages from a Moran process)
		refsimul = kwargs.get('refsimul')
		if refsimul:
			self.reftrees = refsimul.trees
			self.refconbran = refsimul.contempbranches
			self.reftimeslices = refsimul.get_timeslices()
			self.ngen = refsimul.t - 1
			self.times = refsimul.times
		else:
			self.reftrees = kwargs['reftrees']
			self.refconbran = kwargs['contempbranches']
			self.reftimeslices = kwargs['timeslices']
			self.times = kwargs['times']
			self.ngen = kwargs.get('ngen')
			
		# with a proba rootfreq, generate a gene tree for each reference tree;
		# gene trees are like the reference tree set, with an addition of a link of each gene tree node to its reference tree node	
		self.genetrees = [rt.deepcopybelow(keep_lg=True, add_ref_attr=True) for rt in self.reftrees if random.random()<self.rootfreq]
		self.eventsmap = {}
		self.transferrec = {}
		
		if not noTrigger:
			# self.checkdata()
			# if ngen specified, launch simulation for ngen iterations
			if self.ngen: self.evolve(self.ngen)
		
	def get_current_branches(self, l):
		"""retrieves branches reaching the root-to-tip length of l"""
		currbranches = []
		for tree in self.genetrees:
			# beware! after evolution, the node 'tree' stored in the 'self.genetrees' list may not be the root of the tree, due to addition of nodes for transfer or duplication event
			currbranches += tree.go_root().get_comtemporary_branches(l)
		return currbranches
		
	def evolve(self, ngen, verbose=False, stopcondition=(lambda x: (None, None)), **kwargs):
		"""simulation engine, iterates step of the simulation model"""
		evtidgen = kwargs.get('evtidgen', self.eventidgen)
		while self.t < ngen:
			self.t += 1
			print "t=%d"%self.t, 
			timeslice = self.reftimeslices[self.t]
			# get current bene tree branches
			currbranches = self.get_current_branches((timeslice[1]+timeslice[0])/2) # input
			#~ print ": currbranches", [cb.label() for cb in currbranches]
			levents, devents, trec = self.model.stepforward(currbranches, self.refconbran[self.t], timeslice, t=self.t, evtidgen=evtidgen)
			# record what happened
			self.contempbranches.append(currbranches)
			self.eventsrecord.append(levents)
			self.eventsmap.update(devents)
			self.transferrec.update(trec)
			# check if simulation should stop
			s, v = stopcondition(self)
			if s:
				print "Stop at time %d because: %s"%(self.t, v)
				rval = 0
				break
		else:
			rval = 0
		return rval
