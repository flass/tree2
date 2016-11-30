#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Models for simulation of phylogenetic tree(s) forward in time, specifically recording events of duplication, transfer and loss for gene trees relative to a species tree."""

__author__ = "Florent Lassalle <florent.lassalle@imperial.ac.uk>"
__date__ = "27 July 2016"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""

#~ import copy
import tree2
import random
from numpy.random import poisson, geometric, exponential, normal

import sys # for debug

def geometricm1(p, size=None):
	return (geometric(p=p, size=size) - 1)

def normalp1(scale, size=None):
	return normal(loc=1, scale=scale, size=size)

nodelabelprefix = tree2.FwdTreeSim.nodelabelprefix
# RGB components of colours
eventcolcode = {'birth':[0,255,0], 'death':[0,0,255], 'loss':[127,127,127], 'transfer':[200,200,0], 'duplication':[255,0,0]}
# eventcolcode = {'birth':'green', 'death':'blue', 'loss':'grey', 'transfer':'gold', 'duplication':'red'}

######################################
# Birth-Death generative models
# implement the t to t+1 process of tree simulation 
######################################

class BaseModel(object): 
	"""Works with discreet time steps, where branches have to start and stop.
	
	Branches are incrementally elongated, and possibly bifurcated or terminated at each time step, so growth over the tree is homogeneous (let alone for extinct branches) resulting 
	in an ultrametric tree..
	"""
	def __init__(self, **kwargs):
		print 'invoke _BaseModel__init__'
		self.__rseed = kwargs.get('randomseed', None)	# keep track of random seed as instance attribute so simulation can be repeated
		if isinstance(self.__rseed, tuple) and len(self.__rseed)==3: random.setstate(self.rseed)
		elif self.__rseed is None:  self.__rseed = random.getstate()
		else: random.seed(self.__rseed)	# provided seed must be hashable type
		self.tunit = kwargs.get('tunit', 1)
		
	@staticmethod
	def annotateNode(node, eventtype, colorcode=eventcolcode):
		if isinstance(node, tree2.AnnotatedNode):
			node.edit_color(colorcode[eventtype])
		#~ if isinstance(node, tree2.GeneTree):
			#~ pass
		
class SingleTreeModel(BaseModel):
	
	def __init__(self, **kwargs):
		print 'invoke _SingleTreeModel__init__'
		super(SingleTreeModel, self).__init__(**kwargs)
	
	#~ @staticmethod
	#~ def get_extants(currtree, extincts, depthsorted=False):
		#~ sleave = set(currtree.get_leaves())
		#~ extants = list(sleave - set(extincts))
		#~ if depthsorted: extants.sort(key=lambda x: x.depth())
		#~ return extants
		
class MultipleTreeModel(BaseModel):
	
	def __init__(self, **kwargs):
		print 'invoke _MultipleTreeModel__init__'
		super(MultipleTreeModel, self).__init__(**kwargs)
	
	#~ @staticmethod
	#~ def get_extants(currtrees, extincts, depthsorted=False):
		#~ # generate a single list of all extant leaves across all trees in  the currtees set
		#~ sleave = set(sum((currtree.get_leaves() for currtree in currtrees), []))
		#~ extants = list(sleave - set(extincts))
		#~ if depthsorted: extants.sort(key=lambda x: x.depth())
		#~ return extants
		

class UniformDiscreetBirthDeathModel(SingleTreeModel):
	"""Works with discreet time steps, where branches have to start and stop.
	
	Branches are incrementally elongated, and possibly bifurcated or terminated at each time step, so growth over the tree is homogeneous (let alone for extinct branches) resulting 
	in an ultrametric tree.
	Growth of branch is unit over a timeslice; birth or death events occur at a given probability over the time slice.
	"""
	def __init__(self, **kwargs):
		print 'invoke _UniformDiscreetBirthDeathModel__init__'
		super(UniformDiscreetBirthDeathModel, self).__init__(**kwargs)	
		self.bprob = kwargs.get('birthprob', 1)
		self.dprob = kwargs.get('deathprob', 1)
		assert self.bprob >= 0 and self.bprob <= 1
		assert self.dprob >= 0 and self.dprob <= 1
	
	def stepforward(self, simul, t=None, allowdeath=True, evtidgen=None, **kwargs):
		"""Implement the atomic step of simulation iteration. 'extincts' is a list of leaf nodes in the tree"""
		levents = []
		devents = dict(birth=[], death=[])
		extants = simul.get_extants()
		# test for bifurcation first
		for leaf in extants:
			if random.random() <= self.bprob:
				# adds two children, growing their branch at the same time
				for i in range(2):
					newchild = leaf.newnode()
					leaf.link_child(newchild, newlen=self.tunit)
					BDevent('birth', leaf, t, evtidgen, levents, devents)
					self.annotateNode(leaf, 'birth')
			else:
				# just grow the current branch
				leaf += self.tunit
		# then test for extinction, that can thus occur on a newly bifurcated branch
		if allowdeath:
			newextants = simul.get_extants()
			for leaf in newextants:
				if random.random() <= self.dprob:
					simul.extincts.append(leaf)
					BDevent('death', leaf, t, evtidgen, levents, devents)
					self.annotateNode(leaf, 'death')
		return (levents, devents, extants, self.tunit)

class GenericDiscreetBirthDeathModel(SingleTreeModel):
	"""Works with discreet time steps, where branches have to start and stop.
	
	Branches are incrementally elongated, and possibly bifurcated or terminated at each time step, so growth over the tree is homogeneous (let alone for extinct branches) resulting 
	in an ultrametric tree.
	Growth of branch is unit over a timeslice; progeny of a node (including the possibility of zero descendant, i.e. extinction) follows a random process 
	to be specified through 'randprocess' argument (a function) and its _single_ parameter 'randprocparam' (a float).
	"""
	def __init__(self, **kwargs):	
		print 'invoke _GenericDiscreetBirthDeathModel__init__'
		super(GenericDiscreetBirthDeathModel, self).__init__(**kwargs)
		self.randprocess = kwarg['randprocess']
		self.randprocparam = kwarg['randprocparam']	# only support single paprameters
	
	def stepforward(self, simul, allowdeath=True, evtidgen=None, **kwargs):
		"""Implement the atomic step of simulation iteration. 'extincts' is a list of leaf nodes in the tree"""
		# NB: assumes timeslice numbering starts with 1; t=0 would induce null-rate at first step
		levents = []
		devents = dict(birth=[], death=[])
		extants = simul.get_extants()
		t = simul.t
		nb = 0
		nd = 0
		for leaf in extants:
			nprogeny = self.randprocess(self.randprocparam)
			if nprogeny == 0:
				if allowdeath:
					# lineage goes extinct
					simul.extincts.append(leaf)
					BDevent('death', leaf, t, evtidgen, levents, devents)
					leaf.edit_label("%s%d.%d"%(nodelabelprefix['deadtip'], t, nd))
					self.annotateNode(leaf, 'death')
					nd += 1
				else:
					nprogeny = 1
			if nprogeny == 1:
				# just grow the current branch
				leaf += self.tunit
			elif nprogeny > 1:
				# adds nprogeny children, growing their branch at the same time
				BDevent('birth', leaf, t, evtidgen, levents, devents)
				self.annotateNode(leaf, 'birth')
				leaf.edit_label("%s%d.%d"%(nodelabelprefix['node'], t, nb))
				nb += 1
				for i in range(nprogeny):
					newchild = leaf.newnode()
					leaf.link_child(newchild, newlen=self.tunit)
		return (levents, devents, extants, self.tunit)
			
class PoissonBirthDeathModel(GenericDiscreetBirthDeathModel):
	"""The number of daughter lineages is sampled from a Poisson distribution.
	
	Mean number of daughter lineages is mu; default parameter mu=1 leads to distribution with mean of 1"""
	def __init__(self, **kwargs):
		print 'invoke _PoissonBirthDeathModel__init__'
		mu = kwargs.get('mu', 1)
		super(PoissonBirthDeathModel, self).__init__(randprocess=poisson, randprocparam=mu, **kwargs)
			
class GeomBirthDeathModel(GenericDiscreetBirthDeathModel):
	"""The number of daughter lineages is sampled from a shifted geometric distribution.
	
	Mean number of daughter lineages is 1/(p-1) ; default parameter p=0.5 leads to distribution with mean of 1"""
	def __init__(self, **kwargs):
		print 'invoke _GeomBirthDeathModel__init__'
		p = kwargs.get('p', 0.5)
		super(GeomBirthDeathModel, self).__init__(randprocess=geometricm1, randprocparam=p, **kwargs)
		
#~ 
#~ class IndependentBirthDeathModel(SingleTreeModel):
	#~ """Following Nee et al. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. London Ser. B 344:305–11
	#~ 
	#~ Growth of each branch occur at each simulation step, following a continuous random process. 
	#~ Birth or death events then occur independently (as opposed to conjugated birth/death in Hey model, see below) on that branch following a random process given the branch length.
	#~ This can be with or without memory of the other simulation steps, i.e. given the total branch length or that added during the current iteration.
	#~ """
	#~ defaultrates = dict(growth=1, birth=1, death=1)
	#~ defaultprocesses = dict(growth=normalp1, birth=exponential, death=exponential)
	#~ 
	#~ def __init__(self, randprocess=defaultprocesses, rates=defaultrates, randomseed=None):
		#~ super(ContinuousBirthDeathModel, self).__init__(tunit=tunit, randomseed=randomseed)
		#~ self.randprocess = randprocess
		#~ self.rates = rates	# only support single paprameters
	#~ 
	#~ def stepforward(self, currtree, extincts, allowdeath=True, considerWholeBranch=False):
		#~ """Implement the atomic step of simulation iteration. 'extincts' is a list of leaf nodes in the tree"""
		#~ devents = dict(birth=[], death=[])
		#~ sleave = set(currtree.get_leaves())
		#~ # first grow branches
		#~ growbr = lambda: self.randprocess['growth'](self.rates['growth'])
		#~ for leaf in (sleave - set(extincts)):
			#~ newlen = growbr()
			#~ leaf += newlen	# add up on self.__l attribute
			#~ # then trigger births
			#~ if considerWholeBranch: birth = lambda: self.randprocess['birth'](self.rates['birth']/leaf.lg())
			#~ else:  birth = lambda: self.randprocess['birth'](self.rates['birth']/newlen)
			#~ nprogeny = birth()
			

class BaseMoranProcess(BaseModel):
	"""Base class for parenting Moran process classes in diamond with either SingleTreeModel or MultipleTreeModel classes"""

	def __init__(self, **kwargs):
		print 'invoke _BaseMoranProcess__init__'
		super(BaseMoranProcess, self).__init__(**kwargs)
		self.rate = kwargs.get('rate', 1)
		self.popsize = kwargs.get('popsize', 100)
		
	def stepforward(self, simul, allowdeath=True, timeLabelledNodes=True, evtidgen=None, **kwargs):
		"""Implement the atomic step of simulation iteration. 'extincts' is a list of leaf nodes in the tree
		
		'currtrees' argument can be either a single tree object when used within a PartialMoranProcess class instance, 
		or a list of tree objects when used within a  MoranProcess class instance.
		"""
		# NB: assumes timeslice numbering starts with 1; t=0 would induce null-rate at first step
		levents = []
		devents = dict(birth=[], death=[])
		extants = simul.get_extants()
		t = simul.t
		# first grow branches
		newlen = self.newlen(t)
		for leaf in extants:
			leaf += newlen	# add up on self.__l attribute
		# pick a branch among the N (hypothetical) branches for speciation
		ibirth = random.randint(0, self.popsize-1)
		# when modelling the full Moran process (MoranProcess class instance), ibirth should always be within range(len(extants))
		if ibirth in range(len(extants)):
			# birth occurs in one of the current tree's extant lineages
			bleaf = extants[ibirth]
			bleaf.edit_label("%s%d"%(nodelabelprefix['node'], t))
			BDevent('birth', bleaf, t, evtidgen, levents, devents)
			self.annotateNode(bleaf, 'birth')
			for i in range(2):
				newchild = bleaf.newnode()
				bleaf.link_child(newchild, newlen=0)	
				# !!! beware as last speciation will have two 0-length species;
				# should use GSA to deal with that (see Hartman 2012 http://sysbio.oxfordjournals.org/content/59/4/465.full)
		else:
			# should happen only for PartialMoranProcess instances
			BDevent('birth', self.dummynode, t, evtidgen, levents, devents)
		if allowdeath:
			# pick a branch among the N (hypothetical) branches for extinction, other that the one picked for speciation
			ideath = random.randint(0, self.popsize-1)
			while ideath == ibirth: ideath = random.randint(0, self.popsize-1)
			# when modelling the full Moran process (MoranProcess class instance), ideath should always be within range(len(extants))
			if ideath in range(len(extants)):
				# death occurs in ANOTHER one of the current tree's extant lineages
				dleaf = extants[ideath]
				dleaf.edit_label("%s%d"%(nodelabelprefix['deadtip'], t))
				# lineage goes extinct
				simul.extincts.append(dleaf)
				BDevent('death', dleaf, t, evtidgen, levents, devents)
				self.annotateNode(dleaf, 'death')
			else:
				# should happen only for PartialMoranProcess instances
				BDevent('death', self.dummynode, t, evtidgen, levents, devents)
		# print 'BaseMoranProcess.stepforward(): levents:', levents
		# print 'BaseMoranProcess.stepforward(): devents:', devents
		return (levents, devents, extants, newlen)

class MoranProcess(BaseMoranProcess, MultipleTreeModel):
	"""Follownig description of model C in Hey, J. 1992. Using Phylogenetic Trees to Study Speciation and Extinction. Evolution, 46(3), 1992, pp. 627-640
	
	Ideal model is that on the tree, time between two speciation/extinction (i.e. birth/death [B/D]) events is exponentially distibuted length (with rate param B). 
	This can be seens as a branch growing of an extra length (exponentially distibuted with rate param B), at the end of which a Birth event occur; 
	a Death event occurs simultaneously in another lineage. Dificult to implement in a simulation process as, focusing on a branch at each simulation time step, 
	the extinction (Death) event has to happen on a non-yet grown portion of another branch, which could have been speciating in the meantime. Rather do the following:
	
	According to Hey (1992), growth of all branches occur at each simulation step, of a length following a exponential decay process function of the time elapsed 
	(so added length is gradually shorter, accounting for the growing breadth of the tree, if one consider together theextant and extinct lineages). 
	Conjugated birth and death events then occur simultaneously on a randomly selected pair of branches. Assumes a population of size N original species parallely 
	evolving (resuting in several unconected trees, from which only one tree will eventually prevail, or be sampled.
	
	In this class, all trees from the original population are simulated. 
	
	A posteriori, seed species at t0 of the simulation may be considered distantly related (e.g. by a star tree with arbitrarily long branches), so that gene sequences 
	evolved on any of these distantly related species branches can be considered remote homologs within the same gene family; 
	or not at all, so that gene exchange between distantly related species can be the source of gene origination (thus assimilating gene origination process to a transfer process).
	"""
	
	def __init__(self, **kwargs):
		print 'invoke _MoranProcess__init__'
		super(MoranProcess, self).__init__(**kwargs)
		
	#~ def get_extants(self, currtrees, extincts, depthsorted=False):
		#~ """wrapper for method from MultipleTreeModel parent class, with making sure that the size of list is of expected size"""
		#~ extants = super(MoranProcess, MoranProcess).get_extants(currtrees, extincts, depthsorted=depthsorted)
		#~ assert len(extants) == self.popsize
		#~ return extants
		
	def newlen(self, t):
		# place holder; returns constant vaue
		return self.tunit

class PartialMoranProcess(BaseMoranProcess, SingleTreeModel):
	"""Follownig description of model C in Hey, J. 1992. Using Phylogenetic Trees to Study Speciation and Extinction. Evolution, 46(3), 1992, pp. 627-640
	
	Ideal model is that on the tree, time between two speciation/extinction (i.e. birth/death [B/D]) events is exponentially distibuted length (with rate param B). 
	This can be seens as a branch growing of an extra length (exponentially distibuted with rate param B), at the end of which a Birth event occur; 
	a Death event occurs simultaneously in another lineage. Dificult to implement in a simulation process as, focusing on a branch at each simulation time step, 
	the extinction (Death) event has to happen on a non-yet grown portion of another branch, which could have been speciating in the meantime. Rather do the following:
	
	According to Hey (1992), growth of all branches occur at each simulation step, of a length following a exponential decay process function of the time elapsed 
	(so added length is gradually shorter, accounting for the growing breadth of the tree, if one consider together theextant and extinct lineages). 
	Conjugated birth and death events then occur simultaneously on a randomly selected pair of branches. Assumes a population of size N original species parallely 
	evolving (resuting in several unconected trees, from which only one tree will eventually prevail, or be sampled.
	
	In this class, a single tree can be simulated (assuming it will be the one prevailing), allowing the B/D events to occur in lineages out of the tree, 
	i.e. only a fraction ni/N of events will occur on the tree, with ni the number of extant lineages at ti.
	
	!!! While this gives a tree equivalent to one sampled from a tree population from a Moran process, the simulations will differ in that 
	time slices from a PartialMoranProcess will have various (exponentionally distibuted) lengths, whereas time slices from a (full) MoranProcess will have constant length. 
	When using the simulated tree as a reference for a gene tree simulation, e.g. with BirthDeathDTLModel, this will impact the rate of DTL events per reference tree branches, 
	as events have constant rate per time slice.
	"""
	
	def __init__(self, **kwargs):
		print 'invoke _PartialMoranProcess__init__'
		super(PartialMoranProcess, self).__init__(**kwargs)
		
	def newlen(self, t):
		# growth at ti follow expontial law of parameter b*i*(i+1), with compound parameter b = B
		b = float(self.rate)/(self.popsize - 1)
		l = exponential(1/(b * t * (t+1))) * self.tunit
		return l
		
	# generic place holder nodes for filling up event record dictionaries
	dummynode = tree2.Node()
	dummynode.edit_label('out')

class DiscreetDTLModel(GenericDiscreetBirthDeathModel):
	pass

class BirthDeathDTLModel(MultipleTreeModel):
	"""Evolution along a reference tree (collection, as provided by a Moran process), marked by DTL events.
	
	DTL events happen at rates rdup, rtrans, rloss respectively, uniformly over time slices = node intervals.
	T event is emmited at rate rtrans, and recipient is chose among the N contemporaneous branches.
	Preserves extinct branches.
	"""
	def __init__(self, rdup=1e-4, rtrans=1e-4, rloss=1e-4, randomseed=None):
		super(BirthDeathDTLModel, self).__init__(randomseed=randomseed)
		if rloss+rdup+rtrans > 1:
			# assumes that only one kind of event can happen per branch per time slice
			# hence their summed probabilities must be lower than 1 (remainder account for speciation probability)
			raise ValueError, "sum of DTL event probabilities must sum to less than 1"
		self.rdup = rdup
		self.rtrans = rtrans
		self.rloss = rloss
		
	@staticmethod
	def midTimesliceOnBranch(node, timeslice):
		"""Given a reference tree branch and a tuple containing the (start, end) time coordinates of a timeslice, returns the height of the branch if it was to be shortened at the middle time of the timeslice"""
		lgabove = node.distance_root() - node.lg() # time above the branch
		lgonbranch = (sum(timeslice)/2) - lgabove	# height on the branch to match the midtime of the timeslice
		if not lgonbranch >= 0:
			raise ValueError, "lgabove = node.distance_root() - node.lg() = %f - %f\nlgonbranch = sum(%s)/2 - %f = %f"%(node.distance_root(), node.lg(), repr(timeslice), lgabove, lgonbranch)
		return lgonbranch
		
	@classmethod
	def lossEvent(cls, node, timeslice, silent=True):
		#~ if not silent: print 'LOSS event @ node %s'%node.label()
		print 'LOSS event @ node %s'%node.label()
		# get relevant branch lengths
		midt = cls.midTimesliceOnBranch(node, timeslice)
		# remove children
		for i in range(node.nb_children()):
			child = node.get_children()[0]	# operate as with a stack of children
			if not silent: print 'remove child #%d: %s'%(i, child.label())
			node.unlink_child(child, silent=silent)
			del child
		# set partial branch length
		node.set_lg(midt)
		cls.annotateNode(node, 'loss')
		
	@classmethod
	def duplicationEvent(cls, node, timeslice, silent=True):
		#~ if not silent: print 'DUPL event @ node %s'%node.label()
		print 'DUPL event @ node %s'%node.label()
		fat = node.go_father() # can be None if input node is the root
		# get relevant branch lengths
		try:
			midt = cls.midTimesliceOnBranch(node, timeslice)
		except ValueError, e:
			sys.stdout.flush()
			sys.stderr.write('node: %s\n'%node.label())
			sys.stderr.write('fat: %s; is root: %s\n'%(fat.label(), str(fat.is_root())))
			sys.stderr.write(str(node.go_root())+'\n')
			#~ node.go_root().seaview()
			tree2.dump_pickle(node.go_root(), "errortree.pickle")
			raise ValueError, e
		postmidt = node.lg() - midt
		newfat = node.newnode(l=midt)
		# choose to associate the duplicate node to node in reference tree (same as node above in gene tree) ;
		# only for consistency of having a .ref attribute in all gene tree nodes
		newfat.ref = node.ref
		if fat:
			# add duplication node below input node's father
			fat.link_child(newfat, newlen=midt, silent=silent)
			# disconect input node from its father
			fat.unlink_child(node, silent=silent)
		# relocate input node below duplication node and set partial branch length
		newfat.link_child(node, newlen=postmidt, silent=silent)
		# add duplicate sibling and set partial branch length
		dupli = node.deepcopybelow()
		newfat.link_child(dupli, newlen=postmidt, silent=silent)
		# NEW COPY node is annotated as duplicated
		cls.annotateNode(dupli, 'duplication')
		# new parent node is annotated as well
		cls.annotateNode(newfat, 'transfer')
		
	@classmethod
	def transferEvent(cls, donornode, recipientnode, timeslice, silent=True):
		#~ if not silent: print 'TRANS event from node %s to node %s'%(donornode.label(), recipientnode.label())
		print 'TRANS event from node %s to node %s'%(donornode.label(), recipientnode.label())
		fat = donornode.go_father() # can be None if input node (donornode) is the root
		# get relevant branch lengths
		try:
			dmidt = cls.midTimesliceOnBranch(donornode, timeslice)
			rmidt = cls.midTimesliceOnBranch(recipientnode, timeslice)
		except ValueError, e:
			sys.stdout.flush()
			sys.stderr.write('recipientnode: %s\n'%recipientnode.label())
			sys.stderr.write('donornode: %s\n'%donornode.label())
			sys.stderr.write('fat: %s; is root: %s\n'%(fat.label(), str(fat.is_root())))
			sys.stderr.write(str(donornode.go_root())+'\n')
			#~ donornode.go_root().seaview()
			raise ValueError, e
		dpostmidt = donornode.lg() - dmidt
		rpostmidt = recipientnode.lg() - rmidt
		# add transfer node below input node's father
		newfat = donornode.newnode(l=dmidt)
		# choose to associate the transfer node to donor node in reference tree (same as node above in gene tree) ;
		# only for consistency of having a .ref attribute in all gene tree nodes
		newfat.ref = donornode.ref
		if fat:
			# add transfer node below input node's father
			fat.link_child(newfat, newlen=dmidt, silent=silent)
			# disconect input node from its father
			fat.unlink_child(donornode, silent=silent)
		# relocate input node below transfer node and set partial branch length
		newfat.link_child(donornode, newlen=dpostmidt, silent=silent)
		# add transfer recipient and set partial branch length
		trans = recipientnode.deepcopybelow(add_ref_attr=True)
		newfat.link_child(trans, newlen=rpostmidt, silent=silent)
		# RECIPIENT node is annotated as transferred
		cls.annotateNode(trans, 'transfer')
		# new parent node is annotated as well
		cls.annotateNode(newfat, 'transfer')
			
	@staticmethod
	def getUniqueEventId(eventnode, eventtype, recipientnode=None):
		"""generate a code like exODT/ALE annotation of reconciled gene trees, cf. Szollosi et al. 2013, Lateral Gene Transfer from the Dead, Systematic Biology 62(3):386–397."""
		pass
		
	def stepforward(self, currbranches, currrefbranches, timeslice, t, evtidgen=None, **kwargs):
		"""generate event record that refer to the reference tree"""
		# NB: assumes timeslice numbering starts with 1; t=0 would induce null-rate at first step
		levents = []
		devents = {}
		trec = {}
		for cb in currbranches:
			# consider events exclusive, each proba is counted cummulatively so every draw in [0;1] can only point to one event type
			if random.random() <= self.rloss:
				evtype = 'loss'
				self.lossEvent(cb, timeslice)
			elif random.random() <= self.rloss+self.rdup:
				evtype = 'dupl'
				self.duplicationEvent(cb, timeslice)
			elif random.random() <= self.rloss+self.rdup+self.rtrans:
				evtype = 'trans'
				# pick a recipient branch from the current REFERENCE tree branches
				rec = random.choice(currrefbranches)
				self.transferEvent(cb, rec, timeslice)
			else:
				evtype = None
			if evtype: 
				e = DTLevent(evtype, cb, t, evtidgen, levents, devents, trec)
				# annoates the node's subtree labels by appending a string that signifies the event
				if e.recgenenode:
					if evtype=='trans':
						# add a pre-tag to salvage the info about the line of events leading to this node
						# which is not present on the label copied from the reference tree
						prelab = cb.label().split('-', 1)[-1]
					else:
						prelab = ""
					# annotate the transfer/duplication node
					cb.go_father().edit_label("%s%s%d"%(prelab, DTLevent.etshorts[evtype], e.id()), mode='a', sep="-")
					# and the recipient with all its descendants
					e.recgenenode.edit_all_labels("%s%s%d"%(prelab, DTLevent.etshorts[evtype], e.id()), mode='a', sep="-")
				else:
					# annotate the node with the event (only for losses)
					cb.edit_label("%s%d"%(DTLevent.etshorts[evtype], e.id()), mode='a', sep="-")
		return (levents, devents, trec)


###########################################################
## Evolutionary event records to be stored after simulation
###########################################################

# module's unbound generator function
def eventIdGen():
	n = 0
	while True:
		yield n
		n += 1
			
class BaseEvent(object):
	"""descriptor of the DTL event that occurred during the gene tree simulation"""
	# intitiates the class' unique identifier generator object using the module's generator function
	classidgen = eventIdGen()
	
	def __init__(self, eventtype, treenode, t, argidgen, levents=None, devents=None):
		self.eventtype = eventtype
		self.treenode = treenode				# the tree node at which the events occurs
		self.t = t								# the timeslice number == the simulation iteration number
		if argidgen: self.__id = next(argidgen)	# expects a generator object (should be property of the simulator instance)
		else: self.__id = next(classidgen)		# or relies on the class attribute (unbound to instances) 'classidgen'
												# (initialized once when the module is imported for all separate simulations,
												# and not continuable when loading pickled instances)
		if not (levents is None): levents.append(self)
		if not (devents is None): devents.setdefault(treenode.label(), []).append(self)
		# print 'BaseEvent.__init__(): levents:', levents
		# print 'BaseEvent.__init__(): devents:', devents
		
	def id(self):
		return self.__id		
			
class BDevent(BaseEvent):
	"""descriptor of the birth-death event that occurred during the (species) tree simulation"""
	
	evttypes = ['birth', 'death']
	
	def __init__(self, eventtype, treenode, t, argidgen, levents=None, devents=None):
		assert eventtype in self.evttypes
		super(BDevent, self).__init__(eventtype=eventtype, treenode=treenode, t=t, argidgen=argidgen, levents=levents, devents=devents)
			
class DTLevent(BaseEvent):
	"""descriptor of the duplication, transfer or loss event that occurred during the gene tree simulation"""
	
	evttypes = ['dupl', 'trans', 'loss']
	etshorts = {'dupl':'D', 'trans':'T', 'loss':'L'}
	
	def __init__(self, eventtype, dongenenode, t, argidgen, levents=None, devents=None, trec=None):
		assert eventtype in self.evttypes
		self.dongenenode = dongenenode			# the gene tree node from which the events departs (for T) or where it simply occurs (for D and L)
		self.donrefnode = dongenenode.ref		# the reference tree node from which the events departs (for T) or where it simply occurs (for D and L)
		self.recgenenode = dongenenode.go_brother() if eventtype in ['trans', 'dupl'] else None	# the reference tree node into which the events arrives (for T and D only)
		self.recrefnode = dongenenode.go_brother().ref if eventtype in ['trans', 'dupl'] else None	# the reference tree node into which the events arrives (for T and D only)
		# NB avoid dynamic querying of 'donrefnode' and 'recrefnode' as the returned value of dongenenode.ref and dongenenode.go_brother().ref might change later in the simulation
		super(DTLevent, self).__init__(eventtype=eventtype, treenode=dongenenode, t=t, argidgen=argidgen, levents=levents, devents=devents)
		if trec: trec.setdefault(self.recrefnode.label(), []).append(self)
		del self.treenode	# redundant with self.dongenenode
