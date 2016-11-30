#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Enriched phylogenetic tree objects, specifically representing a reconciled gene trees."""

__author__ = "Florent Lassalle <florent.lassalle@ucl.ac.uk>"
__date__ = "14 January 2015"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""

import copy
from math import log

import tree2

class GeneTree(tree2.AnnotatedNode):
	"""Enriched phylogenetic tree object, specifically representing a reconciled gene trees."""	
	def __init__(self, branch_lengths=True, keep_comments=False, score=None, **kw):
		"""Derived from AnnotatedNode (and thus Node); adds the event attribute describing coordinates of events in the reference species tree"""		
		super(GeneTree, self).__init__(branch_lengths=branch_lengths, keep_comments=keep_comments, **kw)	
		self.__event = (None, None, None)		# a 3-tuple containning the type of event and the parameters of this event, as described below, and the event DB identifier
		# if 'transfer', attribute is a 3-tuple: (receptor node_id in species tree, donor node_id in species tree, child node_id in gene tree). for Phylariane/Ancestrome reconcilated gene trees; !!! NB: the node annotations for transfer concern the children branches/nodes.
		# if 'speciation', corresponding node_id in species tree. for Phylariane/Ancestrome reconcilated gene trees.
		# if 'duplication', corresponding node_id in species tree. for Phylariane/Ancestrome reconcilated gene trees.
		# if 'uncertain', None
		self.__unicity = None
#		self.__reftree = reftree
		if score=='unicity':
			self.set_unicities_from_comments()

	def newnode(self, branch_lengths=True, keep_comments=False, **kw):
		"""class-specific instance generator for construction of trees of GeneTree nodes"""
		return GeneTree(branch_lengths=branch_lengths, keep_comments=keep_comments, **kw)
		
	def represented_species(self, excludeSelf=False):
		"""find the set of species represented under the node, exculding the transfered parts"""
		slleaves = set(self.get_leaf_labels())
		dtrans = self.getEvents(eventtype='transfer', returnDict=True)
		for translab in dtrans:
			if excludeSelf and translab==self.label(): continue
			transleaves = set(self.idgetnode(dtrans[translab]['eventlocation'][4]).get_leaf_labels())
			slleaves -= transleaves
		dl2sp = self.dictLeafLabelsToSpecies()
		slspe = set([dl2sp[l] for l in slleaves])
		return slspe
	
	def map_to_reftree_node(self, reftree, force=False):
		"""find the MRCA of species represented under the node, exculding the transfered parts"""
		slspe = self.represented_species()
		mrca = reftree.map_to_node(slspe, force=force)
		return(mrca)
		
	def set_transfer(self, reclab, donlab, childid, eventid=None, lrec=None, ldon=None, caureclab=None, caudonlab=None, reftree=None):
		"""defines a transfer event in a reconciled gene tree in reference to a node in species tree."""
		if reftree:
			llab = reftree.get_children_labels()
			if not (reclab in llab):
				raise IndexError, 'Unvalid receptor node label in species tree: %s'%reclab
			elif not (donlab in llab):
				raise IndexError, 'Unvalid donor node label in species tree: %s'%donlab
		children = self.get_children()
		for child in children:
			if child.nodeid()==childid:
				break
		else:
			raise IndexError, 'Unvalid child node id in gene tree: %d'%childid
		if not lrec:
			if reftree:
				if caureclab==None: crlab = reclab
				else: crlab = caureclab
				mrcrec = reftree[reclab]
				macrec = reftree[crlab]
				lr = mrcrec.path_to(macrec, returnLabels=True)
			else:
				lr = []
		else:
			lr = lrec
		if not ldon:
			if reftree:
				if caudonlab==None: cdlab = donlab	
				else: cdlab = caudonlab
				mrcdon = reftree[donlab]
				macdon = reftree[cdlab]
				largeset = set( macdon.get_children_labels() )
				smallset = set( mrcdon.get_children_labels() )
				donSet = (largeset - smallset) | set([mrcdon.label()])
				ld = list(donSet)			
			else:
				ld = []
		else:
			ld = ldon
		self.__event = ('transfer', (reclab, lr, donlab, ld, childid), eventid)		
		
	def set_singlenodeevent(self, eventtype, eventid=None, locspe=-1, lnodes=None, caulocspe=None, reftree=None, checkfortransfers=False):
		"""defines any event with one coordinate (i.e. those of type duplications, speciaton or uncertain) in a reconciled gene tree in reference to a node in species tree."""
		if locspe!=-1:
			if reftree:
				llab = reftree.get_children_labels()
				for spe in [locspe, caulocspe]:
					if spe!=None and (not (spe in llab)): raise IndexError, 'Unvalid node label in species tree: %s'%spe
				if not lnodes:
					mrca = reftree[locspe]
					if not caulocspe: maca = mrca
					else: maca = reftree[caulocspe]
					l = mrca.path_to(maca, returnLabels=True)
				else: 
					l = lnodes
			else:
				l = lnodes
			loclab = locspe
		elif reftree:
			if not lnodes:
				if checkfortransfers:
					repspe = self.represented_species()
				else:
					repspe = set(self.listSpecies())
				mrca = reftree.map_to_node(repspe)
				if eventtype=='speciation':
					l = [mrca.label()]
				else:
					# find the species represented under the node
					if self.go_father(): 
						# and contrast it with the set of species represented in the full tree to get the MACA
						repfull = set(self.go_root().listSpecies())
						outboth = repfull - repspe
					else:
						outboth = []
					maca = reftree.map_to_ancient_node(repspe, outboth)
					if maca:
						l = mrca.path_to(maca, returnLabels=True)
					else:
						l = [mrca.label()]
			else: 
				l = lnodes
			loclab = mrca.label()
		else:
			# if no information, locates the event at the root of the species tree
			if not lnodes: l = []
			else: l = lnodes
			loclab = 'N1'
		self.__event = (eventtype, (loclab, l), eventid)
		
	def set_duplication(self, locspe=-1, eventid=None, lnodes=None, caulocspe=None, reftree=None, checkfortransfers=False):
		"""defines a duplication event in a reconciled gene tree in reference to a node in species tree."""
		self.set_singlenodeevent('duplication', locspe=locspe, lnodes=lnodes, eventid=eventid, caulocspe=caulocspe, reftree=reftree, checkfortransfers=checkfortransfers)
		
	def set_speciation(self, locspe=-1, eventid=None, lnodes=None, caulocspe=None, reftree=None, checkfortransfers=False):
		"""defines a speciation event in a reconciled gene tree in reference to a node in species tree.
		
		!!! must be done after reconciliation of transfers and/or duplications.
		"""
		self.set_singlenodeevent('speciation', locspe=locspe, lnodes=lnodes, eventid=eventid, caulocspe=caulocspe, reftree=reftree, checkfortransfers=checkfortransfers)
			
	def set_gain(self, locspe=-1, eventid=None, lnodes=None, caulocspe=None, reftree=None, checkfortransfers=False):
		"""defines a speciation event in a reconciled gene tree in reference to a node in species tree.
		
		!!! must be done after reconciliation of transfers and/or duplications.
		"""
		self.set_singlenodeevent('gain', locspe=locspe, lnodes=lnodes, eventid=eventid, caulocspe=caulocspe, reftree=reftree, checkfortransfers=checkfortransfers)
			
	def set_uncertainevent(self, locspe=-1, eventid=None, lnodes=None, caulocspe=None, reftree=None, checkfortransfers=False):
		"""defines a uncertain event (either a duplication or a transfer) in a reconciled gene tree in reference to a node in species tree.	"""
		self.set_singlenodeevent('uncertain', locspe=locspe, lnodes=lnodes, eventid=eventid, caulocspe=caulocspe, reftree=reftree, checkfortransfers=checkfortransfers)
		
	def set_anyevent(self, eventtype, eventloc, eventid):
		self.__event = (eventtype, eventloc, eventid)
		
	def set_eventid(self, eventid):
		eventtype, eventloc, oldeventid = self.__event
		self.__event = (eventtype, eventloc, eventid)
		
	def set_eventtype(self, eventtype):
		oldeventtype, eventloc, eventid = self.__event
		self.__event = (eventtype, eventloc, eventid)
		
	def set_eventloc(self, eventloc, item=-1):
		eventtype, oldeventloc, eventid = self.__event
		if item < 0: neweventloc = eventloc
		else: neweventloc = oldeventloc[0:item]+(eventloc,)+oldeventloc[item+1:]
		self.__event = (eventtype, neweventloc, eventid)
		
	def clean_event(self):
		self.__event = (None, None, None)
		
	def set_unicity(self, unicity):
		self.__unicity = unicity
		
	def set_unicities_from_comments(self):
		for node in self.get_all_children():
			node.set_unicity(float(node.comment()))
	
	def computeUnicity(self, logUnicity=True, dlabspe=None, setValue=True):
		"""compute unicity of species set at leaves under the node."""
		if not dlabspe:
			dlabspe = self.dictLeafLabelsToSpecies()
		lspe = dlabspe.values()
		dspecount = {}
		for spe in lspe:
			dspecount[spe] = dspecount.setdefault(spe, 0) + 1
		unicity = reduce(lambda x,y: x*y, dspecount.values(), long(1))
		#~ unicity = long(1)
		#~ for count in dspecount.values():
			#~ unicity *= count
		if logUnicity: u = log(unicity)
		else: u = unicity
		if setValue: self.__unicity = u
		else: return u
		
	def combineUnicities(self, other, logUnicity=True, dlabspe=None):
		"""compute the unicity of the merged sets of leaves of self node and another node"""
		dnewson = self.dictLeafLabelsToSpecies()
		dother = other.dictLeafLabelsToSpecies()
		dother.update(dnewson)
		unicity = self.computeUnicity(dlabspe=dother, setValue=False, logUnicity=True)
		return unicity
			
	def computeDupliConsistencyScore(self, dlabspe=None):
		"""returns the node Duplication Consistency Score as in TreeBest (Villela et al., Genome Res. 2009)"""
		if not dlabspe:
			dlabspe = self.dictLeafLabelsToSpecies()
		cintersect = None
		cunion = set()
		for c in self.get_children():
			lspe = []
			for llab in c.get_leaf_labels():
				lspe.append(dlabspe[llab])
			if cintersect==None: cintersect = set(lspe)
			else : cintersect = cintersect & set(lspe)
			cunion = cunion | set(lspe)
			#~ print lspe, cintersect, cunion
		return(float(len(cintersect)) / float(len(cunion)) )
		
		
	def computeDuplicationLossRate(self, reftree, mode='pernode', normlen='treelen', verbose=False):
		"""compute the rate of loss necessary to explain a duplication at the node with a Dollo parsimony
		
		can divide the loss count by the node count or total branch length under the duplication event location
		to get a rate per node ('pernode') or per branch length units ('perlength') by 
		or divide only by the sum of ortholog counts at nodes under the event, 
		to get a rate per lineage that had born the gene ('perlineage').
		"""
		if not isinstance(reftree, tree2.ReferenceTree): raise TypeError
		rt = copy.deepcopy(reftree)
		anc = rt.map_to_node(self.listSpecies())
		for c in self.get_children():
			lspe = c.listSpecies(ignoreTransfers=True)
			srt = copy.deepcopy(reftree)
			srt.presenceInTree(lspe)
			srt.DolloPars(inventorlab=anc.label())
			rt += srt
			#~ nloss.append(float(srt.sumEvents('loss')))
		nloss = float(srt.sumEvents('loss'))
		if mode=='pernode':
			nnodes = float(len(anc.get_all_children()))
			if verbose: print nloss, 'loss /', nnodes, 'nodes under', anc.label(), '=',
			return ( nloss / nnodes )
		elif mode=='perlength':
			underlen = anc.treelength()
			if normlen=='treelen': normlen = reftree.treelength()
			else: normlen = float(normlen)
			return ( nloss / (underlen / normlen) )
		elif mode=='perlineage':
			nhomo = float(anc.sumHomologCounts())
			if verbose: print nloss, 'loss /', nhomo, 'homologs under', anc.label(), '=',
			return ( nloss / nhomo )
		else:
			raise ValueError, "mode %s not valid"%mode
			
		
	def increasesUnicity(self, newfat):
		"""States if branching the node 'newson' under the node 'self' increases the conflict = unicicty score at the node 'self'"""
		if self.combineUnicities(newfat) > newfat.unicity(): return True
		else: return False
		
	def event(self):
		return self.__event
	
	def getdicevent(self):
		d = {}
		d['eventtype'] = self.__event[0]
		d['eventlocation'] = self.__event[1]
		d['eventid'] = self.__event[2]
		return d
			
	def eventtype(self):
		return self.__event[0]
		
	def eventloc(self):
		return self.__event[1]
		
	def eventid(self):
		return self.__event[2]
			
	def duplication(self):
		if self.__event[0]=='duplication':
			return self.eventloc()
		
	def speciation(self):
		if self.__event[0]=='speciation':
			return self.eventloc()
			
	def transfer(self):
		if self.__event[0]=='transfer':
			return self.eventloc()
			
	def transferchild(self):
		if self.__event[0]=='transfer':
			return self.eventloc()[4]
			
	def getEvents(self, eventtype=None, returnDict=False, lineage=None, returnKey=None):
		"""search for events in gene tree, and return a dictionary for each, stored as a list or as a dictionary.
		
		[ can be search in the whole tree (if 'lineage' is None) or on the path from node named 'lineage' to the root.
		Can be restricted to the 'maxeventcount' more recent events in each category (if 0, no max), and to a specific category 'eventtype' (if None, no specification).
		"""
		def checkevent(node, events):
			e = node.eventtype()
			if (not eventtype and e) or (eventtype and e==eventtype):
				d = node.getdicevent()
				#~ if returnKey: d = d.get(returnKey)
				if returnDict: events[node.label()] = d
				else: events.append(d)
				
		if returnDict:
			self.complete_internal_labels()
			events = {}
		else:
			events = []	
		if not lineage:
			for node in self.get_sorted_children(order=-1):
				checkevent(node, events)
		else:
			f = self[lineage]
			while f:
				checkevent(f, events)
				f = f.go_father()
		return events
		
	def cleanEvents(self, eventtype=None, lineage=None):
		events = self.getEvents(eventtype=eventtype, returnDict=True, lineage=None)
		for nodelab in events:
			self[nodelab].clean_event()
	
	def unicity(self):
		return self.__unicity
		
	def get_all_ortholog_groups(self, duplabels):
		"""returns the maximum subtrees of orthologs under a duplication node as a list
		
		or as a dictionary which keys are the sequence of duplication node labels from the root to the subtree 
		e.g.: d['.1.3.56.124'] = c if c is the node under the duplication nodes 'N1', 'N3', and 'N56' and is labeled 'N124'.
		"""
		lortho = []
		sdupli = set(duplabels)
		if not (set(self.get_children_labels()) & sdupli):
			# no duplication node at the node: it is a group of orthologs
			lortho.append(self)
		else:
			for c in self.get_children():
				# recursively searches in shalower nodes
				lortho += c.get_all_ortholog_groups(duplabels)
		return lortho
		
	def pruneTransferedSubtrees(self, intactCopy=None, silent=True):
		"""Prune all the subtrees that correspond to transfer recipients. return forest of subtrees"""
		
		def labFromIntactCopy(childid, intactCopy):
			"""when working on a tree already pruned of some leaf compared to the full version where transfered children were annotated"""
			if intactCopy:
				# use an intact copy as a map
				intactchild = intactCopy.idgetnode(childid)
				c = self.map_to_node(intactchild.get_leaf_labels(), force=True)
				return c.label()
			else:
				raise IndexError, "cannot find node with id %d"%childid
		
#		self.seaview(ignoreBS=1)
		dcolapsed = {} 		# d[removed] = staying
		didlab = {}
		for node in self:
			didlab[node.nodeid()] = node.label()
		if not silent: print "\ndidlab", didlab, "\n---"
		ltrans = self.getEvents(eventtype='transfer')
		#~ dtrans = self.getEvents(eventtype='transfer', returnDict=True)
		forest = []
		for trans in ltrans:
			if not silent: print "trans", trans
			#~ trans = dtrans[translab]
			#~ if not silent: print "translab", translab, ":", trans
			childid = trans['eventlocation'][-1]
			child = self.idgetnode(childid)
			if not child:
				#~ childlab = didlab.get(childid, labFromIntactCopy(childid, intactCopy))
				if intactCopy: childlab = didlab.get(childid, labFromIntactCopy(childid, intactCopy))
				else: childlab = didlab[childid]
				if not silent: print "lost childid:", childid, "- retracing child:"
				while childlab in dcolapsed:
					if not silent: print childlab, "->",
					# trace back removed nodes to staying nodes
					childlab = dcolapsed[childlab]
				if not silent: print ""
				child = self[childlab]
				
#			if not child: self.seaview(ignoreBS=1)
			# record track of node to be removed
			if not silent: print "to pop:", child.label()
			collapsednodes = self.pop(child, tellReplacingNode=True)
			dcolapsed[collapsednodes[0]] = collapsednodes[1]
			# actually removes a node by pruning its son
			tst = self.pop(child)	
			if not silent: print "poped:", tst.newick(ignoreBS=1), "\nremains", self.newick(ignoreBS=1), "\ndcolapsed:", dcolapsed, "\n---"
			forest.append(tst)
		return forest
	
	def getOrthoSubtrees(self, reftree=None, silent=True):
		"""retuns a dictionary of gene tree nodes with transfer/duplication events -> their children subtrees of orthologs = orthologous subfamily
		
		if a reference tree is provided, can detect simple gene conversion between close relative, avoiding to split the subfamily
		"""
		genetree = copy.deepcopy(self)
		# dictionary of gene tree nodes with transfer/duplication events -> their children subtrees of orthologs to treat separately
		orthoSubtrees = {}
		# prepare the sequence of nodes to prune: store their leaf list
		levtnodes = []
		devtlabostleaves = {}
		for node in self.get_sorted_children(order=3):
			# chose which child node was the recipient of the duplication/transfer event = the one to pop
			child = None
			if node.duplication():
				# make the hypothesis that duplication recipient is smaller
				child1, child2 = node.get_children()
				if child1.nb_leaves() <= child2.nb_leaves() : child = child1
				else : child = child2
			elif node.transfer():
				child = self.idgetnode(node.transfer()[4])
			else:
				continue
			levtnodes.append(node)
			devtlabostleaves[node.label()] = child.get_leaf_labels()
		
		# sequentially prune the ortholous subtrees indicated by the leaf sets
		for evtnode in levtnodes:
			#~ print "evtnode", evtnode.label()
			#~ print "\tdevtlabostleaves[%s]"%evtnode.label(), devtlabostleaves[evtnode.label()]
			leaf2pop = set(genetree.get_leaf_labels()) & set(devtlabostleaves[evtnode.label()])
			#~ print "\tleaf2pop", leaf2pop
			gtevtnode = genetree.map_to_node(leaf2pop)
			fevtnode = gtevtnode.go_father()
			if reftree and evtnode.transfer() and fevtnode:
				# find the next node above being the root or under a duplication 
				underdup = fevtnode
				funderdup = underdup.go_father()
				while funderdup and funderdup.eventtype()!='duplication':
					underdup = funderdup
					funderdup = underdup.go_father()
				# restrict reftree to sppecies present in this subtree
				rt = reftree.restrictToLeaves(underdup.listSpecies())
				#~ print  "\tgtevtnode", gtevtnode.label()
				# verify that it is not an allelic replacement
				gfevtnode = fevtnode.go_father()
				while gfevtnode and gfevtnode.eventtype()!='duplication':
					#~ print "\tgfevtnode", gfevtnode.label(), gfevtnode.listSpecies()
					gflspe = gfevtnode.listSpecies(asSet=False)
					#~ if rt.is_monophyletic(gflspe):
					if len(set(gflspe))==len(gflspe) and rt.is_monophyletic(gflspe):
						# transfer does not induce duplication nor brings foreign species that would make heterogeneous phylogenetic profile
						# this transfer is rather a simple gene conversion between close relative ; do not separate this subtree from what is above (same subfamily)
						if not silent: print "allelic replacement at", fevtnode.label(), fevtnode.event()
						break # the while gfevtnode loop = ommit to split fam at this event node
					gfevtnode = gfevtnode.go_father()
				else:
					ost = genetree.pop(gtevtnode)
					orthoSubtrees[evtnode] = ost
			else:
				ost = genetree.pop(gtevtnode)
				orthoSubtrees[evtnode] = ost
		else:
			# subtree of the "primary" orthologous lineage initiated by an origination event
			if not self in orthoSubtrees:
				orthoSubtrees[self] = genetree
			else:
				# already an event (must be a duplication) recorded at the root ; must create a double recording the original speciation
				#~ selfspe = copy.deepcopy(self)
				selfspe = copy.copy(self)	# use the Node.__copy__ method
				evloc, evnodes = self.getdicevent()['eventlocation']
				selfspe.set_speciation(locspe=evloc, lnodes=evnodes)
				orthoSubtrees[selfspe] = genetree
		if not silent: print "orthoSubtrees:\t\n", "\t\n".join(["%s: %s"%(eventnode.label(), str(orthoSubtrees[eventnode].get_leaf_labels())) for eventnode in orthoSubtrees])
		return orthoSubtrees
		
			
	def dictUpEventNodes(self, leventnodes):
		"""make a dictionary of each event node to the closest event above in the tree"""
		d = {}
		for i, eventnode in enumerate(leventnodes):
			ort = orthoReftrees[eventnode]
			upeventnode = None
			# searches the event above
			for aftereventnode in leventnodes[i+1:]:
				# if eventnode.is_child(aftereventnode):	# replaced so that the potential shalow copy of the root knows its children
				if eventnode in aftereventnode.get_all_children():
					if aftereventnode.duplication() and aftereventnode==aftereventnode.go_root():
						# special case of the gene tree root when it is a duplication : two aftereventnodes are recorded at the root; must choose the right one
						ost = orthoSubtrees[eventnode]
						upost = orthoSubtrees[aftereventnode]
						# we want the aftereventnode which orthoSubtree is in the same first bi-partition of the gene tree as the orthoSubtree of eventnode (MRCA lower than the root)
						if aftereventnode.map_to_node(ost.get_leaf_labels()+upost.get_leaf_labels())==aftereventnode:
							# avoid the one which MRCA of both leaf set is the root
							continue
					upeventnode = aftereventnode
					break
			d[eventnode] = upeventnode
		return d

		
########################################
### Reference tree mapping based methods:

	def branchingDepth(self, reftree, dlabspe=None):
		"""returns the depth in a reference tree of the MRCA of species represented in the subtree (self) """
		if not dlabspe:
			dlabspe = self.dictLeafLabelsToSpecies()
		lspe = dlabspe.values()
		taxset = set(lspe)
		mrca = reftree.map_to_node(taxset)
		return mrca.depth()

	def combineOrthologs(self, duplabels, onlyOrthologs=False, silent=True):	#, dcombinost={} ; this dictionary can be provided pre-filled to prevent integration of particular subtrees
		"""finds several combinations of unicopy sets of leaves starting from a group of orthologs
		
		return a dictionary 'dcombinost' of explored of orthologous subtree labels to leaf sets
		'onlyOrthologs=True' states that nothing else than strictly orthologous genes can be sampled; when 'False', allow all genes that 
		are not explicitely paralogous (for instance, allow xenologs from transfer events).		
		"""
		if self.unicity() > 0:
			raise IndexError, "the node %s is not unicopy:\n%s"%(self.label(), self.newick(ignoreBS=True))
		dcombinost={}
		if not silent: print "input ost", self.label()
		# dictionary of leaf labels to species identifiers in the full gene tree
		dlabspe = self.go_root().dictLeafLabelsToSpecies()
		# list of the explored duplication nodes out of the input subtree
		exploreddup = ['']
		# input species/leaf label sets
		lspe = self.listSpecies()
		lleaves = self.get_leaf_labels()
		dcombinost[(self.hierachical_path(labels=exploreddup),)] = (lleaves, lspe)
		if not silent: print "dcombinost\n", dcombinost
		n = self
		f = self.go_father()
		while f:
			# identifies the kind of homology relationship
			if not (f.label() in duplabels):
				# no conflict at 'f'
				case='orthologous'
				if not silent: print "f ortho", f.label()
			else:
				if (not onlyOrthologs) and f.duplication()==None:
					# there is some conflict at 'f', but it is not explained by a duplication (more probably by a transfer)
					case='homologous with conflict'
					if not silent: print "f homo wc", f.label()
				else:
					# 'f' is a duplication node above the original group of orthologs: 
					# skips the duplication node, avoiding any kind of paralog
					n = f
					f = n.go_father()
					continue # the while f loop
			# selects additional leaves respecting the unicopy constraint
			for c in f.get_children():
				if c != n:
					if not silent: print "c: %s"%c.label()
					cost = c.get_all_ortholog_groups(duplabels)
					lusedtk = []
					for ost in cost:
						path = ost.path_to(f)
						for p in path:
							if p.label() in duplabels:
								if p.duplication():
									if p.label() not in exploreddup:
										exploreddup.append(p.label())
								else:
									ocase = 'homologous with conflict'
						if not silent: print "exploreddup", exploreddup
						# each ost is tagged with a string concatenate of labels of duplication nodes met since 'f' (turning point in tree)
						kost = ost.hierachical_path(labels=exploreddup)
						if not silent: print "ost: %s, kost: %s"%(ost.label(), kost)
						# higher explored duplication node above ost	
						radical = kost.split('.',2)[1]
						tempd = {}
						for tk in dcombinost:
							# find previous ost combinations where no ost is from this lineage
							for k in tk:
								if not silent: print "k", k
								if k.split('.',2)[1]==radical:
									break
							else:
								if not silent: print "tk", tk
								ll, ls = copy.deepcopy(dcombinost[tk])
								if not silent: print "ll", ll
								# adds leaves from the input orthologous subtree
								templs = []
								templl = []
								for leaf in ost.get_leaf_labels():
									if not dlabspe[leaf] in ls:
										# adds only the leaves of unrepresented species
										templl.append(leaf)
										templs.append(dlabspe[leaf])
									else:
										# do not add leaves from an ost that contains species already present in unicopy combination
										break
									#~ else:
										#~ if case=='orthologous':
											#~ print "\nbad", leaf
											#~ # species set of 'ost' is (should be) disjoint from that of current unicopy subtree
											#~ print "bad inference of orthology for node %s under non-duplication node %s"%(ost.label(), f.label())
											#~ raise ValueError
										# elif case=='homologous with conflict':
											# species set of 'ost' is not disjoint from that of current unicopy subtree, but their leaves are not paralogous: no problem
								else:
									if not silent: print templl
									ll += templl
									ls += templs
									# if ost exploration was not vain, add it to the dictionary
									if not tk in lusedtk: lusedtk.append(tk)
									# index of combined ost
									tkost = tk + (kost,)
									tempd[tkost] = (ll, ls)
						else: # for tk loop, else clause
							if not silent: 
								print "added ust:"
								for tkost in tempd:
									print tkost, tempd[tkost]
							dcombinost.update(tempd)
					else: # for ost loop, else clause
						# purges the the previous combinations that were declined, 
						# in order to combine orthologs gradually distant from input ost 
						# (avoid combining close homologs with distant ones while missing interleaving ones)
						for tk in lusedtk:
							del dcombinost[tk]
			n = f
			f = n.go_father()
			if not silent: print "dcombinost.keys()", dcombinost.keys(), "\n"
		return dcombinost

	def combinationsUnicopySubtrees(self, duplabels, onlyOrthologs=False, silent=True):
		dust = {}
		dco = self.combineOrthologs(duplabels, onlyOrthologs=onlyOrthologs, silent=silent)
		for co in dco:
			lleaves, lse = dco[co]
			genetree = copy.deepcopy(self.go_root())
			ust = genetree.restrictToLeaves(lleaves)
			dust[''.join(co)] = ust
		return dust

			
	def whichMulticopySubtrees(self, focusedSon=None, multi=[], returnLabels=False):
		"""finds the groups of duplicated leaves under a putative duplication node.
		
		return the MRCAs of these groups, with the first being the one under the focused son.
		"""
		def getAnc(para, multi, ancs, returnLabels):
			llab = para.dictSpeciesToLeafLabels(lspe=multi, catValues=True)
			anc = para.map_to_node(llab)	# !!!	cannot distinguish several groups of duplicates : must include combinatorics of multicopy species
			if anc:
				if returnLabels: ancs.append(anc.label())
				else: ancs.append(anc)
				
		if not focusedSon:
			focusedSon = self.get_children()[0]
		elif not focusedSon in self.get_children():
			raise IndexError, "%s is not direct child of %s"%(focusedSon.label(), self.label())
		if not multi:
			# finds the set of species that are multicopy under duplication node
			multi = set(focusedSon.listSpecies())
			for child in self.get_children():
				if child != focusedSon:
					multi &= set(child.listSpecies())
			if not multi:
				raise ValueError, "no duplicated species found under this node: %s"%focusedSon.label()
		#~ print "multi", multi
		# finds common ancestors of duplicated leaves in each children of the node
		ancs = []
		getAnc(focusedSon, multi, ancs, returnLabels)
		for para in self.get_children():
			if para != focusedSon:
				getAnc(para, multi, ancs, returnLabels)
		return(ancs)

		
	def moreParsimoniousScenario(self, child, minbs, originalnodelabs=None, silent=True):
		"""tries to avoid inference of duplications by performing moves (SPRs) of non-supported branches"""
		# determines the sources of conflict
		ancs = self.whichMulticopySubtrees(focusedSon=child)
		## sould try here to compact the groups of paralogs, i.e. put outside of these groups non-supported branches leading to non-duplicated leaves
		# is the position of the focused paralogous group under the duplicated node supported? 
		# finds the smallest supported clade including the focused paralogous group under the duplication node
		focuspara = ancs[0]
		fpath = focuspara.path_to(self)[:-1]	# exculdes the duplication node = 'self'
		for j in range(len(fpath)):
			# goes backward (toward root) to detect support
			if fpath[j].bs() >= minbs:
				focuspara = fpath[j]
		focuspara.computeUnicity()
		focuspara.go_father().computeUnicity()
		if not silent: print "focuspara: %s[%f]"%(focuspara.label(), focuspara.unicity())	#.newick(ignoreBS=True)
		# then, check support on the path from the duplication node to the other paralogous group(s)
		for anc in ancs[1:]:
			lbs = [0] + self.bootstraps_on_path_to(anc) # the duplication node = 'self' will not be crossed: asociates 0 to it to pass the while clause
			path = self.path_to(anc)
			newbro = None
			i = 0	# starts at the duplication node = 'self'
			# goes forward (toward leaves) on path, while not crossing high branch support
			while (i < len(path)-1) and (lbs[i] < minbs):
				## detects other duplication nodes in order to try merge them  (excludes the original duplication node = 'self')
				if (i > 0) & (path[i].findDuplications(minbs, recursive=False, state=True, modifyTopology=False)):
					# find the 'newbro' = the son of the new duplication node ('path[i]') which is not on the 'path' lineage
					underdup = path[i].get_children()
					for udp in underdup:
						if udp != path[i+1]:
							nb = udp
							break	# the for udp loop
					nb.computeUnicity()
					if not silent: print "nb: %s[%f]"%(nb.label(), nb.unicity())	#.newick(ignoreBS=True)
					# verifies that branching 'focuspara' with 'nb' does not increase the conflict under it and will not do a trivial move (graft to its brother)
					if (not focuspara.increasesUnicity(nb)) and (not focuspara.is_brother(nb)):
						# test if the future node above the pair 'focuspara' and 'nb' will have a higher unicity score then the current node above 'focuspara'
						if not silent: print "predicted new node unicity: %f, current unicity %f"%(focuspara.combineUnicities(newbro), focuspara.go_father().unicity())
						if focuspara.combineUnicities(nb) < focuspara.go_father().unicity():
							# displace 'focuspara' to the position next to 'nb'	
							return focuspara.doSPR(nb, originalnodelabs=originalnodelabs)
					#~ # !!! elif clause may never be met as it seem equivalent to if clause; ??(if clause is prefered as it moves the node toward the leaves)
					#~ # or that branching 'nb' with 'focuspara' does not increase the conflict under it
					#~ elif not nb.increasesUnicity(focuspara):
						#~ print 
						#~ # displace 'newbro' to the position next to 'focuspara'
						#~ return nb.doSPR(focuspara, originalnodelabs=originalnodelabs)
				## if no move was performed, goes to the next node
				newbro = path[i+1]
				i += 1
			else:
				if newbro : 
					newbro.computeUnicity()
					if not silent: print "newbro: %s[%f]"%(newbro.label(), newbro.unicity())	#.newick(ignoreBS=True)
					# verifies that branching 'focuspara' next to 'newbro' does not increase the conflict at its father
					if (not focuspara.increasesUnicity(newbro.go_father())) and (not focuspara.is_brother(newbro)):
						# test if the future node above the pair 'focuspara' and 'newbro' will have a higher unicity score then the current node above 'focuspara'
						if not silent: print "predicted new node unicity: %f, current unicity %f"%(focuspara.combineUnicities(newbro), focuspara.go_father().unicity())
						if focuspara.combineUnicities(newbro) < focuspara.go_father().unicity():
							# displace the focuspara to the position next to 'newbro'	
							return focuspara.doSPR(newbro, originalnodelabs=originalnodelabs, silent=silent)
						
				#ajouter : clause sur unicité combinée qui ne doit pas augmenter lors du NNI (prédire la valeur avec combineUnicities() )
			
	def findDuplications(self, minbs, originalnodelabs=None, recursive=True, state=False, modifyTopology=True, verbose=False):
		"""finds unicity leaps from father to sons and returns the list of duplicated nodes and edits the tree by moving branches toward a more parsimonious scenario
		
		'minbs' specifies the maximum (bootstrap) support at a node to allow it being crossed by SPR moves
		'originalnodelabs' is the set of labels used in the tree prior to modification; providing it prevent re-use of labels of deleted nodes
		recursive toward leaves if 'recursive' enabled (default)
		with 'state' option true, just returns boolean stating if the node is a duplication or not (do not intent topology changes)
		"""
#		if not state: print "\nself: %s"%self.label()
		self.computeUnicity()
		if verbose: dbg = True #self.label() in ['N3']
		else: dbg = False
		ldupli = []
		changedTopo = False
		if dbg and not state: print "self: %s\n"%self.label(), self.newick(ignoreBS=True)
		if self.is_leaf():
			if state: return False
			else: return (ldupli, changedTopo)
		else:
			# records original context of the node
			genetree = self.go_root()
			children = self.get_children()
			# computes the intersection of set of species present under each node
			sspe = set(children[0].listSpecies())
			for child in children[1:]:
				sspe &= set(child.listSpecies())
			if not sspe:
				## intersection of child sets is empty: no conflict to resolve at this node
				if dbg: print 'dbg 0'
				if state:
					return False
				elif recursive:
					for childlab in self.children_labels():
						# must keep track of the children to inspect (father-to-child dependencies may change during the inspection of the first son, label of the second son remain stable as the two son subtrees are independent)
						if dbg: print "childlab:", childlab
						l, c = genetree[childlab].findDuplications(minbs, originalnodelabs=originalnodelabs, modifyTopology=modifyTopology, recursive=recursive, verbose=verbose)	
						ldupli += l
						if c: changedTopo = True
			else:	
				## intersection of child sets is not empty: conflict to resolve at this node (Scornavacca et al. 2009, Def. 1)
				if state:
					return True
				else:
					# records original label of the node
					originallab = self.label()
					if not originalnodelabs: originalnodelabs = genetree.get_children_labels()
					if dbg: print 'dbg 1'
					for child in children:
						# tries to improve topology under the node into a more parsimonious scenario where this node is not a duplication
						if modifyTopology:
							delreplnew = self.moreParsimoniousScenario(child, minbs, originalnodelabs=originalnodelabs, silent=(dbg-1))
							if dbg: print delreplnew
						else:
							delreplnew = None
						if delreplnew:
							changedTopo = True
							(deletednode, replacementnode, newnode) = delreplnew
							if dbg: print "delrepnew:", deletednode, replacementnode, newnode
							# updates the 'originalnodelabs' list
							originalnodelabs.append(newnode)
#							print "deletednode: %s, replacementnode: %s, newnode: %s"%(deletednode, replacementnode, newnode)
							# refreshes the list of duplicated nodes (should never happen)
							if genetree[deletednode] in ldupli:
								print "HAZARDOUS\ndeleted dupli:", deletednode
								deli = ldupli.index(genetree[deletednode])
								ldupli[deli] = genetree[replacementnode]
							## new nodes must have labels added: 'self' or 'child' nodes may have been deleted from the gene tree by SPR moves
							if deletednode==originallab: 
								rn = genetree[replacementnode]
								nn = genetree[newnode]
								rn.computeUnicity()
								nn.computeUnicity()
								rnd = rn.depth()
								nnd = nn.depth()
								if rnd <= nnd: node = rn
								else: node = nn							
							else: 
								node = self
#							print "node:\n", node.newick(ignoreBS=True)
							# duplication has been resolved: tests again the node
							if dbg: print 'dbg 1.1'
							if recursive:
								l, c = node.findDuplications(minbs, originalnodelabs=originalnodelabs, modifyTopology=modifyTopology)	
								ldupli += l
								if c: changedTopo = True
							break # the for child loop
					else :
						# no way to improve topology: this node is a duplication
						if not self in ldupli:
							if dbg: print 'dbg 1.2'
							ldupli.append(self)
						if recursive:
							# must keep track of the children to inspect (father-to-child dependencies may change during the inspection of the first son, label of the second son remain stable as the two son subtrees are independent)
							for childlab in self.children_labels():
								if dbg: print childlab
								l, c = genetree[childlab].findDuplications(minbs, originalnodelabs=originalnodelabs, modifyTopology=modifyTopology, verbose=verbose)	
								ldupli += l
								if c: changedTopo = True
		return (ldupli, changedTopo)
						
	def annotateDuplications(self, minbs, maxlossrate=1, reftree=None, duplabels=None, modifyTopology=True, uncertainEvents=True, verbose=False):
		"""reconciliation of gene trees by finding duplications.
		
		Finds unambiguous duplications (not likely to be transfers), i.e. those with unicopy outgroups in the observed set of species
		or those ancestral to the set of species with no pattern of ancestral differential losses between paralogs
		or ambiguous duplications, i.e. those ancestral to the set of species with a pattern of ancestral differential losses between paralogs (potentially a transfer from outside the observed set of species)
		(annotation of ambiguity to be implemented)
		
		returns the list of conflicting nodes (duplications or other events) 
		and annotates the tree with events (!!! no copy is made, GeneTree instance is edited in place).
		"""
					
		if reftree: reftree.complete_node_ids()
		self.complete_internal_labels()
		for node in self:
			node.computeUnicity()
		if verbose: print self.newick(ignoreBS=True)
		# at putative duplication nodes, solves the conflict by moving non-supported branches or stores the dupicated node in
		if duplabels:
			ldupli = [self[lab] for lab in duplabels]
			changedTopo = False
		else:		
			originalnodelabs = self.get_children_labels()
			ldupli, changedTopo = self.findDuplications(minbs, originalnodelabs=originalnodelabs, modifyTopology=modifyTopology)
			# found duplication nodes may have been subsequently deleted
			self.complete_node_ids()
			duplabels = [dup.label() for dup in ldupli]
		## annotation of the duplication event in the tree ; go from leaves to root (post-order traversal)
		ldupli = self.sort(ldupli, order=2)
		for dup in ldupli:
			# evaluation of the confilct : duplication or transfer ?
			#~ cons = dup.computeDupliConsistencyScore()
			#~ if (cons >= mincons):
			if verbose: print dup.label(),
			lossrate = dup.computeDuplicationLossRate(reftree, mode='pernode', verbose=verbose)
			if verbose: print "lossrate", lossrate
			if lossrate <= maxlossrate:
				dup.set_duplication(reftree=reftree)
			elif uncertainEvents:
				dup.set_uncertainevent(reftree=reftree)
			#~ else:
				#~ if verbose: print "\tas transfer" ,
				#~ # conflict is not consistent with a duplication and will be reconciled as a transfer event (DEPRECATED)
				#~ recchild, recspe = dup.findAdditiveTransfer(reftree)
				#~ donchild = recchild.go_brother()
				#~ donspe = donchild.listSpecies()
				#~ ## mapping of transfer receptor and donor nodes on reference tree 
				#~ # simple/reference mapping (most recent common ancestors of leaf sets) 
				#~ # and cautious mapping (takes into account species representation in gene tree as well as topological uncertainties given branch support)
				#~ refreccla, lrec = recchild.possibleReceptors(reftree, recspe, self)
				#~ refdoncla, ldon = donchild.possibleDonors(reftree, donspe, recspe, self, minbs)
				#~ # find transfer node = the father of donor and receptor branches
				#~ genetrans = self.coalesce([recchild, donchild])	
				#~ if verbose: print "@", genetrans.label(), "->", recchild.label()
				#~ genetrans.set_transfer(reclab=refreccla.label(), lrec=lrec, donlab=refdoncla.label(), ldon=ldon, childid=recchild.nodeid(), reftree=reftree)
				#~ if verbose: print genetrans.event()
		# returns the list of duplication node labels and boolean stating of topology has changed
		return (duplabels, changedTopo)
		
	def annotateDuplicationsOrAdditiveTransfers(self, minbs, mincons, reftree, duplabels=None, modifyTopology=True, verbose=None):
		"""alias for annotateDuplications(), in which nodes with unicity conflict not consistent with a duplication are reconciled by seting a transfer event at the node or within the subtree it defines."""
		return self.annotateDuplications(minbs=minbs, mincons=mincons, reftree=reftree, duplabels=duplabels, modifyTopology=modifyTopology, uncertainEvents=False, verbose=verbose)
		
	def findAdditiveTransfer(self, reftree):
		"""knowing this node shows unicity conflict not consistent with a duplication, a transfer event is infered and  mapped at the node or within the subtree it defines."""
		# determines the sources of conflict
		lconfl = self.whichMulticopySubtrees()
		print "lconfl", [n.label() for n in lconfl]
		# mapping of the duplication node in the reference species tree
		mrcadup = reftree.map_to_node(self.listSpecies())
		lspemrca = mrcadup.get_leaf_labels()				
		speratio = 1
		# choose the conflicting node that is thought to be received by transfer
		recspe = []
		for confl in lconfl:
			lspe = confl.listSpecies()
			# criterion for being received by transfer id to have a sparse set of the species under the MRCA of all species present under the duplication node
			ratio = float(len(lspe)) / float(len(lspemrca))	
			if ratio <= speratio:
				speratio = ratio
				tc = confl
				recspe = lspe
		return (tc, recspe)		
		
	def findTaxonomicPerturbation(self, reftree, taxpertcheckednodes={}, minperturbation=1, silent=True):
		"""search for nodes in self that induce an over-estimation of the common ancestor of species represented in the tree
		
		(tentative of) implementation of XD (transfer detection) algorithm from TPMS program (Bigot et al. BMC Bioinfo. 2013?).
		'taxpertcheckednodes' is a dictionary with keys the node labels and values a boolean stating the result (True = there was taxonomic perturbation)
		of nodes that were previously checkded with this function (and the cached results are not expected to have changed)	 
		"""
		ltrans = []			# output list of transfered node objects (avoid later fetching)
		#~ dlabtpbool = {}		# output dictionary recording checked nodes, same structure as 'taxpertcheckednodes'
		dleafspe = self.dictLeafLabelsToSpecies()		
		lspe = self.listSpecies()
		rt = copy.deepcopy(reftree)
		# prune species that are absent in subnode
		for leaflab in rt.get_leaf_labels():
			if leaflab not in lspe:
				rt.pop(leaflab)	
		
		subnodes = self.get_sorted_children(order=2)
		rootnodes = [self]+self.get_children()
		for subnode in subnodes:
			if not subnode in rootnodes:
				# nodes that have a grand-father to compare with
				sublab = subnode.label()
				if sublab in taxpertcheckednodes:
					if taxpertcheckednodes[sublab]:
							ltrans.append(subnode)
				else:
					taxpertcheckednodes[sublab] = False	# records node was checked, default to no transfer
					fatnode = subnode.go_father()
					grandnode = fatnode.go_father()
					#~ print "\nsubnode", subnode.label()
					#~ print "grandnode", grandnode.label()
					grandleaves = set(grandnode.get_leaf_labels())
					#~ print "grandleaves", grandleaves
					ltransleaves = []
					#~ print "ltrans [", 
					for transnode in ltrans:
						if transnode.is_child(grandnode):
							#~ print transnode.label(),
							# this perturbation is already spotted, ignore the leaves
							ltransleaves += transnode.get_leaf_labels()
							# to do as if the transfered subtree was pruned, 
							if transnode.go_father() in [fatnode, grandnode]:
								# if the subnode father or grand-father would bave been collapsed, go one father & grand-father up
								fatnode = grandnode
								grandnode = grandnode.go_father()
					if not grandnode or self.is_child(grandnode):
						# gone too high, no more grand-father node to compare with
						continue # for subnode loop
						
					#~ print "]\ngrandnode", grandnode.label()
					grandleaves = set(grandnode.get_leaf_labels()) - set(ltransleaves)
					#~ print "grandleaves", grandleaves
					grandspe = [dleafspe[leaf] for leaf in grandleaves]
					anc = rt.map_to_node(grandspe)
					refanc = reftree.map_to_node(grandspe)
					# masking the species represented under subnode 
					subleaves = grandleaves - set(subnode.get_leaf_labels())
					#~ print "subleaves", subleaves
					if not subleaves: continue
					subspe = [dleafspe[leaf] for leaf in subleaves]
					subanc = rt.map_to_node(subspe)
					subrefanc = reftree.map_to_node(subspe)
					if subanc != anc:
						# masking may reduce the ancienty of grandnode species ancestor
						pertthresh = minperturbation #+ 1*(len(subnode.listSpecies())==1) # if subnode is a leaf (or a group of duplicated leaves), increase the threshold to avoid trivial detection of the tree reduction
						d = anc.node_distance(subanc)
						if d >= pertthresh :
							if not silent: print sublab, refanc.label(), '%d->'%d, subrefanc.label()
							ltrans.append(subnode)
							taxpertcheckednodes[sublab] = True	# records node was found to be a transfer
					
		return ltrans, taxpertcheckednodes
			
	def possibleReceptors(self, reftree, recspe, inftree, **kw):
		"""mapping of Prunier transfer receptor clade on reference tree (most recent common ancestors of leaf sets)
		
		'inftree' is the gene tree where the transfer inference was made
		'recspe' is the set of species infered to be descendent of the transfer recipient
		does not refer to self, but is let as a GeneTree method for sake of coherence with possibleDonors() function
		"""
		silent = kw.get('silent', True)
		refreccla = reftree.map_to_node(recspe)
		## 'cautious' mapping of Prunier transfer receptor clade on reference tree : 
		# given the set of species represented in the transfered leaf set, 
		# the set of possible receptors in the species tree is the path between their most recent 
		# and their most ancient common ancestor not being a parent of other leaves found in gene tree			
		allspe = set(inftree.listSpecies())
		cautiousreccla = reftree.map_to_ancient_node(recspe, (allspe - set(recspe)) )
		if not silent: print "refreccla", refreccla.label(), "cautiousreccla", cautiousreccla.label()
		lrec = refreccla.path_to(cautiousreccla, returnLabels=True)
		return (refreccla, lrec)
		
	def possibleDonors(self, reftree, donspe, recspe, inftree, **kw):
		"""mapping of Prunier transfer donor clade on reference tree (most recent common ancestors of leaf sets))
		
		'inftree' is the gene tree where the transfer inference was made
		'recspe' is the set of species infered to be descendent of the transfer recipient
		'donspe' is the set of species infered to be descendent of the transfer donor
		"""
		## 'cautious' mapping of Prunier transfer donor clade on reference tree : 
		# the trivial donor in gene tree is the neighbor/brother of receiver clade 'genereccla', but low branch supports around it makes it uncertain,
		# one can define boundaries around the trivial donor in the gene tree ; 
		# those boundaries can be mapped in the species tree as the most recent common ancestor of their leaf sets
		# the set of possible donors is the explorable space between the mapped boundaries 
		# and the the most ancient common ancestor of the upper boundary not being a parent of other leaves found in gene tree
		minbs = kw.get('minbs')
		silent = kw.get('silent', True)
		#~ if not silent: print "minbs", minbs
		#~ if not silent: print "recspe", recspe
		#~ if not silent: print "donspe", donspe
		refdoncla = reftree.map_to_node(donspe)
		genedoncla = self						# 'self' node is considered as the trivial donor in gene tree
		if not silent: print "genedoncla", genedoncla.label()
		allspe = set(inftree.listSpecies())	
		if minbs:
			# explore the possible nodes in gene tree that may be donor as well as 'genedoncla' (given they are not separated by highly supported branches)
			lbound = genedoncla.find_support_boundaries(minbs)
			#~ if not silent: print "lbound", [n.label() for n in lbound]
			# map boundaries to possible nodes of the gene tree to the species tree
			lbref = []
			for bound in lbound:
				# suppose there is no conflict under boundary nodes, i.e. the subtree corresponds to one in the species tree and the mapped node is next to 'refdoncla' in species tree
				# if not, correspondance will be mapped to nodes too high (toward root) in species tree or out of donor's neighborhood
				# this is not a pain however, as the node corresponding to the upper boundary in 'lbound' will be met first when exploring species tree
				blspe = set(bound.listSpecies()) - set(recspe)
				bref = reftree.map_to_node(blspe)
				if bref:
					lbref.append(bref)
			# find the upper/deepest boudary
			upbref = lbref[0]
			for bref in lbref:
				if bref.depth() < upbref.depth():
					upbref = bref
				if bref == reftree:
					break
			# change it into its most ancient ancestor not being a parent of other leaves found in gene tree
			iup = lbref.index(upbref)
			macbref = reftree.coalesce_ancient([bref.label()], (allspe - set(bref.get_leaf_labels())) )
			lbref[iup] = macbref
			if not silent: print "lbref", [n.label() for n in lbref]
			# explore the species tree from 'refdoncla' to mapped boundaries
			ldon = refdoncla.explore_within_boundaries(lbref, returnLabels=True)
		else:
			cautiousdoncla = reftree.map_to_ancient_node(donspe, (allspe - set(donspe)) )
			if not silent: print "refdoncla", refdoncla.label(), "cautiousdoncla", cautiousdoncla.label()
			# possible donor set
			if cautiousdoncla==refdoncla:
				ldon = [refdoncla.label()]
			else:
				# explore the species tree from 'refdoncla' to 'cautiousdoncla' = upper boundary (no downer boundary but 'refdoncla' itself)
				ldon = refdoncla.explore_within_boundaries([refdoncla, cautiousdoncla], returnLabels=True)
		#~ if not silent: if not silent: print "ldon", ldon
		return (refdoncla, ldon)
				
	def addTransferEventFromSubtree(self, reftree, copiedrecnode, copieddonnode=None, copiedsubtree=None, skipRootTransfers=True, skipPreviousTransfers=True, skipUncoherentTransfers=True, **kw):
		"""map on a full gene tree (self) a new transfer event spotted at a node (copiedtransnode) included in a copy of a subtree of self (copiedsubtree)
		
		return the pruned version of transfered node (copiedtransnode), or None if no transfer was annotated.
		if no copiedsubtree is provided, will return a boolean stating if a transfer was annotated.
		"""
		silent = kw.get('silent', True)
		recchild = self.map_to_node(copiedrecnode.get_leaf_labels())
		if copieddonnode: donchild = self.map_to_node(copieddonnode.get_leaf_labels())
		else: donchild = recchild.go_brother()
		## mapping of transfer receptor and donor nodes on reference tree 
		recspe = set(recchild.listSpecies())
		donspe = set(donchild.listSpecies())
		# find transfer node = the father of donor and receptor branches
		genetrans = self.coalesce([recchild, donchild])	
		if genetrans==self:
			if skipRootTransfers:
				# cannot infer a transfer at root of gene tree ; this is rather an origination ('gain')
				if not silent: print "cannot infer a transfer at root of gene tree"
				if copiedsubtree: return None
				else: return False
			#~ else:
				#~ # create an extra node above (!!! experimental !!!)
				#~ self.go_root().create_node_above()
		children = genetrans.get_children()
		if not recchild in children:
			# transfer is higher than the direct father of recchild (because of masked transfer between recchild and donchild)
			for child in children:
				if recchild.is_child(child):
					recchild = child
					donchild = recchild.go_brother()
					recspe = recchild.represented_species()
					donspe = donchild.represented_species()
					break
		if skipUncoherentTransfers and ((donspe <= recspe) or (reftree.coalesce(donspe).is_childorself(reftree.coalesce(recspe)))):
			# donnor clade is included in receptor clade, transfer would be from a clade to its ancestor : impossible ; this is rather duplication+loss
			if not silent: print "donnor clade", donspe, "is included in receptor clade", recspe
			if copiedsubtree: return None
			else: return False
		if skipPreviousTransfers and genetrans.transfer():
			# genetrans is already a (Prunier-reliable) transfer node and changing it (its direction) may mess the reconciliation
			if not silent: print genetrans.label(), "is already a transfer"
			if copiedsubtree: return None
			else: return False
		if not silent: print "recchild", recchild.label(), "recspe", recspe,  "donchild", donchild.label(), "donspe", donspe
		# simple/reference mapping (most recent common ancestors of leaf sets) 
		# and cautious mapping (takes into account species representation in gene tree as well as topological uncertainties given branch support)
		refreccla, lrec = recchild.possibleReceptors(reftree, recspe, self, **kw)
		refdoncla, ldon = donchild.possibleDonors(reftree, donspe, recspe, self, **kw)
		if not silent: print "as transfer @", genetrans.label(), "->", recchild.label()
		#~ print "reclab", refreccla.label(), "lrec", lrec, "donlab", refdoncla.label(), "ldon", ldon, "recchild", recchild.nodeid(), '=', recchild.label()
		genetrans.set_transfer(reclab=refreccla.label(), lrec=lrec, donlab=refdoncla.label(), ldon=ldon, childid=recchild.nodeid(), reftree=reftree)
		if not silent: print genetrans.label(), genetrans.event()
		if copiedsubtree:
			# prune the transfered node from subtree to continue the reconciliation
			prunedcopiedrecnode = copiedsubtree.pop(copiedrecnode)
			return prunedcopiedrecnode
		else:
			return True
		
	def completeReconciliation(self, minbs, reftree, duplications=True, speciations=True, inferTPMStransfers=False, silent=True):
		"""take a gene tree annotated with transfers and complete the reconciliation with the inference of duplications and speciations"""	
		duplabels = []
		if duplications:
			if not silent: print "complete reconciliation at the node in unicity conflict with duplications"+"\nor transfers in presence of topological incongruence between the species and gene tree"*inferTPMStransfers
			# clean previous duplication inferences, as they may not be necessary in face of detected transfers
			self.cleanEvents(eventtype='duplication')
			prunedtree = copy.deepcopy(self)
			# forest of subtrees reconciled for transfers
			forest = prunedtree.pruneTransferedSubtrees(silent=silent) # )
			forest.append(prunedtree)
			duplabels = []
			for tst in forest:
				if not silent: print tst.listSpecies(), "present under", tst.label()
				if not silent: print tst
				# inference of duplication on pruned gene tree
				ldupli, changedTopo = tst.findDuplications(minbs=minbs, modifyTopology=False) #, verbose=(not silent))
				ldupli = tst.sort(ldupli, order=2)
				if inferTPMStransfers:
					dtaxpertcachetst = {} # global cache dict for taxonomic perturbation search for this pruned subtree
					while ldupli:
						# iterative search of transfer explaining the duplication conflict
						dtaxpertcachedupli = copy.deepcopy(dtaxpertcachetst) # cache dict for taxonomic perturbation search at this round of iterative reconciliation
						foundTransfers = False
						for dup in ldupli:
							# check at every node with duplication conflict
							if dup.label() in duplabels:
								# subtree under this node was already checked and the node was defititively considered a duplication event, ignore it
								continue # the for dup loop
							# determines the sources of multicopy conflict
							if not silent: print dup.label()
							multinodes = dup.whichMulticopySubtrees()
							if not silent: print "\tmultinodes", [n.label() for n in multinodes]
							# there may remain undetected transfers that cause the duplication conflict
							# search for sources of potential taxonomic perturbations
							ltaxpertnodes, dtaxpertcheckdup = dup.findTaxonomicPerturbation(reftree, taxpertcheckednodes=dtaxpertcachedupli, silent=silent)
							if not silent: print "\tltaxpertnodes", [n.label() for n in ltaxpertnodes]
							for tpn in ltaxpertnodes:
								for mtn in multinodes:
									if mtn==tpn or (mtn.is_child(tpn) and reftree.is_monophyletic(tpn.listSpecies())):
										# this node is taxonomically out of his place and contains a monophyletic group of taxa causing multi-copy conflict
										# consider the node as an additive transfer
										#~ foundTransfers = self.addTransferEventFromSubtree(reftree, tpn, skipPreviousTransfers=True, skipUncoherentTransfers=True, minbs=minbs)#, silent=silent)
										prunedtpn = self.addTransferEventFromSubtree(reftree, tpn, copiedsubtree=tst, skipPreviousTransfers=True, skipUncoherentTransfers=True, minbs=minbs)#, silent=silent)
										if prunedtpn: foundTransfers = True
							if foundTransfers:
								ldupli, changedTopo = tst.findDuplications(minbs=minbs, modifyTopology=False)
								ldupli = tst.sort(ldupli, order=2)
								break # the for dup loop, continue to next iteration of while dupli loop
							# else: no transfer found, this node is defititively considered a duplication event, will be ignored in next iteration rounds
							duplabels.append(dup.label())
							dtaxpertcachedupli.update(dtaxpertcheckdup)		# store the taxonomic perturbation results gathered at this node as they remain true until the subtree is changed (by pruning a transfer)
						else:
							# no transfer found, list of duplicated nodes is unchanged
							dtaxpertcachetst.update(dtaxpertcachedupli)	# store the taxonomic perturbation results as they will remain true because the subtree will not change anymore
							break # the while dupli loop
				else:
					duplabels += [dup.label() for dup in ldupli]
				if not silent: print tst.label(), "duplabels", duplabels

			#~ # recover identifiers of previously recorded nodes which type of event change from speciation to duplication events and vice-versa
			#~ tspeeventids = ()
			#~ for duplabel in duplabels:
				#~ devent = self[duplabel].getdicevent()
				#~ evtid = devent['eventid']
				#~ evtype = devent['eventtype']
				#~ if evtid and evtype=='speciation': tspeeventids += (evtid,)
			#~ # map duplications on the full gene tree and annotate them as bona fide duplications whatever the post-duplication loss rate (cf. GeneTree.computeDuplicationLossRate() )  
			#~ # (potential transfer-associated conflict was removed at previous step, remaining conflict must be duplication-associated)
			
			#~ self.annotateDuplications(minbs=minbs, maxlossrate=1, duplabels=duplabels, modifyTopology=False, uncertainEvents=False, reftree=reftree) #, verbose=(not silent))

			alldupli = [self[lab] for lab in duplabels]
			for dup in alldupli:
				dup.set_duplication(reftree=reftree, checkfortransfers=True)
				if not silent: print dup.label(), dup.event()
				
		if speciations:				
			if not silent: print "complete reconciliation with speciations"
			# and finally create speciation events where nothing else was mapped
			for node in self:
				if (not node.eventtype()):
					if node==self: node.set_gain(reftree=reftree, checkfortransfers=True)	
					else: node.set_speciation(reftree=reftree, checkfortransfers=True)	
					if not silent: print node.label(), node.event()
				elif node==self and node.eventtype()=='speciation':
					node.set_gain(reftree=reftree, checkfortransfers=True)	
				
		#~ return duplabels, tspeeventids

		
					
###################################
### XML methods

	def reconciliationEventXML(self, ind, indstart):
		xml = ""
		event = ""
		etype = self.eventtype()
		if etype:
			event += "%s<rec:%s>"%(indstart+2*ind, etype)
			if self.eventid():
				event += "%s<rec:event_id>%s</rec:event_id>"%(indstart+3*ind, str(self.eventid()))
			if etype in ['speciation', 'duplication', 'uncertain', 'gain']:
				(locspe, lnodes) = self.eventloc()
				event += "%s<rec:locationSp>%s</rec:locationSp>"%(indstart+3*ind, locspe)
				event += "%s<rec:possibleLocationsSp>%s</rec:possibleLocationsSp>"%(indstart+3*ind, ",".join(lnodes))
			elif etype=='transfer':
				(reclab, lrec, donlab, ldon, childid) = self.eventloc()
				event += "%s<rec:recipientSp>%s</rec:recipientSp>"%(indstart+3*ind, reclab)
				event += "%s<rec:originSp>%s</rec:originSp>"%(indstart+3*ind, donlab)
				event += "%s<rec:transferedChild>%d</rec:transferedChild>"%(indstart+3*ind, childid)
				event += "%s<rec:possibleRecipients>%s</rec:possibleRecipients>"%(indstart+3*ind, ",".join(lrec))
				event += "%s<rec:possibleOrigins>%s</rec:possibleOrigins>"%(indstart+3*ind, ",".join(ldon))
			event += "%s</rec:%s>"%(indstart+2*ind, etype)
			xml += "%s<rec:event>"%(indstart+ind)+event+"%s</rec:event>"%(indstart+ind)
		return xml
		
