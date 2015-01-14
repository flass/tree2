#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Enriched phylogenetic tree objects, specifically representing a species/reference trees."""

__author__ = "Florent Lassalle <florent.lassalle@ucl.ac.uk>"
__date__ = "14 January 2015"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""

import copy
import subprocess

import tree2

class ReferenceTree(tree2.AnnotatedNode):
	"""	"""	
	def __init__(self, branch_lengths=True, keep_comments=False, **kw):
		""" 
		Class of trees representing a species/reference tree on which one can annotate mapped evenents and ancestral gene contents
		
		About the 'events' attribute: it is a dictionary of counts
		of gene gain and loss at this node, duplication and transfer (reception) events and how many time the node was a donor and undergone allelic replacement.
		the first pair of digits (gain,loss) add up to the homolog count difference on between the node and its father ; 
		the second pair give the nature of gene gain(s): by duplication and/or transfer
		the third are complementary annotation
		!!! a transfer received : gain += 1, receptor +=1 ; replacement means simultaneous gain by transfer and loss (see DolloParsFromGeneTree() method)
		so the observed patterns when receiving a transfer (as yielded by newick(comment='events') are:
		transfer brings a new gene -> (1/0-0/1-0/0); transfered gene replaces an previous one -> (1/1-0/1-0/1)
		"""		
		super(ReferenceTree, self).__init__(branch_lengths=branch_lengths, keep_comments=keep_comments, **kw)	
		self.__presence = False		# character (homologous gene family) presence=True/absence=False, or homolog presence count (int), or factor state (str, several alternative states are given as str concatenated by '|')
		self.__donor = None		# if gene was transfered at this Node, indicates the donor Node.
		self.__events = dict(gain=0, loss=0, duplication=0, receptor=0, donor=0, replacement=0)
		self.__misc = {}	# dictionary for miscelaneous info, notably topology of replicon bearing the genes
		
	def newnode(self, branch_lengths=True, keep_comments=False, **kw):
		"""class-specific instance generator for construction of trees of ReferenceTree nodes"""
		return ReferenceTree(branch_lengths=branch_lengths, keep_comments=keep_comments, **kw)
		
#### Reference tree mapping methods
#### to map absence/presence of characters/genes (to make phylogenetic profiles) and for evolutionary events (gene gain/loss, HGT with Receptor/Donor).
	
	### Annotation methods	
	def presenceAtNode(self, state=True):
		"""annotates the Node for the presence state of a character"""
		self.__presence = state
		
	def addCountAtNode(self, count=1):
		"""incremant the annotation at the Node for the count of homologs/characters at this node"""
		self.__presence += count
	
	def presenceInTree(self, lspe, state=True):
		"""annotates the tree below the Node for the presence/absence state of a character given a list of labels of the nodes bearing the character."""
		for spe in lspe:
			if self[spe]:
				self[spe].presenceAtNode(state)
	
	def addCountInTree(self, lspe):
		"""annotates the tree below the Node for the count of homologs/characters given a list of labels of the nodes bearing the character or a dictionary of counts at these nodes."""
		if isinstance(lspe, dict):
			for spe in lspe:
				if self[spe]:
					if isinstance(lspe[spe], list):
						self[spe].addCountAtNode(count=len(lspe[spe]))
					else:
						self[spe].addCountAtNode(count=int(lspe[spe]))
		else:
			for spe in lspe:
				if self[spe]:
					self[spe].addCountAtNode(count=1)
				
	def resetPresence(self, state=False):
		"""annotates homogeneously all the nodes below the Node for the presence state of a character"""
		for child in self.get_all_children():
			child.presenceAtNode(state)

	def transferedFrom(self, donor):
		"""annotates the Node for a transfer event ; donor can be a Node object or a string matching a node label in the Node"""
		if isinstance(donor, Node):
			self.__donor = donor
		elif self.go_root()[str(donor)]:
			self.__donor = self.go_root()[str(donor)]
			
	def addEvent(self, events, count=1):
		"""adds one event to the event count of the node"""
		if isinstance(events, str):
			self.__events[events] += count
		elif isinstance(events, list):
			for event in events:
				self.__events[event] += count
				
	def cleanEvents(self, eventtype=None, reccursive=False):
		if eventtype:
			self.__events[eventtype] = 0
		else:
			for et in self.__events:
				self.__events[et] = 0
		if reccursive:
			for c in self.get_children():
				c.cleanEvents(eventtype=eventtype, reccursive=reccursive)
				
	def misc(self):
		return self.__misc
				
	def __iadd__(self, reftree):
		"""recusively adds event and homolog counts from different instances of the SAME template reference tree"""
		# increment event dict counts
		for eventtype in self.__events:
			self.__events[eventtype] += reftree[self.label()].getEvents()[eventtype]
		# increment presence counts
		self.__presence += reftree[self.label()].homolog_count()
		# increment misc dict counts
		refd = reftree[self.label()].misc()
		selfd = self.__misc
		if selfd or refd:
			for key in refd:
				selfd[key] = selfd.get(key, 0) + refd.get(key, 0)
		# recursion
		for child in self.get_children():
			child += reftree
		return self
				
	def __isub__(self, reftree):
		"""recusively substracts event and homolog counts from different instances of the SAME template reference tree"""
		for eventtype in self.__events:
			self.__events[eventtype] -= reftree[self.label()].getEvents()[eventtype]
		self.__presence -= reftree[self.label()].homolog_count()
		for child in self.get_children():
			child -= reftree
		return self
		
	def __add__(self, reftree):
		sumtree = copy.deepcopy(self)
		sumtree += reftree
		return sumtree
		
	def __sub__(self, reftree):
		diftree = copy.deepcopy(self)
		diftree -= reftree
		return diftree
		
	def diff(self, reftree):
		"""returns a tree with the difference in homolog counts of (self - reftree), detailed as events.
		
		Summary difference of counts is given in the presence attribute (as a relative integer),
		detailed difference of counts is given in the events attribute as gains and losses (always positive integers).
		This allows to sum such trees without loosing the detail of differences.
		!!! differs from stating "self - reftree" (using __sub__()), 
		because the events returned by diff reflect the differences in homolog counts, not the diferences in event counts.
		"""
		
		dif = copy.deepcopy(self)
		dif -= reftree
		for node in self:
			lab = node.label()
			dif[lab].cleanEvents()
			d = node.homolog_count() - reftree[lab].homolog_count()
			if d > 0: dif[lab].addEvent('gain', d)
			if d < 0: dif[lab].addEvent('loss', -d)
		return dif		
				
	def getEvents(self):
		return self.__events
		
	def whichAncEvent(self, eventtype, returnLabels=False, asSet=False):
		l = []
		for anc in self:
			nevt = anc.countEvents(eventtype) 
			if nevt > 0:
				if asSet:
					if returnLabels: l.append(anc.label())
					else: l.append(anc)
				else:
					if returnLabels: l += [anc.label()]*nevt
					else: l += [anc]*nevt
		return l
		
	def countEvents(self, eventtype=None):
		if eventtype: return self.__events[eventtype]
		else: return sum(self.__events.values())
		
	def sumEvents(self, eventtype=None):
		nevt = self.countEvents(eventtype=eventtype)
		for child in self.get_children():
			nevt += child.sumEvents(eventtype=eventtype)
		return nevt
		
	def state(self):
		return self.__presence
		
	def factorStateToDict(self, cleanState=True, alternativeAsFraction=True, recursive=True):
		"""method to transform (dirty) storage of states as strings in __presence attribute to keys in __misc dictionary attribut"""
		states = self.state().split('|')
		if alternativeAsFraction: f = 1.0/len(states)
		else: f = 1
		d = self.misc()
		for state in states: d[state] = f
		if cleanState: self.presenceAtNode(False)
		if recursive:
			for child in self.get_children():
				child.factorStateToDict(alternativeAsFraction=alternativeAsFraction, recursive=recursive)
		
	def children_states(self, excludedNodes=[]):
		states = []
		for child in self.get_children():
			if not (child in excludedNodes):
				states.append(child.state())
		return states

	def homolog_count(self):
		return int(self.__presence)
	
	def is_present(self):
		return bool(self.__presence)
	presence = is_present
	
	def sumHomologCounts(self):
		sumh = self.homolog_count()
		for child in self.get_children():
			sumh += child.sumHomologCounts()
		return sumh
		
	def checkGainLossPresenceCoherence(self, fatcount=0, noError=False, silent=True):
		count = fatcount + self.__events['gain'] - self.__events['loss']
		if count != self.__presence:
			e = "sum of father count + gain - loss != homolog count at %s:\n%d + %d - %d != %d"%(self.label(), fatcount, self.__events['gain'], self.__events['loss'], self.__presence)
			if noError:
				print e
				return False
			else:
				raise ValueError, e
		else:
			check = True
			for child in self.get_children():
				check *= child.checkGainLossPresenceCoherence(fatcount=count, noError=noError)
			if not silent: print "coherent scenario,", ("YES" if check else "NO")
			return bool(check)
		
	def is_present_in_ancestor(self, depth=1):
		f = self
		for n in range(depth):
			if f.__father:
				f = f.go_father()
			else:
				raise IndexError, 'Specified ancestor does not exist'
		return f.is_present()
				
	def donor(self):
		"""if a transfer event reception is recorded at the Node, returns the donor, None otherwise."""
		return self.__donor
				
	### Pylogenetic profiling methods
	def is_unicopy_in_clade(self):
		"""Returns True if all leaves of the Node bear the gene in zero or one copy"""
		leaves = self.get_leaves()
		for leaf in leaves:
			if leaf.homolog_count() > 1:
				return False
		return True
		
	def is_ubiquitous_in_clade(self):
		"""Returns True if all leaves of the Node bear the gene"""
		allp = True
		leaves = self.get_leaves()
		for leaf in leaves:
			allp *= leaf.is_present()
		return allp
		
	def homogeneous_state_in_clade(self, excludedLeaves=[], counts=False):
		"""Tests homogeneity of character states among leaves of the Node, or among a reduced set ; if homogeneous, returns the common state ['True'/'False'], if not, returns 'None'."""
		leaves = set(self.get_leaves()) - set(excludedLeaves)
		if not leaves:
			return None # if excluded leaf set contains the set of leaves under the Node
		states = []
		for leaf in leaves:
			if not counts: states.append(leaf.is_present())
			else: states.append(leaf.homolog_count())
		for state in states:
			if state != states[0]:
				return None
		else:
			return states[0]
		
	def clade_specific_contrast(self, contrast='full', force=False, counts=False, returnforeback=False, excludedLeaves=[]):
		"""Tests existence of contrast between homogeneous character states inside and outside the Node; 
		
		outside is defined either as a node including the Node in the tree (contrast='nodelabel'), 
		as the rest of the full tree (contrast='full'),
		or as the rest of the full tree that have the gene present (contrast='allpres') or absent (contrast='allabs'),
		or as a subtree defined by the depth (in node distance) from the Node (contrast=int());
		if there is a contrast, returns the  common state inside the Node, if not, returns 'None'.
		"""
		if contrast.isdigit():
			depth = int(contrast)
			oldest = self
			while depth and oldest!=self.go_root():
				oldest = oldest.go_father()
				depth -= 1
		else :
			if contrast=='full':
				oldest = self.go_root()
			elif contrast in ['allpres', 'allabs']:
				if contrast=='allpres': lspe = self.presents(returnLabels=True, leavesOnly=True)
				elif contrast=='allabs': lspe = self.absents(returnLabels=True, leavesOnly=True)
				oldest = self.map_to_node(lspe)
				if not oldest:
					if not force: raise ValueError, "wrong contrast definition, '%s' criterion yield species set %s with no matching ancestor"%(contrast, lspe)
					else: return None
			else:
				oldest = self.go_root()[contrast]
				if not oldest: raise ValueError, "wrong contrast definition, %s does not match a node label in tree"%contrast
			if not self.is_child(oldest):
				if not force: raise ValueError, "wrong contrast definition, %s is not a ancestor node of %s"%(contrast, self.label())
				else: return None
				
		foreground = self.homogeneous_state_in_clade(counts=counts)
		background = oldest.homogeneous_state_in_clade(excludedLeaves=self.get_leaves()+excludedLeaves, counts=counts)
		if foreground==background or foreground==None or background==None:
			return None
		else:
			if returnforeback: return (foreground, background)
			else: return foreground
	
	def getPhyloProfile(self, returnLabels=False, leavesOnly=True, counts=True):
		"""returns a dictionary of presence/absence states at leaves of the Node. by default, keys are Node objects ; if returnLabels=True, they are the node labels (strings)."""
		d = {}
		if leavesOnly: nodes = self.get_leaves()
		else: nodes = self.get_all_children()
		for node in nodes:
			if not counts: state = node.state()
			else: state = node.homolog_count()
			if not returnLabels: d[node] = state
			else: d[node.label()] = state
		return d	
		
	def presents(self, returnLabels=False, leavesOnly=True, counts=False):
		lpres = []
		for node in self:
			if leavesOnly and not node.is_leaf():
				continue
			if returnLabels: n = node.label()
			else: n = node
			if counts: lpres += [n]*node.homolog_count()
			else: lpres += [n]*node.presence()
		return lpres
		
	def absents(self, returnLabels=False, leavesOnly=True):
		labs = []
		for node in self:
			if leavesOnly and not node.is_leaf():
				continue
			if returnLabels: n = node.label()
			else: n = node
			if not node.is_present(): labs.append(n)
		return labs
		
	def writePhyloProfileMatrix(self, fmatout, fam="somefamily", reflabs=None, leavesOnly=True, header=False, counts=True):
		"""write the phylogenetic profile of the family into a line of a table file readable by COUNT (Csuros and Miklos, 2008)"""
		if not reflabs:
			if leavesOnly: refnodelabs = self.get_leaf_labels()
			else: refnodelabs = self.get_children_labels()
		else:
			refnodelabs = reflabs
		if header:
			fmatout.write('\t'.join(['family']+refnodelabs)+'\n')
		profile = self.getPhyloProfile(returnLabels=True, leavesOnly=leavesOnly, counts=counts)
		s = fam
		for lab in refnodelabs:
			s +=  '\t'+str(profile[lab])
		s += '\n'
		fmatout.write(s)
		return profile
		
	def writePhyloProfileLine(self, nfmatout, fam="somefamily", reflabs=None, leavesOnly=True, header=False):
		"""wrapper for writePhyloProfileMatrix() for one-line matrix output"""
		fmatout = open(nfmatout, 'w')
		profile = self.writePhyloProfileMatrix(fmatout, fam=fam, reflabs=reflabs, leavesOnly=leavesOnly, header=True)
		fmatout.close()
		return profile
	
	def getTransferEvents(self, returnLabels=False, excludeSelf=False):
		"""scans the whole tree below the Node and returns a dictionary of transfers d[receptor]=donor. by default, keys and values are Node objects ; if returnLabels=True, they are the node labels (strings)."""
		d = {}
		children = self.get_all_children()
		for child in children:
			if child==self and excludeSelf:
				continue
			if child.donor():
				if not returnLabels: d[child] = child.donor()
				else: d[child.label()] = child.donor().label()
		return d
	
	### Ancestral state reconstruction methods:
	
	def FitchPars(self, excludedNodes=[], silent=True, refine=False):
		"""infers ancestral (factor) states by Fitch parsimony (recursive) given phylogenetic profile and annotates the tree with inferred states"""
		# recursive leaf-to-root inference of states
		for child in self.get_children():
			if not silent: print child.label(),
			if not child.is_leaf(): # and not child in excludedNodes:
				child.FitchPars(silent=silent, excludedNodes=excludedNodes, refine=refine)
		if not self in excludedNodes:
			childrenstates = self.children_states(excludedNodes=excludedNodes)
			if not silent: print self.label(), "\tchildrenstates", childrenstates,
			childoptstates = [cs.split('|') for cs in childrenstates]
			if not silent: print "\tchildoptstates", childoptstates,
			optstates = set(childoptstates[0])
			for cos in childoptstates[1:]:
				optstates &= set(cos)
			if not optstates:
				optstates = set(childoptstates[0])
				for cos in childoptstates[1:]:
					optstates |= set(cos)
			if not silent: print ", optstates", optstates
			self.presenceAtNode(tree2.setToStr(optstates))
		if refine and self==self.go_root():
			self.refineFitchPars(excludedNodes=excludedNodes, silent=silent)
			
	def refineFitchPars(self, excludedNodes=[], silent=True):
		"""refine inferred states where several are possible 
		by taking into account the parental and children's state: 
		  when the children have each one a distinct state set, one being the same than the parental state set,
		  refine the node's state to this child-parent common state 
		  (implies one single state change in the other child rather than two changes).
		also considers unknown state as non determining the ancestral states.
		"""
		# recursive leaf-to-root refinement of inference of states
		for child in self.get_children():
			# if (not child.is_leaf()):
			child.refineFitchPars(silent=silent, excludedNodes=excludedNodes)
		states = self.state().split('|')
		if (not self.is_root()) and (not self in excludedNodes):
			son = self
			father = self.go_father()
			line = [self]
			if states==['?']:
				# unknown state, try to know it from ancestors
				# !!! change the states of leaves when "?"
				if not silent: print self.label(), "states", "|".join(states)
				while father and (father not in excludedNodes):
					fatherstates = set(father.state().split('|')) - set(['?'])
					if fatherstates:
						fstates = tree2.setToStr(fatherstates)
						if not silent: print " ", father.label(), "fatherstates", fstates
						for node in line:
							# restrict possible states of the node and intermediary nodes to the informative ancestor
							node.presenceAtNode(fstates)
							if not silent: print '  ->', node.label(), "states", node.state()
						break
					else:
						line.append(father)	# if state of initial node is changed to known state regarding an undirect ancestor, implies that intermediary nodes' states are to be changed too.
						son = father
						father = son.go_father()
			if len(states)>1:
				if not silent: print self.label(), "states", "|".join(states)
				if '?' in states:
					# remove undetermined states if an alternative exists
					states.remove('?')
					self.presenceAtNode("|".join(states))
					if not silent: print '  ->', self.label(), "states", self.state()
				if len(states)>1:
					# if several states are still possible, check ancestors
					sstates = set(states)
					while father and (father not in excludedNodes):
						fatherstates = set(father.state().split('|')) - set(['?'])
						sonstates = set(son.state().split('|')) - set(['?'])
						if not silent: print ' ', son.label(), 'sonstates', son.state(), father.label(), 'fatherstates', father.state()
						fatsonstates = fatherstates & sonstates
						fatselfstates = fatherstates & sstates
						if (fatsonstates and (fatsonstates < sonstates)) and (fatselfstates and (fatselfstates < sstates)):
							# in this ancestor, the state subset inferred from the son is restrictedat least as restrictive as one of the self node's children
							fstates = tree2.setToStr(fatselfstates)
							for node in line:
								# restrict possible states of the node and intermediary nodes to the informative ancestor
								node.presenceAtNode(fstates)
								if not silent: print '  ->', node.label(), "states", node.state()
							if len(fatselfstates)==1:
								# reduced the state set to the minimum, stop here
								break # the while father loop
						elif not fatselfstates:
							# ancestor states are disjoint from the self node states, stop here
							break # the while father loop
						# should look one node above to get possible restricted choice of states in higher ancestors
						line.append(father)	# if state of initial node is restricted regarding an undirect ancestor, implies that intermediary nodes' states are restricted too.
						son = father
						father = son.go_father()							
		
	## to do : treat multiple occurence: unique state is coded A|B|C, multiple state would be A|B|C&D
	#~ def FitchPars2(self, excludedNodes=[], silent=True, refine=False):
		#~ """infers ancestral (factor) states by Fitch parsimony (recursive) given phylogenetic profile and annotates the tree with inferred states"""
		#~ # recursive leaf-to-root inference of states
		#~ for child in self.get_children():
			#~ if not silent: print child.label(),
			#~ if not child.is_leaf(): # and not child in excludedNodes:
				#~ child.FitchPars(silent=silent, excludedNodes=excludedNodes, refine=refine)
		#~ if not self in excludedNodes:
			#~ childrenstates = self.children_states(excludedNodes=excludedNodes)
			#~ if not silent: print self.label(), "\tchildrenstates", childrenstates,
			#~ childoptstates = [cs.split('|') for cs in childrenstates]
			#~ if not silent: print "\tchildoptstates", childoptstates,
			#~ optstates = set(childoptstates[0])
			#~ for cos in childoptstates[1:]:
				#~ optstates &= set(cos)
			#~ if not optstates:
				#~ optstates = set(childoptstates[0])
				#~ for cos in childoptstates[1:]:
					#~ optstates |= set(cos)
			#~ if not silent: print ", optstates", optstates
			#~ self.presenceAtNode(tree2.setToStr(optstates))
		#~ if refine and self==self.go_root():
			#~ self.refineFitchPars(excludedNodes=excludedNodes, silent=silent)
			#~ 
	#~ def refineFitchPars2(self, excludedNodes=[], silent=True):
		#~ """refine inferred states where several are possible 
		#~ by taking into account the parental and children's state: 
		  #~ when the children have each one a distinct state set, one being the same than the parental state set,
		  #~ refine the node's state to this child-parent common state 
		  #~ (implies one single state change in the other child rather than two changes).
		#~ also considers unknown state as non determining the ancestral states.
		#~ """
		#~ # recursive leaf-to-root refinement of inference of states
		#~ for child in self.get_children():
			# if (not child.is_leaf()):
			#~ child.refineFitchPars(silent=silent, excludedNodes=excludedNodes)
		#~ states = self.state().split('|')
		#~ if (not self.is_root()) and (not self in excludedNodes):
			#~ son = self
			#~ father = self.go_father()
			#~ line = [self]
			#~ if states==['?']:
				#~ # unknown state, try to know it from ancestors
				#~ # !!! change the states of leaves when "?"
				#~ if not silent: print self.label(), "states", "|".join(states)
				#~ while father and (father not in excludedNodes):
					#~ fatherstates = set(father.state().split('|')) - set(['?'])
					#~ if fatherstates:
						#~ fstates = tree2.setToStr(fatherstates)
						#~ if not silent: print " ", father.label(), "fatherstates", fstates
						#~ for node in line:
							#~ # restrict possible states of the node and intermediary nodes to the informative ancestor
							#~ node.presenceAtNode(fstates)
							#~ if not silent: print '  ->', node.label(), "states", node.state()
						#~ break
					#~ else:
						#~ line.append(father)	# if state of initial node is changed to known state regarding an undirect ancestor, implies that intermediary nodes' states are to be changed too.
						#~ son = father
						#~ father = son.go_father()
			#~ if len(states)>1:
				#~ if not silent: print self.label(), "states", "|".join(states)
				#~ if '?' in states:
					#~ # remove undetermined states if an alternative exists
					#~ states.remove('?')
					#~ self.presenceAtNode("|".join(states))
					#~ if not silent: print '  ->', self.label(), "states", self.state()
				#~ if len(states)>1:
					#~ # if several states are still possible, check ancestors
					#~ sstates = set(states)
					#~ while father and (father not in excludedNodes):
						#~ fatherstates = set(father.state().split('|')) - set(['?'])
						#~ sonstates = set(son.state().split('|')) - set(['?'])
						#~ if not silent: print ' ', son.label(), 'sonstates', son.state(), father.label(), 'fatherstates', father.state()
						#~ fatsonstates = fatherstates & sonstates
						#~ fatselfstates = fatherstates & sstates
						#~ if (fatsonstates and (fatsonstates < sonstates)) and (fatselfstates and (fatselfstates < sstates)):
							#~ # in this ancestor, the state subset inferred from the son is restrictedat least as restrictive as one of the self node's children
							#~ fstates = tree2.setToStr(fatselfstates)
							#~ for node in line:
								#~ # restrict possible states of the node and intermediary nodes to the informative ancestor
								#~ node.presenceAtNode(fstates)
								#~ if not silent: print '  ->', node.label(), "states", node.state()
							#~ if len(fatselfstates)==1:
								#~ # reduced the state set to the minimum, stop here
								#~ break # the while father loop
						#~ elif not fatselfstates:
							#~ # ancestor states are disjoint from the self node states, stop here
							#~ break # the while father loop
						#~ # should look one node above to get possible restricted choice of states in higher ancestors
						#~ line.append(father)	# if state of initial node is restricted regarding an undirect ancestor, implies that intermediary nodes' states are restricted too.
						#~ son = father
						#~ father = son.go_father()							
		
	def DolloPars(self, excludedNodes=[], inventorlab=None, silent=True):
		"""infers (dircrete numerical) ancestral states by Dollo parsimony given phylogenetic profile and annotates the tree with inferred states"""
#		for node in self.get_all_children():
#			node.presenceAtNode(False)
		excludedLeaves = []
		for node in excludedNodes:
			excludedLeaves += node.get_leaves()
		profile = self.getPhyloProfile(returnLabels=True)
		# restrict to calde where homolog is present
		#~ presents = []
		dpresence = {}
		for leaf in profile:
			#~ if (profile[leaf]==True) and (self[leaf] not in excludedLeaves):
			if (profile[leaf]) and (self[leaf] not in excludedLeaves):
				#~ presents.append(leaf)
				dpresence[leaf] = profile[leaf]
		if inventorlab:
			inventor = self[inventorlab]
		else:
			#~ inventor = self.map_to_node(presents)
			inventor = self.coalesce(dpresence.keys())
			inventorlab = inventor.label()
		#~ inventor.presenceAtNode(True)
		if not inventorlab in dpresence:
			# marks presence at the inventor if not at leaves (because it would be already counted in phylogenetic profile)
			inventorcount = max([dpresence[leaf] for leaf in dpresence])
			inventor.addCountAtNode(count=inventorcount)
		else:
			inventorcount = dpresence[inventorlab]
		if inventor.go_father():
			if inventor.go_father().is_present():
				# precise it is a replacement
				inventor.addEvent('replacement', inventorcount)
		inventor.addEvent('gain', inventorcount) # it is a gain anyway
		# descent from inventor ancestor to leaves
		children = inventor.get_sorted_children(order=0)[1:] # excludes inventor ancestor 
		for child in children:
			if child in excludedNodes: continue	
			fathercount = child.go_father().homolog_count()
			childcount = max([leaf.homolog_count() for leaf in child.get_leaves()])
			child.addEvent('loss', (fathercount - childcount))
			child.presenceAtNode(childcount)
		check = self.checkGainLossPresenceCoherence(noError=True)
		if not check: raise ValueError, "\n%s\n%s"%(self.newick(comment='presence'), self.newick(comment='events'))
				
	def DolloParsWithTransfers(self):
		"""infers ancestral states by custom Dollo parsimony given phylogenetic profile and transfer history and annotates the tree with inferred states
		
		treats independently each tranfered subtree for parsimony inference of presence.
		"""
		dtrans = self.getTransferEvents(excludeSelf=True)
		receptors = dtrans.keys()
		#~ receptors.sort(key=lambda x: x.depth()) 
		self.sort(receptors, order=0)
		self.DolloPars(excludedNodes=receptors)
		for receptor in receptors:
			receptor.DolloParsWithTransfers() 
			receptor.presenceAtNode(True)
			receptor.addEvent('receptor')
			dtrans[receptor].addEvent('donor')
			#~ receptor.transferedFrom(dtrans[receptor])
	
	def readCountAsymmetricWagner(self, nfoutput, **kw):
		silent = kw.get('silent', True)
		#~ silent = False
		fam = kw.get('fam', "somefamily")
		"""read output file from the program AsymmetricWagner from Count (Csuros & Miklos, 2009)"""
		foutput = open(nfoutput, 'r')
		# read presence information ('# FAMILY' lines)
		line = foutput.readline()
		while not line.startswith('# FAMILY'):
			line = foutput.readline()
		header = line.strip('# \n').split('\t')[2:-4]
		line = foutput.readline()
		nfam = 1
		while line:
			if line.startswith('# FAMILY\t%s'%fam):
				break
			else:
				nfam += 1
				line = foutput.readline()
		else:
			raise IndexError, "family %s not found in Count AsymmetricWagner ouput file %s"%(fam, nfoutput)	
		profile = line.strip('# \n').split('\t')[2:-4]
		phyloprofile = dict(zip(header, profile))
		self.resetPresence(0)
		self.addCountInTree(phyloprofile)
		if self.is_present():
			# presence at the root is not associated to an explicit gain event by Count.AsymmetricWagner, must add this implicit gain.
			nhom = self.homolog_count()
			self.addEvent('gain', nhom)
			if not silent: print nhom, 'gain @', self.label(), '(presence at root)'
		if nfam==1:
			# read event information ('# CHANGE' lines)
			dnodeevents = {}
			while not line.startswith('# CHANGE'):
				line = foutput.readline()
			for line in foutput:
				lsp = line.rstrip('\n').split('\t')
				nodelab = lsp[1]
				if nodelab != 'total':
					gene_gain, gene_duplication, gene_loss = [int(lsp[i]) for i in (2,4,6)]
					if gene_gain: 
						self[nodelab].addEvent('gain', gene_gain)
						if not silent: print gene_gain, 'gain @', nodelab
					if gene_duplication:
						self[nodelab].addEvent('duplication', gene_duplication)
						if not silent: print gene_duplication, 'duplication @', nodelab
					if gene_loss:
						self[nodelab].addEvent('loss', gene_loss)
						if not silent: print gene_loss, 'loss @', nodelab
					if (gene_gain + gene_duplication + gene_loss) > 0:
						dnodeevents[nodelab] = self[nodelab].getEvents()
		else:
			print "file %s contains info about several families, cannot extract node-spcific event counts"%nfoutput
			dnodeevents = None
		foutput.close()
		self.checkGainLossPresenceCoherence()
		return dnodeevents
		
	def	readCountPosteriors(self, nfoutput, **kw):
		silent = kw.get('silent', True)
		fam = kw.get('fam', "somefamily")
		"""read output file from the program Posteriors from Count (Csuros & Miklos, 2009)"""
		foutput = open(nfoutput, 'r')
		line = foutput.readline()
		while not line.startswith('Family'):
			line = foutput.readline()
		header = line.strip('\n').split('\t')[1:]
		print line
		while not line.startswith(fam):
			line = foutput.readline()
		print line
		results = [float(p) for p in line.strip('\n').split('\t')[1:]]
		features = ['1', 'm', 'gain', 'loss', 'expansion', 'reduction']
		dresults = dict(zip(features, [{}]*len(features)))
		for i, col in enumerate(header):
			if not ':' in col:
				# family category assignment
				continue
			nodelab, feat = col.split(':')
			dresults[feat][nodelab] = results[i]
		# parses the probability of presence of at least one copy at the node
		dnodeevents = {}
		self.resetPresence(0)
		for feat in features:
			dfeat = dresults[feat]
			for nodelab in dfeat:
				p = dfeat[nodelab]
				if not silent: print nodelab, feat, p
				if round(p): # switch p ~ 1 or 0
					# proba close to 1 (>=0.5 <=> cannot be higher in any other alternative)
					if feat=='1':
						self[nodelab].addCountAtNode(p)
					elif feat=='m':
						# NB1: if p('m')~1, necessarily p('1')=1, but Count will display p('1')=0, so an extra count is added to correct
						# NB2: copy count thus cannot exceed 2 (approximation multiple~2 is OK if one does not work on (largely) multicopy family, 
						#      like subfamilies of orthologs and rare xenologs)
						self[nodelab].addCountAtNode(p+1)
					else:
						if feat in ['gain', 'expansion']:
							self[nodelab].addEvent('gain', p)
						elif feat in ['loss', 'reduction']:
							self[nodelab].addEvent('loss', p)
						if not silent: print p, feat, '@', nodelab
						dnodeevents[nodelab] = self[nodelab].getEvents()
						
		self.checkGainLossPresenceCoherence()
		return dnodeevents
	
	def fixed_scenario_str(self):
		reftree = self.go_root()
		return "fixed scenario:\ngains %s\nlosses %s"%(str(reftree.whichAncEvent('gain', returnLabels=True)), str(reftree.whichAncEvent('loss', returnLabels=True))) 

	def relocate_gain_above(self, gainanc, silent=True):
		"""relocate gain to a more ancient node"""
		if not silent: print "relocate gain above\nfrom Count gain ancestor", gainanc.label(), "to referenced event ancestor", self.label()
		gainanc.addEvent('gain', -1)
		f = gainanc.go_father()
		n = gainanc
		while n!=self:
			f.addCountAtNode(1)
			b = n.go_brother()
			b.addEvent('loss', 1)
			n = f
			f = n.go_father()
		self.addEvent('gain', 1)
		if not silent: print self.fixed_scenario_str()
		return self.go_root().checkGainLossPresenceCoherence(silent=silent)
		
	def relocate_gain_below(self, gainanc, silent=True):
		"""relocate gain to a more recent node"""
		if not silent: print "relocate gain below\nfrom Count gain ancestor", gainanc.label(), "to referenced event ancestor", self.label()
		gainanc.addEvent('gain', -1)
		f = gainanc
		childancs = f.get_children()
		# clean above self
		while self.is_child(f):
			# remove presence on self lineage above self
			f.addCountAtNode(-1)
			for c in childancs:
				if self.is_child(c) or self==c:
					f = c
				else:
					# remove losses out of self lineage above self
					c.addEvent('loss', -1)
			childancs = f.get_children()
		self.addEvent('gain', 1)
		if not silent: print self.fixed_scenario_str()
		return self.go_root().checkGainLossPresenceCoherence(silent=silent)
		
	def correct_genetree_event(self, gainanc, eventnode, ost, modifiedgt, silent=True):
		lspe = ost.listSpecies()
		mrca = self.map_to_node(lspe)
		#~ print lspe, mrca.label()
		if gainanc != mrca:
			# may have to consider transfers in ost to correctly compute the 
			raise IndexError, "uncoherent gain mapped at %s while ancestral location of event is %s"%(gainanc.label(), str(eventnode.event()))
		# gain is indeed found at the MRCA of species represented in orthologous subtree
		# location/receptor mapping was wrong in eventnode, change it in the gene tree
		if not modifiedgt: modifiedgt = copy.deepcopy(eventnode.go_root())	# copy of the full gene tree to be returned
		modeventnode = modifiedgt[eventnode.label()]
		if eventnode.transfer():
			transcoord = modeventnode.eventloc()
			modeventnode.set_transfer(reclab=gainanc.label(), donlab=transcoord[2], ldon=transcoord[3], childid=transcoord[4], reftree=self)
		else:
			modeventnode.set_speciation(locspe=gainanc.label(), reftree=self)
		if not silent: print "corrected gene tree event at",  modeventnode.label(), ":",  modeventnode.event
		
		return modifiedgt
			
	def checkCountGainCoordinates(self, gainanc, reftree, eventnode, ost, modifiedgt, silent=True):
		"""check the coherence of a gain inferred by Count programs and the reconciliation of the gene tree, and fix them if necessary"""
		if gainanc.is_child(self):
			if eventnode.duplication():
				# each subtree under a duplication has equal or fewer represented species that the whole duplicated tree
				# so the MRCA of species where he gain was mapped may be lower for each subtree than fot the duplication ancestor.
				# relocate gain at location of the duplication event
				self.relocate_gain_above(gainanc, silent=silent)
			else:
				modifiedgt = reftree.correct_genetree_event(gainanc, eventnode, ost, modifiedgt, silent=silent)
		elif self.is_child(gainanc):
			# Count locate the gain too high in the tree
			if gainanc==reftree:
				# gain at the root is not counted by Count, but scenario is not the most parsimonious if one counts the root gain
				# relocate gain at location of the node event
				self.relocate_gain_below(gainanc, silent=silent)
		else:
			modifiedgt = reftree.correct_genetree_event(gainanc, eventnode, ost, modifiedgt, silent=silent)
		return modifiedgt
		
	def checkCountInference(self, eventnode, ost, modifiedgt, **kw):
		"""check the coherence of the duplication/gain/loss scenario from Count programs and the reconciliation of the gene tree, and fix them if necessary"""
			
		def checkRefEvent(self, lgainanc, silent=True):
			refgain = 0
			for gainanc in lgainanc:
				if self==gainanc:
					if not silent: print "1 gain @", self.label(), ", justified by the reference event"
					lgainanc.remove(gainanc)
					refgain = 1
			return lgainanc, refgain
			
		def checkExistingTransfers(reftree, lgainanc, forest, silent=True):
			"""filter gains corresponding to already annotated transfers"""
			for tst in forest:
				if tst:
					mrca = reftree.map_to_node(tst.listSpecies())
					if mrca in lgainanc:
						lgainanc.remove(mrca)
						if not silent: print "1 gain @",  mrca.label(), ", justified by a transfer ->", tst.label()
			return lgainanc
		
		# set avriables
		silent = kw.get('silent', True)
		if not silent:
			print "ref event", eventnode.event(), "at", eventnode.label()
			print "self (reference event receptor ancestor)", self.label()
		reftree = self.go_root()
		genetree = eventnode.go_root()
		
		# get rid of duplications
		ldupanc = reftree.whichAncEvent('duplication', asSet=True)
		if ldupanc:
				print "infered duplication in ancestors %s (in gene subtree under %s) : not in model assumption (only unicopy subtrees), will consider it as a transfer"%(str([dupanc.label() for dupanc in ldupanc]), eventnode.label())
				for dupanc in ldupanc:
					dupanc.cleanEvents('duplication')
		lgainanc = reftree.whichAncEvent('gain')
		lgainanc, refgain = checkRefEvent(self, lgainanc, silent=silent)
		if len(lgainanc)==1 and refgain==0:
			# gain inferred at only one ancestor, not at the ancestor referenced in the gene tree event.
			gainanc = lgainanc[0]
			modifiedgt = self.checkCountGainCoordinates(gainanc, reftree, eventnode, ost, modifiedgt, silent=silent)
		else:
			# check if some gains are already justified by a transfer event
			copiedost = copy.deepcopy(ost)
			forest = copiedost.pruneTransferedSubtrees(intactCopy=genetree)
			lgainanc = checkExistingTransfers(reftree, lgainanc, forest, silent=silent)

			# Count infered independent gains (transfers), must map them on the gene tree (modified copy)
			if not modifiedgt: modifiedgt = copy.deepcopy(genetree)	# copy of the full gene tree to be returned
			refgain = 0
			reftree.sort(lgainanc, order=0)
			for gainanc in lgainanc:
				if self==gainanc and not refgain:
					# 1 gain inferred at the ancestor referenced in the gene tree event, OK.
					if not silent: print "1 gain @", self.label(), ", justified by the reference event"
					refgain += 1
					continue
				if not silent: print "check unexpected gain @", gainanc.label()
				#~ if not silent: print "copiedost\n", copiedost.newick(ignoreBS=True)
				if gainanc==reftree:
					# gains infered at the root of reference tree are likely gains with mis-annotated coordinates (for instance when using asymmetric Wagner parsimony extremely biased toward losses (equivalent to Dollo parsimony))
					modifiedgt = self.checkCountGainCoordinates(gainanc, reftree, eventnode, ost, modifiedgt, silent=silent)
				else:
					lspe = gainanc.get_leaf_labels()
					print "lspe", lspe
					copiedrecnode = copiedost.map_to_node(lspe, useSpeDict=True, force=True)
					if copiedrecnode:
						copiedsubtree = copiedost
					else:
						for tst in forest:
							if tst:
								copiedrecnode = tst.map_to_node(lspe, useSpeDict=True, force=True)
								copiedsubtree = tst
								if copiedrecnode: break
					print "copiedrecnode", copiedrecnode.get_leaf_labels()
					if copiedrecnode and copiedrecnode.go_father():
						# infer a transfer only if not at root of gene tree (in which case it is an origination (annoated 'gain' if previously was a 'speciation')
						copieddonnode =  copiedrecnode.go_brother()
						forest = [modifiedgt.addTransferEventFromSubtree(reftree, copiedrecnode, copieddonnode=copieddonnode, copiedsubtree=copiedsubtree, skipRootTransfers=False, skipPreviousTransfers=True, skipUncoherentTransfers=True, **kw)] + forest
					else:
						print "could not find or create a gene tree event explaining the gain @", gainanc.label()
					if not silent: print self.fixed_scenario_str()
		return modifiedgt	
						
	def ancestralStatesFromGeneTree(self, genefamtree, method="Dollo", **kw):
		"""infers ancestral presence states (counts) of a gene family by custom parsimony given the reconciliated gene tree.
		
		Return a copy of the input reference tree (self) annotated with the event summary and number of homologs prensent at ancestral nodes 
		according to the annotation of the input gene tree with speciation, transfer and duplication events and their location in reference species tree.
		Treats independently each duplicated/transfered subtree for parsimony inference of presence.
		
		Parsimony method is either built-in Dollo parsimony [method="Dollo"], or use those implemented in Count (Csuros & Miklos, 2009): 
		Dollo [method='CountDollo'] or assymetric Wagner [method='Count.AsymmetricWagner'] parsimony, or maximum-likelihood estimation [method='Posteriors'].
		Count options are to be passed through the string argument 'CountOptions'.
		"""
								
		def check_replacements(ort, ost):	
			"""check in the orthologous gene subfamily tree (ost) 
			if there are transfers among orthologs (gene conversions among close relatives = allelic replacements) 
			and annotate them in the subfamily history tree (ort)
			"""
			for node in ost:
				if node==ost: continue
				if node.transfer():
					rec = ort[node.eventloc()[0]]
					d = rec.getEvents()
					if d['loss']>=d['gain']:
						don = ort[node.eventloc()[2]]
						rec.addEvent('replacement')
						rec.addEvent('receptor')
						don.addEvent('donor')	
				
						
		if not method in ["Dollo", "Count.AsymmetricWagner", "Count.Posteriors"]:
			raise ValueError, "unproper method name %s"%method
		silent = kw.get('silent', True)
		CountOptions = kw.get('CountOptions', '')
		
		if not silent: print "method: '%s', CountOptions: '%s'"%(method, CountOptions)
		
		modifiedgt = None
		
		orthoSubtrees = genefamtree.getOrthoSubtrees(reftree=self, silent=silent)
		leventnodes = orthoSubtrees.keys()
		leventnodes = genefamtree.sort(leventnodes, order=3)
		if not silent: print "\nleventnodes", [n.label() for n in leventnodes]
		orthoReftrees = {}
		for eventnode in leventnodes:
			# reference tree for the considered orthologous subtree
			reftree = copy.deepcopy(self)
			if not silent: print eventnode.label(),eventnode.getdicevent()
			ost = orthoSubtrees[eventnode]
			ancestorlab = eventnode.eventloc()[0]
			refanc = reftree[ancestorlab]
			lspe = ost.listSpecies(asSet=False)
			if not silent: print "\tlspe:", lspe
			reftree.addCountInTree(lspe)
			#~ if not (set(refanc.get_leaf_labels()) & set(lspe)):
				#~ # there are (few) uncoherent transfers annotated in trees that are inferences coming from different Prunier replicates
				#~ print "!!! uncoherent transfer inference at\n", eventnode.label(),eventnode.getdicevent(), '\nost:', ost.get_leaf_labels(), "\n -> treated as a speciation.\n"
			#~ else:
			p1 = reftree.getPhyloProfile(returnLabels=True, leavesOnly=True)
							
			## ancestral content reconstruction
			if method=="Dollo":
				# check for wrong event mapping stemming from presence of conflicting events in gene trees
				mrca = reftree.map_to_node(lspe)
				if not (mrca==refanc or mrca.is_child(refanc)):
					# wrong gene tree event mapping
					modifiedgt = reftree.correct_genetree_event(mrca, eventnode, ost, modifiedgt, silent=silent)
					ancestorlab = mrca.label()
					refanc = mrca
					# temporarily edit gene tree event before assigning a new one
					eventnode.set_eventloc(mrca.label(), item=0)
					eventnode.set_eventloc([mrca.label()], item=1)
				reftree.DolloPars(inventorlab=ancestorlab, silent=silent)	
				
			elif method.startswith("Count"):
				prog = method.split('.')[1]
				dCount_nodeevents = reftree.Count(prog=prog, options=CountOptions, **kw)
				if dCount_nodeevents.keys()!=[refanc]:
					# Count did not infer a simple gain at the reference ancestor
					modifiedgt = refanc.checkCountInference(eventnode, ost, modifiedgt, **kw)
				if not silent: print "\n"
				
			p2 = reftree.getPhyloProfile(returnLabels=True, leavesOnly=True)
			if p1 != p2:
				dspeleaves = ost.dictSpeciesToLeafLabels()
				for spe in p2:
					p2[spe] -= p1[spe]
					if p2[spe]!=0:
						print spe,'; diff count = ', p[spe], '; ost leaves = ', dspeleaves[spe]
					raise IndexError, "ancestral state reconstruction shoud not modify phylogenetic profiles at leaves:\nbefore:\n%s\nafter:\n%s"%(str(p1), str(p2))
			#~ print eventnode.eventtype(), ancestorlab, '@', refanc.getEvents(), '\n'
			if eventnode.duplication():
				refanc.addEvent('duplication')
			#~ elif eventnode.transfer():
				# moved to check_replacements()
				#~ refanc.addEvent('receptor')
				#~ donorlab = eventnode.getdicevent()['eventlocation'][2]
				#~ reftree[donorlab].addEvent('donor')
				#~ refanc.transferedFrom(reftree[donorlab])
					
			# enrich annotation of the whole family with detection of allelic replacements
			check_replacements(reftree, ost)	
		
			reftree.checkGainLossPresenceCoherence()
			orthoReftrees[eventnode] = reftree
							
		# gather event counts from all orthologous lineages
		famHistory = copy.deepcopy(self)
		for eventnode in orthoReftrees:
			ort = orthoReftrees[eventnode]
			famHistory += ort
		
		famHistory.checkGainLossPresenceCoherence()
		
		if not silent: print famHistory.newick(ignoreBS=True, comment='presence'), "\n"
		return famHistory, orthoReftrees, orthoSubtrees, modifiedgt
			
	def DolloParsFromGeneTree(self, genefamtree, silent=True):
		"""infers ancestral presence states (counts) of a gene family by custom Dollo parsimony given the reconciliated gene tree.
		
		Return a copy of the input reference tree (self) annotated with the event summary and number of homologs prensent at ancestral nodes 
		according to the annotation of the input gene tree with speciation, transfer and duplication events and their location in reference species tree.
		Treats independently each duplicated/transfered subtree for parsimony inference of presence.
		"""
		return ancestralStatesFromGeneTree(self, genefamtree, silent=True, method="Dollo")

####################################
### External program call functions

	def Count(self, tmpdir="/panhome/lassalle/tmp", countdir="/panhome/lassalle/count", mem="2048M", prog="AsymmetricWagner", nfrates="", options="", **kw):
		silent = kw.get('silent', True)
		fam = kw.get('fam', "somefamily")
		"""executes Count program (Miklos and Csuros, Mol. Biol. Evol. 2009) to estimate ancestral contents"""
		if prog=="Posteriors" and not nfrates:
			raise ValueError, "Count Posteriors need a rate file"
		nfreftree = "%s/count_reftree.tmp"%(tmpdir)
		nfprofilematrix = "%s/count_profilematrix.tmp"%(tmpdir)
		nfoutput = "%s/count_output"%(tmpdir)
		self.write_newick(nfreftree, mode="write")
		profile = self.writePhyloProfileLine(nfprofilematrix, leavesOnly=True, fam=fam)
		if not silent: print "Count input profile\n", profile
		countcmd = "java -Xmx%s -cp %s/Count.jar ca.umontreal.iro.evolution.genecontent.%s "%(mem, countdir, prog)
		countcmd += " %s %s %s %s"%(options, nfreftree, nfprofilematrix, nfrates)
		if not silent: print countcmd
		larg = countcmd.split()
		foutput = open(nfoutput, 'w')
		subprocess.check_call(larg, stdout=foutput, stderr=subprocess.PIPE)
		foutput.close()
		#~ foutput.seek(0,0) # replace file’s current position to its start
		# read presence information
		if prog=="AsymmetricWagner":
			dnodeevents = self.readCountAsymmetricWagner(nfoutput, **kw)
		elif prog=="Posteriors":
			dnodeevents = self.readCountPosteriors(nfoutput, **kw)
		#~ return phyloprofile
		return dnodeevents
		

			
###########################################
### XML methods

	def scenarioEventXML(self, ind, indstart):
		if self.label(): label = self.label()
		else: label = 'None'
		xml = ""
		if sum(self.__events.values()) > 0:
			xml += "%s<rec:event>"%(indstart+ind)
			for event in self.__events:
				if self.__events[event] > 0:
					xml += "%s<rec:%s>"%(indstart+2*ind, event)
					if event=='receptor':
						xml += "%s<rec:recipientSp>%s</rec:recipientSp>"%(indstart+3*ind, label)
						if sef.donor():
							xml += "%s<rec:originSp>%s</rec:originSp>"%(indstart+3*ind, self.donor().label())
					else:
						xml += "%s<rec:locationSp>%s</rec:locationSp>"%(indstart+3*ind, label)
					xml += "%s</rec:%s>"%(indstart+2*ind, event)
			xml += "%s</rec:event>"%(indstart+ind)
		return xml

