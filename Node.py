#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Newick format parser module.

Manipulates trees in the Newick format. This format is mainly used to
describe phylogenetic binary trees but can have a much wider use. Some
details on the format can be found here:

http://evolution.genetics.washington.edu/phylip/newicktree.html

Examples of Newick format:

* (Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;

* '((A:0.1,B:0.2,C:0.1)ABCnode:0.2,(D:0.4,E:0.1)99:0.1);'
		
Note that bracked delimited nested and unnested comments are ignored.

Sequences can be assigned to each node for simulating evolution.

Derived from module 'tree.py' from 'alfacinha' package by Leonor Palmeira and Laurent Guéguen (http://pbil.univ-lyon1.fr/software/alfacinha/)

Copyright 2007, Leonor Palmeira <palmeira@biomserv.univ-lyon1.fr>;
Copyright 2013, Florent Lassalle <florent.lassalle@univ-lyon1.fr>;
Copyright 2015, Florent Lassalle  <florent.lassalle@ucl.ac.uk>;
Copyright 2016, Florent Lassalle  <f.lassalle@imperial.ac.uk>.
"""

__author__ = "Leonor Palmeira <palmeira@biomserv.univ-lyon1.fr>, Laurent Guéguen <laurent.gueguen@univ-lyon1.fr>, Florent Lassalle  <f.lassalle@imperial.ac.uk>"
__date__ = "10 January 2016"
__credits__ = """Guido van Rossum, for an excellent programming language; Leonor Palmeira and Laurent Guéguen for initiating the Node class"""

import string
import copy
import subprocess
import re
import shutil, os
#import evol

import tree2
from tree2 import svgNode

import numpy as np

def sumlen(*L):
	"""sum branch lengths (floats or None from list L), with None taking value 0 unless only None are summed in which case None is returned""" 
	l = None
	for k in L:
		if k is not None: 
			if l is not None: l += k
			else: l = k
	return l


#######################################################################
#######################################################################
########  Node class: attributes and instance creation

class Node(object):
	"""Node is the smallest unit of definition for a tree.

	A Node instance is linked to its father-Node and its children-Nodes: 
	it defines a sub-tree for which it is the root. 
	Node instances are reccursively contained in their parents up to the 
	tree root-Node, which is equivalent to the whole tree.
	"""

	def __init__(self, branch_lengths=True, keep_comments=False, combrackets='[]', labquotes=False, namesAsNum=False, leafNamesAsNum=False, bootInComm=False, **kw):
		"""Create a Node.
		
		* Keyword argument to build directly the Node from a Newick string: [newick=string] or [nwk=string]
		or from a file: [file=filepath] or [fic=filepath]
		Examples of Newick format:
		- with branch length and leaf labels:
		'(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268,Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;'
		- with branch length, leaf labels and branch supports:
		'((SpeciesA:0.45,SpeciesB:0.42)0.84:0.75,(SpeciesC:0.58,SpeciesD:0.85)0.95:0.115,SpeciesE:0.50);'
		- with branch length, leaf labels, branch supports and bracketted comments:
		'((SpeciesA[Phenotype1]:0.45,SpeciesB[Phenotype1]:0.42)0.84[Phenotype1]:0.75,(SpeciesC[Phenotype1]:0.58,SpeciesD[Phenotype2]:0.85)0.95[Phenotype2]:0.115,SpeciesE[Phenotype3]:0.50);'
		- with branch length, leaf labels, bracketted comments and internal labels instead of branch supports.
		'((SpeciesA[Phenotype1]:0.45,SpeciesB[Phenotype1]:0.42)Clade1[Phenotype1]:0.75,(SpeciesC[Phenotype1]:0.58,SpeciesD[Phenotype2]:0.85)Clade2[Phenotype2]:0.115,SpeciesE[Phenotype3]:0.50);'
		
		* Default behaviour is to discard comments ; can be switched by specifying [keep_comments=True]. 
		Characters used for comment brackets can be changed using [combrackets=charset] (avoid parenthesis !!!).
		Bracketted comments can be located on either end of the label (or support for internal nodes), and these can be mixed:
		'([Phenotype1]SpeciesC:0.58,[Phenotype2]SpeciesD:0.85);' <=> '(SpeciesC[Phenotype1]:0.58,SpeciesD[Phenotype2]:0.85);' <=> '(SpeciesC[Phenotype1]:0.58,[Phenotype2]SpeciesD:0.85);'
		
		* Numerals after a semicolon indicating branch lengths are expected for proper parsing, 
		but in the absence of branch length information (cladogram), 
		this default behavour can be changed by specifying [branch_lengths=False].
		
		* Leaf or internal node labels sometimes happen to be integers, as is the case of ClonalFrame trees (http://www.xavierdidelot.xtreemhost.com/clonalframe.htm).
		Specify [namesAsNum=True] or [leafNamesAsNum=True] for adequate parsing of numeric labels at all nodes or leaf nodes only, respectively.
		
		* A new tree with a multifurcated ('star') toppology can be created using the keyword 'lleaves', specifying the names of the leaves directly attached to the instance (root node).
		
		* Names and values of arbitrary attributes can be passed as a dict (or iterable coercible to a dict) via keyword argument 'setattr'; 
		if these attributes have same name as default ones, the values passed via 'settatr' override the default values. 
		"""
	
		self.__l=kw.get('l') 		# length of the branch to the father-Node
		self.__lab=kw.get('lab')	# Node label
		self.__seq=None 			# Node sequence
		self.__boot=None 			# bootstrap value at the Node
		self.__comment=None 		# nested comment attached to Node
		self.__father=None			# father Node instance
		self.__children=[]			# list of child Node instances
		
		for kwnwk in ['newick', 'nwk']:
			if kwnwk in kw:
				# read from a Newick string
				self.parser(s=kw[kwnwk], branch_lengths=branch_lengths, combrackets=combrackets, labquotes=labquotes, keep_comments=keep_comments, namesAsNum=namesAsNum, leafNamesAsNum=leafNamesAsNum, bootInComm=bootInComm)
				break
		else:
			for kwfic in ['file', 'fic']:
				if kwfic in kw:
					# read from file containing a Newick string
					self.read_nf(a=kw[kwfic], branch_lengths=branch_lengths, combrackets=combrackets, labquotes=labquotes, keep_comments=keep_comments, namesAsNum=namesAsNum, leafNamesAsNum=leafNamesAsNum, bootInComm=bootInComm)
					break
			else:
				if 'lleaves' in kw:
					# generate multifurcated tree with all leaves with names in 'lleaves' attached to the root (self) node
					for leaf in kw['lleaves']:
						c = self.newnode()
						c.add_label(leaf)
						self.link_child(c, newlen=0, newboot=0)
		if 'setattr' in kw:
			try:
				# try to coerce the object to a dict
				dattr = dict(kw['setattr'])
			except ValueError, e:
				raise ValueError, "encountered when setting attribute names and values for Node instance, expecting a dict or object coercible to a dict: '%s'"%str(e)
			self.__dict__.update(dattr)

			
#####################################################
############## methods for access to object atributes:

	def newnode(self, branch_lengths=True, keep_comments=False, **kw):
		"""class-specific instance generator for construction of trees of Nodes"""
		return Node(branch_lengths=branch_lengths, keep_comments=keep_comments, **kw)

	def labelgetnode(self, n, mustMatch=False):
		"""Return the first child of Node on a pre-order traversal labelled with intput string 'n'.
		
		If mustMatch=True, raises an IndexError when no node if found; default is False, as used in __getitem__().
		Assume type of input is controlled ahead of call (as in __getitem__()) so does not check for efficiency saving. 
		"""
		if self.__lab==n: return self
		for c in self.__children:
			if c.label()==n:
				return c
			else:
				t=c.labelgetnode(n, mustMatch=mustMatch)  # (recursive function)
				if t:
					return t
		if not mustMatch: return None
		else: return IndexError, "no node labelled %s in tree"%n

	def __getitem__(self,n):
		"""<==> self[n]. Return the Node with input label. If input is an iterable returning labels or node objects, find their MRCA"""
		if type(n) is str:
			return self.labelgetnode(n)
		#~ elif type(n) in (tuple, list):
		elif hasattr(type(n), '__iter__'):
			return self.coalesce(n)
		else:
			raise TypeError, "unexpected type %s for key: %s"%(str(type(n)), repr(n))
		return None
		
	def __getattr__(self, attr):
		if attr=='boot': return self.__boot
		elif attr=='l': return self.__l
		elif attr in ['lab']: return self.__lab
		elif attr in ['fat', 'father']: return self.__father
		elif attr=='children': return self.__children
		else: raise AttributeError
	
	def __str__(self):
		"""Return printable string in Newick format of the sub-tree defined by the Node."""
		return self.newick()
		
	def __repr__(self):
		"""Return printable string naming the class, number of children nodes and label the Node."""
		return "<%s object: %d nodes; label:%s>"%(str(type(self)).split("'")[1], len(self), repr(self.label()))
		
	def __iter__(self):
		#~ return self.generator()
		return self.preordertraversal_generator()
	
     # Python 3: 
	def __bool__(self): 
		"""Boolean value of an instance of this class (True). 
	 
		NB: If this method is not defined, but ``__len__``  is, then the object 
		is considered true if the result of ``__len__()`` is nonzero. We want 
		Node instances to always be considered True. 
		(stolen from BioPython)
		__len__() points to get_all_children() which should always return a
		non-empty list with at least [self]; this saves time though.
		"""
		return True
	
	# Python 2: 
	__nonzero__ = __bool__ 
		
	def __len__(self):
		return len(self.get_all_children())
		
	def __copy__(self):
		"""produces a shalow copy of the node, with new deep copies of all attributes but shalow copies of the links to father and children
		
		the new copy keep knowing where it is located in the tree, but the other nodes in the tree will not know it exists (they know the original [self] instance)
		"""
		newinstance = self.newnode()
		for attr in self.__dict__:
			if not attr in ['_Node__children', '_Node__father']:
				newinstance.__dict__[attr] = copy.deepcopy(self.__dict__[attr])
			elif attr == '_Node__children':
				newinstance.__dict__[attr] = copy.copy(self.__dict__[attr])	# create a new list container, but with same element pointers in it
			elif attr == '_Node__father':
				newinstance.__dict__[attr] = self.__dict__[attr]	# just use the same pointer
		return newinstance
		
	def deepcopybelow(self, keep_lg=False, add_ref_attr=False, shallow_copy_attr=[]):
		"""return similar subtree object than pop(self) result, but as a copy (not the same object, and the original instance remains untouched attached to the tree).
		
		compared to regular __copy__(), avoids duplicating the tree above. On a root node and with keep_lg=True, equivalent to copy.deepcopy(self).
		"""
		newinstance = self.newnode()
		for attr in self.__dict__:
			if (attr == '_Node__father') or (attr == '_Node__l' and not keep_lg):
				newinstance.__dict__[attr] = None
			elif attr in shallow_copy_attr:
				newinstance.__dict__[attr] = self.__dict__[attr]
			else:
				newinstance.__dict__[attr] = copy.deepcopy(self.__dict__[attr])
		for child in newinstance.get_children():
			# associate the copied children with the newinstance, rather than an unlinked deep copy of their father
			child.change_father(newinstance, newlen=child.lg(), newboot=child.bs())
		if add_ref_attr:
			# tag copy nodes with reference to their original counterpart (add an attribute)
			newnodes = newinstance.get_all_children()
			oldnodes = self.get_all_children()
			# rely on the fact that object are structured equally so that ordering of nodes in lists generated by pre-order traversal will be the same
			for i in range(len(newnodes)):
				newnodes[i].ref = oldnodes[i]
		return newinstance
		
	def deepcopyabove(self, shallow_copy_attr=[]):
		"""return similar root-side tree object than self after calling pop(), but as a copy (not the same object, and the original instance remains untouched with subtree attached).
		
		compared to regular __copy__(), avoids duplicating the tree below.
		"""
		newinstance = self.newnode()
		for attr in self.__dict__:
			if attr == '_Node__children':
				newinstance.__dict__[attr] = []
			elif attr in shallow_copy_attr:
				newinstance.__dict__[attr] = self.__dict__[attr]
			else:
				newinstance.__dict__[attr] = copy.deepcopy(self.__dict__[attr])
		# associate the copied father with the newinstance, rather than with an unlinked deep copy of its child
		oldfat = newinstance.go_father()
		if oldfat:
			# find the apprropriate child toremove
			children = oldfat.get_children()
			appro = []
			for i in range(len(children)):
				child = children[i]
				if child.label()==newinstance.label() and child.lg()==newinstance.lg() and child.bs()==newinstance.bs(): appro.append(i)
			if len(appro)<1: raise IndexError, "could find no appropriate deep copy of the node to replace with the new instance"
			elif len(appro)>1: raise IndexError, "found several candidate deep copies of the node to replace with the new instance, cannot choose"
			else:
				if children[i] is newinstance:
					print "deep copy of the child is already the newinstance; surprising, but that's fine: pass"
				else:
					oldfat.rm_child(children[i])
					oldfat.add_child(newinstance)
		return newinstance
		
		
		
	# defining a cmp method is dangerous as many tests involve testing identity of nodes using '==' operator 
	# (and by default comparing values yielded by hash(object)), instead of clean 'is' operator
	#~ def __cmp__(self, other):
		#~ """to sort sequences with deep nodes first (should prefer the use of sort(key=fun()))"""
		#~ sd = self.depth()
		#~ so = other.depth()
		#~ if sd > so: return -1
		#~ elif sd == so: return 0
		#~ else: return 1
		
	#~ def generator(self):
		#~ children = self.get_all_children()
		#~ for child in children:
			#~ yield child

	def sequence(self):
		"""Return the Sequence of the Node."""
		return self.__seq

	def lg(self):
		"""Return the length of the edge to the father."""
		return self.__l

	def bs(self):
		"""Return the bootstrap value at Node."""
		return self.__boot

	def label(self, comments=False):
		"""Return the label."""
		#~ return self.__lab
		text = self.__lab
		if comments==True and self.__comment: text += '[%s]'%self.__comment
		return text

	def comment(self):
		"""Return the comment."""
		return self.__comment
		
	def getattr_down_n_nodes(self, attr, n, ommitself=False, stopatnodes=[], leafval='same'):
		"""return a list of the chosen attribute value for all branches/nodes down the tree from the focal node.
		
		Operate reccursively over n children (maximum, stops at leaves) if n >= 0,
		or until exploration of all leaves if n < 0.
		"""
		lval = []
		if not ommitself:
			if leafval=='same': val =  getattr(self, attr)
			else: val = leafval if self.is_leaf() else getattr(self, attr)
			if callable(val): lval.append(val())
			else: lval.append(val)
		if n!=0:
			for sn in self.__children:
				if not (sn in stopatnodes):
					lval += sn.getattr_down_n_nodes(attr, n-1, stopatnodes=stopatnodes)
		return lval
			
	def preordertraversal_generator(self):
		"""recursively yield all the nodes below the Node, including itself, in a pre-order traversal"""
		yield self
		for child in self.__children:
			for rechild in child.preordertraversal_generator():
				yield rechild
				
	def get_preordertraversal_children(self):
		"""Return the list of all nodes below the Node, including itself, in a pre-order traversal"""
		return list(self.preordertraversal_generator())
				
	#~ def get_all_children(self):
		#~ """Return the list of all nodes below the Node, including itself, in a pre-order traversal"""
		#~ a=[self]
		#~ for i in self.__children:
			#~ a+=i.get_all_children()
		#~ return a
	# rather use an alias for code consistency
	get_all_children = get_preordertraversal_children
				
	def get_all_children_cache(self, cachedict={}):
		"""Return the list of all nodes below the Node, including itself, in a pre-order traversal"""
		a=[self]
		for i in self.__children:
			if i in cachedict:
				a+=cachedict[i]
			else:
				a+=i.get_all_children_cache(cachedict)
		return a
	
	def postordertraversal_generator(self, deepfirst=False, deeplast=False):
		"""recursively yield all the nodes below the Node, including itself in a post-order traversal"""
		children = self.__children
		if deepfirst: children.sort(key=lambda x: x.max_leaf_depth())
		elif deeplast: children.sort(key=lambda x: -1*x.max_leaf_depth())
		for child in children:
			for rechild in child.postordertraversal_generator():
				yield rechild
		yield self
		
	def get_postordertraversal_children(self, deepfirst=False, deeplast=False):
		"""Return the list of all nodes below the Node, including itself, in a post-order traversal"""
		#~ a=[]
		#~ children = self.__children
		#~ if deepfirst: children.sort(key=lambda x: x.max_leaf_depth())
		#~ elif deeplast: children.sort(key=lambda x: -1*x.max_leaf_depth())
		#~ if children!=[]:
			#~ for i in children:
				#~ a+=i.get_postordertraversal_children()
		#~ a += [self]
		#~ return a
		return list(self.postordertraversal_generator(deepfirst=deepfirst, deeplast=deeplast))
	
	def midordertraversal_generator(self, righttfirst=False):
		"""recursively yield all the nodes below the Node, including itself in a mid-order traversal; 
		
		i.e. the left child first, then the Node, then the other children from left to right
		"""
		children = self.__children
		ordch = range(len(children))
		if righttfirst: ordch.reverse() # use list of index to not change order of self.__children list!
		if children:
			firstch = children[ordch[0]]
			for rechild in firstch.midordertraversal_generator(righttfirst=righttfirst):
				yield rechild
		yield self
		for i in ordch[1:]:
			for rechild in children[i].midordertraversal_generator(righttfirst=righttfirst):
				yield rechild
	
	def get_midordertraversal_children(self, righttfirst=False):
		"""Return the list of all nodes below the Node, including itself, in a mid-order traversal"""
		#~ a=[]
		#~ children = self.__children
		#~ ordch = range(len(children))
		#~ if righttfirst: ordch.reverse()
		#~ if children:
			#~ firstch = children[ordch[0]]
			#~ a+=firstch.get_midordertraversal_children(righttfirst=righttfirst)
		#~ a += [self]
		#~ for i in ordch[1:]:
			#~ a+=children[i].get_midordertraversal_children(righttfirst=righttfirst)
		#~ return a
		return list(self.midordertraversal_generator(righttfirst=righttfirst))
		
	def sorted_generator(self, order=1):
		"""yield all the nodes below the Node, including itself, based on a specified order
		
		order=-1: ordered by increasing depth (i.e. decreasing node distance from root).
		order= 0: pre-oder traversal (classic root-to-leaves exploration)
		order= 1: ordered by decreasing depth (i.e. increasing node distance from root).
		order= 2: post-oder traversal (exploration of each group of leaves, then the nodes above ; compatible with Count's rate files node enumeration)
		order= 3: post-oder traversal (exploration of each group of leaves, then the nodes above ; compatible with Count's rate files node enumeration) with always exploring the deepest node first
		order= 4: mid-order traversal (first the left child(ren), then the father, then the right child(ren))
		order= 5: mid-order traversal (first the right child(ren), then the father, then the left child(ren))
		order= 6: post-oder traversal (exploration of each group of leaves, then the nodes above ; compatible with Count's rate files node enumeration) with always exploring the deepest node last
		"""
		if order in [-1, 1]:
			a = self.get_all_children()
			a.sort(key=lambda x: x.depth()*order)
			for n in a: yield n
		elif order == 0:
			for n in self.preordertraversal_generator(): yield n
		elif order == 2:
			for n in self.postordertraversal_generator(): yield n
		elif order == 3:
			for n in self.postordertraversal_generator(deepfirst=True): yield n
		elif order == 4:
			for n in self.midordertraversal_generator(righttfirst=False): yield n
		elif order == 5:
			for n in self.midordertraversal_generator(righttfirst=True): yield n
		elif order == 6:
			for n in self.postordertraversal_generator(deeplast=True): yield n

	def get_sorted_children(self, order=1):
		"""Return the list of all nodes below the Node, including itself
		
		order=-1: ordered by increasing depth (i.e. decreasing node distance from root).
		order= 0: pre-oder traversal (classic root-to-leaves exploration)
		order= 1: ordered by decreasing depth (i.e. increasing node distance from root).
		order= 2: post-oder traversal (exploration of each group of leaves, then the nodes above ; compatible with Count's rate files node enumeration)
		order= 3: post-oder traversal (exploration of each group of leaves, then the nodes above ; compatible with Count's rate files node enumeration) with always exploring the deepest node first
		order= 4: mid-order traversal (first the left child(ren), then the father, then the right child(ren))
		order= 5: mid-order traversal (first the right child(ren), then the father, then the left child(ren))
		order= 6: post-oder traversal (exploration of each group of leaves, then the nodes above ; compatible with Count's rate files node enumeration) with always exploring the deepest node last
		"""
		
		if order in [-1, 0, 1]:
			a = self.get_all_children()
			a.sort(key=lambda x: x.depth()*order)
		elif order == 2:
			a = self.get_postordertraversal_children()
		elif order == 3:
			a = self.get_postordertraversal_children(deepfirst=True)
		elif order == 4:
			a = self.get_midordertraversal_children(righttfirst=False)
		elif order == 5:
			a = self.get_midordertraversal_children(righttfirst=True)
		elif order == 6:
			a = self.get_postordertraversal_children(deeplast=True)
		return a
		
	def get_comtemporary_branches(self, t, dp_drootdist=None, tagwithlabels=True):
		"""retrieves all the branches of the tree encompassing the time t from the root"""
		if not dp_drootdist: drootdist={}
		else: drootdist=dp_drootdist
		cb = []
		if tagwithlabels:
			stag = self.label()
			if self.__father: ftag = self.__father.label()
			else: ftag = None
		else:
			stag = self
			ftag = self.__father
		# dynamic programing to obtain node-to-root distances as the tree is explored; avoids individual parcours of each node's lineage
		if stag in drootdist: raise IndexError, "stag '%s' already recorded in drootdist"%stag # should not happen, test could be removed to improve performance
		d = drootdist.setdefault(stag, drootdist.get(ftag, 0)+self.__l)
		if d > t:
			cb += [self]
		else:
			for c in self.__children:
				# recursive filling of the list
				cb += c.get_comtemporary_branches(t, dp_drootdist=drootdist, tagwithlabels=tagwithlabels)
		return cb
		
	def get_branch_chronology(self, ltimes, drootdist={}, dtimeindexes={}, tagwithlabels=True):
		"""for each time t in list 'ltimes', retrieves all the branches of the tree encompassing the time t from the root; return a dictionary of times to lists of branches"""
		dtimebranches = {}
		if tagwithlabels:
			stag = self.label()
			if self.__father: ftag = self.__father.label()
			else: ftag = None
		else:
			stag = self
			ftag = self.__father
		# dynamic programing to obtain node-to-root distances as the tree is explored; avoids individual parcours of each node's lineage
		d = drootdist.setdefault(stag, drootdist.get(ftag, 0)+self.__l)
		# similar for indexes of the time list; avoids exploring times that are known to be younger than the node's father (and hence than the node)
		itime = dtimeindexes.get(ftag, 0)
		
		while itime<len(ltimes):
			if d > ltimes[itime]:
				dtimebranches.setdefault(ltimes[itime], []).append(self)
				itime += 1
			else:
				break
		dtimeindexes[stag] = itime
		for c in self.__children:
			# recursive filling of the dict
			tree2.updateappend(dtimebranches, c.get_branch_chronology(ltimes, drootdist=drootdist, dtimeindexes=dtimeindexes, tagwithlabels=tagwithlabels))
		return dtimebranches
		
	def get_all_parents(self):
		"""Return the list of all nodes above the Node, excluding itself, ordered by increasing depth (i.e. decreasing node distance from root)."""
		lparents = []
		f = self.go_father()
		while f:
			lparents.append(f)
			f = f.go_father()
		return lparents
		
	def sort(self, lnodes, order=1, apendOutSelf=True, labels=False):
		"""sorts a list of nodes contained in self in different orders (cf. Node.get_sorted_children())"""
		orderedchildren = self.get_sorted_children(order=order)
		ln = copy.copy(lnodes)
		l = []
		for oc in orderedchildren:
			if not labels: o = oc
			else: o = oc.label()
			while o in ln :
				l.append(ln.pop(ln.index(o))) 
		if apendOutSelf:
			l += ln
		elif lnodes:
			raise IndexError, "nodes in %s are not in %s (self)"%(str([n.label() for n in lnodes]), self.label())
		return l
		
				
				
##############################################################
### HOGENOM Species/Protein identifiers correspondency methods

	def get_prot_labels_from_aln(self, aln):
		outtree = copy.deepcopy(self)
		## retrieving of full protein labels from alignments
		lleaves = self.get_leaf_labels()	
		dspe_prot = {}
		faln = open(aln, 'r')
		aln = faln.read()
		faln.close()		
		for leaf in lleaves:
			s = re.compile(leaf+"_._.+?[ \n]") 
			ss = s.search(aln)
			sprot = ss.group().rstrip(" \n")
			dspe_prot[leaf] = sprot
		for leaf in lleaves:
			outtree[leaf].edit_label(dspe_prot[leaf])
		return outtree
		
	def listSpecies(self, llab=None, ignoreTransfers=False, asSet=True, splitparam=('_',1), **kw):
		"""Return list of identifiers of species present in leaves 
		
		by default, split is done at 1st occurence of '_'.
		I one assume labels are under the HOGENOM model SPECIES_NUMREPLICON_PROTID -> return SPECIES.
		Separator string and number in string can be specified otherwise. 
		
		can be restricted to a particular species set.
		"""
		lspe = []
		lleaves = self.get_leaves()
		if isinstance(self, tree2.GeneTree) and ignoreTransfers:
			ltransleaves = []
			for node in self:
				if node.transfer():
					ltransleaves += self.idgetnode(node.transferchild()).get_leaves()
			lleaves = list(set(lleaves) - set(ltransleaves))
		for leaf in lleaves:
			if (not llab) or (leaf.label() in llab):
				spe = leaf.label().split(splitparam[0], splitparam[1])[0]
				lspe.append(spe)
		if asSet: lspe = list(set(lspe))
		return lspe
		
	def dictLeafLabelsToSpecies(self, llab=None, splitparam=('_',1)):
		"""Return dictionary of leaf labels (under the HOGENOM model SPECIES_NUMREPLICON_PROTID) to species identifiers
			
		can be restricted to a particular species set	
		Needs SPECIES[_ANYTHING] type (HOGENOM-type) labels at leaves."""
		d = {}
		for leaf in self.get_leaves():
			if (not llab) or (leaf.label() in llab):
				spe = leaf.label().split(splitparam[0], splitparam[1])[0]
				d[leaf.label()] = spe
		return d
		 		
	def dictSpeciesToLeafLabels(self, lspe=None, catValues=False):
		"""Return dictionary of species identifiers (under the HOGENOM model SPECIES_NUMREPLICON_PROTID) to leaf labels
			
		can be restricted to a particular leaf label set
		Needs SPECIES[_ANYTHING] type (HOGENOM-type) labels at leaves."""
		if not catValues: d = {}
		else: d = []
		for leaf in self.get_leaves():
			spe = leaf.label().split('_',1)[0]
			if (not lspe) or (spe in lspe):
				if not catValues: d[spe] = d.setdefault(spe, []) + [leaf.label()]
				else: d.append(leaf.label())
		return d

	def dictLeavesToSpecies(self, llab=None):
		"""Return dictionary of Node objects (leaf) (under the HOGENOM model SPECIES_NUMREPLICON_PROTID) to species identifiers
			
		can be restricted to a particular species set	
		Needs SPECIES[_ANYTHING] type (HOGENOM-type) labels at leaves."""
		d = {}
		for leaf in self.get_leaves():
			if (not llab) or (leaf.label() in llab):
				spe = leaf.label().split('_',1)[0]
				d[leaf] = spe
		return d
		 		
	def dictSpeciesToLeaves(self, lspe=None, catValues=False):
		"""Return dictionary of species identifiers (under the HOGENOM model SPECIES_NUMREPLICON_PROTID) to Node objects (leaf)
			
		can be restricted to a particular leaf label set
		Needs SPECIES[_ANYTHING] type (HOGENOM-type) labels at leaves."""
		if not catValues: d = {}
		else: d = []
		for leaf in self.get_leaves():
			spe = leaf.label().split('_',1)[0]
			if (not lspe) or (spe in lspe):
				if not catValues: d[spe] = d.setdefault(spe, []) + [leaf]
				else: d.append(leaf)
		return d
		
#####################################################
############## methods for editing object atributes:

	def __setattr__(self, name, value):
		self.__dict__[name] = value
		
	def set_lg(self,l):
		"""Set the length of the edge to the father to $2 if it is >=0."""
		try:
			lg = float(l)
		except TypeError, e:
			if l==None:
				lg = None
			else:
				raise TypeError, e
		self.__l = lg
			
	def set_bs(self, bs):
		try:
			boot = float(bs)
		except TypeError, e:
			if bs==None:
				boot = None
			else:
				raise TypeError, e
		self.__boot = boot

	def add_label(self, lab):
		try:
			label = str(lab)
		except TypeError, e:			
			if lab==None:
				label = None
			else:
				raise TypeError, e
		self.__lab = label
			
	def edit_label(self, label, mode="write", sep=", "):
		"""edits the label attached to the node"""
		if mode in ["w", "write"] or (self.__lab in ["", None]):
			self.__lab = str(label)
		elif mode in ["a", "append"]:
			self.__lab += sep+str(label)
		else:
			raise ValueError
	
	def edit_all_labels(self, label, mode="append", sep="-"):
		"""edits the label attached to all the nodes below the node, including itself; defaults to appending a string"""
		for node in self:
			node.edit_label(label, mode=mode, sep=sep)
		
	def edit_comment(self, comment, mode="write"):
		"""edits the comment attached to the node"""
		if mode in ["w", "write"] or not self.__comment:
			self.__comment = str(comment)
		elif mode in ["a", "append"]:
			self.__comment += ', '+str(comment)
		else:
			raise ValueError


	@staticmethod
	def collapse(f, c, gf=None, tellReplacingNode=False, keepRefToFather=False, silent=True):
		"""how to operate fusion of two given nodes, either by removing the intermediate node (newpop) or mutating the intermediate into the top one (oldpop)"""
		def oldpop(f, c):				
			""" old implementation: `mutation' of the father node into its non-poped child"""
			if not silent: print "oldpop(%s, %s)"%(str(f.label()), str(c.label()))
			df = f.__dict__
			dc = c.__dict__
			for attr in df:
				if attr=='_Node__l':
					f.set_lg(sumlen(f.lg(),c.lg()))
				elif attr=='_Node__boot':
					f.set_bs(max(f.bs(), c.bs()))
				elif attr=='_Node__children':
					f._Node__children = []
					for gc in c.get_children():
						f.link_child(gc, newlen=gc.lg(), newboot=gc.bs(), silent=silent)
						if gc.go_father() is c: raise IndexError
					c._Node__children = []
				elif attr=='_Node__father':
					pass
				else:
					df[attr] = dc[attr]
			del c
					
		
		def newpop(gf, f, c):
			""" new implementation : no `mutation' of the node object, the father node is disconnected from above and below, and the child is reconnected to the grand-father"""
			if not silent: print "newpop(%s, %s, %s)"%(str(gf.label()), str(f.label()), str(c.label()))
			gf.add_child(c, silent=silent)
			c.change_father(gf, newlen=sumlen(f.lg(),c.lg()), newboot=max(f.bs(), c.bs()), silent=silent)
			gf.rm_child(f, silent=silent)
			f.rm_child(c, silent=silent)
		
		if keepRefToFather or (not gf):
			# poping under the root or another node to keep referenced
			if tellReplacingNode:
				return (c.label(), f.label())
			else:
				oldpop(f, c)				
		else:
			if tellReplacingNode:
				return (f.label(), c.label())
			else:
				newpop(gf, f, c)
		return None

	def pop(self, name, keepRefToFather=False, tellReplacingNode=False, noCollapse=False, verbose=False):
		"""Return the Node with this name, and removes it from the tree.
		
		New implementation: father node is completely disconected from the tree, its other child are grafted to the grand-father
							if the father is the root (no grand-father), uses old implementation
		Old implementation: If node removal leaves only one node under the father node, the father node takes all the attributes from the brother node.
							> used for poping node under the root or under a node the user would keep referenced
		In both implementations, removing a node leads to merge two branches: their lengths are summed and the highest support is kept.
		Caution : function used in Node.reRoot() (coded while using old implementation).
		
		if tellReplacingNode is True, will not change the tree but will return the pair of node to be collapsed if a node is poped as a tuple (removed_node, staying_node).
		
		if noCollapse is True, the  father node of the node to pop is only deleted if it ends up with no child leading to an original leaf.
		Otherwise, the node is kept with up to one child, leading to branches not being fused / nodes not being collapsed into a single one
		but be made of several segment delimited by nodes. Useful in the xODT/ALE representation of speciation-loss events in reconciled gen trees.
		"""
		
		# find node to pop 'np'
		if isinstance(name, type(self)):
			np = name
		elif isinstance(name, str):
			np = self[name]
		else:
			raise TypeError, "element to be poped must be a %s object or a string (label) corresponding to a node in self"%type(self)
		# NB: if types are right but no node from self is matched, None is returned
		if np:
			f=np.go_father()
			if f:
				if not tellReplacingNode: 
					f.unlink_child(np)
				gf = f.go_father()
				if noCollapse:
					if verbose: print "noCollapse on node %s"%str(f.label())
					if f.nb_children()==0:
						# a lineage leading to false leaf (i.e. a terminal node that is not a true tip) has been created;
						# go up the lineage and cut off fals leaf nodes until meeting a node with more than 1 child (i.e. with another descendant lineage than the false leaf)
						while gf and gf.nb_children()==1:
							f = gf
							gf = f.go_father()
						if gf:
							gf.unlink_child(f) # (great)-grandfather node 'gf' gets rid of the leafless lineage
						else:
							pass # reached the root; proceed as usual (collapse nodes with oldpop)
					else:
						return np # stop here
				else:
					# find 'c', the brother node of 'np', to collapse it with f
					if f.nb_children()==1 or tellReplacingNode: # edit rule (tellReplacingNode and f.nb_children()==2) into just (tellReplacingNode) to restore support for multifurcated trees... must test if that does not lead to bugs
						for child in f.get_children():
							if child is not np:
								c = child
								break
						else:
							raise IndexError, "where is the brother of node to pop 'np'?"
						collapsed = self.collapse(f, c, gf, tellReplacingNode=tellReplacingNode, silent=(not verbose))
					if tellReplacingNode: return collapsed
			else:
				if tellReplacingNode: return (np.label(), np.label())
		return np
		
	def dictCollapsed(self, prunedleaves, everyStep=False, new2old=True, silent=True):
		"""check the pattern of node collapsing if one would prune those leaves"""
		dcollapsed = {} 		# d[staying] = removed
		for leaf in prunedleaves:
			collapsednodes = self.pop(leaf, tellReplacingNode=True)
			# follow the series of pruning steps
			while collapsednodes[0] in dcollapsed.values():
				collapsednodes = self.pop(collapsednodes[0], tellReplacingNode=True)
			dcollapsed[collapsednodes[1]] = collapsednodes[0]
		return dcollapsed
		
	def restrictToLeaves(self, lleaves, useSpeDict=False, force=False, returnCopy=True):
		"""returns a copy of the tree restricted to the input leaf set"""
		if not returnCopy:
			st = self
		else:
			c = copy.deepcopy(self)
			# maps to the common ancestor of all leaves
			mrca = c.map_to_node(lleaves, useSpeDict=useSpeDict, force=force)
			# extract this subtree
			if mrca is c: st = c
			elif mrca:    st = c.pop(mrca)
			else:         return None
		# remove the other leaves
		for leaflab in st.get_leaf_labels():
			if not useSpeDict: leaf = leaflab
			else: leaf = leaflab.split('_')[0]
			if leaf not in lleaves: st.pop(leaflab)
		return st

	def complete_label(self, prefix='N', labels=None, force=True):
		"""gives label to an internal node given a set of pre-existing labels in the tree"""
		if not labels: labels = []
		children = self.get_sorted_children()
		for c in children:	  # first builds a list of existing labels to avoid redundancy of new names
			cl = c.label()
			if (cl!="") and (not cl in labels):
				labels.append(cl)
		n = 1
		if self.label()=="" or force:
			lab = "%s%d"%(prefix, n)
			while lab in labels:
				n += 1
				lab = "%s%d"%(prefix, n)
			self.add_label(lab)

	def complete_internal_labels(self, prefix='N', labels=None, ffel=False, force=False, exclude=[], excludeLeaves=False, onlyLeaves=False, silent=True, order=1, fast=False):
		"""give a numeric label following a specified order to ALL nodes that lack a label; naturally skips editing the leaf that should be already labelled.
		
		a list of nodes of self to not edit can be passed as 'exclude'.
		a list of labels can be passed as 'labels' to extend the list of labels to avoid (in addition to the labels already present in the tree).
		
		if 'force' is True, all nodes are renamed, including leaves (unless 'excludeLeaves' is set to True).
		if 'fast' is True, no lookup of what is already present in the tree is done; 
		safer used in combination with 'force=True' (and 'excludeLeaves=True'), as provided by combo option 'ffel'.
		"""
		if isinstance(order, list):
			# orderred node list is given
			children = order
		else:
			children = self.sorted_generator(order=order)
		if not (fast or ffel):
			# first builds a list of existing labels to avoid redundancy of new names
			if not labels: labs = []
			else: labs = labels
			for c in children:
				cl = c.label()
				if (cl not in ["", None]) and (not cl in labs):
					labs.append(cl)
		n = 1
		lab = "%s%d"%(prefix, n)
		for c in children:
			if c in exclude: continue
			if (excludeLeaves or ffel) and c.is_leaf(): continue
			if onlyLeaves and not c.is_leaf(): continue
			if (force or ffel) or (c.label() in ["", None]):
				if (fast or ffel):
					lab = "%s%d"%(prefix, n)
					n += 1
				else:
					while lab in labs:
						n += 1
						lab = "%s%d"%(prefix, n)
				c.add_label(lab)
				if not (fast or ffel): labs.append(lab)
				n += 1
				if not silent: print lab
				
	def check_unique_labelling(self):
		nodelabs = self.get_children()
		assert len(set(nodelabs))==len(nodelabs)
		
#####################################################
############## methods for access to object's child atributes:

	def get_children(self):
		"""Return the list of direct children of the Node"""
		return self.__children
		
	def get_target_child(self, node=None, label=None, nodeid=None):
		if node:
			n = node
		elif label:
			n = self[label]
		elif nodeid:
			if isinstance(self, tree2.AnnotatedNode):
				n = self.idgetnode(nodeid)
			else:
				raise ValueError, "need an AnnotatedNode instance to fetch a node from an nodeid"
		else:
			raise ValueError, "need a clue (Nde instance, label or nodeid) to find the target child"
		for c in self.get_children():
			if n in c:
				return c

	def nb_children(self):
		"""Return the number of direct children."""
		return len(self.__children)

	def nb_all_children(self):
		"""Return the number of descendants."""
		return len(self.get_all_children())

	def children_labels(self):
		"""Return the list of direct child labels."""
		llc = []
		for c in self.__children:
			llc.append(c.label())
		return llc

	def children_comments(self):
		"""Return the list of child comments."""
		lcc = []
		for c in self.__children:
			lcc.append(c.comment())
		return lcc
		
	def get_children_labels(self, sorted=1):
		"""Return the list of labels of all nodes below the Node, including itself"""
		if sorted:
			children = self.get_sorted_children(order=sorted)
		else :
			children = self.get_all_children()
		l = []
		for c in children:
			l.append(c.label())
		return l
		
	def project_to_children_labels(self, lspe):
		l = []
		for spe in lspe:
			l += self[spe].get_children_labels(sorted=False)
		l += lspe
		return set(l)
		
	def get_children_node_ids(self):
		lid = []
		self.complete_node_ids()
		for node in self:
			lid.append(node.nodeid())
		return lid
			
#####################################################
############## methods for access to object's leaf atributes:
		
	def nb_leaves(self):
		"""Return the number of leaves under this node."""
		return len(self.get_leaves())	
			
	def iter_leaf_labels(self, comments=False):
		"""Return an iterator of labels of the leaves under the Node, recursing down the tree in a post-order traversal."""
		if self.__children:
			for c in self.__children:
				for leaflab in c.iter_leaf_labels(comments=comments):
					yield leaflab
		else:
			yield self.label(comments=comments)
		

	def get_leaf_labels(self, comments=False):
		"""Return the list of labels of the leaves defined under the Node, following a post-order traversal."""
		#~ a=[]
		#~ if self.__children!=[]:
			#~ for c in self.__children:
				#~ a += c.get_leaf_labels(comments=comments)
		#~ else:
			#~ a.append(self.label(comments=comments))
		#~ return a
		return [leaflab for leaflab in self.iter_leaf_labels()]
		
	def get_sorted_leaf_labels(self, comments=False, order=1):
		return [node.label(comments=comments) for node in self.get_sorted_children(order=order) if node.is_leaf()]
	
	def get_leaf_sequences(self):
		"""Return the list of sequences of the leaves defined by the Node."""
		a=[]
		if self.__children!=[]:
			for i in self.__children:
				a+=i.get_leaf_sequences()
		else:
			a+=[self.__seq]
		return a

	def get_leaves(self):
		"""Return the list of leaves composing the clade defined by the Node."""
		a=[]
		if self.__children!=[]:
			for i in self.__children:
				a+=i.get_leaves()
		else:
			a+=[self]
		return a
	
	def get_one_leaf(self):
		a = self
		while a.children:
			a = a.get_children()[0]
		else:
			return a
	
	def iter_leaves(self):
		"""Return iterator of leaves composing the clade defined by the Node."""
		if self.__children!=[]:
			for i in self.__children:
				#~ for leaf in i.get_leaves(): yield leaf
				for leaf in i.iter_leaves(): yield leaf
		else:
			yield self	
	
#####################################
############## file input methods:		
	
	def read_nf(self, a, branch_lengths=True, keep_comments=False, combrackets='[]', labquotes=False, namesAsNum=False, leafNamesAsNum=False, bootInComm=False):
		"""Read the $1 file containing a unique tree in Newick format and builds the Node from it."""

		if type(a)==str:
			f=open(a,'r')
			l=f.readline()
			f.close()
			self.parser(l, branch_lengths, keep_comments, combrackets=combrackets, labquotes=labquotes, namesAsNum=namesAsNum, leafNamesAsNum=leafNamesAsNum, bootInComm=bootInComm)
		else:
			raise ValueError, "Invalid file name."
		return


#####################################
############## file output methods:
		
	def write_newick(self, nf, **kw): #, ignoreBS=False, comment='comment', branch_lengths=True, unrooted=False
		"""Write the Node in Newick format in $1 file.

		Keyword argument to determine whether 'append' or 'write'
		(overwrites) mode is turned on: [mode=string]."""
		tree2.write_to_file(self.newick(**kw), nf, **kw)
	
	def write_phylip(self, nf, **kw):
		"""Write the Node in Phylip format in $1 file.

		Keyword argument to determine whether 'append' or 'write'
		(overwrites) mode is turned on: [mode=string]."""
		tree2.write_to_file(self.phylip(), nf, **kw)

	def write_matrix(self, nf, leavesOnly=False, header=True, rownames=False, distMethod='distance', **kw):
		"""Write the Node in simple format in $1 file.
	
		Keyword argument to determine whether 'append' or 'write'
		(overwrites) mode is turned on: [mode=string]."""
		tree2.write_to_file(self.matrix(leavesOnly=leavesOnly, header=header, rownames=rownames, distMethod=distMethod), nf, **kw)
	
#####################################
############## screen output methods:

#### Newick format-based methods

	def commentAsString(self, comment='comment'):
		"""returns a string representation of a node property to be writen in comment field '[xxx]' of a newick tree."""
		# dict of alternative representations of attributes with valid class info
		extras = {'gainlosscount':tree2.ReferenceTree, 'location':tree2.ReferenceTree, 'locationcount':tree2.ReferenceTree, \
		'locationgainloss':tree2.ReferenceTree, 'locationgainlosscount':tree2.ReferenceTree, 'translocationcount':tree2.ReferenceTree, 'locationreplacement':tree2.ReferenceTree}
		#~ drepliabbrev = {'primary':'I', 'secondary':'II', 'pAt':'pAt', 'pTi':'pTi', 'plasmid':'p', 'absent':'-', 'unknown':'?'}
		locstates = ['I', 'II', 'pAt', 'pTi', 'p', '?']
		if comment==None:
			return None
		elif callable(comment):
			com = comment(self)
		elif comment in extras:
			if not isinstance(self, extras[comment]):
				raise AttributeError, "unvalide extra attribute %s for class %s"%(comment, self.__class__.__name__)
			com = {}
		elif comment in self.__dict__:
			com = self.__dict__[comment]
		else:
			for attr in self.__dict__:
				if comment in attr:		# string is included in attribute name, ex: 'comment' in '_Node__comment'
					com = self.__dict__[attr]
					break
			else:
				raise AttributeError, "unvalide attribute %s for class %s"%(comment, self.__class__.__name__)
		if com==None:
			return None
		if isinstance(com, str):
			comm = com
		elif isinstance(com, bool):
			if com: comm = '+'
			else : comm = '-'
		elif isinstance(com, list):
			try:
				comm='/'.join(com)
			except TypeError:
				comm='/'.join(str(com)[1:-1].split(', '))
		elif isinstance(com, dict):
			#~ if isinstance(self, tree2.ReferenceTree) and comment=='gainlosscount':
			if comment=='gainlosscount':
				e = self.getEvents()
				comm="+%d-%d=%d"%(e['gain'], e['loss'],  self.homolog_count())
			elif 'location' in comment:
				d = self.misc()
				comm=" "
				if comment=='location':
					for key in d:
						if isinstance(key, str):
							comm += "%s "%key
				elif comment=='locationcount':
					for key in locstates:
						comm +="%s=%d "%(key, d.get(key, 0))
				elif comment=='translocationcount':
					for key0 in locstates:
						if key0=='?': continue
						for key1 in locstates:
							if key0=='?': continue
							if (key0, key1) in d:
								comm += "%s->%s=%d "%(key0, key1, d[(key0, key1)])
				elif comment=='locationgainloss':
					for key in locstates:
						comm +="%s=+%d-%d "%(key, d.get(('-', key), 0), d.get((key, '-'), 0))
				elif comment=='locationgainlosscount':
					for key in locstates:
						comm +="%s=+%d-%d=%d "%(key, d.get(('-', key), 0), d.get((key, '-'), 0), d.get(key, 0))
				elif comment=='locationreplacement':
					for key in locstates:
						comm +="%s=%d "%(key, d.get(('replacement', key), 0))
			else:
				if isinstance(self, tree2.ReferenceTree) and comment=='events':
					llkeys=[['gain','loss'],['duplication','receptor'],['donor','replacement']]	# yield [G/L-D/r-d/R] pattern
				else:
					llkeys=[com.keys()]
				lc = []
				for lkeys in llkeys:
					l=[str(com[key]) for key in lkeys]
					lc.append('/'.join(l))
				comm='-'.join(lc)
		elif isinstance(com, tuple):
			if isinstance(self, tree2.GeneTree) and comment=='event':
				lt = [com[0], com[1][0]]
				lt.append( ','.join(com[1][1]))
				if len(lt) > 2:
					lt.append(com[1][2])
					lt.append( ','.join(com[1][3]))
					lt.append(com[1][4])
				comm='/'.join(lt)
		else:
			comm=str(com)
		return comm
	
	def newick(self, **kw):
		"""Returns (extented) Newick string of the sub-tree defined by the Node.

		Follows the (extended) Newick/New Hampshire (NHX) standard for coding trees. Supports bracketed comments, filled with tree.Node attribute as specified by 'comment'.
		Newick format example:
		   '((A:0.1,B:0.2,C:0.1)ABCnode:0.2,(D:0.4,E:0.1)99:0.1);'
		"""
		comment = kw.get('comment','comment')
		ignoreBS = kw.get('ignoreBS',False)
		branch_lengths = kw.get('branch_lengths',True)
		unrooted = kw.get('unrooted', getattr(self, 'unrooted', False))
		ignore_trivial_nodes = kw.get('ignore_trivial_nodes',False)
		nullBranchesAsZeroLength = kw.get('nullBranchesAsZeroLength',False)
		
		comm=self.commentAsString(comment=comment)
		s=''
		if not self.__father: end=';'
		else: end=''
		children = self.__children
		n = len(children)-1		# last index in children
		i = 0					# counter of children
		if ignore_trivial_nodes and n==0:
			# to avoid printing something like ((a:0.321)95:0.15)92 which is not readable by some representation programs (including seaview),
			# ignore the node to print (a:0.471)95
			ccopy = copy.deepcopy(children[n])
			if branch_lengths and ccopy.lg()!=None and self.lg()!=None:
				# sum the branch length
				ccopy.set_lg(ccopy.lg()+self.lg())
			if not ignoreBS:
				ccopy.set_bs(max(ccopy.bs(),self.bs()))
			# directly write the child
			s+=ccopy.newick(**kw) # comment=comment, ignoreBS=ignoreBS, branch_lengths=branch_lengths)
		else:
			if (unrooted and self==self.go_root() and self.nb_children() < 3):
				# for unrooted rendering in newick output, join (at least) three clades by commas : (A,B,C) and do not write any branch attribute at the root
				j = 0 					# counter of picked additional clades
				outlen = 0				# stack for length of clade to be split before redistribution
				def enough(i, j):
					return ((i+j >= 3) and (i > n))
				if children[0].nb_children() < children[1].nb_children():
					# put the largest child in first place, so that lone clade is writen last and length can be added
					children.reverse()
				s+='('
				while not enough(i, j):
					grandchildren = children[i].get_children()
					if (j <= 1) and grandchildren:
						# still need to pick grandchildren
						k = 0
						while (k < len(grandchildren)):
	#						print i, j
							# add the grandchildren
							s+=grandchildren[k].newick(**kw) # comment=comment, ignoreBS=ignoreBS, branch_lengths=branch_lengths)
							k += 1
							j += 1
							if k < len(grandchildren): s+=','
						outlen += children[i].lg()
						i += 1
						if not enough(i, j): s+=','
	#					print s
					else:
	#					print i, j
						# temporarily add outgroup branch length
						children[i].set_lg(children[i].lg() + outlen)
						s+=children[i].newick(**kw) # comment=comment, ignoreBS=ignoreBS, branch_lengths=branch_lengths)
						children[i].set_lg(children[i].lg() - outlen)
						outlen = 0
						i += 1
						if not enough(i, j): s+=','
	#					print s
				s+=')'
				
			else:
				# normal recursive writing of clades
				while i <= n:
					if i == 0: s+='('
					s+=children[i].newick(**kw) # comment=comment, ignoreBS=ignoreBS, branch_lengths=branch_lengths)
					if i < n: s+=','
					else: s+=')'
					i += 1
			if not unrooted:
				# write branch attributes
				if (not self.is_leaf()) and (self.__boot!=None and ignoreBS==False):
					s+=str(self.__boot)
					if comm!=None: s+='[%s]'%comm
				elif self.__lab!=None:
					s+=self.__lab
					if comm!=None: s+='[%s]'%comm
				if branch_lengths:
					if self.__l!=None:
						s+=':'+str(self.__l)
					elif nullBranchesAsZeroLength:
						s+=':0'

		s+=end
		return s

	def phylip(self):
		"""Return output of sequences of the leaves defined by the
		Node in Phylip interleaved format."""
		
		s=''
		l=self.get_leaf_sequences()
		s+=str(len(l))+'\t'+str(len(l[0]))+'\n'
		a=[]
		for i in l:
			this=''
			name=string.join(i.name().split(),'_')
			length=len(name)
			if length>10:
				this+=i.name()[:10]
			else:
				this+=i.name()+' '*(10-length)
			seq=i.seq()
			n=len(seq)/10
			for j in range(n):
				this+=seq[j*10:j*10+10]+' '
			this+=seq[n*10:]
			a+=[this]
		for i in a:
			s+=i[0:31]+'\n'
		m=len(a[0][31:])/21 # all lines are of the same length
		j=11
		while j<=m*21:
			s+='\n'
			for i in a:
				s+='		  '+i[j+21:j+21*2]+'\n' # exactly 10 white spaces needed...
			j+=22
		s+='\n'
		for i in a:
			s+='		  '+i[j+21:]+'\n'
		return s

	def matrix(self, leavesOnly=False, header=True, rownames=False, distMethod='distance'):
		"""Return additive distances matrix of the leaves."""
	
		d=self._matrix(leavesOnly=leavesOnly, distMethod=distMethod)
		n=d.keys()
		n.sort() #order lines and columns
		s=''
		if header:
			if rownames: s += '\t'
			s += '\t'.join(n) + '\n'
		for i in n:
			if rownames:
				s += i+'\t'
			for j in n:
				s+=str(d[i][j])+'\t'
			s+='\n'
		s=string.strip(s)
		return s
		
	def dichotomic_list(self, sort=None):
		"""returns an ordrer list of node labels in the tree, with father node loacted between its two children.
		
		by default order the nodes by decrerasing size; specifying 0 or 1 makes them appear as read in the tree.
		"""
		l = []
		children = self.__children
		if children :
			if len(children)==2:
				if not sort:
					if children[0].nb_leaves() <= children[1].nb_leaves():
						sort = 1
					else:
						sort = 0				
				l += children[0-sort].dichotomic_list(sort=sort)
				l.append(self.__lab)
				l += children[1-sort].dichotomic_list(sort=sort)
		else:
			l.append(self.__lab)
		return l

############################
### graphical output methods
		
	def arborescence_ASCII(self, depth=0, uncleBranches=None, lastSon=False, out="print"):
		"""prints a multi-line ASCII drawing of the tree, in a hirearchical arborescence form.
		
		If 'out' argument is "print" (default), the ouput is written to STDOUT ; else, the output is APPENDED to file named 'out'.
		"""
		if out=="print":
			from sys import stdout
			Out = stdout
		else :
			Out = open(out, 'a')
		if not uncleBranches: uncleBranches = []
		branchlength = self.go_root().max_leaf_depth()
		lenpat = max(uncleBranches+[0])
		#~ print uncleBranches, lenpat, depth, branchlength, lastSon
		pattern = ["    "]*(lenpat+2)
		for branch in uncleBranches:
			pattern[branch+1] = "|   "
		spattern = "".join(pattern)+'\n'
		if self.bs()!=None and not self.is_root():
			spattern = spattern[:-3]+"%.2g\n"%self.bs()
		Out.write( spattern )
		if self.is_leaf():
			Out.write( "".join(pattern[:-1])+"+----"+"----"*(branchlength-lenpat)+" "+self.label()+'\n' )
			if not out=="print":
				Out.close()
		else:
			Out.write( "".join(pattern[:-1])+"+"+"---+"*int(not self.is_root())+" - -"*(branchlength-lenpat)*int(bool(self.label()))+" "+self.label()+'\n' )
			if out!="print":
				Out.close()
			children = self.get_children()
			if lastSon:
				del uncleBranches[-1]
			for child in children[:-1]:
				child.arborescence_ASCII(depth=depth+1, uncleBranches=uncleBranches+[depth], out=out)
			children[-1].arborescence_ASCII(depth=depth+1, uncleBranches=uncleBranches+[depth], lastSon=True, out=out)
	
	def writeSvgTree(self, nfout, mode='w', **kw):
		"""make use of tree2.svgNode module"""
		fout = open(nfout, mode)
		fout.write(svgNode.svgTree(self, **kw))
		fout.close()
	
####################################
### External program call functions

	homedir = os.environ['HOME']
	def seaview(self, seaview_exec_path='seaview', tmpfiledir='%s/tmp'%homedir, filename='tmptree2seaview.phb', bg=True, **kw): #, ignoreBS=False, branch_lengths=True, unrooted=False, comment='comment'):
		ntmpfile = "%s/%s_%s"%(tmpfiledir, filename, hex(hash(self)))	# adds a unique (object-specific) string to the file name
		self.write_newick(ntmpfile, mode='write', ignore_trivial_nodes=True, **kw) #, ignoreBS=ignoreBS, branch_lengths=branch_lengths, unrooted=unrooted, comment=comment)
		callline = "%s %s"%(seaview_exec_path, ntmpfile)
		if bg: callline += ' &'
		print callline
		subprocess.call(callline, shell=True)
		#~ os.remove(ntmpfile)

	def phyml(self, alnfilepath, outtreefilepath, PhyMLpath='/usr/bin/phyml', tmpfiledir='%s/tmp'%homedir, quiet=True, **kw):
		"""perform PhyML command with the node as an user input tree
		
		'alnfile' must be the path of an input sequence alignment in Phylip format 
		(default interleaved format, '-q' option for sequential format).
		
		PhyML options are provided either as a command line argument string or through **kw arguments. 
		e.g.: for '-c 8' command line option, add 'c=8' keyword argument.
		'-u' and '-i' options are ignored as they are defined by node, alnfilepath
		default options make the program to recompute branch lengths and SH-like supports 
		from the input topology (Node object) and nucleic alignment under a GTR+G8+I model.
		"""
		options = ""
		if 'options' in kw:
			options = kw['options']
			if ('-i' in options) or ('-u' in options):
				raise ValueError, "Option string must not contain '-i' or '-u'. Input alignment must be provided through 'alnfilepath' argument, input tree is the Node object"
		else:
			for option in kw:
				if option in ['d','q','n','p','b','m','f','t','v','c','a','s','o']:
					if kw[option] == '': options += " -%s"%option
					else: options += " -%s %s"%(option, kw[option])
		if not options: options="-d nt -q -o lr -m GTR -c 8 -a e -v e -b -3"
		if quiet: options += " --quiet"
		cwd = os.getcwd()
		tmpintreefile = "%s/phyml_tmpintree.phb"%tmpfiledir
		self.write_newick(tmpintreefile, comment=None, mode='write')
		alnfile = alnfilepath.split('/')[-1]
		alntmp = "%s/%s"%(tmpfiledir, alnfile)
		shutil.copy(alnfilepath, alntmp)
		command = PhyMLpath + " -u %s -i %s %s"%(tmpintreefile, alntmp, options)
		larg = command.split()
		subprocess.check_call(larg)	# raises CallProcessError if RAxML does not return with status 0.
		# PhyML call went well, erase the temporary input tree and load the output tree
		os.remove(tmpintreefile)
		phymltreefilepath = "%s/%s_phyml_tree.txt"%(cwd, alnfile)
		shutil.move(phymltreefilepath, outtreefilepath)
		
	def recomputeBranchLengthsAndSupports(self, alnfilepath, outtreefilepath, o='lr', b='-3', **kw):
		"""default options make the program to recompute branch lengths and SH-like supports."""
		self.phyml(alnfilepath=alnfilepath, outtreefilepath=outtreefilepath, o=o, b=b, **kw)
		optgenetree = tree2.GeneTree(fic=outtreefilepath)
		# updates the branch lengths/supports from RAxML/PhyML output to the current GeneTree object
		for optnode in optgenetree:
			lleaves = optnode.get_leaf_labels()
			node = self.map_to_node(lleaves)
			node.set_lg(optnode.lg())
			node.set_bs(optnode.bs())
		
	
#####################################
########### distance measure methods:

	def depth(self):
		"""Return the depth of this Node as its 'distance' (in number of branches) to the root.
	
		WARNING: Branch lengths are not used.
		"""
	
		d=0
		if self.__father:
			d=self.__father.depth()+1
		return d
		
	def max_leaf_depth(self):
		"""returns the maximum of all leaf depths (from the tree root) under the Node."""
		mld = 0
		for leaf in self.get_leaves():
			if leaf.depth() > mld :
				mld = leaf.depth()
		return mld

	def distance_root(self, nullBranchesAsZeroLength=False):
		"""Return the distance (sum of branch lengths) from this Node to the root."""
		if self.__l==None:
			if nullBranchesAsZeroLength:
				l = float(0)
			else:
				raise ValueError, "No branch length at node %s"%(repr(self))
		else:
			l = self.__l
		if self.__father:
			l += self.__father.distance_root(nullBranchesAsZeroLength=nullBranchesAsZeroLength)
		return l
			
	def node_distance(self,n):
		"""Return the node distance (number of nodes) on the path separating the Node from node $1."""
		if self.go_root() is not n.go_root():
			raise IndexError, "nodes are not parts of the same tree"	
		elif self is n:
			return 0
		else:
			d1 = self.depth()
			d2 = n.depth()
			selfch = self.get_all_children()	
			if n in selfch:
				return (d2-d1)
			elif self in n.get_all_children():
				return (d1-d2)
			else:
				d = 1
				f = self.go_father()
				fch = f.get_all_children_cache({self:selfch})
				while (n not in fch):
					d += 1
					f = f.go_father()
					fch = f.get_all_children_cache({f:fch})
				d3 = f.depth()
				return (d+(d2-d3))
		
	def distance(self,n):
		"""Return the distance (sum of branch lengths) on the path separating the Node from node $1."""
		if self.go_root() is not n.go_root():
			raise IndexError, "nodes are not parts of the same tree"		
		elif self is n:
			return float(0)
		else:
			d1 = self.distance_root(nullBranchesAsZeroLength=True)
			d2 = n.distance_root(nullBranchesAsZeroLength=True)
			selfch = self.get_all_children()
			if n in selfch:
				return (d2-d1)
			elif self in n.get_all_children():
				return (d1-d2)
			else:
				d = self.__l
				f = self.go_father()
				fch = f.get_all_children_cache({self:selfch})
				while (n not in fch):
					d += f.lg()
					f = f.go_father()
					fch = f.get_all_children_cache({f:fch})
				d3 = f.distance_root(nullBranchesAsZeroLength=True)
				return (d+(d2-d3))
			
	def leaf_distances(self, excludedLeaves=None):
		"""Return the list of distances from the node to the leaves below it.
		
		Ignore leaves in excludedLeaves for distance calculation.
		"""
		if not excludedLeaves: excludedLeaves = []
		ld=[]
		for leaf in self.get_leaves():
			if (leaf.label() in excludedLeaves) or (leaf in excludedLeaves):
				if self is leaf:
					ld.append(0)
				else:
					continue
			else:
				ld.append(self.distance(leaf))
		return ld
			
	def mean_leaf_distance(self, excludedLeaves=None):
		"""Return the mean distance from the node to the leaves below it in a R-friendly manner.
		
		Ignore leaves in excludedLeaves for distance calculation.
		"""
		ld = self.leaf_distances(excludedLeaves=excludedLeaves)
		sld = sumlen(*ld)
		if sld is None:
			return 'NA'
		else:
			return sld/len(ld)
	
	def max_leaf_distance(self, excludedLeaves=None):
		"""Return the maximum distance from the node to the leaves below it in a R-friendly manner.
		
		Ignore leaves in excludedLeaves for distance calculation.
		"""
		ld = self.leaf_distances(excludedLeaves=excludedLeaves)
		mld = max(ld)
		if mld is None:
			return 'NA'
		else:
			return mld
			
	def lineage_age(self, excludedLeaves=None):
		"""Similar to Node.mean_leaf_distance(), except that it adds half of its own branch length to give an approximate age of the lineage"""
		age = self.mean_leaf_distance( excludedLeaves=excludedLeaves)
		if self.__l and not age=='NA':
			age += (self.__l)/2
		return age
		
	def lineage_distance(self, n):
		"""Similar to Node.distance(), except that it substracts half of the nodes' branch lengths to give an approximate distance on the lineage path"""
		d = self.distance(n)
		if d:
			for node in [self, n]:
				if node.lg():
					d -= (node.lg())/2
		return d
			
	def treelength(self, excludeSelf=False):
		"""return the total length of the tree under the node"""
		lg = 0
		for node in self:
			if node is self and excludeSelf: continue
			if node.lg() is not None:
				lg += node.lg()
		return lg
		
	def treemedian(self, lignored=None):
		if not lignored: lignored = []
		l = []
		for node in self:
			if (not node in lignored) and (node.lg() is not None):
				l.append(node.lg())
		l.sort()
		if len(l)%2==1:
			return l[(len(l)+1)/2]
		else:
			return (l[len(l)/2] + l[(len(l)/2)+1])/2
				
	def _matrix(self, leavesOnly=False, distMethod='distance'):
		"""Return the distance matrix of all leaves, starting at the root."""
		self=self.go_root()
		if leavesOnly:
			l=self.get_leaves()
		else:
			l=self.get_all_children()
		d={}
		for i in range(len(l)):
			d[l[i].label()]={}
			for j in range(len(l)):
				if distMethod=='distance':
					d[l[i].label()][l[j].label()]=l[i].distance(l[j])
				elif distMethod=='lineage_distance':
					d[l[i].label()][l[j].label()]=l[i].lineage_distance(l[j])
				else:
					raise ValueError, "bad distance method descriptor"
		return d	

#####################################
############## miscellaneous methods:
	def get_parents(self):
		lparents = []
		f = self.go_father()
		while f:
			lparents.append(f)
			f = f.go_father()
		return lparents				
		
	def get_parent_labels(self):
		"""Return the list of labels of all nodes above the Node, excluding itself, ordered by increasing depth (i.e. decreasing node distance from root)."""
		lparents = []
		f = self.go_father()
		while f:
			lparents.append(f.label())
			f = f.go_father()
		return lparents	

	def get_dicCladeToLeaves(self, prefix):
		"""Returns a dictionary of all leaves below each node under Node"""
		self.complete_internal_labels(prefix)
		dic = {}
		children = self.get_sorted_children()
		for c in children:
			dic[c.label()] = c.get_leaf_labels()
		return dic
	
	def __imul__(self, x):
		""" Recursively multiplies the length of the branches that are under this node by factor x (>0)."""
		if x>=0:
			for i in self.__children:
				i.__l*=x
				i*=x
		return self
		
	def __idiv__(self, x):
		""" Recursively divides the length of the branches that are under this node by factor x (>0)."""
		if x>0:
			for i in self.__children:
				i.__l/=x
				i/=x
		return self
		
	def __iadd__(self, x):
		"""adds the given length to the node's branch, i.e. the branch ABOVE the node, unlike for __idiv__ and __imul__."""
		if self.__l is None: self.__l = 0
		self.__l += x
		return self
	
		
########################################
### Topology / taxonomic comparison methods:
	
	def is_bifurcated(self):
		if len(self.__children) == 2:
			return True
		else:
			return False

	def is_leaf(self):
		"""returns bollean stating if Node is a terminal node/leaf"""
		return (not bool(self.__children))
			
	def is_cherry(self):
		"""returns bollean stating if Node has only leaves below"""
		for c in self.__children:
			if not c.is_leaf():
				return False
		return True
			
	def is_root(self):
		if self.__father: return False
		else: return True
		
	def is_subroot(self):
		if not self.__father:
			return False
		else:
			if self.__father.is_root():
				return True
			
	def is_child(self, node1):
		"""returns bollean stating if the Node is below node1 in the tree"""
		#~ root = self.go_root()
		#~ f = self.__father
		#~ while f:
			#~ if f is node1:
				#~ return True
			#~ elif f is root:
				#~ return False
			#~ else:
				#~ f=f.go_father()	
		f = self.__father
		while f:
			if f is node1:
				return True
			else:
				f=f.go_father()
		else:
			return False	
				
	def is_parent_of_any(self, lnodes, returnList=False):
		lchild = []
		for node in lnodes:
			if node.is_child(self):
				if returnList: lchild.append(node)
				else: return True
		else:
			if returnList: return lchild
			else: return False
				
	def is_childorself(self, node1):
		return ((self is node1) or self.is_child(node1))
				
	def is_brother(self, node1):
		"""returns bollean stating if the Node is neighbor of node1 in the tree"""
		try:
			bro = self.go_brother()
		except ValueError:
			raise IndexError, "Node %s has no brother"%self.__lab
		if node1 is bro:
			return True
		else:
			return False
	
	def is_NNI(self, cautiousrec, cautiousdon, lleaves):
		"""Tests if a transfer with given cautious donor and receptor clades (objects or labels) in reference tree (self) is a nearest-neighbor interchange (NNI).
		
		'cautious' donors and receptors are nodes returned by function 'map_to_ancient_node'.
		reference tree is pruned to display only leaves in provided in 'lleaves'.
		"""
		reftree = self
		if isinstance(cautiousdon, Node) and isinstance(cautiousrec, Node):
			caudon = reftree[cautiousdon.label()]
			caurec = reftree[cautiousrec.label()]
		elif isinstance(cautiousdon, str) and isinstance(cautiousrec, str):
			caudon = reftree[cautiousdon]
			caurec = reftree[cautiousrec]
		else:
			raise ValueError, "'cautiousrec' and 'cautiousdon' must be either both tree.Node objects or both strings (labels of tree.Node objects)"
		if caurec.depth() >= 2 and caudon is caurec.go_uncle():
			return True
		elif caudon.depth() >= 2 and caurec is caudon.go_uncle():
			return True
		else:
			return False		
		
	def map_to_node(self, lleaves, force=False, useSpeDict=False, **kw):
		"""finds deepest node/clade in tree that includes all leaves in input label list = MRCA = most recent common ancestor of the leaf set
		
		uses a desending algorithm
		"""
		def getSpeSet(node, useSpeDict):
			if not useSpeDict: s = set(node.get_leaf_labels())
			else: s = set(node.listSpecies(**kw))
			return s
			
		ll = getSpeSet(self, useSpeDict)
		si = set(lleaves)
		if not si <= ll:
			if force:
				si = si & ll
			else:
				raise IndexError, "Input leaf set is larger than leaves present in reference tree:\nsi: %r\nll: %r"%(si, ll)
		if not si:
			return None
		clade = self
		scl = ll
		#~ print "si", si
		#~ print "scl", scl
		while si <= scl:
			children = clade.get_children()
			for child in children:
				scl = getSpeSet(child, useSpeDict)
				if si <= scl:
					#~ print "scl", scl
					clade = child
					break
			else:
				break
		return clade
		
	def mrca(self, lleaves, force=False, useSpeDict=False):
		"""function alias to Node.map_to_node"""
		return self.map_to_node(lleaves, force=force, useSpeDict=useSpeDict)
		
	def map_to_ancient_node(self, lleaves, excludedLeaves=[], useSpeDict=False):
		"""finds most ancient node/clade in tree that includes all leaves in $1 input label list 'lleaves' but is not a parent of any leaf present in $2 input label list 'excludedLeaves'."""
		clade = self.map_to_node(lleaves, useSpeDict=useSpeDict)
		se = set(excludedLeaves)
		f = clade.go_father()
		while f and not (se & set(f.get_leaf_labels())):
			clade = f
			f = f.go_father()
		return clade		
	
	def coalesce(self, lnodes, useSpeDict=False):
		"""finds deepest (closest to leaves) node/clade in tree that includes all nodes in input label/node list = MRCA = most recent common ancestor of the node set"""
		lleaves = []
		for node in lnodes:
			try:
				lleaves += node.get_leaf_labels()
			except AttributeError:
				lleaves += self[node].get_leaf_labels()
		return self.map_to_node(lleaves, useSpeDict=useSpeDict)
		
	def coalesce_ancient(self, lnodes, excludedLeaves=[]):
		"""finds shalowest (closest to root) node/clade in tree that includes all nodes in input label/node list = MACA = most ancient common ancestor of the node set but is not a parent of any leaf present in $2 input label list 'excludedLeaves'."""
		lleaves = []
		for node in lnodes:
			try:
				lleaves += node.get_leaf_labels()
			except AttributeError:
				lleaves += self[node].get_leaf_labels()
		return self.map_to_ancient_node(lleaves, excludedLeaves=excludedLeaves)
		
	def map_to_nr_nodes(self, reftree, seed=None, lnodes=None, useSpeDict=True, returnLabels=False, silent=True):
		"""searches in a gene tree (self) the non-redundant nodes corresponding to nodes in a reference tree
		
		set of reference tree nodes to match can be restricted providing a node label list (lnodes;
		target clade can be narrowed by precsing a seed leaf close to which it is sought.
		"""
		if not silent: print "seed", seed
		d = {}
		if not lnodes:l = reftree.get_children_labels()
		else:l = lnodes
		for node in l:
			refleaves = reftree[node].get_leaf_labels()
			if not silent: print "refleaves", refleaves
			clade = self.map_to_node(refleaves, force=True, useSpeDict=useSpeDict)
			# clade=None if the interscection of refleaves and gene tree species set is null
			# if so, attempts to find the closest node in reference tree mapable to a node in gene tree
			father = reftree[node]
			while (not clade) and father:
				refleaves = father.get_leaf_labels()
				# searches the intersection of 'refleaves'
				clade = self.map_to_node(refleaves, force=True, useSpeDict=useSpeDict)
				# enlarges leaf set to leaves of brother nodes
				father = father.go_father()
			if not silent: print "refleaves", refleaves
			if clade:
				if seed and not clade.is_leaf():
					# search if one can narrow the matching clade to a subtree above seed
					if not seed in self.get_leaf_labels(): raise IndexError, "seed %s must be the label of a leave in tree\n%"%(seed, self.newick())
					nextclade = clade
					nextcladespe = cladespe = set(clade.listSpecies()) & set(refleaves)
					st = self.get_target_child(label=seed)
					while st:
						if not silent: print "clade", clade.label()
						nextclade = st.map_to_node(refleaves, force=True, useSpeDict=useSpeDict)
						if not nextclade: break
						nextcladespe = set(nextclade.listSpecies()) & set(refleaves)
						if nextclade.is_child(clade):
							# shows amelioriation (narrowing of the mapped clade close to the target), go straight ot it
							pivot = st.map_to_node(nextclade.get_leaf_labels()+[seed])
							if not silent: print "pivot", pivot.label()
							if not (nextcladespe < cladespe):
								st = pivot
								clade = nextclade
							else:
								### envisager clause alternative ou nextcladespe puisse etre plus reduit que cladespe
								### si cela permet de vraiment rapprocher le clade de seed: utiliser la distance au pivot? le fait de passer par des duplications?
								### accepter de perdre les espèces qui seraient absente dans le clade correspondant de la famille de la prot de ref ?						
								break
						else:
							st = st.get_target_child(label=seed)			
				
				if not returnLabels: d.setdefault(clade, []).append(node)
				else: d.setdefault(clade.label(), []).append(node)
				
		return d
		
	def map_from_collapsed_node(self, leavesunder, representedleaves):
		"""find in a tree the one or several nodes matching a node of a pruned tree (same tree with pruned branches)
		
		the node of the partial tree is described by the list of leaves under it (leavesunder),
		and the the partial tree is described by the list of leaves represented in it (representedleaves).
		"""
		# refernce node is the MRCA of leaves in full tree
		refnode = self.map_to_node(leavesunder)
		srl = set(representedleaves)
		srn = set(self.get_leaf_labels())
		rl = refnode.label()
		lredundant = [rl]
		if slu == srn:
			# full tree and partial tree are identical
			return lredundant
		elif slu <= srn:
			# some leaves from the full tree are missing in the partial tree
			prunedleaves = srn - slu
			# check the pattern of node collpasing if one would prune those leaves
			dcollapsed = self.dictCollapsed(prunedleaves) 		# d[staying] = removed
		# adds potential collased nodes leading to the reference node
		while rl in dcollapsed:
			rl = dcollapsed[rl]
			lredundant.append(rl)
		return lredundant
		
	def map_collapsed_nodes(self, partialtree):
		"""make correspondancy between nodes in a tree to the one or several nodes matching it in a pruned tree (same tree with pruned branches)
		
		the node of the partial tree is described by the list of leaves under it (leavesunder),
		and the the partial tree is described by the list of leaves represented in it (representedleaves).
		"""
		dnodecollapsed = {}
		slu = set(partialtree.get_leaf_labels())
		srn = set(self.get_leaf_labels())
		#~ if slu == srn:
			#~ # full tree and partial tree are identical
			#~ for node in partialtree:
				#~ dnodecollapsed[node] = [node.label()]
			#~ return dnodecollapsed
		#~ el
		if slu <= srn:
			# some leaves from the full tree are missing in the partial tree
			prunedleaves = srn - slu
			# check the pattern of node collpasing if one would prune those leaves
			dcollapsed = self.dictCollapsed(prunedleaves) 		# d[staying] = removed
			for node in partialtree:
				leavesunder = node.get_leaf_labels()
				# refernce node is the MRCA of leaves in full tree
				refnode = self.map_to_node(leavesunder)
				rl = refnode.label()
				lredundant = [rl]
				# adds potential collased nodes leading to the reference node
				while rl in dcollapsed:
#					print 'rl', rl
					rl = dcollapsed[rl]
					lredundant.append(rl)
				dnodecollapsed[node] = lredundant
		else:
			raise ValueError, "leaf set of partialtree is not included within that of the node"
		return dnodecollapsed		
		
	def is_monophyletic(self, lleaves, ltrans=[], useSpeDict=False):
		"""tests is a clade defined by a set of leaf labels is monophyletic in Node, excluding transfered leaves from comparison"""
		for leaf in ltrans:
			if leaf in lleaves:
				lleaves.remove(leaf)
		ancestor = self.map_to_node(lleaves, useSpeDict=useSpeDict)
		if not useSpeDict:
			clade = ancestor.get_leaf_labels()
		else:
			clade = ancestor.dictLeafLabelsToSpecies().values()
		for leaf in ltrans:
			if leaf in clade:
				clade.remove(leaf)
		if set(clade) == set(lleaves):
			return True
		else:
			return False
			
	def whichClade(self, lspe, returnLabels=True, init=True, holes=0):
		"""identifies the monophyletic groups present from a list of species names, possibly present in several copies
		
		if 'holes' is specified, will allow the recognition of a clade which is not complete, i.e.
		missing 'holes' leaves (if 'holes' is an int) or
		missing a fraction of 'holes' leaves (if 'holes' is a float).
		"""
		lcla = []
		def get_spe(core, sspe, i):
			core.append(sspe.pop(i))
			
		if lspe:
			if init: sspe = self.sort(lspe, order=0, labels=True)
			else: sspe = lspe
			core = [sspe.pop(0)]
			i = 0
			while i < len(sspe):
				#~ print "sspe", sspe
				if sspe[i] in core:
					i += 1
					continue
				cat = core+sspe[i:i+1]
				if self.is_monophyletic(cat):
					get_spe(core, sspe, i)
				else:
					mrca = self.map_to_node(cat)
					allspe = mrca.get_leaf_labels()
					otherspe = set(allspe) - set(core)
					if otherspe:
						if otherspe <= set(sspe):
							for ospe in otherspe:
								o = sspe.index(ospe)
								get_spe(core, sspe, o)
							i = 0
							continue
						elif holes:
							missing = otherspe - set(sspe)
							if isinstance(holes, float): h = max(len(allspe)*holes, 1)
							else: h = holes
							#~ print "len(missing)", len(missing),"h", h
							if (len(missing) <= h):
								for ospe in (otherspe - missing):
									o = sspe.index(ospe)
									get_spe(core, sspe, o)
								i = 0
								continue
							#~ print "missing in otherspe", missing, "for matching", mrca.label()
					i += 1
			cla = self.map_to_node(core)
			if returnLabels: lcla.append(cla.label())
			else: lcla.append(cla)
			lcla += self.whichClade(sspe, returnLabels=returnLabels, init=False, holes=holes)
		return lcla
		
	def youngest(self, llabels):
		y = llabels[0]
		for lab in llabels[1:]:
			if self[lab].depth > self[y].depth():
				y = lab
		return self[y]
		
	def lowBSleaves(self, minbs, useSpeDict=False):
		lowbsleaves = []
		dleafspe = self.dictLeavesToSpecies()
		for leaf in self.get_leaves():
			lbs = leaf.bootstraps_on_path_to(self)
			for bs in lbs:
				if bs >= minbs:
					break
				else:
					# no branch support is high enough on the path from the leaf to the root, do not consider this leaf
					if useSpeDict: lowbsleaves.append(dleafspe[leaf])
					else: lowbsleaves.append(leaf.label())
		return lowbsleaves
		
		
	def hasSameRoot(self, other, useSpeDict=False):	#, ignoreLowBS=None):
		"""compare roots on the basis of distribution of leaves in subroot nodes 
		
		with useSpeDict=False, both 'self' and 'other' trees must share the same leaf label set; 
		with useSpeDict=True, they only need to have matching species prefix for their leaf label, 
		  and for this reason both trees must be unicopy (not more than one leaf per species).
		with ignoreLowBS=float, leaves under nodes with low support in 'other'
		  and their counterpart in 'self' will not be considered.
		"""
		if not isinstance(other, Node):
			raise TypeError, "expected a tree.Node instance"
		subRoots1 = self.get_children()
		subRoots2 = other.get_children()
		sameSR = []
#		ignoreleaves = other.lowBSleaves(ignoreLowBS)
		for sr1 in subRoots1:
			for sr2 in subRoots2:
				if useSpeDict: 
					s1 = set(sr1.listSpecies()) #- set(ignoreleaves)
					s2 = set(sr2.listSpecies()) #- set(ignoreleaves)
				else:
					s1 = set(sr1.get_leaf_labels()) #- set(ignoreleaves)
					s2 = set(sr2.get_leaf_labels()) #- set(ignoreleaves)
				if s1 == s2:
					sameSR.append(sr2)
					break
		if set(subRoots2)==set(sameSR):
			return True
		else:
			return False
			
	def hasSameTopology(self, other, checkInternalLabels=False):
		"""recursively test if al nodes in $1 (other) are monophyletic in self (reference tree)"""
		sleaflab = set(self.get_leaf_labels())
		oleaflab = set(other.get_leaf_labels())
		if sleaflab==oleaflab:
			if checkInternalLabels:
				sintlab = set(self.get_children_labels()) - sleaflab
				ointlab = set(other.get_children_labels()) - oleaflab
				if sintlab!=ointlab: return False
			# recursion
			ochildren = other.get_children()
			schildren = self.__children
			scrange = range(len(schildren))
			for ochild in ochildren:
				for j in scrange:
					schild = schildren[j]
					# explore every child to match clades
					if schild.hasSameTopology(ochild, checkInternalLabels=checkInternalLabels):
						# remove value from range, so this combination is not evaluated in further children matching search
						scrange.remove(j)
						break # for schild, escape the return False
				else:
					return False
			else:
				return True
		else:
			return False
				
	def countCommonBiparts(self, other, useSpeDict=False):
		n = 0
		for node in other:
			if useSpeDict: lleaves = node.listSpecies()
			else: lleaves = node.get_leaf_labels()
			if self.is_monophyletic(lleaves):
				n += 1
		return n

	def compareToCollection(self, ncolfile, exactMatch=False, branch_lengths=True):
		"""tracks nodes with topology equivalent to those of subtrees of self in a replicate tree collection
		
		if exactMatch=True, tests the exact concordance of topologies; else, just tests the monophily of the node.
		return a tuple with a dictionary of observed branch length in matching nodes and the total number of replicates.
		"""
		colfile = open(ncolfile, 'r')
		dlab_length = {}
		repcount = 0
		for line in colfile:
			repcount += 1
			replicate = self.newnode(newick=line.rstrip('\n'), branch_lengths=branch_lengths)
			replicate.cureTree(reftree=self, branch_lengths=branch_lengths)
			for node in replicate:
				lleaves = node.get_leaf_labels()
				refnode = self.map_to_node(lleaves)
				if (not exactMatch and refnode.is_monophyletic(lleaves)) or (exactMatch and refnode.hasSameTopology(node)):
					if branch_lengths: ndlg = node.lg()
					else: ndlg = 1
					dlab_length[refnode.label()] = dlab_length.setdefault(refnode.label(), []) + [ndlg]
		colfile.close()
		return (dlab_length, float(repcount))			
							
	def getAverageLengthsFromCollection(self, ncolfile, exactMatch=False):
		"""set the length of nodes in self to the average of those observed in nodes with equivalent topology in a replicate tree collection"""
		(dlab_length, repcount) = self.compareToCollection(ncolfile, exactMatch=exactMatch)
		for node in self:
			if node.label() in dlab_length:
				lengths = dlab_length[node.label()]
				node.set_lg( sum(lengths)/float(len(lengths)) )

	def getBootstrapFromCollection(self, ncolfile, exactMatch=False, branch_lengths=True):
		"""set the bs of nodes in self to the number of observed nodes with equivalent topology in a replicate tree collection"""
		(dlab_length, repcount) = self.compareToCollection(ncolfile, exactMatch=exactMatch, branch_lengths=branch_lengths)
		root = self.go_root()
		for node in self:
			if not node.is_leaf() or node is root:
				if node.label() in dlab_length:
					lengths = dlab_length[node.label()]
					node.set_bs( float(len(lengths))/float(repcount) * 100)
					
	def getBsAndLenFromCollection(self, ncolfile, exactMatch=False):
		"""mixture of getAverageLengths... and getBootstrap..."""
		(dlab_length, repcount) = self.compareToCollection(ncolfile, exactMatch=exactMatch)
		treeroot = self.go_root()
		for node in self:
			if node.label() in dlab_length:
				lengths = dlab_length[node.label()]
				node.set_lg( sum(lengths)/float(len(lengths)) )
				if not node.is_leaf() or node is treeroot:
					node.set_bs( float(len(lengths))/float(repcount) * 100)		
	
#####################################	   
### 'walking along the tree' methods:

	def go_root(self):
		"""Return Root Node."""	
		if self.__father:
			self=self.__father.go_root()
		return self

	def go_father(self, n=1):
		"""Return ancestor Node of this Node. Go n generations above the node, 1 by default (the father)"""
		if n==1:
			return self.__father
		elif n>1:
			if self.__father:
				return self.__father.go_father(n-1)
			else:
				raise IndexError, "Already at root when n=%d. No more ancestors to go to."%n
		else:
			raise ValueError, "please provide ancestor rank n > 0"
		

	def go_child(self,n):
		"""Return $1 (number) child of this Node."""
		if self.__children and n<len(self.__children):
			self=self.__children[n]
		else:
			raise IndexError, 'Already at leaf. No more children to go to.'
		return self
	
	def go_brother(self):
		"""return brother node (other child from biffurcated father node)"""
		fat = self.go_father()
		if fat:
			c = fat.get_children()
#			print fat.get_leaf_labels()
			if len(c)==2:
				if self is c[0]:
					return c[1]
				elif self is c[1]:
					return c[0]
				else:
					raise ValueError, "Node %s not in his father's child list"%self.__lab
			else:
				raise ValueError, "Node %s's father has %s child(ren): go_brother() is not determined"%(self.label(),len(c))
		else:
			raise ValueError, "Node %s has no father"%self.__lab
			
	def go_uncle(self):
		"""return uncle node (other child from biffurcated grandfather node)"""
		fat = self.go_father()		
		if fat:
			uncle = fat.go_brother()
			return uncle	
			
	def lineage(self, value="node"):
		"""return the list of node objects/ids/labels located above the node (excludinfg itself)"""
		line = []
		fat = self.go_father()
		while fat:
			if value=="node":
				line.append(fat)
			elif value=="id":
				line.append(fat.nodeid())
			elif value=="label":
				line.append(fat.label())
			fat = fat.go_father()
		return line	
		
	def path_to(self, n, returnLabels=False):
		"""returns list of nodes on the path from the Node to node $1 (including both)"""
		if self.go_root() != n.go_root():
			print "roots (self, n):", [self.go_root(), n.go_root()], self.go_root().label(), n.go_root().label()
			print "self\n", self
			print "n\n", n
			raise IndexError, "nodes are not parts of the same tree"
		elif self is n:
			if returnLabels: return [self.label()]
			else: return [self]
		else:
			if n in self.get_all_children():
				path = n.path_to(self, returnLabels=returnLabels)
				path.reverse()
				return path
			elif self in n.get_all_children():
				f = self.go_father()
				if returnLabels: return [self.label()] + f.path_to(n, returnLabels=returnLabels)
				else: return [self] + f.path_to(n, returnLabels=returnLabels)
			else:
				f = self
				path = []
				while (n not in f.get_all_children()):
					if returnLabels: path.append(f.label())
					else: path.append(f)	
					f = f.go_father()
				return path + f.path_to(n, returnLabels=returnLabels)
				
	def bootstraps_on_path_to(self, n, excludeLeaves=False, excludeTransfers=False):
		"""returns list of branch bootstraps on the path from the Node to node $1"""
		if excludeTransfers and not isinstance(self, tree2.GeneTree): raise TypeError, "cannot test transfer event type on %s"%(str(type(self)))
		treePath = self.path_to(n)			
#		for node in treePath:
#			print node.label(),
		lbs = []
		for i in range(1, len(treePath)):
			if treePath[i-1].depth() > treePath[i].depth():
				if not ((excludeLeaves & treePath[i-1].is_leaf()) or (excludeTransfers and treePath[i-1].transfer())):
					lbs.append(treePath[i-1].bs())
			elif treePath[i].depth() > treePath[i-1].depth():
				if not ((excludeLeaves & treePath[i].is_leaf()) or (excludeTransfers and treePath[i].transfer())):
					lbs.append(treePath[i].bs())
			else:
				raise IndexError, 'neighbour nodes %s and %s should be at different depth'%(treePath[i-1].label(), treePath[i].label())
		return lbs
		
	def lengths_on_path_to(self, n, excludeLeaves=False):
		"""returns list of branch bootstraps on the path from the Node to node $1"""
		treePath = self.path_to(n)
#		for node in treePath:
#			print node.label(),
		llg = []
		for i in range(1, len(treePath)):
			if treePath[i-1].depth() > treePath[i].depth():
				if not (excludeLeaves & treePath[i-1].is_leaf()):
					llg.append(treePath[i-1].lg())
			elif treePath[i].depth() > treePath[i-1].depth():
				if not (excludeLeaves & treePath[i].is_leaf()):
					llg.append(treePath[i].lg())
			else:
				raise IndexError, 'neighbour nodes %s and %s should be at different depth'%(treePath[i-1].label(), treePath[i].label())
		return llg		
	
	#~ def findboundariesbeneath(self, minbs):
		#~ boundbeneath = []
		#~ boundabove  = []
		#~ if (self.bs() < minbs):
			#~ children = self.get_children()
			#~ for child in children:
				#~ tc = child.findboundariesbeneath(minbs)
				#~ boundbeneath += tc[0]
				#~ boundabove += tc[1]
		#~ else:
			#~ boundbeneath += [self]
		#~ return (boundbeneath, boundabove)
		#~ 
	#~ def findboundariesabove(self, minbs):
		#~ boundbeneath = []
		#~ boundabove  = []
		#~ fat = self.go_father()
		#~ if fat:
			#~ if (fat.bs() < minbs):
				#~ bro = self.go_brother()
				#~ tb = bro.explorebeneath(minbs)
				#~ tf = fat.exploreabove(minbs)
				#~ boundbeneath += tb[0] + tf[0]
				#~ boundabove += tb[1] + tf[1]
			#~ else:
				#~ boundabove += [fat]
		#~ return (boundbeneath, boundabove)
		
	def explorebeneath(self, minbs, findBoundaries=False, lbound=None):
		"""explore non-supported path from starting node and below it (recursive)
		
		return either the list of accessible (non supported) nodes 
		or the list of boundary (supported) nodes + tree leaves (natural boundaries).
		"""
		if not lbound: lbound = []
		nodes = []
		if (self.bs() < minbs) and not (self in lbound):
			children = self.get_children()
			if not findBoundaries or not children: nodes += [self]
			for child in children:
				nodes += child.explorebeneath(minbs, findBoundaries=findBoundaries, lbound=lbound)
		else:
			nodes += [self]
		return nodes
		
	def exploreabove(self, minbs, findBoundaries=False, lbound=None):
		"""explore non-supported path above the starting node, and beneath the found nodes (recursive)
		
		return either the list of accessible (non supported) nodes 
		or the list of boundary (supported) nodes + tree leaves (natural boundaries).
		"""
		nodes = []
		if not lbound: lbound = []
		if not self in lbound:
			fat = self.go_father()
			if fat:
				if fat.bs() < minbs:
					if not findBoundaries: nodes += [fat]
					bro = self.go_brother()
					nodes += bro.explorebeneath(minbs, findBoundaries=findBoundaries, lbound=lbound)
					nodes += fat.exploreabove(minbs, findBoundaries=findBoundaries, lbound=lbound)
				else:
					nodes += [fat]
		else:
			if findBoundaries: nodes += [self]
		return nodes		
		
	#~ def explorebeneath(self, minbs, boundbeneath=[], boundabove=[]):
		#~ """explore non-supported path from starting node and below it (recursive)
		#~ 
		#~ return either the list of accessible (non supported) nodes or the list of boundary (supported) nodes.
		#~ """
		#~ nodes = []
		#~ if (self.bs() < minbs) and not (self in lbound):
			#~ nodes += [self]
			#~ children = self.get_children()
			#~ for child in children:
				#~ nodes += child.explorebeneath(minbs, lbound=lbound)
		#~ else:
			#~ nodes += [self]
		#~ return nodes
		#~ 
	#~ def exploreabove(self, minbs, boundbeneath=[], boundabove=[]):
		#~ """explore non-supported path above the starting node, and beneath the found nodes (recursive)
		#~ 
		#~ return either the list of accessible (non supported) nodes or the list of boundary (supported) nodes.
		#~ """	
		#~ nodes = []
		#~ fat = self.go_father()
		#~ if fat:
			#~ if (fat.bs() < minbs) and not (fat in lbound):
				#~ nodes += [fat]
				#~ bro = self.go_brother()
				#~ nodes += bro.explorebeneath(minbs, lbound=lbound)
				#~ nodes += fat.exploreabove(minbs, lbound=lbound)
			#~ else:
				#~ nodes += [fat]
		#~ return nodes
		
	def explore_nonsupported_paths(self, minbs):
		"""explore all non-supported paths starting from the node, finding all non supported nodes around"""
		nodes = self.explorebeneath(minbs)
		nodes += self.exploreabove(minbs)
		return nodes
		
	#~ def find_support_boundaries(self, minbs):
		#~ """explore all non-supported paths starting from the node, finding the surrounding supported nodes"""
		#~ boundbeneath, boundabove = self.findboundariesbeneath(minbs)
		#~ t = self.findboundariesabove(minbs)
		#~ boundbeneath += t[0]
		#~ boundabove += t[1]
		#~ return (boundbeneath, boundabove)
		
	def find_support_boundaries(self, minbs, returnLabels=False):
		"""explore all non-supported paths starting from the node, finding the surrounding supported nodes"""
		nodes = self.explorebeneath(minbs, findBoundaries=True)
		nodes += self.exploreabove(minbs, findBoundaries=True)
		if not returnLabels: return nodes
		else: return tree2.getlabels(nodes)

	def explore_within_boundaries(self, lbound, maxbs=101, returnLabels=False):
		"""explore all paths starting from the node while not crossing one of the listed boudaries"""
		nodes = []
		if self in lbound:
			# the start point is a boundary itself ; if alone, it serve as both upper and lower boundaries,
			# otherwise, must determine if it constitues an upper or lower boundary, or both
			lup = []
			ldown = []
			for bound in lbound:
				if bound.is_child(self):
					ldown += bound
				elif self.is_child(bound):
					lup += bound
			if ldown:
				# if there are other listed boundaries beneath, this should not be a lower boundary, rather an upper boundary
				nodes += self.explorebeneath(maxbs, lbound=ldown)
				nodes += self.exploreabove(maxbs, lbound=lbound)
			elif lup:
				# if there are other listed boundaries above and not beneath, this should not be an upper boundary, rather an lower boundary
				nodes += self.explorebeneath(maxbs, lbound=lbound)
				nodes += self.exploreabove(maxbs, lbound=lup)
			else:
				# the start point is alone, serve as both upper and lower boundaries
				nodes += self.explorebeneath(maxbs, lbound=lbound)
				nodes += self.exploreabove(maxbs, lbound=lbound)
		else:	
			nodes += self.explorebeneath(maxbs, lbound=lbound)
			nodes += self.exploreabove(maxbs, lbound=lbound)
		if not returnLabels: return nodes
		else: return tree2.getlabels(nodes)
		
	#~ def explore_within_boundaries(self, boundbeneath, boundabove, maxbs=101, returnLabels=False):
		#~ """explore all paths starting from the node while not crossing one of the listed boudaries"""
		#~ nodes = self.explorebeneath(maxbs, boundbeneath=boundbeneath, boundabove=boundabove)
		#~ nodes += self.exploreabove(maxbs, boundbeneath=boundbeneath, boundabove=boundabove)
		#~ if not returnLabels:return nodes
		#~ else: return tree2.getlabels(nodes)
		
	def hierachical_path(self, labels=None, labprefix='N'):
		root = self.go_root()
		if not labels: labels = root.get_children_labels()
		else: labels = []
		sufix = ""
		path = root.path_to(self)
		for node in path:
			if (node.label() in labels) or node is self:
				sufix += '.'+node.label().lstrip(labprefix)
		return sufix

#####################################
####### EvolSequence related methods:

	def assign_seq(self,seq):
		"""Assign a sequence to the Node.
	
		WARNING: Intended to be used with the EvolSequence class to
		simulate evolution.
		"""
		
		self.__seq=seq
		return

	def evolve_seq(self,seq,mod, **kw): #verifier que les arguments marchent!!!
		"""Assign a sequence $1 to the Node, and evolve it along the sub-tree defined by the Node, according to model $2.
		
			WARNING: Must be used with the EvolSequence class in order to
			be able to simulate evolution.
	
		Keyword argument to determine how the number of substitutions
		on the sequence is approximated: [approx=string] 'Rounded' is
		the default. Keyword argument to determine the algorithm used
		to simulate neighbor-dependent substitutions: [algo=string]
		'Berard' is the default.
	
			Optional argument 'positions' is a list of all allowed
			positions. This argument defaults to all positions of the
			EvolSequence.
			
			WARNING: When an empty list is given, all positions are
			considered allowed.
		"""
	
		self.__seq=seq.copy()
		self.__seq.g_name(self.__lab)
		self.__seq.evolve(mod,self.__l, **kw)
			
		for c in self.__children:
			c.evolve_seq(self.__seq,mod, **kw)
			
	def evolve_seg_seq(self,seq,dmod,**kw):
		"""Assign a sequence $1 to the Node, and evolve it along the
		sub-tree defined by the Node, according to models in
		dictionary $2.
	
			The dictionary items are (deb,fin):mod where the modele mod is
			applied to the positions in range [deb:fin]. If ranges
			overlap, a random model is applied on each substitution that
			occurs in the overlap.
			
			WARNING: Must be used with the EvolSequence class in order to
			be able to simulate evolution.
	
		Keyword argument to determine how the number of substitutions
		on the sequence is approximated: [approx=string] 'Rounded' is
		the default. Keyword argument to determine the algorithm used
		to simulate neighbor-dependent substitutions: [algo=string]
		'Berard' is the default.
		"""

		self.__seq=seq.copy()
		self.__seq.g_name(self.__lab)
		self.__seq.evolve_seg(dmod,self.__l, **kw)
			
		for c in self.__children:
			c.evolve_seg_seq(self.__seq,dmod, **kw)
			
#####################################
#################### parsing methods:
	
	def clean(self,s):
		"""Clean the string in Newick format.
	
		Bracket-delimited nested and unnested comments are eliminated.
		"""
		s=s.strip()
		if s[len(s)-1]==';':
			if s.count('(')!=s.count(')'):
				raise ValueError, "Opening parenthesis do not match closing parenthesis."
			else:
				brackopen=s.count('[')
				brackclose=s.count(']')
				if brackopen!=0 or brackclose!=0:
					if brackopen==brackclose:
						open=s.find('[')
						x=1
						for i in range(open+1,len(s)):
							if x>0:
								if s[i]=='[':
									x+=1
								elif s[i]==']':
									x-=1
							else:
								break
						s=s[0:open]+s[i:len(s)]
						s=self.clean(s)
					else:
						raise ValueError, "Opening brackets do not match closing brackets."
		else:
			raise ValueError, "Missing ';' at end."
		return s

	def read_commented_lab(self, s, combrackets='[]', labquotes=False, namesAsNum=False, leafNamesAsNum=False):
		"""Parse node labels and bootstraps; deals with nested comments located next to labels (usually within brackets)"""
		if labquotes:
			lquote = s.find(labquotes)
			rquote = s.find(labquotes)
			self.__lab = s[lquote+1:rquote]
			self.read_commented_lab(s[lquote:]+s[:rquote+1], namesAsNum=namesAsNum, combrackets=combrackets, labquotes=False, leafNamesAsNum=leafNamesAsNum)
		else:
			lbrack = s.find(combrackets[0])
			rbrack = s.find(combrackets[1])
			if lbrack==-1 and rbrack ==-1:
				# clean case, no bracketed comments
				laboot = s
			elif lbrack == 0:
				# comment on the left side of string
				self.__comment = s[lbrack+1:rbrack]
				laboot = s[rbrack+1:]
			elif rbrack == len(s)-1:
				# comment on the right side of string
				self.__comment = s[lbrack+1:rbrack]
				laboot = s[0:lbrack]
			else:
				# comment in the middle, wrong syntax
				raise ValueError, "Comments brackets in the middle of the label, wrong syntax."
			if laboot:
				if namesAsNum or (self.is_leaf() and leafNamesAsNum):
					self.__lab = str(laboot)
				else:
					try:
						self.__boot=float(laboot)
					except ValueError, e:
						self.__lab = str(laboot)
						
	def read_commented_branch(self, s, combrackets='[]', bootInComm=False):
		"""Parse branch annotation (usually length); deals with nested comments located next to annotation (usually within brackets)
		
		This is found in RAxML's 'RAxML_bootstrapBranchLabels*' and derived 'RAxML_rootedTrees*' output files.
		See Czech L. et al. (2017) "A Critical Review on the Use of Support Values in Tree Viewers and Bioinformatics Toolkits." Mol Biol Evol. Jun; 34(6):1535–1542 doi:10.1093/molbev/msx055"
		"""
		lbrack = s.find(combrackets[0])
		rbrack = s.find(combrackets[1])
		if lbrack==-1 and rbrack ==-1:
			# clean case, no bracketed comments
			sbl = s
			c = None
		elif lbrack == 0:
			# comment on the left side of string
			c = s[lbrack+1:rbrack]
			sbl = s[rbrack+1:]
		elif rbrack == len(s)-1:
			# comment on the right side of string
			c = s[lbrack+1:rbrack]
			sbl = s[0:lbrack]
		else:
			# comment in the middle, wrong syntax
			raise ValueError, "Comments brackets in the middle of the branch annotation, wrong syntax."
		if sbl:
			try:
				self.__l = float(sbl)
			except ValueError, e:
				raise ValueError, "Incorrect branch length value: '%s' -> must be numerical."%(sbl)
		if c:
			if bootInComm:
				try:
					self.__boot=float(c)
				except ValueError, e:
					raise ValueError, "Incorrect branch support value: '%s' -> must be numerical."%(c)
			else:
				self.__comment = c

	def _parser(self, s, branch_lengths=True, keep_comments=False, combrackets='[]', labquotes=False, namesAsNum=False, leafNamesAsNum=False, bootInComm=False):
		"""Should not be directly used. Use parser() instead."""
		parenth=s.rfind(')')
		# deal with tree structure
		t=s[1:parenth] #eliminate external parenthesis
		x=0
		cuts=[0] #cutting points for nodes of the same depth
		incomment = 0
		for i in range(len(t)):
			if t[i]=='[': incomment += 1
			if t[i]==']': incomment -= 1
			if keep_comments and incomment: continue
			if t[i]=='(':
				x+=1
			elif t[i]==')':
				x-=1
			elif x==0 and t[i]==',':
				cuts+=[i+1]
		if cuts!=[0]:
			cuts+=[len(t)]
			i=0
			while i<len(cuts)-1:
				child=t[cuts[i]:cuts[i+1]].strip(',')
				self.__children+=[self.newnode()]
				for j in self.__children:
					j.__father=self
				self.__children[len(self.__children)-1]._parser(child, branch_lengths=branch_lengths, keep_comments=keep_comments, combrackets=combrackets, labquotes=labquotes, namesAsNum=namesAsNum, leafNamesAsNum=leafNamesAsNum)
				#print "***********Going to the children*************"
				i+=1
		# annotate node
		if branch_lengths:
			semicol=s.rfind(':')
			if semicol>parenth:
				# deal with branch annotations
				self.read_commented_branch(s[semicol+1:len(s)], combrackets=combrackets, bootInComm=bootInComm)
				# deal with node annotations
				self.read_commented_lab(s[parenth+1:semicol], namesAsNum=namesAsNum, combrackets=combrackets, labquotes=labquotes, leafNamesAsNum=leafNamesAsNum)
			elif (semicol==-1 and parenth!=-1) or (semicol < parenth): # added (semicol < parenth) condition to include case when root node has no ':brlength' terminal tag (before ultimate ';')
				self.read_commented_lab(s[parenth+1:], namesAsNum=namesAsNum, combrackets=combrackets, labquotes=labquotes, leafNamesAsNum=leafNamesAsNum)
			elif semicol==-1:
				raise ValueError, "Incorrect syntax."
		else:
			if not namesAsNum:
				try:
					self.__boot=float(s[parenth+1:])
				except ValueError, e:
					self.read_commented_lab(s[parenth+1:], namesAsNum=namesAsNum, combrackets=combrackets, labquotes=labquotes, leafNamesAsNum=leafNamesAsNum)
			else:
				self.read_commented_lab(s[parenth+1:], namesAsNum=namesAsNum, combrackets=combrackets, labquotes=labquotes, leafNamesAsNum=leafNamesAsNum)
		
		return
		
		# old version:
		#~ """Should not be directly used. Use parser() instead."""
		#~ parenth=s.rfind(')')
		#~ if branch_lengths:
			#~ semicol=s.rfind(':')
			#~ if semicol>parenth: # deal with node values
				#~ try:
					#~ self.__l=float(s[semicol+1:len(s)])
				#~ except ValueError, e:
					#~ raise ValueError, "Incorrect branch value -> must be numerical."
				#~ comlaboot = s[parenth+1:semicol]
				#~ if not self.is_leaf():
					#~ try:
						#~ self.__boot=float(comlaboot)
					#~ except ValueError, e:
						#~ self.read_commented_lab(comlaboot, namesAsNum=namesAsNum, combrackets=combrackets, leafNamesAsNum=leafNamesAsNum)
				#~ else:
					#~ self.read_commented_lab(comlaboot, namesAsNum=namesAsNum, combrackets=combrackets, leafNamesAsNum=leafNamesAsNum)
			#~ elif (semicol==-1 and parenth!=-1):
				#~ self.read_commented_lab(s[parenth+1:], namesAsNum=namesAsNum, combrackets=combrackets, leafNamesAsNum=leafNamesAsNum)
			#~ elif semicol==-1:
				#~ raise ValueError, "Incorrect syntax."
		#~ else:
			#~ if not namesAsNum:
				#~ try:
					#~ self.__boot=float(s[parenth+1:])
				#~ except ValueError, e:
					#~ self.read_commented_lab(s[parenth+1:], namesAsNum=namesAsNum, combrackets=combrackets, leafNamesAsNum=leafNamesAsNum)
			#~ else:
				#~ self.read_commented_lab(s[parenth+1:], namesAsNum=namesAsNum, combrackets=combrackets, leafNamesAsNum=leafNamesAsNum)
				#~ 
		#~ s=s[1:parenth] #eliminate external parenthesis
		#~ x=0
		#~ cuts=[0] #cutting points for nodes of the same depth
		#~ for i in range(len(s)):
			#~ if s[i]=='(':
				#~ x+=1
			#~ elif s[i]==')':
				#~ x-=1
			#~ elif x==0 and s[i]==',':
				#~ cuts+=[i+1]
		#~ if cuts!=[0]:
			#~ cuts+=[len(s)]
			#~ i=0
			#~ while i<len(cuts)-1:
				#~ child=s[cuts[i]:cuts[i+1]].strip(',')
				#~ self.__children+=[self.newnode()]
				#~ for j in self.__children:
					#~ j.__father=self
				#~ self.__children[len(self.__children)-1]._parser(child, branch_lengths=branch_lengths, combrackets=combrackets, namesAsNum=namesAsNum, leafNamesAsNum=leafNamesAsNum)
				#~ #print "***********Going to the children*************"
				#~ i+=1
#~ #		print self.bs(), self.comment(), self.lg()
		#~ return
	
	def parser(self, s, branch_lengths=True, keep_comments=False, combrackets='[]', labquotes=False, namesAsNum=False, leafNamesAsNum=False, bootInComm=False):
		"""Fill the Node's attributes from parsing $1 string.
	
		Follows the Newick standard for coding trees.
		Bracked delimited nested and unnested comments are ignored.
			
		Newick format example:
			   '((A:0.1,B:0.2,C:0.1)ABCnode:0.2,(D:0.4,E:0.1)99:0.1);'
	
		WARNING: all interior nodes (if given) should be strings ->
		numerical values will be interpreted as the node's bootstrap
		value.
		"""
	
		if not keep_comments:
			s=self.clean(s)
		s="".join(s.split()).strip(';') # remplace s=s.strip(';')
	
		#~ if branch_lengths:
			#~ # if no branch length defined at root, insert zero branch length value
			#~ parent = s.rfind(')')
			#~ if not ':' in s[parent+1:]:
				#~ s += ":0" 
			
		if s:
			self._parser(s, branch_lengths=branch_lengths, keep_comments=keep_comments, combrackets=combrackets, labquotes=labquotes, namesAsNum=namesAsNum, leafNamesAsNum=leafNamesAsNum, bootInComm=bootInComm)
			if branch_lengths and self==self.go_root():
				self.__l = 0 # sets distance above root at 0.
		else:
			raise ValueError, "This should not have happened! Check behind your back."
			

######################################################################
############################# DEV

	def doNNI(self, originalnodelabs=None, silent=True):
		uncle = self.go_uncle()
		if not uncle: raise IndexError, "node must have an uncle node to be attached to"
		return self.doSPR(newbro, originalnodelabs=originalnodelabs, silent=silent)
		
	def doSPR(self, newbro, originalnodelabs=None, silent=True):
		"""Subtree Prune and Regraft : prunes the node and regrafts it next to destination node newbro
		
		returns a tuple containing (the node that has been detelted from the gene tree, the node that will replace it, and the new node)
		NB: lengths and supports of new branches are set to None since modifying the topology would make previous estimates false according to the new topology
		"""
		genetree = self.go_root()
		if not originalnodelabs: originalnodelabs = genetree.get_children_labels()
		else: originalnodelabs = []
		oldfatlab = self.go_father().label()
		oldbrolab = self.go_brother().label()
		if not silent: print "SPR on %s: (%s,%s)%s"%(self.label(), oldbrolab, self.label(), oldfatlab),
		subtree = genetree.pop(self, keepRefToFather=False)
		(newfat, grandfat) = subtree.graftTo(newbro, originalnodelabs=originalnodelabs, silent=silent)
		newfatlab = newfat.label()
		if not silent:
			newcladestr = "(%s,%s)%s"%(newbro.label(), subtree.label(), newfatlab)
			if grandfat: print "-> (XXX,%s)%s"%(newcladestr, grandfat.label())
			else: print "-> %s"%newcladestr
		return (oldfatlab, oldbrolab, newfatlab)
		
	def graftTo(self, newbro, originalnodelabs=None, silent=True):
		"""Grafts the node as the neighbour of newbro 
		
		return a tuple containing (the new node (newfat) created above newbro and self, and its potential new father (grandfat))
		NB: lengths and supports of new branches are set to None since modifying the topology would make previous estimates false according to the new topology
		"""
		subtree = self
		grandfat = newbro.go_father()	# None if newbro is the root
		if grandfat: grandfat.unlink_child(newbro)
		newfat = self.newnode()
		newfat.complete_label(labels=originalnodelabs)
		newfat.link_child(subtree, newlen=subtree.lg(), newboot=subtree.bs())
		newfat.link_child(newbro, newlen=newbro.lg(), newboot=newbro.bs())
		if grandfat: grandfat.link_child(newfat)
		return (newfat, grandfat)
			
	def as_leaf(self, newlabel=None, silent=True):
		"""unlinks node from its potential children, making it a leaf"""
		if not silent: print "(before as_leaf) self.children_labels:", self.children_labels()
		while self.__children:
			c = self.__children[0]
			self.unlink_child(c, silent=silent)
		if not silent: print "(after as_leaf) self.children_labels:", self.children_labels()
		if newlabel: self.edit_label(newlabel)

	def link_child(self, newchild, newlen=None, newboot=None, silent=True):
		"""attach a child node to the node; falcutatively sets the properties of the interleaving branch in the child's attributes."""
		self.add_child(newchild, silent=silent)
		newchild.change_father(self, newlen=newlen, newboot=newboot, silent=silent)
		
	def unlink_child(self, child, silent=True):
		"""detach a child node from the node and erase the interleaving branch properties from the child's attributes.
		
		the properties of the former interleaving branch are return as a tuple (length, support).
		"""
		if not silent:
			print "self.__children:", self.__children
			print "child.__father:", repr(child.go_father())
		lbs = (child.lg(), child.bs())
		self.rm_child(child, silent=silent)
		child.change_father(None, silent=silent)
		if not silent:
			print "self.__children:", self.__children
			print "child.__father:", repr(child.go_father())
		return lbs
	
	def change_father(self, newfat, newlen=None, newboot=None, silent=True):
		self.__father = newfat
		self.__l = newlen
		self.__boot = newboot
		if not silent:
			newfatlab = newfat.label() if newfat else None
			print "%s <- %s %s:%s"%(newfatlab, self.label(), str(self.lg()), str(self.bs()))
	
	def add_child(self, newchild, silent=True):
		self.__children.append(newchild)
		if not silent: print "%s -< %s"%(self.label(), newchild.label())
	
	def rm_child(self, child, silent=True):
		self.__children.remove(child)
		if not silent: print "%s -/- %s"%(self.label(), child.label())
	
	def create_node_above(self, label=None, nodeid=None, labprefix='N', newlen=0, newboot=1, silent=True):
		"""add a new node betwwen the node and its father (if it has one)
		
		The branch length attribute value of the bottom node (self) is transferred the new node (thus the old branch length becomes the top branch length).
		The branch support attribute value of the bottom node (self) is kept as is (thus the old branch support becomes the bottom branch support). 
		This preserves the absence of (meaningful) support at leaf nodes.
		The new bottom branch (between self and the new node) will have a length of 0, unless specified otherwise through 'newlen'. 
		The new top branch (between the potential father and the new node) will have a support of 1, unless specified otherwise through 'newboot'. 
		returns the new node.
		"""
		newnode = self.newnode()
		if label:
			# decrement node label -- deprecated, do not use
			if label=='decrement': newnode.add_label(labprefix+str(int(self.label().strip(labprefix))-1))
			else: newnode.add_label(label)
		if nodeid:
			# decrement node_id -- deprecated, do not use
			if nodeid=='decrement': newnode.set_node_id(self.nodeid()-1)
			else: newnode.set_node_id(nodeid)
		fat = self.father
		if fat:
			brlen, bs = fat.unlink_child(self)
			fat.link_child(newnode, brlen, newboot)
		newnode.link_child(self, newlen=newlen, newboot=bs, silent=silent)
		return newnode
	
	def create_leafnode_below(self, label=None, nodeid=None, labprefix='N', newlen=0, newboot=None, silent=True):
		"""add a new leaf node below the node
		
		the new branch (between self and the new node) will have a length of 0 and a null support, unless specified otherwise. 
		returns the new node.		
		"""
		newnode = self.newnode()
		self.link_child(newnode, newlen=newlen, newboot=newboot, silent=silent)
		if label: newnode.add_label(label)
		if nodeid: newnode.set_node_id(nodeid)
		return newnode
		
	def reverse_fathertochild(self, silent=True):
		"""changes father of the Node into its child ; Node becomes orphan : no father, no branch length"""
		oldfat = self.go_father()
		oldfat.rm_child(self, silent)
		oldfat.change_father(self, self.__l, self.__boot, silent)
		self.add_child(oldfat, silent)
		self.__l=0
		self.__father=None
		self.__boot=None
	
	def path_from(self, fat):
		"""returns ordered list of nodes from $1 to Node, excluding both ; $1 must be an ancestor of Node
		
		DEPRECATED : better use self.path_to(n) where n and self have no constrained relative positions.
		"""
#		raise Warning, "DEPRECATED : better use self.path_to(n) where n and self have no constrained relative positions."
		if self not in fat.get_all_children():
			raise ValueError, "%s does not descend from %s"%(self.label(), fat.label())
		path = []
		while self not in fat.get_children():
			for c in fat.get_children():
				if self in c.get_all_children():
					fat = c
					path.append(fat)
					break
		return path

	def resolveNode(self, outgroups, recursive=False, trifurcateRoot=False, maxfurcate=2):
		"""turns unrooted tree/trifurcated node into explicitely rooted tree/bifurcated node
		
		Must create node(s) for that : branch length(s) from new node(s) to root/original node is set to 0.
		Acts recursively if there are more than 3 children at the node ; for that reason, outgroup is an ordered list of sub-root nodes, sequentially used to root unresolved nodes.
		By default, resolve recursively all children nodes until total resolution of the node or exhaustion of outgroups list. 
		By default, recursion is limited to the resolution of the focal node, but it can be propagated down the tree with recursive=True. 
		If outgroups='subroot', outgroup is always chosen as the first child ; with recursive=True, resolution will proceed that way in the whole tree, with arbitrary branching of new nodes.
		If maxfurcate=k, resolves down to a k-furcated node (can still be less than k-furcated if so at first).
		"""
		k = maxfurcate - 1
		if trifurcateRoot and self.is_root(): k = max(k, 2)
		if isinstance(outgroups, str):
			if outgroups=="subroot":
				og = self.get_children()[0:k]
				nog = "subroot"
			else:
				raise ValueError, "unproper outgroup definition"
		elif isinstance(outgroups, list) or isinstance(outgroups, tuple):
			og = outgroups[0:k]
			nog = outgroups[k:]
		else:
			raise ValueError, "unproper outgroup definition"
		children = self.get_children()
		if len(children) > k+1 and og:
			if not set(og) <= set(children):
				raise IndexError, "please provide sub-root node(s) as outgroup(s)"
			newSubRoot = self.newnode()
			for child in children:
				if child not in og:
					child.change_father(newSubRoot, child.lg(), child.bs())
					newSubRoot.add_child(child)
			newSubRoot.change_father(self, 0)
			self.__children = og + [newSubRoot]
			newSubRoot.resolveNode(nog, recursive=recursive, trifurcateRoot=trifurcateRoot, maxfurcate=maxfurcate)
		if recursive:
			for child in self.get_children():
				child.resolveNode(nog, recursive=recursive, trifurcateRoot=trifurcateRoot, maxfurcate=maxfurcate)
		
	# deprecated name
	defineNode = resolveNode

	def reRoot(self, node1, node2, silent=True, branch_lengths=True):	
		"""reroots the Node between children node1 and node2"""
		if not (isinstance(node1, Node) and isinstance(node2, Node)):
			raise TypeError, "expected a pair of tree.Node instances"
		if self!=self.go_root():
			raise ValueError, "Reroots only full tree/root node"		
		if node1 is self or node2 is self:
			raise ValueError, "Cannot reroot between an internal node and present root, please provide internal nodes or leaves"
		if node1.go_father() is node2:
			uncle = node2
			nephew = node1
		elif node2.go_father() is node1:
			uncle = node1
			nephew = node2
		elif node1.go_father()==node2.go_father()==self:
			if len(self.__children) > 2: 	# unresolved root
				raise ValueError, "Unrooted tree! please provide explicitely rooted tree"
		else:
			raise ValueError, "Unprecise location in tree, please provide adjacent nodes"

		subrootNodes = self.get_children()
		for i in range(len(subrootNodes)):
			if (node1 in subrootNodes[i].get_all_children()) and (node2 in subrootNodes[i].get_all_children()):
				newsr = subrootNodes.pop(i)	 # NB: this is list.pop() function. newsr = newsubroot = root of the subtree to change
				if not silent: 'newsr', newsr.label(), newsr.lg(), newsr.bs()
				break
		else:
			raise ValueError, "Topological problem, node1 and node2 should be in the same sub-root clade"
		
		if len(subrootNodes) > 1:
			# <=> root had > 2 children, after poping one from the list above
			supsr = self.newnode()					  # if subtree not needing to be changed is a multifurcation, a node is created in place of old root
			for sr in subrootNodes:
				sr.change_father( supsr, sr.lg(), sr.bs(), silent=silent )
				supsr.add_child(sr, silent)
		else:
			supsr = subrootNodes[0]
		if not silent:
			print 'supsr', supsr.label(), supsr.lg(), supsr.bs()
		if branch_lengths: nl = sumlen(newsr.distance_root(nullBranchesAsZeroLength=True), supsr.lg())
		else: nl = None
		supsr.change_father( newsr, newlen=nl, newboot=supsr.bs(), silent=silent ) # subtree not needing to be changed are rooted to newsr ; old root is skipped in topology
		newsr.add_child(supsr, silent)
	
		pathToReverse = nephew.path_from(newsr)		# !!! use of _DEPRECATED_ FUNCTION path_from() ; must be replaced by path_to().
		for interNode in pathToReverse:
			if not silent:
				print 'interNode', interNode.label()
			interNode.reverse_fathertochild(silent=silent)	 # path from old newsr to new root where father-child links must be reversed

		# sets new root
		if branch_lengths: 
			if isinstance(branch_lengths, list):
				if len(branch_lengths)!=2: raise ValueError, "'branch_lengths' parameter should be either a bool or a list of float of length 2"
				brlenneph, brlenuncl = branch_lengths
			else:
				brlenneph, brlenuncl = [nephew.lg()/2]*2
		else: 
			brlenneph, brlenuncl = [None]*2
		self.__children = [uncle, nephew]
		uncle.rm_child(nephew, silent)
		nephew.change_father(self, newlen=brlenneph, silent=silent)
		uncle.change_father(self, newlen=brlenuncl, silent=silent)
		return self

	def findRootBranch(self, other):
		"""in other, finds the branch where Node tree (self) is rooted
		
		(requires both trees have the same unrooted topology) ; returns a tuple of the two nodes of other that define the root branch.
		"""
		if self.hasSameRoot(other):
			return 0
		else:
			subRoots = self.get_children()
			mappedSubRoots = []
			for sr in subRoots:
				ll = sr.get_leaf_labels()
				msr = other.map_to_node(ll)
				if msr and msr != other and not msr.is_leaf(): # excludes mapping hits that match the entire other (means the searched leaf set forms a paraphyletic group in it) and unique leaves directly branching at root
					mappedSubRoots.append(msr)
	
			if len(mappedSubRoots)>=1:
				for msr in mappedSubRoots:
	#			   print msr
					if msr.go_father()!=other: # should not be the case anyway if Node tree and other are not rooted similarly AND ARE TOPOLOGICALLY EQUIVALENT
						# general case where one subRoot node exists in other as an internal monophyletic group and not the other (is paraphyletic and comprizes the root of other)
						return (msr, msr.go_father())	   
				raise ValueError, "Trees are not topologically equivalent"
				# subRoots of Node mapped in other were all subRoots too, what means they are rooted the same... but that subroot nodes are not topologically equivalent, ex: multifurcation at root in one, biffurcated node in the other
			else :
				raise ValueError, "No mapped subRoots"  
			
	def newOutgroup(self, outgroup, branch_lengths=True, silent=True):
		"""reroots the tree so that the specified node is the outgroup"""
		if outgroup is self:
			raise ValueError, "New outgroup cannot be the present root"
		outfat = outgroup.go_father()
		if not outfat.is_root():
			self.reRoot(outgroup, outfat, branch_lengths=branch_lengths, silent=silent)
			
	def outgroupInClade(self, outgroup):
		"""Returns the node where to root the tree so that 'outgroup' (a Node object being a leaf in 'self') is contained by the largest possible clade. 
		
		Finds the root where the outgroup is within the largest clade
		Used for rooting properly a gene tree knowing the most basal species in a reference species tree present in the gene tree.
		!!! should be used only on unicopy trees (trees where each species is represented at most once)!!!
		"""
		ll = self.get_leaves()
		if not outgroup in ll:
			raise TypeError, "'outgroup' must be a Node object being a leaf in 'self'"
		d_score_root = {}
		for node in self.get_all_children():
			if node != self:
				self.newOutgroup(node)
				subRootNodes = self.get_children()
				for subrootchild in subRootNodes:
					if outgroup.label() in subrootchild.get_leaf_labels():
						nl = outgroup.depth()
						d_score_root[nl] = d_score_root.setdefault(nl, []) + [node.label()]
		maxOutNode = self[ d_score_root[ max(d_score_root.keys()) ][0] ]
		return maxOutNode
		
	def rootGeneTreeAsRefTree(self, reftree, ltrans, silent=True, branch_lengths=True, ignoreLowBS=None, findARoot=False):
		"""Re-roots the gene tree (self) like the reference tree backbone (pruned of transfered leaves)"""
		if not silent:
			print "initial genetree"
			print self.newick()#ignoreBS=True)
			print "ltrans"
			print ltrans
		lleaves = self.get_leaf_labels()	
		### loads copy of gene tree pruned of transfered leaves
		prunedGenetree = copy.deepcopy(self)
		lgeneleaves = prunedGenetree.get_leaf_labels()
		dictleafspe = self.dictLeafLabelsToSpecies()
		dictspeleaf = self.dictSpeciesToLeafLabels()
		for leaf in lgeneleaves:
			if dictleafspe[leaf] in ltrans:
				prunedGenetree.pop(leaf)	
		if not silent:
			print "prunedGenetree"	
			print prunedGenetree.newick()				
		### loads copy of reference tree with the same set of leaf as in genetree pruned of transfered leaves
		prunedReftree = copy.deepcopy(reftree)
		lrefleaves = prunedReftree.get_leaf_labels()
		for refleaf in lrefleaves:
			if (not refleaf in dictspeleaf) or (refleaf in ltrans):
				prunedReftree.pop(refleaf)
		if not silent:
			print "prunedReftree"			
			print prunedReftree
		### reroots pruned gene tree consistently with pruned reference tree
		if not prunedGenetree.hasSameRoot(prunedReftree, useSpeDict=True):	#, ignoreLowBS=ignoreLowBS):
			srref = prunedReftree.get_children()
			# for each of the two groups under the reftree root
			for srnode in srref:
				if not silent: print "srnode", srnode.get_leaf_labels()
				target = prunedGenetree.map_to_node(srnode.get_leaf_labels(), useSpeDict=True)
				if not silent: print "target\n", target
				if target.is_monophyletic(srnode.get_leaf_labels(), useSpeDict=True):
					prunedGenetree.newOutgroup(target, branch_lengths=branch_lengths)
					break
			else:
				if findARoot:
					# find the biggest monophyletic group in gene tree that is included in one of the reference tree outgroups
					lgroups = []
					for srnode in srref:
						ssrl = set(srnode.get_leaf_labels())
						for srl in srnode.get_leaf_labels():
							gtl = prunedGenetree[srl]
							while gtl.go_father():
								if set(gtl.go_father().get_leaf_labels()) <= ssrl:
									gtl = gtl.go_father()
								else:
									if not gtl in lgroups:
										lgroups.append(gtl)
									break
					maxgroup = lgroups[0]
					if not silent:
						print "lgroups"
						print maxgroup.get_leaf_labels()
					for group in lgroups[1:]:
						if not silent: print group.get_leaf_labels()
						if group.nb_leaves() > maxgroup.nb_leaves():
							maxgroup = group
					if not silent: print "new outgroup:", maxgroup.get_leaf_labels()
					prunedGenetree.newOutgroup(group, branch_lengths=branch_lengths)
				else:
					raise IndexError, "Could not find similar root for 'prunedGenetree' and 'prunedReftree'"
		if not silent:
			print "rerooted prunedGenetree"
			print prunedGenetree.newick()#ignoreBS=True)
		### reroots (full) gene tree consistently with pruned gene tree
		sr = prunedGenetree.get_children()
		sr0 = self.map_to_node(sr[0].get_leaf_labels())
		sr1 = self.map_to_node(sr[1].get_leaf_labels())
		if not silent: print "subroot nodes:\n", sr0, "\n", sr1
		#~ if sr0 in sr1.get_children() or sr1 in sr0.get_children():
			#~ if not silent: print "root between:\n", sr0, "\n", sr1
			#~ self.reRoot(sr0, sr1, branch_lengths=branch_lengths)
		#~ elif sr[0].depth > sr[1].depth():
			#~ if not silent: print "new outgroup: %s\n"%sr0.label(), sr0
			#~ self.newOutgroup(sr0, branch_lengths=branch_lengths)
		#~ elif sr[1].depth > sr[0].depth():
			#~ if not silent: print "new outgroup: %s\n"%sr1.label(), sr1
			#~ self.newOutgroup(sr1, branch_lengths=branch_lengths)
		if sr0 in sr1:
			if not silent: print "new outgroup: %s\n"%sr0.label(), sr0
			self.newOutgroup(sr0, branch_lengths=branch_lengths)
		elif sr1 in sr0:
			if not silent: print "new outgroup: %s\n"%sr1.label(), sr1
			self.newOutgroup(sr1, branch_lengths=branch_lengths)
		else:
			if not silent: print "root arbitrarily with: %s\n"%sr0.label(), sr0
			self.newOutgroup(sr0, branch_lengths=branch_lengths)
		#~ else:
			#~ raise IndexError, "Unable to root gene tree consistently"
		if not silent:
			print "rerooted genetree"
			print self.newick(ignoreBS=True)
			
	def reimplantSubtree(self, outgroup):
		"""change local rooting of a subtree inside a full tree, given an outgroup for this subtree"""
		stf = self.go_father()
		if not stf:
			self.newOutgroup(outgroup)
		else:
			# uncouples the gene tree and subtree
			self.change_father(None, newlen=self.lg(), newboot=self.bs())
			stf.rm_child(self)
			# reroot the subtree
			self.newOutgroup(outgroup)
			# re-connects the two parts
			self.change_father(stf, newlen=self.lg(), newboot=self.bs())
			stf.add_child(self)
			
	def implantForGroupMonophyly(self, subtree, recleaftoreroot):
		"""change local rooting of a subtree inside a full tree, so that the input leaf set is momophyletic
		
		'subtree' can be a part of 'self' (same Node object, full gene tree) or a distinct Node object matching the topology of a subtree of 'self'.
		"""
		# we must find a "safe" place to reroot, i.e. where new bipartitions and incidental new branch supports will not contradict the preceding Prunier analysis of the subtree.
		outgroupnode = None
		# such a "safe" place is found on the branch leading to the clade of 'recleaves' that is not monophyletic (the 'recleaftoreroot'); 
		# as 'recleaftoreroot' is not monophyletic according to the current rooting, the outgroup used to root is the remaining part of 'subtree' (the other side of the branch where Prunier has cut)
		for initleaf in recleaftoreroot:
			initnode = subtree[initleaf]
			upnode = initnode.go_father()	
			# 'recleaftoreroot' is split in two parts, one is monophyletic, the other is paraphyletic, being the brother node of the sought outgroup
			# starting from one pruned leaf, find the upper node that contain only pruned leaves (last value of 'initnode')
			while set(upnode.get_leaf_labels()) <= set(recleaftoreroot):
				initnode = upnode
				upnode = initnode.go_father()
			if upnode is subtree:
				# initleaf is located in the monophyletic part of 'recleaftoreroot'
				continue # the for initleaf loop, to find one under the other part of 'recleaftoreroot'
			# brother of 'initnode' will serve as outgroup
			outgroupnode = initnode.go_brother()
			if not subtree in self:
				# reroot the subtree
				subtree.newOutgroup(outgroupnode, branch_lengths=False)
			break # the for initleaf loop
		else:
			raise IndexError, "could not find the branch where to root"
		if not subtree in self:
			# must reroot the correspoding subtree in-place in the full gene tree
			gtst = genetree.map_to_node(subtree.get_leaf_labels())
			gtog = genetree.map_to_node(outgroupnode.get_leaf_labels())
		else:
			gtst = subtree
			gtog = outgroupnode
		gtst.reimplantSubtree(gtog)
		
	def cureTree(self, reftree, ltrans=[], branch_lengths=True, tryRooting=True, justTryrooting=False):
		""" !!!! DEPRECATED FONCTION, USING reRoot() (hazardous)"""
		subRootNodes = self.__children
		self.resolveNode(subRootNodes[0:1])
		self.complete_internal_labels()	
		if tryRooting:
			try:
				self.rootGeneTreeAsRefTree(reftree, ltrans=ltrans, silent=True, branch_lengths=branch_lengths)
			except IndexError, e:
				if not justTryrooting:
					raise IndexError, e
				
	def getMRP(self, boot_thresh=0.5, taxset=None, matrixformat='dictoflist', writeto=None, order=1, **kw):
		"""return the Matrix Representation with Parsimony of the tree as in Daubin et al. (2004), Genome Res. 12:1080-1090
		
		writeto is None (default): matrix is returned as a dict of lists 
		writeto is a str : write PHYLIP format matrix to this file path; matrix still is returned as a dict of lists 
		writeto is False : matrix is returned as a 2-D nested lists;
		"""
		def addtodictoflist(mat, key, val, addrow=False):
			if addrow: mat[key] = []
			mat[key].append(val)
		def addtolistoflist(mat, key, val, addrow=False):
			if addrow: mat.append([])
			mat[-1].append(val)
		
		addfun = eval('addto'+matrixformat)
		sleaves = self.get_leaf_labels()
		if taxset: lleaves = list(taxset)
		else: lleaves = sleaves
		# dictionary of taxa (rows) to bipartitions (colunms)
		mrp = {}
		n = 0
		for node in self.get_sorted_children(order=order):
			# ignore root node (trivial bipartition and nodes with low branch support)
			if node.is_root() or node.bs()<boot_thresh: continue
			ll = node.get_leaf_labels()
			for leaflab in lleaves:
				if leaflab in ll: addfun(mrp, leaflab, True, (not bool(n)))
				elif leaflab in sleaves: addfun(mrp, leaflab, False, (not bool(n)))
				else: addfun(mrp, leaflab, None, (not bool(n)))
			n += 1
		if writeto:
			if matrixformat!='dictoflist': raise ValueEror, "requires matrixformat='dictoflist' for writing out matrix"
			tree2.writeMRP(mrp, writeto, **kw)
		return mrp
		
	### automated classification of sequences into coherent clades/genotypes 

	def prune_genotypes(self, relvarthresh=5, minseqwithin=5, minvarwithin=0, minbs=0, fixnbcut=0, returnLabels=True, annotateOnTree=False, multipass=True, silent=True):
		"""reccursively search for branches that increase significantly the variance of branch lengths compared to the remainder of the branches in the clade below, and prune the outliers"""
		t = copy.deepcopy(self)
		lst = []
		nbcut= 0
		for node in t.get_sorted_children(order=3):
			if node.nb_all_children()>=minseqwithin:
				if not silent: print node.label(), node.nb_children(), 'children'
				brlens = dict([(n, n.lg()) for n in node.get_sorted_children(order=3) if n!=node])
				brboots = dict([(n, n.bs()) for n in node.get_sorted_children(order=3) if n!=node])
				varwithin = np.array(brlens.values()).var()
				# test the difference of variance between the branch lengths within the subtree with all branches or when removing one of the branches
				trhough = False
				while (not trhough) and len(brlens)>=minseqwithin:
					for n in brlens:
						varexcl = np.array([brlens[k] for k in brlens if k!=n]).var()
						meanexcl = np.array([brlens[k] for k in brlens if k!=n]).mean()
						if not silent: print 'varwithin', varwithin, 'varexcl', varexcl
						if brboots[n]>minbs and brlens[n]>meanexcl and varexcl>=minvarwithin and (varwithin/varexcl)>relvarthresh:
							# inner branch is significantly longer than the rest of the subtree, record the branch for cutting
							if (fixnbcut and nbcut<fixnbcut):
								nbcut += 1
								if not silent: print '\tcut', n.label()
								st = t.pop(n, keepRefToFather=True)
								lst.append(st)
								# update the distribution of branch lengths
								brlens = dict([(n, n.lg()) for n in node.get_sorted_children(order=3) if n!=node])
								varwithin = np.array(brlens.values()).var()
								break
					else:
						trhough = True
				if (not node.is_root()) and len(brlens)>=minseqwithin and varwithin>=minvarwithin:
					# test the difference of variance between the branch lengths within the subtree or when including the parent branch
					meanwithin = np.array(brlens.values()).mean()
					if node.is_subroot() and node.go_father().is_bifurcated():
						parlen = node.lg() + node.go_brother().lg()
						parboot = max([0]+[n.bs() for n in (node, node.go_brother())])
					else:
						parlen = node.lg()
						parboot = node.bs() if node.bs() else 0
					varwithfat = np.array(brlens.values()+[parlen]).var()
					if not silent: print 'varwithfat', varwithfat, 'varwithin', varwithin
					if parboot>minbs and parlen>meanwithin and varwithin>=minvarwithin and (varwithfat/varwithin)>relvarthresh:
						# parent branch is significantly longer than the subtre below, cut the subtree
						if (fixnbcut and nbcut<fixnbcut):
							nbcut += 1
							if not silent: print '\tcut', node.label()
							st = t.pop(node, keepRefToFather=True)
							lst.append(st)
		lst.append(t)
		if (fixnbcut and nbcut>=fixnbcut) or (not multipass):
			llst = lst
		else:
			if nbcut>=1:
				if not silent: print 'next pass (recursive)'
				print len(lst), fixnbcut
				llst = []
				for st in lst:
					if (fixnbcut and nbcut>=fixnbcut):
						# already reached the maximum number of cuts, just transfer the subtree
						llst.append(st)
					else:
						# !!! reccursive call !!!
						rlst = st.prune_genotypes(relvarthresh=relvarthresh, minseqwithin=minseqwithin, minvarwithin=minvarwithin, fixnbcut=(fixnbcut-nbcut), returnLabels=False, annotateOnTree=False, multipass=multipass, silent=silent)
						nbcut += len(rlst)-1
						llst += rlst
			elif fixnbcut:
				if not silent: print 'next pass (recursive)'
				print len(lst), fixnbcut
				llst = []
				for st in lst:
					if (fixnbcut and nbcut>=fixnbcut):
						# already reached the maximum number of cuts, just transfer the subtree
						llst.append(st)
					else:
						# nothing was found, continue but with all the thresholds lowerred by 2
						# !!! reccursive call !!!
						rlst = st.prune_genotypes(relvarthresh=relvarthresh/2, minseqwithin=minseqwithin/2, minvarwithin=minvarwithin/2, fixnbcut=(fixnbcut-nbcut), returnLabels=False, annotateOnTree=False, multipass=False, silent=silent)
						nbcut += len(rlst)-1
						llst += rlst
			else:
				llst = lst
		genolabs = [n.get_leaf_labels() for n in llst]
		if annotateOnTree:
			tt = copy.deepcopy(self)
			n = 0
			for gl in genolabs:
				n += 1
				for lab in gl:
					tt[lab].edit_label(lab+'_#%d'%n)
			return tt
		if returnLabels: return genolabs
		else: return llst

