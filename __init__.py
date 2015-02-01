#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Object-oriented phylogenetic tree manipulation package.

Copyright 2007, Leonor Palmeira <palmeira@biomserv.univ-lyon1.fr>;
Copyright 2013, Florent Lassalle <florent.lassalle@univ-lyon1.fr> or <florent.lassalle@ucl.ac.uk>.
"""

__all__ = ["Node", "AnnotatedNode", "ReferenceTree", "GeneTree"]


__author__ = "Florent Lassalle <florent.lassalle@ucl.ac.uk>"
__date__ = "14 January 2015"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""


import copy
import os
import cPickle

#~ from tree2.Utils import *

from tree2.Node import Node
from tree2.AnnotatedNode import AnnotatedNode
from tree2.ReferenceTree import ReferenceTree
from tree2.GeneTree import GeneTree

#######################################################################
########## misc functions

def getlabels(lnodes):
	llab = []
	for node in lnodes:
		llab.append(node.label())
	return llab
	
def setToStr(s, sep="|"):
	l = list(s)
	l.sort()
	return sep.join(l)
	
def setToStr2(ls, andsep="&"):
	l = []
	for s in ls:
		l.append(setToStr(s))
	return andsep.join(l)
	

	
#######################################################################
########## file I/O functions

def read_multiple_newick(nf, treeclass="Node", **kw):
	"""reads a multiple Newick file and returns a list of trees"""
	tc = eval(treeclass)
	if not issubclass(tc, Node): raise ValueError, "wrong tree class specification"
	f = open(nf, 'r')
	lt = []
	for line in f:
		t = tc(newick=line, **kw)
		#~ if treeclass=="Node": t = Node(newick=line, **kw)
		#~ elif treeclass=="AnnotatedNode": t = AnnotatedNode(newick=line, **kw)
		#~ elif treeclass=="ReferenceTree": t = ReferenceTree(newick=line, **kw)
		#~ elif treeclass=="GeneTree": t = GeneTree(newick=line, **kw)
		#~ else: raise ValueError, "wrong tree class specification"
		lt.append(t)
	return lt	
	
def read_nexus(nf, treeclass="Node", returnDict=True, translateLabels=True, getTaxLabels=False, allLower=True, **kw):
	#~ tc = getattr(tree2, treeclass)
	tc = eval(treeclass)
	if not issubclass(tc, Node): raise ValueError, "wrong tree class specification"
	datadim = None
	ltax = []
	dtax = {}
	ltrees = []
	ltreenames = []
	block = None
	labelrows = False
	f = open(nf, 'r')
	for line in f:
		lst = line.strip(' \t\n\r')
		if allLower: lst = lst.lower()
		if lst.startswith('#') or lst=='': continue
		if lst.startswith('begin') and lst.endswith(';'): block = lst.rstrip(';').split('begin ')[1]
		if lst.startswith('end') and lst.endswith(';'): block = None
		if (not labelrows) and ltax and getTaxLabels=='list': return ltax
		if (not labelrows) and dtax and getTaxLabels=='dict': return dtax
		if block:
			if block=='taxa':
				if lst.startswith('dimensions'): datadim = int(lst.split('dimensions ntax=')[1].strip(';'))
				elif lst.startswith('taxlabels'): labelrows = True
				elif labelrows and lst==';': labelrows = False ; assert datadim in [len(ltax), None]
				elif labelrows: ltax.append(lst)
			if block=='trees':
				n = 1
				multilinecomment = False
				tt = None
				#~ print labelrows, 1
				if lst.startswith('['): multilinecomment = True
				if multilinecomment:
					if ']' in lst: multilinecomment = False
					continue
				elif lst.startswith('translate'): labelrows = True
				elif labelrows and lst==';': labelrows = False ; assert datadim in [len(dtax), None] # ; print dtax
				elif labelrows: tnum = lst.rstrip(',').split() ; dtax[tnum[0].strip("'")] = tnum[1].strip("'")
				elif lst.startswith('tree '): tt = lst.split('tree', 1)[1].strip(' ').partition(' = ')
				if tt:
					treename, sep, nwk = tt
					if nwk.lower().startswith('[&r]'): nwk = nwk.lstrip(' [&rR]')
					t = tc(newick=nwk, leafNamesAsNum=(len(dtax)>0), **kw)
					#~ if treeclass=="Node": t = Node(newick=nwk, leafNamesAsNum=(len(dtax)>0), **kw)
					#~ elif treeclass=="AnnotatedNode": t = AnnotatedNode(newick=nwk, leafNamesAsNum=(len(dtax)>0), **kw)
					#~ elif treeclass=="ReferenceTree": t = ReferenceTree(newick=nwk, leafNamesAsNum=(len(dtax)>0), **kw)
					#~ elif treeclass=="GeneTree": t = GeneTree(newick=nwk, leafNamesAsNum=(len(dtax)>0), **kw)
					#~ else: raise ValueError, "wrong tree class specification"
					if translateLabels and dtax:
						for leaf in t.get_leaves():
							leaf.edit_label(dtax[leaf.label()])
						#~ for num in dtax:
							#~ t[num].edit_label(dtax[num])
					ltrees.append(t)
					if treename:
						if treename in ltreenames: ltreenames.append("%s-%d"%(treename, n))
						else: ltreenames.append(treename)
					else:
						ltreenames.append(n)
					n += 1
	f.close()
	if returnDict:
		dtree = dict(zip(ltreenames, ltrees))
		if dtax: dtree['translate'] = dtax
		dnexus = {'tree':dtree}
		if ltax: dnexus['taxlabels'] = ltax
		if datadim: dnexus['dimensions'] = datadim
		return dnexus
	else:
		return ltrees

def write_to_file(s, nf, **kw):
	"""Write a string (e.g. representing the Node in several formats) in $1 file.

	Keyword argument to determine whether 'append' or 'write'
	(overwrites) mode is turned on: [mode=string]."""
	mode = kw.get('mode', 'write')
	if type(nf)==str:
		if mode=='append':
			f=open(nf, 'a')
		elif mode=='write':
			f=open(nf, 'w')
		f.write(s+'\n')
		f.close()
	else:
		raise ValueError, "Unvalid file name. Unable to write."

def load_pickle(nfile):
	"""loads a pickled tree2 object"""
	fpickle = open(nfile, 'r')
	gtpickle = cPickle.Unpickler(fpickle)
	genetree = gtpickle.load()
	fpickle.close()
	return genetree

def dump_pickle(genetree, nfile):
	"""dumps a tree2 object as pickle"""
	fpickle = open(nfile, 'w')
	gtpickle = cPickle.Pickler(fpickle, 2)
	gtpickle.dump(genetree)
	fpickle.close()

####################################
### External program call functions

### Visualization functions
homedir = os.environ['HOME']
def seaview(node, seaview_exec_path='seaview', tmpfiledir='%s/tmp'%homedir, ignoreBS=False, branch_lengths=True, unrooted=False, bg=True, comment='comment'):
	if node.seaview:
		node.seaview(seaview_exec_path=seaview_exec_path, tmpfiledir=tmpfiledir, ignoreBS=ignoreBS, branch_lengths=branch_lengths, unrooted=unrooted, bg=bg, comment=comment)
	else:
		raise TypeError, "no appropriate method for %s object"%str(type(node))

def figtree(node, figtree_exec_path='figtree', tmpfiledir='%s/tmp'%homedir, ignoreBS=False, branch_lengths=True, unrooted=False, bg=True):
	if node.figtree:
		node.figtree(figtree_exec_path=figtree_exec_path, tmpfiledir=tmpfiledir, ignoreBS=ignoreBS, branch_lengths=branch_lengths, unrooted=unrooted, bg=bg)
	else:
		raise TypeError, "no appropriate method for %s object"%str(type(node))

def archaeopteryx(node, treename="", arch_exec_path='archaeopteryx', tmpfiledir='%s/tmp'%homedir, bg=True):
	if node.archaeopteryx:
		node.archaeopteryx(treename=treename, arch_exec_path=arch_exec_path, tmpfiledir=tmpfiledir, bg=bg)
	else:
		raise TypeError, "no appropriate method for %s object"%str(type(node))
		

### Phylogenetic computation functions

def phyml(node, alnfilepath, outtreefilepath, PhyMLpath='/usr/bin/phyml', quiet=True, **kw):
	"""perform PhyML command with the node as an user input tree
	
	'alnfile' must be the path of an input sequence alignment in Phylip format.
	PhyML options are provided through **kw arguments. ex: for '-c 8' command line option, add 'c=8' keyword argument.
	'-u' and '-i' options are ignored as they are efined by node, alnfilepath
	default options make the program to recompute branch lengths and SH-like supports 
	from a given topology (input Node object) and a nucleic alignment under a GTR+G8+I model.
	"""
	node.phyml(alnfilepath=alnfilepath, outtreefilepath=outtreefilepath, PhyMLpath=PhyMLpath, quiet=quiet, **kw)



###################################
### distance-associated	functions

def dNoverdS(dNtree, dStree, factor=3, ignoreThreshold=0.05, medianForIgnored=True):
	"""from two trees bearing as branch lengths the counts of non-synonimous or synonimous substitutions, respectively, makes a tree with dN/dS ratio as branch lengths
	
	can take as input trees outputed by substitution mapping analysis with bppML (Pouyet, Guéguen and Bailly-Béchet, 2012).
	ratio is dN/(k*dS), with k being the ratio of non-synonimous-to-synonimous possible substitutions given the observed codons in the alignment ; 
	k is generally close to 3 (if the full genetic code is represented), the default.
	root is ignored (length set to 0).
	if dS tree has branches of length 0 (or under a length threshold), the information is not sufficient and
	the output tree will have corresponding branches with ratio set to the median of other ratio in the tree.
	"""	
	dNtree.complete_internal_labels()
	ratiotree = copy.deepcopy(dNtree)
	for rationode in ratiotree:
		rationode.set_lg(None)
	# define length thresholds under wich a branch is ignored because information is too small
	# i.e. if length is < ignoreThreshold * mean branch length in the tree
	toosmalldN = ignoreThreshold * (dNtree.treelength()/len(dNtree.get_all_children()))
	toosmalldS = ignoreThreshold * (dStree.treelength()/len(dStree.get_all_children()))
	lignored = []
	for dNnode in dNtree:
		dSnode = dStree.map_to_node(dNnode.get_leaf_labels())
		if (not dSnode.lg()) or (dSnode.lg() <= toosmalldN) or (dNnode.lg() <= toosmalldS):
			lignored.append(ratiotree[dNnode.label()])
		else:
			ratio = dNnode.lg()/(factor*dSnode.lg())
			ratiotree[dNnode.label()].set_lg(ratio)
	if medianForIgnored:
		medianratio = ratiotree.treemedian(lignored=lignored)
		for rationode in lignored:
			rationode.set_lg(medianratio)
	return ratiotree
	
#########################################
#### functions dealing with multiple trees

def concatMRP(treelist, boot_thresh=0.5, taxset=None, writeto=None, **kw):
	if taxset:
		lleaves = list(taxset)
	else:
		stax = set(treelist[0].get_leaf_labels())
		for tree in treelist[1:]:
			stax |= set(tree.get_leaf_labels())
		lleaves = list(stax)
	lmrp = [tree.getMRP(writeto=None, boot_thresh=boot_thresh, taxset=lleaves) for tree in treelist]
	cmrp = {}
	for leaf in lleaves:
		for i in range(len(lmrp)):
			if lmrp[i]:
				row = cmrp.setdefault(leaf, [])
				row += lmrp[i][leaf]
	if writeto!=None: writeMRP(cmrp, writeto, **kw)
	return cmrp
	
def writeMRP(mrp, nfout, outfmt='phylip', codeout={True:'1', False:'0', None:'?'}):
	fout = open(nfout, 'w')
	if outfmt.startswith('phylip'):  fout.write("    %d    %d\n"%(len(mrp), len(mrp.values()[0])))
	for leaf in cmrp: 
		if outfmt=='phylip': fout.write(leaf.ljust(10))
		elif outfmt=='phylip.relax': fout.write("%s\t"%leaf)
		elif outfmt=='fasta':  fout.write(">%s\n"%leaf)
		else: raise ValueError, "unsuported output format: %s"%outfmt
		fout.write(''.join([codeout[bi] for bi in cmrp[leaf]])+"\n")
	fout.close()
