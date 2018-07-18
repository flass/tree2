#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Object-oriented phylogenetic tree manipulation package.

Copyright 2007, Leonor Palmeira <palmeira@biomserv.univ-lyon1.fr>;
Copyright 2013, Florent Lassalle <florent.lassalle@univ-lyon1.fr>;
Copyright 2015, Florent Lassalle  <florent.lassalle@ucl.ac.uk>;
Copyright 2016, Florent Lassalle  <f.lassalle@imperial.ac.uk>.
"""

__all__ = ["Node", "AnnotatedNode", "ReferenceTree", "GeneTree", "FwdTreeSim", "Exceptions", "svgNode"]


__author__ = "Florent Lassalle <florent.lassalle@ucl.ac.uk>"
__date__ = "14 January 2015"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the Node class."""


import copy
import os
import cPickle

#~ from tree2.Utils import *

from tree2.Node import Node, sumlen
from tree2.AnnotatedNode import AnnotatedNode
from tree2.ReferenceTree import ReferenceTree
from tree2.GeneTree import GeneTree
#~ from tree2.Exceptions import *

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
	
def updateappend(d1, d2):
	"""update dictionary d1 with items from dictionary d2; similar to dict.update(), but appends LIST object (value from d2) to existing values (in d1) rather than overwriting it"""
	for key in d2:
		if key in d1: d1[key] += d2[key]
		else: d1[key] = d2[key]
	
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
				elif lst.startswith('translate'):
					labelrows = True
				elif labelrows:
					tnum = lst.rstrip(';,').split()
					if tnum: dtax[tnum[0].strip("'")] = tnum[1].strip("'")
					if lst.endswith(';'):
						labelrows = False ; assert datadim in [len(dtax), None] # ; print dtax
				elif lst.startswith('tree '): tt = lst.split('tree', 1)[1].strip(' ').partition(' = ')
				if tt:
					treename, sep, nwk = tt
					if nwk.lower().startswith('[&'): 
						rootstatus, sep, nwk = nwk.partition('] ')
						rootstatus = rootstatus.lstrip('[&').lower()
					else:
						rootstatus = None
					t = tc(newick=nwk, leafNamesAsNum=(len(dtax)>0), **kw)
					# tag the tree for being originally unrooted or not (all Node objects are intrinsically rooted) for later writing back into the right form
					t.unrooted = True if rootstatus=='u' else False
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
		
def write_nexus(ltrees, nfout, ltax, dtranslate={}, ltreenames=[], mode='w', onlytrees=False, **kw):
	"""combine in one nexus file a collection of trees with the same set of taxa represented"""
	if ltreenames: assert len(ltreenames)==len(ltrees)
	if dtranslate:
		assert len(dtranslate)==len(ltax)
		assert isinstance(dtranslate.keys()[0], str)
	if isinstance(nfout, file):
		if nfout.mode==mode and not nfout.closed:
			fout = nfout
	elif os.path.exists(os.path.dirname(nfout)):
		fout = open(nfout, mode)
	else:
		raise ValueError, "output file argument is nor a valid path nor a file writeable in %s mode"%mode
	if not onlytrees:
		fout.write('#NEXUS\n')
		## taxa block
		fout.write('begin taxa;\n\tdimensions ntax=%d;\n\t\ttaxlabels\n'%len(ltax))
		for tax in ltax:
			fout.write('\t\t%s\n'%tax)
		else:
			fout.write('\t\t;\nend;\n')
	## trees block
	fout.write('begin trees;\n')
	if dtranslate:
		fout.write('\ttranslate\n')
		nums = [int(num) for num in dtranslate.keys()]
		nums.sort()
		for num in nums:
			fout.write('\t\t%d\t%s,\n'%(num, dtranslate[str(num)]))
		else:
			fout.write('\t\t;\n')
	# tree blocks
	for k in range(len(ltrees)):
		if ltreenames: treename = ltreenames[k]
		else: treename = 'tree_%d'%k
		fout.write('   tree %s = %s\n'%(treename, ltrees[k].newick(**kw)))
	else:
		fout.write('end;')
	fout.close()
			

def write_to_file(s, nf, **kw):
	"""Write a string (e.g. representing the Node in several formats) in $1 file.

	Keyword argument to determine whether 'append' or 'write'
	(overwrites) mode is turned on: [mode=string]."""
	mode = kw.get('mode', 'write')
	if type(nf)==str:
		if mode in ['a', 'append']:
			f=open(nf, 'a')
		elif mode in ['w', 'write']:
			f=open(nf, 'w')
		else:
			raise ValueError, "unknown file mode: %s"%mode
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
	
#######################################################################
########  XML Tree Format related functions

def extract_XMLmarker(s, mark, ind="", mtype=""):
	"""extracts the information nested in a XML marker from a XML string."""
	if mtype:
		marker = "%s type='%s'"%(mark, mtype)
	else:
		marker = mark
	if s.startswith('%s<%s>'%(ind, marker)):
		return s.split('<%s>'%marker)[-1].split('</%s>'%mark)[0]
	else:
		return None

def concat_phyloxml(lnftree, nfout):
	"""Concatenates phyloXML trees from list of file names 'lnftree' into file named 'nfout' """
	lines = []
	head = ""
	tail = ""	
	for nftree in lnftree:
		ftree = open(nftree, 'r')
		tree = ftree.readlines()
		ftree.close()
		if not (tree[0].startswith('<?') and tree[1].startswith('<phyloxml')):
			raise SyntaxError, "Wrong syntax for phyloXML header:\n%s"%'\n'.join(tree[0:2])
		head = tree[0:2]
		if tree[-1]!='</phyloxml>':
			raise SyntaxError, "Wrong syntax for phyloXML tail:\n%s"%tree[-1]
		tail = tree[-1]
		lines += tree[2:-1]	
	fout = open(nfout, 'w')
	fout.writelines(head)
	fout.writelines(lines)
	fout.write(tail)
	fout.close()
	
def read_phyloXML(nfin, branch_lengths=True, keep_comments=False, ind="  ", indstart="\n"):
	"""Reads the $1 file containing (multiple) tree(s) in PhyloXML format and builds the Node(s) from it.
	
	Returns a dictionary with keys as tree name strings and values as Nodes. 
	"""
	if type(nfin)==str:
		f=open(nfin,'r')
		s=f.read()
		f.close()
		l = s.split(indstart)
	else:
		raise ValueError, "Invalid file name."
		
	d = {}
	start = 0
	end = len(l)-1
	while l[start].startswith('<?'):
		start += 1
	if l[start].startswith('<phyloxml'):
		start += 1
	else:
		raise SyntaxError, "Invalid phyloXML syntax in file %s at line %d."%(nfin, start+1)
	if l[end]=='</phyloxml>':
		end -= 1
	else:
		raise SyntaxError, "Invalid phyloXML syntax in file %s at line %d."%(nfin, end+1)
	
	ntree = l.count('</phylogeny>')
	for k in range(ntree):
		treestart = start
		treestop = l[start+1:].index('</phylogeny>') + start + 1
		start = treestop+1
		if l[treestart].startswith('<phylogeny'):
			treestart += 1
		else:
			raise SyntaxError, "Invalid phyloXML syntax in file %s at line %d."%(nfin, treestart+1)	
		name = extract_XMLmarker( l[treestart], 'name', ind)
		if name:
			treestart += 1
		else:
			name = str(k+1)
		node = AnnotatedNode()
		s = indstart.join(l[treestart:treestop])
		node.parse_phyloXML(s=s, branch_lengths=branch_lengths, ind=ind, indstart=indstart+ind)
		d[name] = node
	return d
	
def write_phyloXML(dtrees, nfout, ind="  ", indstart="\n", normparams=None):
	""" writes the tree in PhyloXML format, readable by Achaeopteryx (http://www.phylosoft.org/archaeopteryx/)"""
	fout = open(nfout, 'w')
	fout.write("<?xml version='1.0' encoding='UTF-8'?>%s<phyloxml xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' xsi:schemaLocation='http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd' xmlns='http://www.phyloxml.org'>"%(indstart))
	for name in dtrees:
		fout.write("%s<phylogeny rooted='true'>"%indstart+ind)
		fout.write("%s<name>%s</name>"%(indstart+ind,name))
		if normparams:
			fout.write(dtrees[name].phyloXML_norm(ind, indstart+ind*2, normparams))
		else:
			fout.write(dtrees[name].phyloXML(ind, indstart+ind*2))
		fout.write("%s</phylogeny>"%indstart+ind)
	fout.write("%s</phyloxml>"%indstart)
	fout.close()	
		
def read_obo(nfobo):
	"""reads in OBO XML file and parse OBO structure to make a tree of term hierarchy ([to be] developped for reading Gene Ontology [GO] terms hierarchy)"""


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
		

### Phylogenetic inference functions

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
### selection detection functions

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
	
def getnodefromtreelist(n, treelist):
	for t in treelist:
		k = t[n]
		if k: return t[n]
	else:
		return None

def concatMRP(treelist, boot_thresh=0.5, taxset=None, writeto=None, **kw):
	"""generate a matrix representation of the trees, one line per bipartition.
	
	default coding: 0=absent from the tree, 1=present in the tree.
	for use in Matrix Representation with Parsimony (MRP) inferences.
	"""
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
	
def writeMRP(mrp, nfout=None, outfmt='phylip', codeout={True:'1', False:'0', None:'?'}):
	fout = open(nfout, 'w')
	if outfmt.startswith('phylip'):  fout.write("    %d    %d\n"%(len(mrp), len(mrp.values()[0])))
	for leaf in mrp:
		rowout = [codeout[bi] for bi in mrp[leaf]]
		if outfmt=='phylip': fout.write(leaf.ljust(10))
		elif outfmt=='phylip.relax': fout.write("%s\t"%leaf)
		elif outfmt=='fasta':  fout.write(">%s\n"%leaf)
		else: raise ValueError, "unsuported output format: %s"%outfmt
		fout.write(''.join(rowout)+"\n")
	fout.close()

###################################################
### tree-scoring functions

def unicity(node, lspecies=[], dlabspe={}, useSpeDict=False, **kw):
	"""compute unicity score of species set at leaves under the node.
	
	Unicity score is defined at the node of a genetree as the product of 
	non-zero counts of gene copies from a single species/organism in the subtree.
	cf. Bigot, Lassalle, Daubin and Perriere, 2011. BMC Bioinformatics.
	"""
	if lspecies: lspe = lspecies
	elif dlabspe: lspe = [dlabspe[lab] for lab in node.iter_leaf_labels()]
	elif useSpeDict: lspe = node.listSpecies(**kw)
	else: lspe = node.get_leaf_labels()
	return reduce(lambda x,y: x*y, [lspe.count(spe) for spe in set(lspe)], long(1))
