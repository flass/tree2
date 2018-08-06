#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Enriched phylogenetic tree objects, can be exported to densely annotated formats like PhyloXML and Nexus."""

__author__ = "Florent Lassalle <florent.lassalle@ucl.ac.uk>"
__date__ = "14 January 2015"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""

import copy
import subprocess

import tree2
from os import environ
homedir = environ['HOME']

class AnnotatedNode(tree2.Node):
	"""	"""	
	def __init__(self, branch_lengths=True, keep_comments=False, addAttr=[], **kw):
		""" """		
		super(AnnotatedNode, self).__init__(branch_lengths=branch_lengths, keep_comments=keep_comments, **kw)		
		self.__color=[] 			# list of RGB component for coloring (ex: [0, 255, 78])		
		self.__nodeid = -1			# for Phylariane/Ancestrome trees.
		self.__taxid = None	
		for attr in set(addAttr):
			# create extra public attribute slot, with value None by default
			setattr(self, attr, None)
		if kw.has_key('xml'):
			ind='  '
			indstart='\n'
			if kw.has_key('ind'): ind=kw['ind']
			if kw.has_key('indstart'): indstart=kw['indstart']
			self.parse_phyloXML(s=kw['xml'], branch_lengths=branch_lengths, ind=ind, indstart=indstart)
			
	def newnode(self, branch_lengths=True, keep_comments=False, **kw):
		"""class-specific instance generator for construction of trees of AnnotatedNodes"""
		an = AnnotatedNode(branch_lengths=branch_lengths, keep_comments=keep_comments, **kw)
		for attr in self.__dict__:
			if not attr.startswith('_'):
				# copy the extra public attribute slot, with value None by default
				setattr(an, attr, None)
		return an
	
	def color(self):
		"""Return the color of the Node."""
		return self.__color
	
	def hexrgbcolor(self):
		"""Return the RGB string representation (as in SVG standard) of the color of the Node."""
		return '#'+''.join(["%02x"%val for val in self.__color])
			
	def edit_color(self, col):
		"""edits the color attached to the node"""
		self.__color = col
		
	def taxid(self):
		return self.__taxid
		
	def set_taxid(self, taxid):
		self.__taxid = taxid
		
	def nodeid(self):
		return self.__nodeid
			
	def set_node_id(self, nodeid):
		self.__nodeid = nodeid
		
	def idgetnode(self, nodeid):
		for node in self:
			if node.nodeid() == nodeid:
				return node
				
	def idgetlab(self, nodeid):
		for node in self:
			if node.nodeid() == nodeid:
				return node.label()
				
	def labgetid(self, lab):
		return self[lab].nodeid()
	
	def children_nodeids(self):
		"""Return the list of direct child node ids."""
		llc = []
		for c in self.__children:
			llc.append(c.nodeid())
		return llc
		
	def father_nodeid(self):
		f = self.father
		if f: return f.nodeid()
		else: return None
				
	def complete_node_ids(self, order=1, force=False):
		"""numbering of nodes; default by increasing depth = from root to leaves, but not following the tre structure"""
		nodeid=0
		if isinstance(order, list):
			# orderred node list is given
			children = order
		else:
			children = self.get_sorted_children(order=order)
		for node in children:
			if force or (not node.nodeid()>=0):
				node.set_node_id(nodeid)
				nodeid += 1
			
	def get_leftright_index(self, i=-1, d=None):
		if i<0: self.complete_node_ids()
		if not d: dnodeid_leftright = {}
		else: dnodeid_leftright = d
		i += 1
		ileft = i
		for child in self.get_children():
			j, dc = child.get_leftright_index(i=i, d=dnodeid_leftright)
			i = j
			dnodeid_leftright.update(dc)
		i += 1
		iright = i
		dnodeid_leftright[self.nodeid()] = (ileft, iright)
		return (i,  dnodeid_leftright)
			
###########################
#### File input methods
					
	def parse_phyloXML(self, s, branch_lengths=True, ind="  ", indstart="\n"):		
		"""Parses an indented multi-line PhyloXML string and fills the Node with the described tree."""
		l = s.split(indstart)
		l[0]=l[0].lstrip(ind)
		if l[0]!='<clade>' or l[-1]!='</clade>':
			raise SyntaxError, "Invalid phyloXML clade delimitation syntax at:\n%s"%indstart.join(l)

		n = 0
		while n < len(l):
			nind = l[n].count(ind)
			if nind == 0:
				if l[n]=='<clade>' or l[n]=='</clade>':
					n += 1
				else:
					u = max([0, n-3])
					v = min(len(l)-1, n+3)
					raise SyntaxError, "Invalid phyloXML clade delimitation syntax at:\n%s"%indstart.join(l[u:v])
			elif nind == 1:
				if l[n]=='%s<clade>'%ind:
					cladestart = n
					while not l[n]=='%s</clade>'%ind:
						n += 1
					n += 1
					clade = indstart.join(l[cladestart:n])
					node = AnnotatedNode(xml=clade, branch_lengths=branch_lengths, ind=ind, indstart=indstart+ind)
					self.link_child(node, newlen=node.lg())					
				else:
					name = extract_XMLmarker(l[n], 'name', ind=ind)
					if name:
						self.add_label(name)
					else:
						length = extract_XMLmarker(l[n], 'branch_length', ind=ind)
						if length:
							if branch_lengths:
								self.set_lg(length)
						else:
							bs = extract_XMLmarker(l[n], 'confidence', ind=ind, mtype='bootstrap')
							if bs:
								self.set_bs(bs)
					n += 1
			else:
				# misses deeper indentations (like within <color> markers) if not surounded by <clade> markers
				n += 1
			
###########################
### screen output methods
			
	def nexus(self, treename="", figtree=True, ignoreBS=False, branch_lengths=True, unrooted=False, mini=False, **kw):
		"""Return NEXUS string of the sub-tree defined by the Node."""
		def hexstr(rgb):
			hexa='0123456789abcdef'
			s = ''
			for col in rgb:
				s += hexa[int(col/16)]+hexa[int(col%16)]
			return s
		outtree = copy.deepcopy(self)
		if treename =="":
			treename = 'tree_1'
		if figtree==True:
			for node in outtree.get_all_children():
				if node.bs() or node.color() or isinstance(node, tree2.GeneTree):
					node.edit_comment('&')
					if node.bs():
						node.edit_comment('bootstrap = %f'%node.bs(), mode='append' )
					if node.color():
						node.edit_comment('!color = #%s'%hexstr(node.color()), mode='append' )
					if isinstance(node, tree2.GeneTree):
						node.edit_comment('!events = %s'%str(node.getEvents()), mode='append' )
		nex0 = '[&R] '+outtree.newick(ignoreBS=ignoreBS, branch_lengths=branch_lengths, unrooted=unrooted)
		if mini: return nex0
		nex1 = '#NEXUS\nbegin taxa;\n\tdimensions ntax=%d;\n\ttaxlabels\n\t%s\n;\nend;\n'%( self.nb_leaves(), '\n\t'.join(self.get_leaf_labels(comments=True)) )
		nex2 = '\nbegin trees;\n\ttree %s = %s\nend;\n'%(treename, nex0)
		nex3 = tree2.figtreeblock if figtree else ''
		return nex1+nex2+nex3

	def phyloXML(self, ind="  ", indstart="\n", treetype=None):
		"""Returns PhyloXML string of the sub-tree defined by the Node.
		
		"events" specifiy if nodes have annotated events
		"""
		xml = "%s<clade>"%indstart
		if self.label()!="":
			xml += "%s<name>%s</name>"%(indstart+ind, self.label())
		if self.lg != None:
			xml += "%s<branch_length>%f</branch_length>"%(indstart+ind, self.lg())
		if self.bs() != None:
			xml += "%s<confidence type='bootstrap'>%s</confidence>"%(indstart+ind, self.bs())
		if self.color():
			col = ['red', 'green', 'blue']
			xml += "%s<color>"%(indstart+ind)
			for i in range(3):
				xml += "%s<%s>%d</%s>"%(indstart+ind*2, col[i], self.color()[i], col[i])
			xml += "%s</color>"%(indstart+ind)
		if self.nodeid() >= 0:
			xml += "%s<node_id>%d</node_id>"%(indstart+ind, self.nodeid())					
		if self.get_children():
			for c in self.get_children():
				xml += c.phyloXML(ind, indstart+ind, treetype=treetype)
		if self.taxid():
			xml += "%s<taxonomy provider='NCBI'>"%(indstart+ind)
			xml += "%s<id>%s</id>"%(indstart+ind*2, str(self.taxid()))
			xml += "%s</taxonomy>"%(indstart+ind)
		if treetype:
			if treetype=='scenario':
				# for gene family evolutionary scenarios
				xml += self.scenarioEventXML(ind, indstart)
			elif treetype=='reconciliation':
				# for Phylariane/Ancestrome reconcilated gene trees. !!! NB: the node annotations for transfer concern the children branches/nodes
				xml += self.reconciliationEventXML(ind, indstart)
			else:
				raise ValueError, "Non handled tree type: %s"%treetype
		xml += "%s</clade>"%indstart
		return xml
		
	def phyloXML_norm(self, ind="  ", indstart="\n", params='default'):
		"""Returns PhyloXML string of the sub-tree defined by the Node, following Phylariane/Ancestrome/Agrogenom norma.
		
		"params" are 2-tuples of strings expliciting the parameters of the reconciliation program.
		"""
		self.complete_node_ids()
		xml = ''
		if not isinstance(params, list):
			# Prunier standard
			params = [('fwd.depth','2'), ('boot.thresh.conflict','0.9'), ('max.bp', '1.0'), ('multi_root.bool','false'), ('prog', 'Prunier'), ('prog-version', 'v.2.0'), ('prog-arguments', '....')]
		for param in params:
			fm = (indstart,)+param
			xml += '%s<rec:param name="%s" value="%s"/>'%fm
		if isinstance(self, tree2.GeneTree):
			xml += self.phyloXML(ind=ind, indstart=indstart, treetype='reconciliation')
		if isinstance(self, tree2.ReferenceTree):
			xml += self.phyloXML(ind=ind, indstart=indstart, treetype='scenario')
		return xml

################################
### File output methods
		
	def write_phyloXML(self, nfout, treename="", ind="  ", indstart="\n", normparams=None):
		""" writes the tree in PhyloXML format, readable by Achaeopteryx (http://www.phylosoft.org/archaeopteryx/)"""
		header =  "<?xml version='1.0' encoding='UTF-8'?>"
		header += "%s<phyloxml xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' xsi:schemaLocation='http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd' xmlns='http://www.phyloxml.org'"%(indstart)
		# addition of rec module specification to header is non compliant with late version of achaeopteryx:
		# if isinstance(self, tree2.GeneTree) or isinstance(self, tree2.ReferenceTree): header += " xmlns:rec='http://www.phyloxml_rec.org' rec:schemaLocation='http://www.phyloxml_rec.org http://www.phyloxml_rec.org/1.10/phyloxml.xsd'>"
		# else:  header += ">"
		header += ">"
		header += "%s<phylogeny rooted='true'>"%indstart
		header += "%s<name>%s</name>"%(indstart+ind,treename)
		fout = open(nfout, 'w')
		fout.write( header )
		if normparams is not None:
			# Phylariane website/database-readable version of phyloXML
			fout.write(self.phyloXML_norm(ind, indstart+ind, normparams))
		else:
			# standard/Archaeopteryx-readable version of phyloXML
			fout.write(self.phyloXML(ind, indstart+ind))
		fout.write("%s</phylogeny>"%indstart)
		fout.write("%s</phyloxml>"%indstart)
		fout.close()
		
	def write_nexus(self, nfout, treename="", figtree=True, ignoreBS=False, branch_lengths=True, unrooted=False, mode='w'):
		""" writes the tree in NEXUS format, readable by FigTree (http://tree.bio.ed.ac.uk/software/figtree/)""" 
		fouttree = open(nfout, mode)
		fouttree.write(self.nexus(treename=treename, figtree=figtree, ignoreBS=ignoreBS, branch_lengths=branch_lengths, unrooted=unrooted) )
		fouttree.close()	
		
	
####################################
### External program call functions

	def archaeopteryx(self, treename="", arch_exec_path='archaeopteryx', tmpfiledir='%s/tmp'%homedir, bg=True):
		ntmpfile = "%s/tmptree2archaeopteryx.xml"%tmpfiledir
		self.write_phyloXML(ntmpfile, treename=treename)
		callline = "%s %s"%(arch_exec_path, ntmpfile)
		if bg: callline += ' &'
		print callline
		subprocess.call(callline, shell=True)	
			
	def figtree(self, figtree_exec_path='figtree', tmpfiledir='%s/tmp'%homedir, ignoreBS=False, branch_lengths=True, unrooted=False, bg=True):
		ntmpfile = "%s/tmptree2figtree.phb"%tmpfiledir
		self.write_nexus(ntmpfile, ignoreBS=ignoreBS, branch_lengths=branch_lengths, unrooted=unrooted)
		callline = "%s %s"%(figtree_exec_path, ntmpfile)
		if bg: callline += ' &'
		print callline
		subprocess.call(callline, shell=True)			

	
