#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Models for simulation of phylogenetic tree(s) forward in time, specifically recording events of duplication, transfer and loss for gene trees relative to a species tree."""

__author__ = "Florent Lassalle <florent.lassalle@imperial.ac.uk>"
__date__ = "27 July 2016"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""


### not functional ; still under velopment

class Genome(object):
	""""""
	def __init__(self, dreplicons):
		self.dreplicons = dreplicons
		
	def get_all_genes(self):
		lgenes = []
		for r in self.dreplicons.values():
			lgenes += r.lgenes
		return lgenes
		
class Replicon(object):
	""""""
	replicontypes = ['chromosome', 'plasmid', 'chromid']
	
	def __init__(self, lgenes, circular=True, **kwargs):
		self.lgenes = lgenes
		self.circular = circular
		self.name = kwargs.get('name')
		self.type = kwargs.get('type', "chromosome")
		
	def insert_gene(self, gene, locus=None):
		"""add gene to the ordered list of genes; by defautls at its end, or if specified at its index 'locus'"""
		if not locus is None: self.lgenes.insert(locus, gene)
		else: self.lgenes.append(gene)
		gene.attach_to_replicon(self)
		
	def insert_gene(self, gene, locus=None):
		"""remove gene from the ordered list of genes; by defautls use the object as its own identifier, but instead use its 'locus' index if specified"""
		if locus is None:
			g = gene
			self.lgenes.remove(g)
		else:
			g = self.lgenes[locus]
		g.attach_to_replicon(None)
	
class Gene(object):
	""""""
	genetypes = ['core', 'heterogeneous', 'ORFan']
	
	def __init__(self, treenode, **kwargs):
		self.treenode = treenode
		self.name = kwargs.get('name')
		self.type = kwargs.get('type', "core")
		self.replicon = None
		
	def attach_to_replicon(self, repli):
		self.replicon = repli
		#replicon.insert_gene(self, locus=locus)
