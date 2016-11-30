#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Errors related to phylogenetic tree manipulation"""

__author__ = "Florent Lassalle <florent.lassalle@imperial.ac.uk>"
__date__ = "27 July 2016"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""

import tree2

class TreeReferenceError(LookupError):
	def __init__(self, node, tree, verbose=True):
		self.node = node
		self.tree = tree
		self.verbose = verbose
	def __str__(self):
		if self.verbose: return "node %s:\n%s\nnot in tree %s:\n%s"%(repr(self.node), str(self.node), repr(self.tree), str(self.tree))
		else: return "node %s not in tree %s"%(repr(self.node), repr(self.tree))
		
class AggregateTreeReferenceError(TreeReferenceError):
	"""allow to report several reference error at once, e.g. after systematic exploration of possible references"""
	def __init__(self, nodes, trees, verbose=True):
		self.nodes = nodes
		self.trees = trees if (trees is isinstance(tree, list)) else [trees]
		self.verbose = verbose
	def __str__(self):
		if self.verbose: return "nodes:\n%s\nnot in tree %s:\n%s"%("\n".join(["%s:\n%s"%(repr(node), str(node)) for node in self.nodes]), \
		 "\n".join(["%s:\n%s"%(repr(tree), str(tree)) for tree in self.trees]))
		else: return "nodes [%s] not in any of those trees %s"%(", ".join([repr(node) for node in self.nodes]), ", ".join([repr(tree) for tree in self.trees]))

class DuplicateNodeLabel(Warning):
	def __init__(self, label, nodes, verbose=True):
		self.label = label
		self.nodes = nodes
		self.verbose = verbose
	def __str__(self):
		if self.verbose: return "nodes:\n%s\nshare the same label '%s'"%("\n".join(["%s:\n%s"%(repr(node), str(node)) for node in self.nodes]), self.label)
		else: return "nodes [%s] share the same label '%s'"%(", ".join([repr(node) for node in self.nodes]), self.label)
		
