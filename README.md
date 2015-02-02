# tree2
Manipulating phylogenetic trees: topologies, annotations, gene/species tree reconciliations

Description
-----------

*tree2* is a source-package of object-oriented Python modules for dealing with phylogenetic trees. It contains several classes based on the basic class [*Node*](https://github.com/flass/tree2/blob/master/Node.py), that represents a node of a tree. Internal nodes contain their children nodes as attributes, so by reccursion the root node contains all its descendent and is equated to the complete tree. It follows that tree objects are intrisincally rooted, even though the input unrooted trees are supported and represented as multifurcated at the root. Other attributes are attached to each *Node* instance, like branch length, label and comment, so it is sufficient to represent everything that could be found in the [Newick format](http://evolution.genetics.washington.edu/phylip/newicktree.html) plus bracketted comments like this:

((SpeciesA[Phenotype1]:0.45,SpeciesB[Phenotype1]:0.42)0.84:0.75,(SpeciesC[Phenotype1]:0.58,SpeciesD[Phenotype2]:0.85)0.95:0.115, SpeciesE[Phenotype3]:0.50);

It accomodates labelling of internal nodes (that can come handy in exploring species trees), with labels replacing support values:

((SpeciesA[Phenotype1]:0.45,SpeciesB[Phenotype1]:0.42)Clade1:0.75,(SpeciesC[Phenotype1]:0.58,SpeciesD[Phenotype2]:0.85)Clade2:0.115, SpeciesE[Phenotype3]:0.50);

*Node* class methods are inherited by descendant classes and include most of the methods for topology manipulation, like pruning, re-rooting, etc. It also provides a method to vizualize instantly (from within the Python interpretter) a tree object using [SEAVIEW](http://doua.prabi.fr/software/seaview)'s basic graphic engine, an option that comes quite handy during code development to evaluate the properties of manipulated trees.

Further annotations and handling of more complicated tree formats (NEXUS and phyloXML) is dealt be the [*AnnotatedNode*](https://github.com/flass/tree2/blob/master/AnnotatedNode.py) class. It notably includes extended attributes such as branch color, node ID and taxonomic ID. Export of an *AnnoatedNode* tree instance of information is supported without annotation loss by export in phyloXML format. This allows smarter graphic representation of a tree through an external call to [Figtree](http://tree.bio.ed.ac.uk/software/figtree/) or [Archaeopteryx](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx) programs.

Based on *AnnoatedNode* are [*GeneTree*](https://github.com/flass/tree2/blob/master/GeneTree.py) and [*ReferenceTree*](https://github.com/flass/tree2/blob/master/ReferenceTree.py) classes, that are mostly enriched in methods for gene/species tree reconciliation procedures. These methods otfen require to involve both types of objects that represent the gene and species trees, respectively, and perform annotations on them given the output of the reconciliations, notably evolutionary events like gene duplication, horizontal gene transfer or gene loss.


Requirements 
------------

This package was developped under LINUX but as it is pure Python (except for certain shell calls to external) it should run OK on other operating systems. It runs on Python version 2.(>=5).
When installing, do not forget to update yout $PYTHONPATH environment variable.
Optional programs can be installed for use through the API proveded by *tree2*:
- **PhyML**, phylogenetic reconstruction - used for re-estimating branch length and supports after topological manipulations (download [here](http://www.atgc-montpellier.fr/phyml/binaries.php) or available as a standard Debian package)
- **SeaView**, graphic platform for sequence alignment, tree reconstruction and vizualisation - used here only fore tree vizualisation (download [here](http://doua.prabi.fr/software/seaview) or available as a standard Debian package)
- **Figtree**,  graphic representation of annotated phylogenetic trees, using NEXUS format (download [here](http://tree.bio.ed.ac.uk/software/figtree/) or available as a standard Debian package)
- **Archaeopteryx**, graphic representation of richly annotated phylogenetic trees, using phyloXML format (download [here](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx))
- **Count**, parsimonious and ML estimation of gene gain/loss scenrios along a species tree given a phylogenetic profile (download [here](http://www.iro.umontreal.ca/~csuros/gene_content/count.html))

Credit
------

*tree2* package is derived from the *tree* module from [*alfacinha* package](https://github.com/leonorpalmeira/alfacinha) by Leonor Palmeira and Laurent Guéguen (doc here: http://pbil.univ-lyon1.fr/software/alfacinha/)


Usage
-----

To get started, you have to create a tree from a Newick string:

```python
import tree2
t=tree2.Node(newick="(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927)0.89:0.08386)0.94:0.06124)0.94:0.15057)0.90:0.54939,Mouse:1.21460)0.86:0.10;") # typped in
```
or
```python
import tree
t=tree2.Node(fic="/path/to/tree.newick") # from a file
```

Then you can print out your tree to verify it:
```python
print t			# standard Newick string representation
t.arborescence_ASCII()	# hierarchical arborescence representation in text mode
t.seaview()		# graphic representation
```

You can now access its several attributes, globally or for specific nodes
```python
t.get_leaf_labels()	# list of leaf labels
o = t['Orang']		# access a node through its indexed label (must be unique)
print o.label()		# node label
print o.lg()		# length of the branch leading to the node (above the node)
print o.bs()		# support of the branch
```

You may want to find the clade that include all great ape species (the node representing their last common acestor) and label it accordingly:
```python
ga = t.map_to_node(['Orang', 'Gorilla', 'Chimp', 'Human'])
ga.edit_label('GreatApes')
print t.newick(ignoreBS=True)	# ignore branch supports to display internal node labels
t.seaview(ignoreBS=True)
```
and then you can access the node instance through its label:
```python
print t['GreatApes'].newick(ignoreBS=True)
```

You can iterate on the nodes of the tree:
```python
# using iterator
for n in t:
  print n

You may want to remove some species from your dataset while keeping the paroperties of the rest of the tree (notably consistent branch lengths and supports)

```python
g = t.pop('Gorilla') # prune the Gorilla branch identified by its label
print t
print g
ho = t.map_to_node(['Gorilla', 'Chimp', 'Human']) # identify the node of the Homininae clade
h = t.pop(ho) # prune the Homininae branch identified by its obect reference
print t
print h
```

Many other features are accessible and editable, but for these please refer to the specific documentation of each function and class methods through *help()* function.

Need help?
-----------

If you have any problems or comments about *tree2*, please contact me.
