# tree2
Manipulating phylogenetic trees: topologies, annotations, gene/species tree reconciliations

Description
-----------

*tree2* is a source-package of object-oriented Python modules for dealing with phylogenetic trees. It contains several classes based on the basic class [*Node*](https://github.com/flass/tree2/blob/master/Node.py), that represents a node of a tree. Internal nodes contain their children nodes as attributes, so by reccursion the root node contains all its descendent and is equated to the complete tree. It follows that tree objects are intrisincally rooted, even though the input unrooted trees are supported and represented as multifurcated at the root. Other attributes are attached to each *Node* instance, like branch length, label and comment, so it is sufficient to represent everything that could be found in the [Newick format](http://evolution.genetics.washington.edu/phylip/newicktree.html) plus bracketted comments like this:

((SpeciesA[Phenotype1]:0.45,SpeciesB[Phenotype1]:0.42)0.84:0.75,(SpeciesC[Phenotype1]:0.58,SpeciesD[Phenotype2]:0.85)0.95:0.115, SpeciesE[Phenotype3]:0.50);

It accomodates labelling of internal nodes (that can come handy in exploring species trees), with labels replacing support values:

((SpeciesA[Phenotype1]:0.45,SpeciesB[Phenotype1]:0.42)Clade1:0.75,(SpeciesC[Phenotype1]:0.58,SpeciesD[Phenotype2]:0.85)Clade1:0.115, SpeciesE[Phenotype3]:0.50);

*Node* class include most of the methods for topology manipulation, like pruning, re-rooting, etc. that are inherited by descendant classes. It also includes the basic graphic representation of a tree through an external call to [SEAVIEW](http://doua.prabi.fr/software/seaview) program (exists as a standard Debian package).

Further annotations and handling of more complicated tree formats (NEXUS and phyloXML) is dealt be the [*AnnotatedNode*](https://github.com/flass/tree2/blob/master/AnnotatedNode.py) class. It notably includes extended attributes such as branch color, node ID and taxonomic ID. Export of an *AnnoatedNode* tree instance of information is supported without annotation loss by export in phyloXML format. it follows that it includes smarter graphic representation of a tree through an external call to [Archaeopterix](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx) program.

Based on *AnnoatedNode* are [*GeneTree*](https://github.com/flass/tree2/blob/master/GeneTree.py) and [*ReferenceTree*](https://github.com/flass/tree2/blob/master/ReferenceTree.py) classes, that are mostly enriched in methods for gene/species tree reconciliation procedures. These methods otfen require to involve both types of objects that represent the gene and species trees, respectively, and perform annotations on them given the output of the reconciliations, notably evolutionary events like gene duplication, horizontal gene transfer or gene loss.


Requirements 
------------

This package was developped under LINUX but as it is pure Python (except for certain shell calls to external) it should run OK on other operating systems. It runs on Python version 2.(>=5).
When installing, do not forget to update yout $PYTHONPATH environment variable.

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
t=tree2.Node(fic="tree.newick") # from a file
```

Then you can verify your tree:
```python
t.get_leaf_labels()
print t
print t['Orang']
t.seaview()
```

Now you may want to find the clade that include all great ape species (the node representing their last common acestor) and label it accordingly:
```python
ga = t.map_to_node(['Orang', 'Gorilla', 'Chimp', 'Human'])
ga.edit_label('GreatApes')
print t.newick(ignoreBS=True)
t.seaview(ignoreBS=True)
```
and then you can access the node instance through its label:
```python
print t['GreatApes'].newick(ignoreBS=True)
```

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

Many other features are accessible and editable, but for these please refer to the spefici docuentation of each function and class methods though *help()* function.

Need help?
-----------

If you have any problems or comments about *tree2*, please contact me.
