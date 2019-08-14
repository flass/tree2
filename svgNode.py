#!/usr/bin/python
# -*- coding: utf-8 -*-

"""module providing functions for SVG representation of phylogenetic trees from tree2.Node module

Copyright 2016, Florent Lassalle  <f.lassalle@imperial.ac.uk>.
"""

__author__ = "Florent Lassalle  <f.lassalle@imperial.ac.uk>"
__date__ = "14 January 2015"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""

from random import uniform

# variable to be uptdated during the calls
svgmarker = {'arrow':[6, 9]}

eventcolors = {'D':'#dd0000', 'T':'#33ee00', 'L':'#0000ee', '#':'#dd8800'}
for i in range(10): eventcolors[str(i)] = eventcolors['#']

# generic SVG formating functions
def groomStyle(style):
	lstel = style.split(';')
	dqual = {}
	lqual = ['stroke', 'stroke-width'] # enforce presence of those qualif
	for stel in lstel:
		if not stel.strip(' '): continue
		qualif, val = [s.strip(' ') for s in stel.split(':')]
		dqual[qualif] = val
		if not qualif in lqual: lqual.append(qualif) # keep appearance order
	return '; '.join(["%s:%s"%(qualif, dqual[qualif]) for qualif in lqual])+'; '

def hexrgb(rgbtup):
	return '#'+''.join(["%02x"%val for val in rgbtup])

def randRGBColGen(ranges=None):
	colranges = {'R':(.3, .9), 'G':(.3, .9), 'B':(.3, .9)}
	colranges.update(ranges)
	rgb = [int(255*uniform(a, b)) for a,b in [colranges[col] for col in 'RGB']]
	return hexrgb(rgb)

def svgAddArrowMarkerDef(name='arrow', widthheight=None, color='red'):
	# set class variables for svgTransfer() calls
	if widthheight: svgmarker[name] = widthheight
	path = "M0,0 L0,%d L%d,%d z"%(svgmarker[name][0], svgmarker[name][1], float(svgmarker[name][0])/2)
	s =  '<defs>\n'
	s += '<marker id="%s" markerWidth="10" markerHeight="10" refX="0" refY="3" orient="auto" markerUnits="strokeWidth">\n'%name
	s += '<path d="%s" fill="%s" />\n'%(path, color)
	s += '</marker>\n'
	s += '</defs>\n'
	return s
	
# functions to draw the tree
def getDictNodeCoords(tree, **kw):
	"""generates coordinates of nodes in the tree for further drawing
	
	**kw (optional) arguments:
	'interleaves' define spacing (y axis, in pixels) between leaves
	'phylofact'   define the scaling factor (y axis) for the whole tree (from the tree distance metric into pixels)
	'cladofact'   define the scaling factor (x axis) for the whole tree (from the tree 'interleaves' parameter)
	'nodeoverlap' specify if the function must try to print branches in a non-overlaping fashion (default, True)
	'offset'      tupple of x,y vector by which the position of all points is offset
	'treeorder'   order of tree traversal defining the representation (see tree2.Node.get_sorted_children) 
	              NB: noderowoverlap=True is only compatible with treeorder={2,3}
	"""
	cladogram = kw.get('cladogram', False)
	interleaves = kw.get('interleaves', 20)
	cladofact = kw.get('cladofact', 5)
	phylofact = kw.get('phylofact', 200)
	noderowoverlap = kw.get('noderowoverlap', True)
	treeorder = kw.get('treeorder', 2)
	offset = kw.get('offset')
	if noderowoverlap and (treeorder not in [2,3]):
		raise ValueError, "treeorder=%d; noderowoverlap=True is only compatible with treeorder={2,3}"%treeorder
	dxy = {}
	yl = 0
	internodes = interleaves*cladofact
	nodes = tree.get_sorted_children(order=treeorder)
	# fill Y-axis coordinates
	for node in nodes:
		if node.is_leaf() or (not noderowoverlap):
			h = yl
			yl += interleaves
		else:
			# get mean height of children
			h = int(float(sum([dxy[c.label()][1] for c in node.get_children()]))/node.nb_children())
		dxy[node.label()] = (None, h)
	# fill X-axis coordinates
	for node in tree.get_sorted_children(order=2):
		# has to be done in post-order traversal; does not impact the leaf order in tree layout
		if cladogram:
			if node.is_leaf():
				x = 0
			else:
				# get max depth of children
				x = max([dxy[c.label()][0] for c in node.get_children()])+internodes
		else:
			x = node.distance(tree)*phylofact
		dxy[node.label()] = (x, dxy[node.label()][1])
	if cladogram:
		# correct the depth of nodes to get the root at x=0
		totdepth = dxy[tree.label()][0]
		for node in tree.get_sorted_children(order=0):
			x, y = dxy[node.label()]
			dxy[node.label()] = (totdepth - x, y)
	if offset:
		for nodelab in dxy:
			dxy[nodelab] = (dxy[nodelab][0]+offset[0], dxy[nodelab][1]+offset[1])
	return dxy

def svgPathAB(a, b, dxy=None, squared=True, bezier=False, style=None, marker=None, aorbit=(0,0), borbit=(0,0), half=None):
	"""writes a line from point a to point b; these can be given directly as coordinate tupples, or as key to a coordinate dictionary dxy."""
	if not style: st = "stroke:black; stroke-width:1; fill:none; "
	else: st = style
	if marker:
		if isinstance(marker, str): st += 'marker-end:url(#%s) '%marker
		else: st += 'marker-end:url(#arrow) '
	if dxy:
		ax, ay = dxy[a]
		bx, by = dxy[b]
	else:
		ax, ay = a
		bx, by = b
	pathdata = "M %d %d"%(ax, ay)
	correctend = svgmarker.get(marker, [6, 9])[1] if marker else 0
	if bezier:
		# define control point
		concavex = 1 if bezier > 0 else -1	
		# distances
		dx = bx - ax
		dy = by - ay
		# midpoint with offset
		cx = ax + dx/2 - dy/5
		cy = ay + dy/2 - concavex*dx/5
		pathdata += " Q %d %d %d %d"%(cx, cy, bx, by - correctend)
	elif squared:
		pathdata += " L %d %d %d %d"%(bx, ay, bx, by - correctend)
	else:
		pathdata += " L %d %d"%(bx, by - correctend)   ##### correction for marker coords should account for angle
	return '<path style="%s" d="%s" />\n'%(st, pathdata)

def svgPathToFather(node, dxy, squared=True, style=None, half=None):
	"""path data for SVG representatioon of the branch from the node to its father"""
	f = node.go_father()
	if f:
		return svgPathAB(node.label(), f.label(), dxy=dxy, squared=squared, style=style, half=half)
	else:
		return ''

def svgTransfer(don, rec, dxy, tree, donbrheight=.5, recbrheight=.5, style=None, marker='arrow', orbit=(0,0), arc=True, **kw):
	color = kw.get('color', 'black')
	width = float(kw.get('width', 2.0))
	if not style: st = "stroke:%s; stroke-width:%d; fill:none; "%(color, width)
	else: st = style
	donx, dony = dxy[don]
	recx, recy = dxy[rec]
	recfatx, recfaty = dxy[tree[rec].go_father().label()]
	donfatx, donfaty = dxy[tree[don].go_father().label()]
	depxy = (donfatx+(donx-donfatx)*donbrheight, dony) # startpoint of transfer
	arrxy = (recfatx+(recx-recfatx)*recbrheight, recy) # endpoint of transfer
	xml  = svgPathAB(depxy, arrxy, style=st, marker=marker, bezier=((arrxy[0]-depxy[0]) if arc else False), borbit=orbit)
	xml += svgEventLab('T', x=arrxy[0], y=arrxy[1], color=color, inCircle=True, orbit=orbit)
	return xml
		
def svgPathPadded(node, dxy, leaves=False, xend=None, xshift=None, style="stroke:black; stroke-width:1; stroke-dasharray:5,5; ", **kw):
	if leaves or (not node.is_leaf()):
		if not xshift: xs = (55 if leaves else 35)
		else: xs = xshift
		if not xend:
			xe, ye = svgPlotSize(dxy)
			xe += (100 if leaves else 50)
		else:
			xe = xend
		sxy = dxy[node.label()]
		pathdata = "M %d %d"%(sxy[0]+xs, sxy[1])
		pathdata += " L %d %d"%(xe, sxy[1])
		return '<path style="%s" d="%s" />\n'%(style, pathdata)
	else:
		return ''

def svgPlotSize(dxy):
	xend = max([xy[0] for xy in dxy.values()])
	yend = max([xy[1] for xy in dxy.values()])
	return xend, yend

def svgEventLab(text, nodelab=None, dxy=None, x=None, y=None, orbit=(0, 0), style=None, fontsize=10, color=None, inCircle=False, circlecolor=None):
	"""add a label at a particular coordinate; particular handling of DTL event-type labels for colouring"""
	xlab, ylab = (x, y)
	if not xlab: xlab = dxy[nodelab][0]
	if not ylab: ylab = dxy[nodelab][1]
	xml = ''
	col = color if color else 'black'
	# draw circular frame
	if inCircle:
		circlecol = circlecolor if circlecolor else eventcolors.get(text[0], ('black' if (col not in ['black', '#000000']) else 'white'))
		rad = fontsize*.5+(fontsize*len(text)*.1)
		#~ rad = fontsize*len(text)*.3
		xml += '<circle cx="%d" cy="%d" r="%d" fill="%s" stroke="%s"/>'%(xlab+orbit[0], ylab+orbit[1], rad, circlecol, col)
	# draw text label
	if not style: st = 'fill="%s" font-size="%dpx" '%(col, fontsize)
	xml += '<text x="%d" y="%d" %s >%s</text>\n'%(xlab+orbit[0]-(fontsize*len(text)*.3), ylab+orbit[1]+(fontsize*.4), st, text)
	return xml

def svgTree(tree, labels=True, supports=True, comment=None, fontsize=10, textorbit=5, lgscale=0.05, padleaves=False, padinternalnodes=False, svgwrap=True, dnodecoord=None, **kw):
	"""XML string describing a SVG representatioon of the tree
	
	In style descriptors, only the last occurence of a qualifier is considered by the groomStyle() parser ;
	hence 'defaultstyle' is overridden by 'modstyle', which is overridden by 'padstyle'.
	
	Keyword arguments 'duplications', 'transfers', and 'losses' allow to add tree decorators indicating the corresponding events;
	they can be passed as lists of 2-tuples containing the node name and the frequency of the event for duplications and losses, and as 3-tuples for transfers (donor, receiver, freq).
	Keyword argument 'counts' allow to add tree decorators indicating the copy number at the relevant tree node.
	
	Other **kw (optional) arguments passed on to svgNode.getDictNodeCoords:
	
	general graphic options:
	'interleaves' define spacing (y axis, in pixels) between leaves
	'phylofact'   define the scaling factor (y axis) for the whole tree (from the tree distance metric into pixels)
	'cladofact'   define the scaling factor (x axis) for the whole tree (from the tree 'interleaves' parameter)
	'nodeoverlap' specify if the function must try to print branches in a non-overlaping fashion (default, True)
	'offset'      tupple of x,y vector by which the position of all points is offset
	'treeorder'   order of tree traversal defining the representation (see tree2.Node.get_sorted_children) 
	              NB: noderowoverlap=True is only compatible with treeorder={2,3}
	              
	 HGT-specific graphic options:
	'transfercolor'      colour of transfer edges 
	'transferpathtype'   'arc' or 'line' (default)
	'transferwidth'      width of transfer edges
	"""
	brwidthfac = 10
	
	xml = ''
	dnodestyle = kw.get('dnodestyle', {})
	defaultstyle = kw.get('defaultstyle', "stroke:black; stroke-width:1; fill:none; ")
	modstyle = kw.get('modstyle', "")
	padstyle = kw.get('padstyle', "stroke-dasharray: 5, 5; ")
	treeorder = kw.get('treeorder', 2)
	bw = kw.get('branchwidths')
	dbranchwidth = bw if isinstance(bw, dict) else None
	branchwidthattr = bw if isinstance(bw, str) else None
	if 'transfers' in kw:
		# prepare for tracing of transfer arrows
		xml += svgAddArrowMarkerDef(name='arrow')
		transfercol = kw.get('transfercolor', (0,0,0))
		transfercol = transfercol if isinstance(transfercol, str) else hexrgb(transfercol)
		transpathtype = kw.get('transferpathtype', 'arc')
		transferwidth = kw.get('transferwidth')
	if dnodecoord: dxy = dnodecoord
	else: dxy = getDictNodeCoords(tree, offset=(100, 100), **kw)
	xend, yend = svgPlotSize(dxy)
	### draw tree
	# (optionally) first a tree with branch thickness proprtional of presence probability
	if branchwidthattr or dbranchwidth:
		for node in tree.get_sorted_children(order=treeorder):
			if branchwidthattr: nodebw = getattr(node, branchwidthattr, 0)
			elif dbranchwidth: nodebw = dbranchwidth.get(node.label(), 0)
			if nodebw:
				# only draw branches with some thichness
				nodecol = getattr(node, '_AnnotatedNode__color', 'grey')
				nodecols = nodecol if isinstance(nodecol, str) else hexrgb(nodecol)
				widstyle = "stroke:%s; stroke-width:%f; fill:none; stroke-linejoin:round; stroke-linecap:round; "%(nodecols, nodebw*brwidthfac)
				xml += svgPathToFather(node, dxy, style=widstyle)
	# (then) a tree with simple solid branches
	for node in tree.get_sorted_children(order=treeorder):
		style = groomStyle(dnodestyle.get(node, defaultstyle+modstyle))
		pstyle = groomStyle(style+padstyle)
		xml += svgPathToFather(node, dxy, style=style)
		x, y = dxy[node.label()]
		# draw padding
		if node.is_leaf():
			if padleaves:  xml += svgPathPadded(node, dxy, style=pstyle, leaves=padleaves, **kw)
		else:
			if padinternalnodes: xml += svgPathPadded(node, dxy, style=pstyle, **kw)
		if labels:
			if padleaves: xlab = xend + 100
			else: xlab = x
			if node.is_leaf():	
				# write on the right side of the node
				xml += '<text x="%d" y="%d"  font-size="%dpx" >%s'%(xlab+abs(textorbit), y+(fontsize*.4), fontsize, node.label())
			else:
				# write on the left-top side of the node
				xml += '<text x="%d" y="%d" font-size="%dpx" text-anchor="end" >%s'%(x-textorbit, y+(fontsize*.4)+(textorbit if textorbit>0 else 0), fontsize, node.label())
			if comment: xml += ' [%s]'%node.commentAsString(comment=comment)
			xml += '</text>\n'
		if supports and node.bs():
			if padinternalnodes:
				# avoid printing over the padding line
				xml += '<text x="%d" y="%d"  font-size="%dpx" >%.2f</text>\n'%(x+textorbit*2, y+(fontsize*.4)-textorbit, fontsize, node.bs())
			else:
				xml += '<text x="%d" y="%d"  font-size="%dpx" >%.2f</text>\n'%(x+textorbit*2, y+(fontsize*.4), fontsize, node.bs())
				
	if not kw.get('cladogram'):			
		# get the coordinates for bottom-middle placement of scale bar
		xr = dxy[tree.label()][0] # root 
		xs = (xend + xr)*.5
		ys = yend + 200
		# get the distance scale
		n = tree.get_children()[0]
		xn = dxy[n.label()][0]
		scale = (xn - xr)*(float(lgscale)/n.lg())
		# draw scale bar
		xml += '<path style="%s" d="M %d %d L %d %d" />\n'%(defaultstyle, xs-(scale*.5), ys, xs+(scale*.5), ys)
		xml += '<text x="%d" y="%d"  >%s</text>\n'%(xs, ys+50, str(lgscale))
	
	# plot DTL events
	if kw.get('treetype')=='species':
		# count number of events per node; allow to increment height of event logo representation on the branch (orbit on x axis) so they do not overlap.
		# not ideal as events are arbitrarily ordred instead of time-ordered.
		# instead of this increment, the height can be set from an additional member (positioned last) of the event descriptor tupple 
		nnodeevts = {} 
		for te in kw.get('transfers', []):
			print te
			don, rec, freq = te[:3]
			if len(te)>3: hev = te[4]
			else: hev = nnodeevts.setdefault(rec, 0)*5
			transcol = transfercol if transfercol else randRGBColGen({'G':(0.0, 0.6)})
			if not transferwidth:
				transwid = 2.0*brwidthfac
			elif transferwidth=='freq':
				transwid = freq*brwidthfac
			else:
				transwid = transferwidth*brwidthfac
			xml += svgTransfer(don, rec, dxy, tree, marker=None, color=transcol, width=transwid, orbit=(hev, 0), arc=(transpathtype=='arc'))
			nnodeevts[rec] += 1
			#~ xml += svgEventLab('%.2f'%freq, x=dxy[rec][0], y=dxy[rec][1], color=tcol)
		for de in kw.get('duplications', []):
			nodelab, freq = de[:2]
			if len(de)>2: hev = de[3]
			else: hev = nnodeevts.setdefault(nodelab, 0)*5
			xml += svgEventLab('D', nodelab=nodelab, dxy=dxy, inCircle=True, orbit=(hev, 0))
			#~ xml += svgEventLab('%.2f'%freq, nodelab=nodelab, dxy=dxy, textorbit=10)
			nnodeevts[nodelab] += 1
		for le in kw.get('losses', []):
			nodelab, freq = le[:2]
			if len(le)>2: hev = le[3]
			else: hev = nnodeevts.setdefault(nodelab, 0)*5
			xml += svgEventLab('L', nodelab=nodelab, dxy=dxy, inCircle=True, orbit=(hev, 0))
			#~ xml += svgEventLab('%.2f'%freq, nodelab=nodelab, dxy=dxy, textorbit=10)
			nnodeevts[nodelab] += 1
		for nodelab, freq in kw.get('counts', []):
			#~ xml += svgEventLab('%d'%int(round(freq)), nodelab=nodelab, x=xlab+50, dxy=dxy, inCircle=True)
			for k in range(int(round(freq))):
				xml += svgEventLab('#', nodelab=nodelab, x=xlab+70+(k*15), dxy=dxy, inCircle=True)
	if kw.get('treetype')=='reconciliation':
		for nodelab, freq in kw.get('transfers', []):
			xml += svgEventLab('T', nodelab=nodelab, dxy=dxy, inCircle=True)
		for nodelab, freq in kw.get('duplications', []):
			xml += svgEventLab('D', nodelab=nodelab, dxy=dxy, inCircle=True)
		for nodelab, freq in kw.get('losses', []):
			xml += svgEventLab('L', nodelab=nodelab, dxy=dxy, inCircle=True)
		
	
	if dbranchwidth:
		# key for branch width coding
		maxcount = int(round(max(dbranchwidth.values())))+1
		keyunit = 200
		xkey = xend*2/3
		ykey = yend+150
		xml += svgEventLab('Gene copy number', x=xkey+(maxcount*keyunit/2), y=ykey-40, fontsize=20)
		for k in range(maxcount):
			widstyle = "stroke:grey; stroke-width:%d; fill:none; stroke-linejoin:round; stroke-linecap:round; "%(k*brwidthfac)
			xml += svgPathAB((xkey+(k*keyunit), ykey), (xkey+((k+1)*keyunit), ykey), style=widstyle)
			xml += svgPathAB((xkey+(k*keyunit), ykey), (xkey+((k+1)*keyunit), ykey))
			xml += svgEventLab(str(k), x=xkey+((k+.5)*keyunit), y=ykey+((k+1)*(brwidthfac/2))+10, fontsize=20)
	
	if svgwrap:
		xmax = xend + (100+55 if padleaves else (35 if padinternalnodes else 0)) + 400
		ymax = yend + 300
		svgxml =  '<?xml version="1.0" encoding="ISO-8859-1" standalone="no"?>\n'
		svgxml += '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" xmlns:xlink="http://www.w3.org/1999/xlink"\n'
		svgxml += 'width="%f"\nheight="%f" >\n'%(xmax, ymax)
		svgxml += xml
		svgxml += '</svg>'
		return svgxml
	else:
		return xml
