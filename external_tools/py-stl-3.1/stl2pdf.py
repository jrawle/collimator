#! /usr/bin/env python
# -*- python coding: utf-8 -*-
# Copyright © 2012 R.F. Smith <rsmith@xs4all.nl>. All rights reserved.
# $Date: 2012-06-04 21:17:53 +0200 $
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY AUTHOR AND CONTRIBUTORS ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL AUTHOR OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.

'''Program for converting a view of an STL file into a PDF file.'''

import sys
import os
import stl
import xform
import cairo

name = ('stl2pdf [ver. ' + '$Revision: 3.1 $'[11:-2] + 
       '] ('+'$Date: 2012-06-04 21:17:53 +0200 $'[7:-2]+')')

def usage():
    print name
    print "Usage: stl2pdf infile [outfile] [transform [transform ...]]"
    print "where [transform] is [x number|y number|z number]"

## This is the main program ##
# Process the command-line arguments
validargs = ['x', 'y', 'z', 'X', 'Y', 'Z']
if len(sys.argv) == 1:
    usage()
    sys.exit(0)
infile = sys.argv[1]
if len(sys.argv) < 3 or sys.argv[2] in validargs:
    outfile = None
    del sys.argv[:2]
else:
    outfile = sys.argv[2]
    del sys.argv[:3]
tr = xform.Xform()
while len(sys.argv) > 1:
    if not sys.argv[0] in validargs:
        print "Unknown argument '{}' ignored.".format(sys.argv[0])
        del sys.argv[0]
        continue
    try:
        ang = float(sys.argv[1])
        if sys.argv[0] in ['x','X']:
            tr.rotx(ang)
        elif sys.argv[0] in ['y','Y']:
            tr.roty(ang)
        else:
            tr.rotz(ang)
        del sys.argv[:2]
    except:
        print "Argument '{}' is not a number, ignored.".format(sys.argv[1])
        continue
# Open the file
try:
    stlobj = stl.Surface(infile)
except:
    print "The file '{}' cannot be read or parsed. Exiting.".format(sys.argv[1])
    sys.exit(1)
# Process the file
for result in stlobj.processfacets:
    print result
# Remove spaces from name
stlobj.name = stlobj.name.strip()
# Apply transformations
if tr.unity == False:
    stlobj.xform(tr)
# Calculate viewport and transformation
xmin, xmax, ymin, ymax, zmin, zmax = stlobj.extents()
pr = xform.Zpar(xmin, xmax, ymin, ymax)
# Prepare output.
if outfile == None:
    outbase = os.path.basename(infile)
    if outbase.endswith((".stl", ".STL")):
        outbase = outbase[:-4]
    outfile = outbase+".pdf"
#out = canvas.Canvas(outfile, (pr.w, pr.h), pageCompression=1)
out = cairo.PDFSurface(outfile, pr.w, pr.h)
ctx = cairo.Context(out)
ctx.set_line_cap(cairo.LINE_CAP_ROUND)
ctx.set_line_join(cairo.LINE_JOIN_ROUND)
ctx.set_line_width(0.25)

# Calculate the visible facets
vizfacets = [f for f in stlobj.facets if pr.visible(f.n)]
# Next, depth-sort the facets using the largest z-value of the three vertices.
vizfacets.sort(None, lambda f: max([f.v[0].z, f.v[1].z, f.v[2].z]))
# Project and illuminate the facets
pf = (stl.ProjectedFacet(f, pr) for f in vizfacets)
# Draw the triangles
for f in pf:
    path = ctx.new_path()
    ctx.move_to(f.x1, f.y1)
    ctx.line_to(f.x2, f.y2)
    ctx.line_to(f.x3, f.y3)
    ctx.close_path()
    ctx.set_source_rgb(f.gray, f.gray, f.gray)
    ctx.fill_preserve()
    ctx.stroke()
# Send output.
out.show_page()
out.finish()
