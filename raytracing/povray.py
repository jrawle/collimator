# -*- coding: utf-8 -*-
#
#
# Copyright (c) 2020 Linus Pithan linus.pithan@esrf.fr
# Distributed under the GNU LGPLv3. See LICENSE for more info.

"""
generates ring pattern in front of collimator to test aceptance volume

call this script with yml file as first parameter. E.g.

$ python raytracing/rings_only_run_povray.py ring_raytracing.yml
"""

from vapory import *
from subprocess import Popen,PIPE
import numpy as np
from moviepy.editor import ImageSequenceClip
from silx.math.histogram import Histogramnd
import yaml
import sys
import os
import re
import matplotlib.pyplot as plt
import cairo

detdist = 500

# load configuration
#with open(sys.argv[-1], "r") as stream:
#    config = yaml.load(stream, Loader=yaml.SafeLoader)
with open("rings.yml", "r") as stream:
    config = yaml.load(stream, Loader=yaml.SafeLoader)

dis = config["ringpattern"]["pinhole_distance"]
# (*distance lightsource pinhole in mm*)
widthInRadian =  np.deg2rad(config["ringpattern"]["reflection_width"]) #(*natual width of the reflections*)
ringpos = np.deg2rad(
    config["ringpattern"]["ringpositions"]
)  # (*angular position of diffraction features*)

def readSoller(soller):
    #catch input object name
    object_name=None
    with open(soller, "r") as file:
        for line in file:
          if re.search("declare", line):
             object_name=line.rsplit("=")[0].split("declare")[-1].strip(" ")
    assert object_name is not None
    if "-" in object_name:
        newname=object_name.replace("-","_")
        Popen(["sed", "-i", "-e", 's/'+object_name+'/'+newname+'/g' , soller]).communicate()
        object_name=newname
    #print("using ",object_name )

    ### make sure that soller does not show with povray but only produces shadow
    if "no_image" not in Popen(["tail", "-n", "1" ,soller], stdout=PIPE).communicate()[0].decode():
        Popen(["sed", "-i" ,'$s/}/no_image}/',soller]).communicate()
    return object_name

def mkgif(a):
    frames=[]
    frames.append(a)
    clip = ImageSequenceClip(frames, fps=5)
    clip.write_gif("shadow.gif")


def CalcMinMaxRad(r):
    return np.transpose(
        np.array(
            [
                np.tan(r - widthInRadian / 2) * dis,
                np.tan(r + widthInRadian / 2) * dis,
            ]
        )
    )


rings = CalcMinMaxRad(ringpos)
rings = np.array(rings).flatten()

##position in the simulation
sim_pos = np.arange(
    config["simulation_range"]["min"],
    config["simulation_range"]["max"],
    config["simulation_range"]["step"],
)

# povray output parameters
pov_width = config["povray_output"]["width"]
pov_height = config["povray_output"]["width"]
pov_antialiasing = config["povray_output"]["antialiasing"]
soller = config["input_file"]
object_name = readSoller(soller)

## create output dir if needed
if not os.path.exists(config["output_dir"]):
    os.makedirs(config["output_dir"])

# parse camera configuration
cam_list = list()
for x in config["camera"]:
    if isinstance(x, dict):
        xx = list(x.items())[0]
        if isinstance(xx[1], list):
            cam_list.extend(xx)
        else:
            cam_list.append(str(xx[0]) + " " + str(xx[1]))
    else:
        cam_list.append(x)

fieldofview = 300

# rotation matrix around Y axis (as Y and Z swapped in stl2pov)
#camera = Camera('orthographic','angle', 45,'location', [450, 129, 0], 'direction', [1, 0, 0], 'up', [0, 1, 0], 'right', [0, 0, 1])
#camera = Camera( 'orthographic', 'location', [0, 129, 0], 'direction', [1,0,0], 'up', [0,fieldofview,0], 'right', [0,0,fieldofview], 'rotate', [0, 0, delta], 'rotate', [0, gamma, 0])
#camera = Camera( 'orthographic', 'location', [450, 129, 0], 'direction', [1,0,0], 'up', [0,1,0], 'right', [0,0,1], 'angle', camangle)
#camera = Camera( 'orthographic','angle 35', 'up', [0,100,0],'right', [100,0,0] ,'location', [120, 100, 120], 'look_at', [120,600, 120])
#camera = Camera(*cam_list)
#print(cam_list)

def plot(image, colorbar=False):
    plt.imshow(image)
    if colorbar:
        plt.colorbar()
    plt.show()

def plotclr(image, colorbar=True):
    plt.imshow(image.astype(int).sum(2))
    if colorbar:
        plt.colorbar()
    plt.show()



def circle(d, twotheta, delta2t=0.025):
    global detdist
#    pixscale = pov_width/(2 * 50 * np.tan(np.radians(camangle)/2))
#    pixscale = 1.552
    pixscale = 601/300

    radius = (detdist * np.tan(np.radians(twotheta)))*pixscale
    deltarad = (np.tan(np.radians(twotheta+delta2t)) - np.tan(np.radians(twotheta-delta2t))) * detdist * pixscale

    data = np.zeros((pov_width, pov_height, 4), dtype=np.uint8)
    surface = cairo.ImageSurface.create_for_data(
        data, cairo.FORMAT_ARGB32, pov_width, pov_height)
    cr = cairo.Context(surface)

    # draw white ring in RGB (centre, radius, start end angles)
    cr.arc(pov_width/2, pov_width/2 + (288/2-15)*pixscale, radius, -np.pi, 0)
    cr.set_line_width(deltarad)
    cr.set_source_rgb(1.0, 1.0, 1.0)
    cr.stroke()

    # make into mask then multiply data by it
    mask = data[:,:,0:3]/255
    #mask = mask + 0.5 # for debugging
    return (d * mask).astype(int)

def pattern(scene):
    # no soller for normalising fall-off from source
    #a = scene(0, y=0, useSoller=False, showMask=False).render(width=pov_width, height=pov_height, antialiasing=pov_antialiasing, remove_temp=False, auto_camera_angle=False)

    # with soller
    #b = scene(0, y=0, useSoller=True, showMask=False).render(width=pov_width, height=pov_height, antialiasing=pov_antialiasing, remove_temp=False, auto_camera_angle=False)

    p = []
    for i in np.arange(0, 31, 0.1):
        c = circle(scene, i, 0.1)
        p.append(np.sum(c))
    return p

def scene(x, y=0, z=0, useSoller=False, showMask=True, soller=soller, delta=0, gamma=0):
    """ Returns the scene at time 't' (in seconds) """
    global detdist, fieldofview

    camera = Camera( 'orthographic', 'location', [0, 129, 0], 'direction', [1,0,0], 'up', [0,fieldofview,0], 'right', [0,0,fieldofview], 'rotate', [0, 0, delta], 'rotate', [0, gamma, 0])

    box = Box([x + dis + 0.001, z+1.5, y+1.5], [x + dis - 0.001, z-1.5, y-1.5])

    cylinders_thin = []
    cylinders_thick = []
    for i in range(0, len(rings)):
        cylinders_thin.append(
            Cylinder([x + dis - 0.001, z, y], [x + dis + 0.001, z, y], rings[i])
        )
        cylinders_thick.append(
            Cylinder([x + dis - 0.002, z, y], [x + dis + 0.002, z, y], rings[i])
        )

    difs = []
    for i in range(2, len(rings), 2):
        difs.append(Difference(cylinders_thin[i], cylinders_thick[i - 1]))

    myunion = cylinders_thin[0]
    for i in range(0, len(difs)):
        myunion = Union(myunion, difs[i])

    myunion = Union(myunion, Difference(box, cylinders_thick[-1]))

    # Pilatus 2M dimensions 253.7 * 288.8
    detector = Box([detdist, -15, -253.7/2], [detdist+1, 288.8-15, 253.7/2], Texture(Pigment("color", [1, 1, 1])))

    s = Scene(
        camera,
        [
            LightSource([x, z, y], "color", [1, 1, 1]),
            Object(detector, 'rotate', [0, 0, delta], 'rotate', [0, gamma, 0])
            #Plane([1, 0, 0], detdist+1, Texture(Pigment("color", [1, 0, 0]))),
        ],
    )
    if showMask:
        s.objects.append(myunion)

    if useSoller:
        object_name = readSoller(soller)
        ob = Object(object_name, 'rotate', [0, 0, delta], 'rotate', [0, gamma, 0])
        s.objects.append(ob)
        s.included = [soller]

    return s

def run():
    global pblank, ra, rb, xaxis
    #blank = scene(10, y=0, useSoller=False, showMask=False).render(width=pov_width, height=pov_height, antialiasing=pov_antialiasing, remove_temp=False, auto_camera_angle=False)
    #pblank = pattern(blank)

    ra = 0
    #rb = 0

    for i in np.arange(-10,10.5,0.5):

        print("i = " + str(i))
    
        a = scene(i, y=0, useSoller=True, showMask=False, soller="../models/c11.inc").render(width=pov_width, height=pov_height, antialiasing=pov_antialiasing, remove_temp=False, auto_camera_angle=False)

        #pa = pattern(a)
        if type(ra) == int:
            ra = a
        else:
            ra = np.dstack((ra, a))

        ## b = scene(i, y=0, useSoller=True, showMask=False, soller="../models/c2.inc").render(width=pov_width, height=pov_height, antialiasing=pov_antialiasing, remove_temp=False, auto_camera_angle=False)
        ## pb = pattern(b)
        ## if type(rb) == int:
        ##     rb = pb
        ## else:
        ##     rb = np.dstack((rb, pb))

    np.save("c11.npy", ra)
##    np.save("c2.npy", rb)    
    
#    xaxis=np.arange(0,31,0.1)
#    plt.plot(xaxis, np.array(np.sum(ra, axis=2)[0])/np.array(pblank)/ra.shape[2], label="Variable channels" )
#    plt.plot(xaxis, np.array(np.sum(rb, axis=2)[0])/np.array(pblank)/rb.shape[2], label="Fixed channels" )
#    plt.legend()
#    plt.show()

def run2(delta=0, gamma=0):
    global na
    aa = []
    print("delta = " + str(delta) + ", gamma = " + str(gamma), end="", flush=True)
    for i in np.arange(-5,5.5,0.5):
        print(" " + str(i), end="", flush=True)
        a = scene(i, y=0, useSoller=True, showMask=False, soller="../models/c9.inc", delta=delta, gamma=gamma).render(width=pov_width, height=pov_height, antialiasing=pov_antialiasing, remove_temp=False, auto_camera_angle=False)
        aa.append(a)
    na = aa[0]
    for i in range(1, len(aa)):
        na = np.dstack((na, aa[i]))
    np.save("/scratch/soller/test_0.7_d" + str(delta) + "g" + str(gamma)+ "_stack.npy", na)

    #plt.imshow(np.sum(na[:,:,30:92], axis=2))
    #plt.colorbar()
    #plt.show()

def rot_threshold(t):
    global na, nb
    #na=np.load("c6_stack.npy")
    #nb=np.load("c7_stack.npy")
    o1=np.sum(na[:,:,30:92], axis=2)/np.sum(bg[:,:,30:92], axis=2)
    o2=np.sum(nb[:,:,30:92], axis=2)/np.sum(bg[:,:,30:92], axis=2)
    p1=np.where(o1>o2, o1, o2)
    p1=np.where(p1>=t, p1, np.nan)
    plot(p1, True)
#    p2=np.where(o2>t, o2, np.nan)
#    plot((p1+p2)/2, True)
        
def arbrot(r, t):
    global na
    o=np.sum(na[:,:,30:92], axis=2)/np.sum(bg[:,:,30:92], axis=2)
    pivot=[560,300]
    padX = [o.shape[1] - pivot[0], pivot[0]]
    padY = [o.shape[0] - pivot[1], pivot[1]]
    op = np.pad(o, [padX, padY], 'constant')
    op1=np.nan_to_num(op,nan=0)
    op2=np.nan_to_num(op,nan=0)
    op2r=ndimage.interpolation.rotate(op2, r, reshape=False)
    out = (np.where(op1>t, op1, op2r) + np.where(op2r>t, op2r, op1))/2
    print(out.max())
    plot(out[0:680, 300:900], True)
        
def add(delta, gamma):
    global p, q, bg
    nb = np.load("/scratch/soller/c2d" + str(delta) + "g" + str(gamma) + "_stack.npy")
    o = np.sum(nb[:,:,21:42], axis=2)/np.sum(bg[:,:,30+21:93-21], axis=2)
    p1 = np.where(o>=1, o, 0)
    try:
        if p.shape == (601, 601):
            p = p + p1
        else:
            p = p1
    except AttributeError:
        p = p1
    q=np.where(p==0, 1, 0)

####### render 

#frames=[]

#a = scene(0, useSoller=False, showMask=False).render(width=pov_width, height=pov_height, antialiasing=pov_antialiasing, remove_temp=False)
#frames.append(a)

#clip = ImageSequenceClip(frames, fps=5)
#clip.write_gif("shadow.gif")
