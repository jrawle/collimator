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

# load configuration
with open(sys.argv[-1], "r") as stream:
    config = yaml.load(stream, Loader=yaml.SafeLoader)


dis = config["ringpattern"]["pinhole_distance"]
# (*distance lightsource pinhole in mm*)
widthInRadian =  np.deg2rad(config["ringpattern"]["reflection_width"]) #(*natual width of the reflections*)
ringpos = np.deg2rad(
    config["ringpattern"]["ringpositions"]
)  # (*angular position of diffraction features*)


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
print("using ",object_name )

### make sure that soller does not show with povray but only produces shadow
if "no_image" not in Popen(["tail", "-n", "1" ,soller], stdout=PIPE).communicate()[0].decode():
    Popen(["sed", "-i" ,'$s/}/no_image}/',soller]).communicate()

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


# camera = Camera( 'orthographic','angle 57', 'up', [0,100,0],'right', [100,0,0] ,'location', [0, 100, 0], 'look_at', [0,400, 0])
#camera = Camera( 'orthographic','angle 35', 'up', [0,100,0],'right', [100,0,0] ,'location', [120, 100, 120], 'look_at', [120,600, 120])
camera = Camera(*cam_list)
print(cam_list)


def scene(p, useSoller=False, showMask=True):
    """ Returns the scene at time 't' (in seconds) """
    # p=t #*10-20
    box = Box([1.5, p + dis + 0.001, 1.5], [-1.5, p + dis - 0.001, -1.5])

    cylinders_thin = []
    cylinders_thick = []
    for i in range(0, len(rings)):
        cylinders_thin.append(
            Cylinder([0, p + dis - 0.001, 0], [0, p + dis + 0.001, 0], rings[i])
        )
        cylinders_thick.append(
            Cylinder([0, p + dis - 0.002, 0], [0, p + dis + 0.002, 0], rings[i])
        )

    difs = []
    for i in range(2, len(rings), 2):
        difs.append(Difference(cylinders_thin[i], cylinders_thick[i - 1]))

    myunion = cylinders_thin[0]
    for i in range(0, len(difs)):
        myunion = Union(myunion, difs[i])

    myunion = Union(myunion, Difference(box, cylinders_thick[-1]))

    s = Scene(
        camera,
        [
            LightSource([0, p, 0], "color", [1, 1, 1]),
            Plane([0, 1, 0], 800, Texture(Pigment("color", [1, 1, 1]))),
        ],
    )
    if showMask:
        s.objects.append(myunion)

    if useSoller:
        s.objects.append(object_name)
        s.included = [soller]

    return s


####### render without soller

frames = []
radial_intensities = []

for i in sim_pos:
    a = scene(i).render(
        width=pov_width,
        height=pov_height,
        antialiasing=pov_antialiasing,
        remove_temp=False,
    )
    b = np.sum(a, axis=2)
    frames.append(a)
    print(i)

clip = ImageSequenceClip(frames, fps=5)
clip.write_gif("no_soller.gif")

##### render with soller

frames = []
radial_intensities = []

for i in sim_pos:
    a = scene(i, useSoller=True).render(
        width=pov_width,
        height=pov_height,
        antialiasing=pov_antialiasing,
        remove_temp=False,
    )
    b = np.sum(a, axis=2)
    frames.append(a)
    print(i)


clip = ImageSequenceClip(frames, fps=5)
clip.write_gif("soller.gif")

##### render without collimator and without pinhole

frames = []
radial_intensities = []
a = scene(0, useSoller=False, showMask=False).render(
    width=pov_width, height=pov_height, antialiasing=pov_antialiasing, remove_temp=False
)
frames.append(a)

a = scene(0, useSoller=True, showMask=False).render(
    width=pov_width, height=pov_height, antialiasing=pov_antialiasing, remove_temp=False
)
frames.append(a)

clip = ImageSequenceClip(frames, fps=5)
clip.write_gif("shadow.gif")
