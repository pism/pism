#!/usr/bin/env python
import pylab as plt
import numpy as np
from argparse import ArgumentParser
import importlib
import sys
import os.path

""" Produce pie charts using PISM's profiling output produced using
the -profile option. """

parser = ArgumentParser()
parser.add_argument("FILE", nargs=1)
options = parser.parse_args()

filename = options.FILE[0]

dirname, basename = os.path.split(filename)
sys.path.insert(0, dirname)

modulename = os.path.splitext(basename)[0]
r = importlib.import_module(modulename)

colors = [(141, 211, 199), (255, 255, 179), (190, 186, 218), (251, 128, 114),
          (128, 177, 211), (253, 180, 98), (179, 222, 105), (252, 205, 229),
          (217, 217, 217), (188, 128, 189), (204, 235, 197), (255, 237, 111)]
colors = np.array(colors) / 255.0

n_procs = r.size
s = r.Stages["time-stepping loop"]

big_events = ["basal_yield_stress",
              "stress_balance",
              "surface",
              "ocean",
              "age",
              "energy",
              "basal_hydrology",
              "fracture_density",
              "mass_transport",
              "calving",
              "bed_deformation",
              "io"]

small_events = {}
small_events["energy"] = ["ice_energy", "btu"]
small_events["stress_balance"] = ["stress_balance.shallow", "stress_balance.modifier",
                                  "stress_balance.strain_heat", "stress_balance.vertical_velocity"]
small_events["stress_balance.modifier"] = ["sia.bed_smoother",
                                           "sia.gradient", "sia.flux", "sia.3d_velocity"]
small_events["io"] = ["io.backup", "io.extra_file", "io.model_state"]

better_names = {"stress_balance.shallow": "SSA",
                "stress_balance.modifier": "SIA",
                "stress_balance.strain_heat": "Strain heating",
                "stress_balance.vertical_velocity": "Vertical velocity"}


def get_event_times(event, n_procs):
    result = [s[event][j]["time"] for j in range(n_procs)]

    max = np.max(result)
    min = np.min(result)

    if max > 0:
        return max, min, min / max
    else:
        return max, min, 0.0


total_time = np.max([s["summary"][j]["time"] for j in range(n_procs)])


def get_data(event_list):
    "Get event data from the time-stepping loop stage."
    return {e: get_event_times(e, n_procs) for e in event_list if e in list(s.keys())}


def aggregate(data, total_time):
    "Combine small events."
    d = data.copy()
    other = [0, 0, 0]
    other_label = ""
    for event in data:
        if data[event][0] / float(total_time) < 0.01:
            print("Lumping '%s' (%f%%) with others..." % (event,
                                                          100.0 * data[event][0] / total_time))
            del d[event]
            other[0] += data[event][0]
            other[1] += data[event][1]
            if other[0] > 0:
                other[2] = other[1] / other[0]
            else:
                other[2] = 0.0
            other_label += "\n{}".format(event)

    d["other"] = other
    return d


def plot(data, total, grand_total):

    events = [(e, data[e][0]) for e in data]
    events.sort(key=lambda x: x[1])

    def better_name(n):
        if n in list(better_names.keys()):
            return better_names[n]
        else:
            return n

    times = [e[1] for e in events]
    times_percent = [100.0 * t / float(total) for t in times]

    if grand_total is not None:
        labels = ["\n(%3.1f s, %3.1f%%)" % (time, 100.0 * time / grand_total) for time in times]
    else:
        labels = ["\n(%3.1f s)" % time for time in times]

    labels = [name + comment for name, comment in zip(names, labels)]

    explode = [0.05]*len(times)
    plt.pie(times_percent, autopct="%3.1f%%", labels=labels, colors=colors, startangle=0.0, explode=explode)
    plt.axis('equal')


def figure(title, event_list, total, grand_total=None):
    plt.figure()
    plt.title("%s (%s)" % (title, filename))
    data = get_data(event_list)
    plot(aggregate(data, total), total, grand_total)
    # plot(data, total, grand_total)
    return data


big = figure("Time-stepping loop",
             big_events,
             total_time)

energy = figure("Energy step",
                small_events["energy"],
                big["energy"][0], total_time)

stressbalance = figure("Stress balance",
                       small_events["stress_balance"],
                       big["stress_balance"][0], total_time)

sia = figure("SIA",
             small_events["stress_balance.modifier"],
             stressbalance["stress_balance.modifier"][0], total_time)

io = figure("I/O during run",
            small_events["io"],
            big["io"][0], total_time)

plt.show()
