from bokeh.plotting import figure, output_notebook, output_file, show
from bokeh.models import Band, ColumnDataSource
from bokeh.layouts import gridplot

import numpy as np
from scipy.interpolate import interp2d
import netCDF4

FS = ['aas1','aas2','cma1','fpa2','ghg1','jvj1','mmr1','oga1','rhi1','rhi3','spr1','ssu1','yko1']
BP = ['ahu1', 'ahu2', 'bds1', 'cma2', 'fpa1', 'fsa1', 'mbr1', 'rhi2', 'tpa1']
models = {"FS" : FS, "BP" : BP}

# Set to "True" to remove models that appear to be obviously wrong or produce poor-quality
# results (oscillations).
remove_outliers = True

outliers = {"a" : [],
            "b" : [],
            "c" : ["mbr1"],
            "d" : ["rhi1", "rhi2", "rhi3"]}

def plot_experiment(fig, experiment, length_scale, model_type, N_samples=501, color="blue", plot_models=True):

    filename = "ismip-hom-{exp}-{length}.npz".format(exp=experiment, length=length_scale)

    raw_data = np.load(filename)

    participating_models = [model for model in models[model_type] if model in raw_data]

    if remove_outliers:
        for model in outliers[experiment]:
            if model in participating_models:
                participating_models.remove(model)

    data = np.array([raw_data[model] for model in participating_models])

    xs = raw_data["x"]

    if plot_models:
        for model in participating_models:
            v = raw_data[model]
            fig.line(xs, v,
                     line_color=color,
                     name="{} ({})".format(model, model_type))

    vx_mean = np.mean(data, axis=0)
    vx_std = np.std(data, axis=0)

    fig.line(xs, vx_mean,
             legend_label="{} mean".format(model_type),
             line_width=4,
             name="{} mean".format(model_type),
             line_color=color)

def plot(experiment, length_scales):
    plots = []
    for length_scale in length_scales:
        p = figure(tooltips="$name, $x, $y")
        p.title.text = "Experiment {}, {} km".format(experiment.upper(), int(length_scale))
        p.xaxis.axis_label = 'x (normalized)'
        p.yaxis.axis_label = 'vx (m / year)'

        models = True
        plot_experiment(p, experiment, length_scale, "BP", color="green", plot_models=models)
        plot_experiment(p, experiment, length_scale, "FS", color="orange", plot_models=models)

        p.legend.location = "top_left"

        plot_pism(p, experiment, length_scale)

        plots.append(p)

    return plots

def plot_pism(fig, experiment, length_scale):
    "Plot PISM's ISMIP-HOM results"
    filename = pism_prefix + "pism-hom-{}-{}.nc".format(experiment.upper(), length_scale)

    L = 1e3 * int(length_scale)

    try:
        f = netCDF4.Dataset(filename)
        x = f.variables["x"][:] / L
        if experiment in "bd":
            v = f.variables["uvelsurf"][0, 1, :]
            fig.line(x, v, name="PISM", line_color="red", line_width=4, legend_label="PISM")
        else:
            V = f.variables["uvelsurf"][0, :, :]
            y = f.variables["y"][:] / L
            v = interp2d(x, y, V)(x, 0.25)
            fig.line(x, v, name="PISM", line_color="red", line_width=4, legend_label="PISM")

    finally:
        f.close()

def grid_plot(experiment_name):
    output_file("ismip-{}.html".format(experiment_name),
                title="ISMIP HOM Experiment {}".format(experiment_name.upper()),
                mode="cdn")

    row1 = plot(experiment_name, ["005", "010", "020"])
    row2 = plot(experiment_name, ["040", "080", "160"])

    grid_plot = gridplot([row1, row2], plot_width=400, plot_height=400)

    show(grid_plot)

if __name__ == "__main__":
    pism_prefix = "./"

    for ex in "abcd":
        grid_plot(ex)
