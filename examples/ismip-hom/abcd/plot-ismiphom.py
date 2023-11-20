import matplotlib
import matplotlib.pylab as plt
import netCDF4
import numpy as np
from scipy.interpolate import interp2d

FS = [
    "aas1",
    "aas2",
    "cma1",
    "fpa2",
    "ghg1",
    "jvj1",
    "mmr1",
    "oga1",
    "rhi1",
    "rhi3",
    "spr1",
    "ssu1",
    "yko1",
]
BP = ["ahu1", "ahu2", "bds1", "cma2", "fpa1", "fsa1", "mbr1", "rhi2", "tpa1"]
models = {"FS": FS, "BP": BP}

# Set to "True" to remove models that appear to be obviously wrong or produce poor-quality
# results (oscillations).
remove_outliers = True

outliers = {"a": [], "b": [], "c": ["mbr1"], "d": ["rhi1", "rhi2", "rhi3"]}


def plot_experiment(
    ax,
    experiment,
    length_scale,
    model_type,
    N_samples=501,
    color="blue",
    plot_models=True,
):
    filename = "ismip-hom-{exp}-{length}.npz".format(
        exp=experiment, length=length_scale
    )

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
            ax.plot(xs, v, "-", lw=1, color=color, alpha=0.5)

    vx_mean = np.mean(data, axis=0)
    vx_std = np.std(data, axis=0)

    ax.plot(xs, vx_mean, label="{} mean".format(model_type), lw=2, color=color)


def plot(experiment, length_scales, axs):
    for length_scale, ax in zip(length_scales, axs):
        ax.set_title("{} km".format(int(length_scale)))

        ax.set_xlabel("x (normalized)")
        ax.set_ylabel("vx (m / year)")

        models = True
        plot_experiment(
            ax, experiment, length_scale, "BP", color="green", plot_models=models
        )
        plot_experiment(
            ax, experiment, length_scale, "FS", color="orange", plot_models=models
        )

        plot_pism(ax, experiment, length_scale)

        ax.legend()


def plot_pism(ax, experiment, length_scale):
    "Plot PISM's ISMIP-HOM results"
    filename = pism_prefix + "pism-hom-{}-{}.nc".format(
        experiment.upper(), length_scale
    )

    L = 1e3 * int(length_scale)

    try:
        f = netCDF4.Dataset(filename)
        x = f.variables["x"][:] / L
        if experiment in "bd":
            v = f.variables["uvelsurf"][0, 1, :]
            ax.plot(x, v, color="red", lw=2, label="PISM")
        else:
            V = f.variables["uvelsurf"][0, :, :]
            y = f.variables["y"][:] / L
            v = interp2d(x, y, V)(x, 0.25)
            ax.plot(x, v, color="red", lw=2, label="PISM")

    finally:
        f.close()


def grid_plot(experiment_name):
    fig, axs = plt.subplots(2, 3)
    fig.dpi = 100
    fig.set_size_inches(12, 8)
    # fig.suptitle("ISMIP HOM Experiment {}".format(experiment_name.upper()))
    fig.tight_layout(h_pad=4)
    fig.subplots_adjust(top=0.9, bottom=0.1)

    row1 = plot(experiment_name, ["005", "010", "020"], axs[0])
    row2 = plot(experiment_name, ["040", "080", "160"], axs[1])

    fig.savefig("ismiphom-{}.png".format(experiment_name))


if __name__ == "__main__":
    pism_prefix = "./"

    for ex in "abcd":
        grid_plot(ex)
