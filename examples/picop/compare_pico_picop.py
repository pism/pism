import matplotlib as mpl
import matplotlib.pylab as plt
import matplotlib.colors as mcolors
import xarray as xr

pico = xr.open_dataset("spatial_pico_1s.nc").squeeze()
picop = xr.open_dataset("spatial_picop_1s.nc").squeeze()

rc_params = {
    "font.size": 6,
}

with mpl.rc_context(rc=rc_params):
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 3.2), sharey=True)
    pico["pico_basal_melt_rate"].plot(ax=axs[0], norm=mcolors.LogNorm(), vmin=0.1, add_colorbar=False)
    picop["picop_basal_melt_rate"].plot(ax=axs[1], norm=mcolors.LogNorm(), vmin=0.1)
    axs[0].set_title("PICO")
    axs[1].set_title("PICOP")
    for ax in axs:
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    fig.savefig("pico_vs_picop.png", dpi=300)
