import matplotlib as mpl
import matplotlib.pylab as plt
import matplotlib.colors as mcolors
import xarray as xr

pico = xr.open_dataset("spatial_pico_1s.nc").squeeze()
picop = xr.open_dataset("spatial_picop_1s.nc").squeeze()

rc_params = {
    "font.size": 6,
}


# SymLogNorm: symmetric log scale with linear region around zero
norm = mcolors.SymLogNorm(
      linthresh=0.1,   # linear threshold - values between -0.1 and 0.1 are linear
      linscale=1.0,    # scaling factor for the linear region
      vmin=-100,       # minimum value (negative)
      vmax=100,        # maximum value (positive)
      base=10          # log base
)                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                              
with mpl.rc_context(rc=rc_params):
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 3.2), sharey=True)
    pico["pico_basal_melt_rate"].plot(ax=axs[0], norm=norm, add_colorbar=False, cmap="twilight_shifted")
    picop["picop_basal_melt_rate"].plot(ax=axs[1], norm=norm, cmap="twilight_shifted")
    axs[0].set_title("PICO")
    axs[1].set_title("PICOP")
    for ax in axs:
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    fig.savefig("pico_vs_picop.png", dpi=300)
