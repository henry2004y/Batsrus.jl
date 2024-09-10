# ---
# title: Subplot Line Animation
# id: demo_lines_animation
# date: 2024-07-15
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.10.4
# description: Subplot 1D line animation using pyplot
# ---

This example shows how to create multi-panel 1D space line animation from a series of FLEKS outputs.

* Twin axes sharing the same x-axis for showing two quantities in one frame.
* Time-dependent title.
* Fixed value ranges.

Since currently we don't have true 1D outputs, in this demo we show the way to extract the first column from 2D cuts.

```julia
using Batsrus, Printf, PyPlot

function create_figure()
   fig, axs = subplots(5, 1, figsize=(12, 8),
      sharex=true, sharey=false, constrained_layout=true)

   lw1 = 1.0
   lw2 = 1.0
   lw3 = 1.0

   xmin, xmax = -6300.0, -1300.0

   ρemin, ρemax = 0.0, 0.05 # mi/me = 400
   ρimin, ρimax = 0.0, 20.0
   vmin, vmax = -600.0, 300.0
   pmin, pmax = 1e-3, 1.0
   Bmin, Bmax = -3.0, 30.0
   Emin, Emax = -5.0, 5.0

   l11, = axs[1].plot([], [], lw=lw1)
   l211, = axs[2].plot([], [], lw=lw1, label="Uex")
   l212, = axs[2].plot([], [], lw=lw1, label="Uey")
   l213, = axs[2].plot([], [], lw=lw1, label="Uez")
   l311, = axs[3].plot([0.0, 1.0], [pmin, pmax], lw=lw1, alpha=0.8, label=L"$P_{exx}$")
   l312, = axs[3].plot([0.0, 1.0], [pmin, pmax], lw=lw1, alpha=0.8, label=L"$P_{eyy}$")
   l313, = axs[3].plot([0.0, 1.0], [pmin, pmax], lw=lw1, alpha=0.8, label=L"$P_{ezz}$")
   l41, = axs[4].plot([], [], lw=lw2, label="x")
   l42, = axs[4].plot([], [], lw=lw2, label="y")
   l43, = axs[4].plot([], [], lw=lw2, label="z")
   l51, = axs[5].plot([], [], lw=lw2, label="x")
   l52, = axs[5].plot([], [], lw=lw2, label="y")
   l53, = axs[5].plot([], [], lw=lw2, label="z")

   axs[1].set_xlim(xmin, xmax)
   axs[1].set_ylim(ρemin, ρemax)
   axs[2].set_ylim(vmin, vmax)
   axs[3].set_yscale("log")
   axs[3].set_ylim(pmin, pmax)
   axs[4].set_ylim(Bmin, Bmax) # B
   axs[5].set_ylim(Emin, Emax) # E

   for ax in axs
      ax.grid(true)
   end

   axs[4].set_ylabel("B [nT]"; fontsize)
   axs[5].set_ylabel("E [mV/m]"; fontsize)

   axs[end].set_xlabel("x [km]"; fontsize)

   axs[1].set_ylabel(L"$\rho_e$ [amu/cc]"; fontsize, color="tab:blue")
   axs[1].tick_params(axis="y", labelcolor="tab:blue")

   ax12 = axs[1].twinx()
   ax12.set_ylim(ρimin, ρimax)

   ax12.set_ylabel(L"$\rho_i$ [amu/cc]"; fontsize, color="tab:red")

   l12, = ax12.plot([], []; lw=lw1, color="tab:red", alpha=0.8)
   ax12.tick_params(axis="y", labelcolor="tab:red")

   axs[2].set_ylabel(L"$U_e$ [km/s]"; fontsize)

   ax22 = axs[2].twinx()
   ax22.set_ylim(vmin, vmax)

   ax22.set_ylabel(L"$U_i$ [km/s]"; fontsize)

   l221, = ax22.plot([], [], lw=lw3, label="Uix", color="tab:red")
   l222, = ax22.plot([], [], lw=lw3, label="Uiy", color="tab:purple")
   l223, = ax22.plot([], [], lw=lw3, label="Uiz", color="tab:brown")

   axs[3].set_ylabel(L"$P_e$ [nT]"; fontsize)

   ax32 = axs[3].twinx()
   ax32.set_ylabel(L"$P_i$ [nT]"; fontsize)

   fake_range = [0.0, 1.0]
   p_range = [pmin, pmax]
   p_alpha = 0.8
   l321, = ax32.plot(fake_range, p_range, lw=lw3, alpha=p_alpha, label=L"$P_{ixx}$",
      color="tab:red")
   l322, = ax32.plot(fake_range, p_range, lw=lw3, alpha=p_alpha, label=L"$P_{iyy}$",
      color="tab:purple")
   l323, = ax32.plot(fake_range, p_range, lw=lw3, alpha=p_alpha, label=L"$P_{izz}$",
      color="tab:brown")

   ax32.set_yscale("log")
   ax32.set_ylim(p_range...)

   leg21 = axs[2].legend(;loc=(0.34, 0.05), ncols=3, frameon=false)
   leg22 = ax22.legend(;loc=(0.0, 0.05), ncols=3, frameon=false)
   leg31 = axs[3].legend(;loc=(0.34, -0.05), ncols=3, frameon=false)
   leg32 = ax32.legend(;loc=(0.0, -0.05), ncols=3, frameon=false)
   leg4 = axs[4].legend(;loc="upper left", ncols=3, frameon=false)
   leg5 = axs[5].legend(;loc="upper left", ncols=3, frameon=false)

   # set the linewidth of each legend object
   legs = (leg21, leg22, leg31, leg32, leg4, leg5)
   for leg in legs
      for legobj in leg.legend_handles
         legobj.set_linewidth(1.5)
      end
   end

   ls = (l11, l12, l211, l212, l213, l221, l222, l223, l311, l312, l313, l321, l322, l323,
      l41, l42, l43, l51, l52, l53)

   axs, ls
end

function slice1d_avg(bd, var, dir::Int=2)
   mean(bd[var], dims=dir) |> vec
end

function animate(files::Vector{String}, axs, ls; outdir="out/", overwrite::Bool=false,
   icut::Int=1, nbox::Int=1, dir=2, doAverage::Bool=false)
   l11, l12, l211, l212, l213, l221, l222, l223, l311, l312, l313, l321, l322, l323,
      l41, l42, l43, l51, l52, l53 = ls

   for (i, file) in enumerate(files)
      @info "$i in $(length(files))"
      outname = outdir*lpad(i, 4, '0')*".png"
      if !overwrite
         isfile(outname) && continue
      end

      bd = load(joinpath(filedir, file))

      x = @views bd.x[:,1,1]

      if doAverage
         slice = let bd = bd, dir = dir, nbox = nbox
            var -> moving_average(slice1d_avg(bd, var, dir), nbox)
         end
      else
         slice = let bd = bd, dir = dir, nbox = nbox, icut = icut
            var -> moving_average(slice1d(bd, var, icut, dir), nbox)
         end
      end

      d = slice("rhos0")
      l11.set_data(x, d)
      d = slice("rhos1")
      l12.set_data(x, d)
      d = slice("Uxs0")
      l211.set_data(x, d)
      d = slice("Uys0")
      l212.set_data(x, d)
      d = slice("Uzs0")
      l213.set_data(x, d)
      d = slice("Uxs1")
      l221.set_data(x, d)
      d = slice("Uys1")
      l222.set_data(x, d)
      d = slice("Uzs1")
      l223.set_data(x, d)
      d = slice("PXXS0")
      l311.set_data(x, d)
      d = slice("PYYS0")
      l312.set_data(x, d)
      d = slice("PZZS0")
      l313.set_data(x, d)
      d = slice("PXXS1")
      l321.set_data(x, d)
      d = slice("PYYS1")
      l322.set_data(x, d)
      d = slice("PZZS1")
      l323.set_data(x, d)

      d = slice_data("Bx")
      l41.set_data(x, d)
      d = slice_data("By")
      l42.set_data(x, d)
      d = slice_data("Bz")
      l43.set_data(x, d)
      d = slice_data("Ex") ./ 1000
      l51.set_data(x, d)
      d = slice_data("Ey") ./ 1000
      l52.set_data(x, d)
      d = slice_data("Ez") ./ 1000
      l53.set_data(x, d)

      title_str = @sprintf "t = %4.1f s" bd.head.time
      axs[1].set_title(title_str; fontsize)

      savefig(outname, bbox_inches="tight", dpi=200)
   end

   return
end

#################

# Data directory
filedir = "./"

const fontsize = 16

pick_file = file -> startswith(file, "z") && endswith(file, ".out")

files = filter(pick_file, readdir(filedir))

axs, ls = create_figure()
animate(files, axs, ls; icol=1, outdir="figures/", overwrite=true)

close()
```
