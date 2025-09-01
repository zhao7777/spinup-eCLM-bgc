# Integrated spinup-stability analysis SpinupStability_BGC_v11.ncl 
# https://github.com/ESCOMP/CTSM/blob/master/tools/contrib/SpinupStability_BGC_v11.ncl
# With geocat.viz-style plots

import xarray as xr
import numpy as np         
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geocat.viz.util as gvutil
import os
os.environ["MPLBACKEND"] = "Agg" 
from glob import glob

# -----------------------------
# Helpers
# -----------------------------

def _open_cycle_dataset(cdir, years):
    files = []
    for y in years:
        files.extend(sorted(glob(f"{cdir}/{caseid}.clm2.{hist_ext}.{y:04d}*.nc")))
    if not files:
        raise FileNotFoundError(f"No files for years {years} found in {cdir}")
    # Decode CF time with cftime objects so noleap calendars work reliably
    return xr.open_mfdataset(
        files,
        combine="by_coords",
        decode_times=True,
        use_cftime=True,
    )
    
def _annual_means_by_year(ds):
    """Return dict {year: dataset_of_annual_means_for_that_year}.
    Handles monthly CLM outputs with CF calendars (e.g., noleap) by averaging 12 timesteps.
    If dataset has no time coord, treat it as already annual for start_y.
    """
    if "time" not in ds.dims and "time" not in ds.coords:
        return {start_y: ds}
    grouped = ds.groupby("time.year").mean("time")
    groups = {int(y): grouped.sel({"year": y}).drop_vars("year") for y in grouped["year"].values}
    return groups

def _get_last_year_annual(groups, end_year):
    """From {year: ds_year}, return ds_year at end_year. If missing, fall back to max available."""
    if end_year in groups:
        return groups[end_year]
    last_key = sorted(groups.keys())[-1]
    return groups[last_key]

def _as_float(x):
    """Turn a 0-D DataArray / dask array / numpy scalar into a Python float."""
    try:
        # xarray.DataArray
        if hasattr(x, "compute"):
            x = x.compute()       # trigger Dask
        if hasattr(x, "values"):
            x = x.values          # get ndarray / scalar
        return float(x)            # Python float
    except Exception:
        return float(np.asarray(x))

def _global_aggregate(varname, annual_da, landarea):
    """Return scalar global aggregate for the annual field, with units aligned to NCL logic.
    - Carbon pools (TOTECOSYSC/TOTSOMC/TOTVEGC): sum over land area, convert g->Pg.
    - GPP: convert gC s^-1 -> gC yr^-1, area-sum, then to Pg yr^-1.
    - TLAI: area-weighted mean (unit m2/m2).
    - TWS: area-weighted mean (mm), then convert to meters.
    annual_da dims: (lat, lon) or with extra singleton dims.
    """
    # align dims
    da = annual_da
    if varname in ("TOTECOSYSC", "TOTSOMC", "TOTVEGC"):
        res = (da * landarea).sum() * 1e-15                 # Pg
        return _as_float(res)
    if varname == "GPP":
        res = (da * secinyr * landarea).sum() * 1e-15       # Pg/yr
        return _as_float(res)
    if varname == "TLAI":
        res = (da * landarea).sum() / landarea.sum()        # m2/m2 (weighted mean)
        return _as_float(res)
    if varname == "TWS":
        res = ((da * landarea).sum() / landarea.sum()) / 1e3  # m (mm→m)
        return _as_float(res)
    raise ValueError(f"Unsupported variable {varname}")

def _compute_deltas_over_subper(series, subper):
    series = np.asarray(series)
    return (series[1:] - series[:-1]) / float(subper)

def _equilibrium_cycle(series, thresh, subper):
    """Replicate NCL's sustained-to-end equilibrium logic.
    Returns the cycle index (0-based) after which all subsequent deltas are within threshold,
    or None if not equilibrated by the end.
    """
    d = _compute_deltas_over_subper(series, subper)
    if d.size == 0:
        return None
    ok = np.abs(d) < thresh
    if not ok.any():
        return None
    if ok[-1]:
        bad = np.where(~ok)[0]
        return (bad[-1] + 1) if bad.size else 0
    return None

def _squeeze_2d(da):
    # Pick a single slice on 'time' if present, and drop any leftover singleton dims
    if "time" in da.dims:
        da = da.isel(time=0)
    return da.squeeze()
    
def _pct_landarea_disequil_totecosysc(fields_per_cycle, landarea2d, cell_thresh, subper):
    """Compute % land area in disequilibrium for TOTECOSYSC only.
    fields_per_cycle: list of annual TOTECOSYSC fields (xarray DataArray) at the end of each cycle (units gC m^-2)
    Returns numpy array of % per (cycle-1).
    """
    # concat along a new 'cycle' dim
    da = xr.concat(fields_per_cycle, dim="cycle")
    diffs = da.diff("cycle") / float(subper)       
    mask = np.abs(diffs) > float(cell_thresh)
    
    bad_area = (mask.astype("float64") * landarea2d).sum(dim=("lat", "lon"))
    total_area = landarea2d.sum(dim=("lat", "lon"))
    pct = 100.0 * (bad_area / total_area)
    return pct.values

def _ncl_axes_style(ax, ylabel=None, xlabel=None, title=None):
    ax.grid(False, linewidth=0.6)
    for side in ("top", "right", "left", "bottom"):  
        ax.spines[side].set_visible(True)            
        ax.spines[side].set_linewidth(1.0)  
    ax.tick_params(which="both", direction="in", top=True, right=True) 
    ax.grid(False, linewidth=0.6)
    ax.tick_params(labelsize=8)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=10)
    if title:
        ax.set_title(title, fontsize=10, weight="bold")

def _infer_extent(da, pad_deg=1.0):
    lon = da["lon"].values
    lat = da["lat"].values
    lon_min, lon_max = float(np.nanmin(lon)), float(np.nanmax(lon))
    lat_min, lat_max = float(np.nanmin(lat)), float(np.nanmax(lat))
    # If your lon is 0..360, keep as-is; Cartopy PlateCarree handles both.
    return (lon_min - pad_deg, lon_max + pad_deg, lat_min - pad_deg, lat_max + pad_deg)

def make_total_panel(cycles, series, deltas, eq_cycle, pct_area_noeq, glob_thresh, glob_thresh_area, subper, lons, lats):
    from matplotlib.gridspec import GridSpec
    import matplotlib as mpl

    mpl.rcParams.update({
        "font.family": "sans-serif",  
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],  
        "mathtext.fontset": "stixsans",  
        "axes.unicode_minus": False,
    })

    plt.rcParams.update({"font.size": 8})

    order = [
        ("TOTECOSYSC", "Pg C"), ("ΔTOTECOSYSC", "Pg C yr$^{-1}$"),
        ("TOTSOMC", "Pg C"),    ("ΔTOTSOMC", "Pg C yr$^{-1}$"),
        ("TOTVEGC", "Pg C"),    ("ΔTOTVEGC", "Pg C yr$^{-1}$"),
        ("TLAI", "m$^2$ m$^{-2}$"), ("ΔTLAI", "m$^2$ m$^{-2}$ yr$^{-1}$"),
        ("GPP", "Pg C yr$^{-1}$"), ("ΔGPP", "Pg C yr$^{-2}$"),
        ("TWS", "m"),           ("ΔTWS", "m yr$^{-1}$"),
        ("PCTAREA", "%"),
    ]

    fig = plt.figure(figsize=(12, 8), constrained_layout=True)
    gs = GridSpec(4, 4, figure=fig, wspace=0.10, hspace=0.15)

    def add_thresh(ax, thr, color="gray"):
        ax.axhline(thr, ls="--", lw=0.6, color=color)
        ax.axhline(-thr, ls="--", lw=0.6, color=color)
        ax.axhline(0.0, color="black", lw=0.6)
        
    positions = [
        (0,0), (0,1),
        (0,2), (0,3),
        (1,0), (1,1),
        (1,2), (1,3),
        (2,0), (2,1),
        (2,2), (2,3),
        (3,0)
    ]

    varpairs = [
        ("TOTECOSYSC","TOTECOSYSC"), ("TOTECOSYSC","Δ"),
        ("TOTSOMC","TOTSOMC"),       ("TOTSOMC","Δ"),
        ("TOTVEGC","TOTVEGC"),       ("TOTVEGC","Δ"),
        ("TLAI","TLAI"),             ("TLAI","Δ"),
        ("GPP","GPP"),               ("GPP","Δ"),
        ("TWS","TWS"),               ("TWS","Δ"),
        ("PCTAREA", None)
    ]

    years = (cycles - cycles[0]) * subper  # years since spinup start
    
    for (row,col), (v, which) in zip(positions, varpairs):
        ax = fig.add_subplot(gs[row, col])
        if v == "PCTAREA":
            ax.plot(years[1:], pct_area_noeq, marker="o", markersize=3, lw=0.6, color="black")
            ax.axhline(glob_thresh_area, ls="--", lw=0.6, color="gray")
            _ncl_axes_style(ax, ylabel="%",title="% Land Area in TOTECOSYSC Disequil")
            ax.set_xlim(years[0], years[-1])
            print(pct_area_noeq)
            continue
        if which != "Δ":
            y = series[v]
            ax.plot(years, y, marker="o", markersize=3, lw=0.6, color="black")
            if eq_cycle.get(v) is not None:
                ax.axvline(eq_cycle[v], ls=":", lw=0.8, color="gray")
                ax.text(0.98, 0.02, f"EqCycle={eq_cycle[v]}", transform=ax.transAxes,va="bottom", ha="right", fontsize=8)
            _ncl_axes_style(ax, ylabel=dict(
                TOTECOSYSC="Pg C", TOTSOMC="Pg C", TOTVEGC="Pg C",
                TLAI="m$^2$ m$^{-2}$", GPP="Pg C yr$^{-1}$", TWS="m"
            )[v], title=v)
            ax.set_xlim(years[0], years[-1])
        else:
            y = deltas[v]
            ax.plot(years[1:], y, marker="o", markersize=3, lw=0.6, color="black")
            add_thresh(ax, glob_thresh[v], color="gray")
            _ncl_axes_style(ax, ylabel=dict(
                TOTECOSYSC="Pg C", TOTSOMC="Pg C", TOTVEGC="Pg C",
                TLAI="m$^2$ m$^{-2}$", GPP="Pg C yr$^{-1}$", TWS="m"
            )[v], xlabel="Cycle", title=f"Δ {v}")
            if v in ("TOTECOSYSC","TOTSOMC","TOTVEGC","TLAI","GPP"):
                ax.set_ylim(-0.8, 0.8)
            elif v in ("TWS","TLAI",):
                ax.set_ylim(-0.04, 0.04)
            ax.set_xlim(years[0], years[-1])


    last = (fields_totecosysc[-1] - fields_totecosysc[-2]) / float(subper)
    prev = (fields_totecosysc[-2] - fields_totecosysc[-3]) / float(subper)

    rotated_pole_kwargs = dict(
    pole_longitude=-162.0,           
    pole_latitude=39.25,              
    central_rotated_longitude=0.0, globe=None
)
    proj_map = ccrs.RotatedPole(**rotated_pole_kwargs)

    # vmax = float(np.nanmax(np.abs([last.values, prev.values])))
    # vmin = -vmax
    vmax = 10
    vmin = -vmax

    ax_map1 = fig.add_subplot(gs[3, 1], projection=proj_map) 
    ax_map2 = fig.add_subplot(gs[3, 2], projection=proj_map)  
    for ax, data, ttl in [
        (ax_map1, last, f"ΔTOTECOSYSC Yr: {cycles[-2]*subper+1}–{cycles[-1]*subper+1}"),
        (ax_map2, prev, f"ΔTOTECOSYSC Yr: {cycles[-3]*subper+1}–{cycles[-2]*subper+1}")
    ]:
        # ax.set_extent(extent, crs=ccrs.PlateCarree()) 
        ax.add_feature(cfeature.COASTLINE, linewidth=0.4)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle=":")
        mesh = ax.pcolormesh(
            lons, lats, data.values,
            transform=ccrs.PlateCarree(),
            shading="auto", vmin=vmin, vmax=vmax, cmap="RdBu_r"
        )
        ax.set_title(ttl + " (gC m$^{-2}$ yr$^{-1}$)", fontsize=8)
    
    cbar = fig.colorbar(
        mesh, ax=[ax_map1, ax_map2], orientation="vertical",
        pad=0.02, shrink=0.85, extend='both'
    )
    cbar.ax.tick_params(labelsize=8, direction="in")
    fig.savefig("figs/spinup_TOTAL_PANEL.png", dpi=200, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)

# -----------------------------
# Configuration
# -----------------------------
caseid   = "clmoas"
hist_ext = "h0"
start_y = 1960
n_years = 10
subper  = n_years
do_plot = True

# NOTE:
# This spinup iteratively runs from 1950 to 1960.
# Each cycle (e.g., cycle00, cycle01, ...) is a repetition of that period.
# Files are located in directories like:
#  - 19600101_19700101_cycle00
#  - 19600101_19700101_cycle01
#  - ... up to cycle10

# Paths
# Spinup directory pattern (adjust to your system)
base_spinup_dir = f"/p/scratch/cjibg36/jibg3674/eCLM_BGC_SPINUP/AD"
cycle_dirs = []
for d in sorted(glob(f"{base_spinup_dir}/{start_y:04d}0101_{start_y+n_years:04d}0101_cycle*/")):
    # require that restart file exists in that cycle directory, otherwise skip
    restart_file = os.path.join(d, f"{caseid}.clm2.r.{start_y+n_years:04d}-01-01-00000.nc")
    if os.path.exists(restart_file):
        cycle_dirs.append(d)

# Path to grid/domain file
path_grid = "/p/project1/cjibg36/jibg3674/shared_DA/setup_eclm_cordex_444x432_v9/input_clm"
file_grid = "domain.lnd.EUR-11_EUR-11.230216_mask.nc"
domain = xr.open_dataset(os.path.join(path_grid, file_grid))

lat = domain['yc']
lon = domain['xc']

# Variables & thresholds
variables = ["TOTECOSYSC", "TOTSOMC", "TOTVEGC", "TLAI", "GPP", "TWS"]
glob_thresh = {
    "TOTECOSYSC": 0.02,   # Pg C / yr
    "TOTSOMC":    0.02,   # Pg C / yr
    "TOTVEGC":    0.02,   # Pg C / yr
    "TLAI":       0.02,   # m2/m2 per yr (change in global mean)
    "GPP":        0.02,   # Pg C / yr^2 (change in PgC/yr per yr)
    "TWS":        0.001,  # m / yr
}

# For % land area in disequilibrium (TOTECOSYSC only)
glob_thresh_area = 3.0   # %
cell_thresh      = 1.0   # gC m^-2 yr^-1 (per cell absolute delta threshold)

secinyr = 60*60*24*365

# -----------------------------
# Load cycles, compute annual fields for end-of-cycle year, and aggregate
# -----------------------------
print(f"Found {len(cycle_dirs)} cycle directories.")

# Per-variable global series across cycles (scalars)
series = {v: [] for v in variables}
# For TOTECOSYSC maps and area-% we also keep the annual field at the end of each cycle
fields_totecosysc = []
# Cache a representative landarea for area-% calc
landarea_da = None

for ic, cdir in enumerate(cycle_dirs):
    years = list(range(start_y, start_y + n_years))
    ds = _open_cycle_dataset(cdir, years)
    if "area" not in ds.variables or "landfrac" not in ds.variables:
        raise KeyError("Dataset must contain 'area' (m2 gridcell) and 'landfrac' (0-1)")
    area = ds["area"] * 1e6  # convert km^2 -> m^2 if needed; adjust if already in m^2
    landfrac = ds["landfrac"]
    landarea = area * landfrac
    if landarea_da is None:
        landarea_da = landarea.load()

    groups = _annual_means_by_year(ds)
    end_year = years[-1]
    ds_last = _get_last_year_annual(groups, end_year)

    # Compute scalars per variable
    for v in variables:
        if v not in ds_last.variables:
            if v == "TWS":
                raise KeyError("TWS is missing; add a reconstruction if needed or drop TWS from variables")
            else:
                raise KeyError(f"Variable {v} missing in dataset")
        val = _global_aggregate(v, ds_last[v], landarea)
        series[v].append(val)

    # Save field for TOTECOSYSC area-% and maps (keep in gC m^-2)
    fields_totecosysc.append(ds_last["TOTECOSYSC"].load())

    ds.close()

cycles = np.arange(len(cycle_dirs))

# -----------------------------
# Compute deltas and equilibrium cycles
# -----------------------------
print("Computing deltas and equilibrium cycles (sustained-to-end)...")

deltas = {v: _compute_deltas_over_subper(series[v], subper) for v in variables}

eq_cycle = {}
for v in variables:
    eq_cycle[v] = _equilibrium_cycle(series[v], glob_thresh[v], subper)
    
area2d     = _squeeze_2d(area)
landfrac2d = _squeeze_2d(landfrac)
landarea2d = (area2d * landfrac2d)
if ic == 0:
    landarea_da = landarea2d

# % land area in disequilibrium for TOTECOSYSC only
pct_area_noeq = _pct_landarea_disequil_totecosysc(fields_totecosysc, landarea2d, cell_thresh, subper)

# For threshold-vs-sustained equilibrium year on area %: we want last sustained-to-end under thresh
area_eq_cycle = _equilibrium_cycle(100.0 - pct_area_noeq, 100.0 - glob_thresh_area, 1)  # same logic: we want >= (100-thresh) to be OK

file_lsm = '/p/project/cjibg36/jibg3674/shared_DA/EUR-11_TSMP_FZJ-IBG3_444x432_LAND-LAKE-SEA-MASK.nc'
grid_centre = xr.open_dataset(file_lsm, decode_times=False)
lats = grid_centre['lat'].values
lons = grid_centre['lon'].values
grid_centre.close()

make_total_panel(cycles, series, deltas, eq_cycle, pct_area_noeq, glob_thresh, glob_thresh_area, subper, lons, lats)
