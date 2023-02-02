### OPTI_SIF_CAMS41
### V.Bastrikov (vladislav.bastrikov@lsce.ipsl.fr), P.Peylin (philippe.peylin@lsce.ipsl.fr)
### Laboratoire des Sciences du Climat et de l'Environnement, Gif-sur-Yvette, France

# Import main libraries
import numpy as np, netCDF4, os, sys
from optparse import OptionParser
import  numpy
# Boolean 'TRUE' flag possible notations
ok_flags = ["y", "yes", "true"]

# CTESSEL vegetation classes names, zonal band names, month names
pftnames = { 1: "Crops, mixed farming",  2: "Short grass",            3: "Evergreen needleleaf",  4: "Deciduous needleleaf",  5: "Deciduous broadleaf",
             6: "Evergreen broadleaf",   7: "Tall grass",             8: "Desert",                9: "Tundra",               10: "Irrigated crops",
            11: "Semidesert",           12: "Ice caps and glaciers", 13: "Bogs and marshes",     14: "Inland water",         15: "Ocean",
            16: "Evergreen shrubs",     17: "Deciduous shrubs",      18: "Mixed forest - Wood",  19: "Interrupted forest",   20: "Water-land mixtures" }
bandnames = { "N" : "North", "T" : "Tropics", "S" : "South" }
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# List of keywords required in task file (run.def)
required_items = ["sif_path", "sif_var", "obs_path", "obs_var", "sim_path", "sim_var",
                  "cvh_path", "cvh_var", "cvl_path", "cvl_var", "tvh_path", "tvh_var", "tvl_path", "tvl_var",
                  "a_prior", "a_min", "a_max", "a_error", "b_prior", "b_min", "b_max", "b_error", "veg_threshold"]

# Get task filename from command line arguments (or use "run.def" by default)
parser = OptionParser()
(options, args) = parser.parse_args()
task_path = "run.def" if len(args) == 0 else args[0]
if not os.path.exists(task_path):
    print "Can't find task file :", task_path
    sys.exit()

# Load and interpret task file
task = dict()
with open(task_path) as stream:
    for line in stream:
        # skip lines starting from '#'
        if line.strip().startswith("#"): continue
        # split each line on '=' considerling the left block to be the keyword and the right block to be the value
        if "=" in line:
            key, value = map(str.strip, line.split("#")[0].split("="))
            if key == "" or value == "": continue
            task[key] = value

# Check presence of required items in task file
for item in required_items:
    if item not in task:
        print "No data for the field '%s' in task file" % item
        sys.exit()

# Import additional libraries if needed
if "opt_method" in task and task["opt_method"] == "bfgs":
    from scipy.optimize import fmin_l_bfgs_b
if "ok_plot" in task and task["ok_plot"].lower() in ok_flags:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap

# Prepare smoothing window value (3 for default setting, 0 for turning smoothing off)
if "ok_smooth" in task and task["ok_smooth"].lower() in ok_flags:
    window = 3
elif "ok_smooth" in task and task["ok_smooth"].isdigit():
    window = int(task["ok_smooth"])
else:
    window = 0

# Load data for SIF, OBS-GPP, SIM-GPP and PFT maps, mask invalid data entries, apply scaling coefficients (if present)
data = dict()
for item in ["sif", "obs", "sim", "cvh", "cvl", "tvh", "tvl"]:

    # Check that file exists and open it
    if not os.path.exists(task["%s_path" % item]):
        print "Can't find %s file : %s" % (item.upper(), task["%s_path" % item])
        sys.exit()
    nc = netCDF4.Dataset(task["%s_path" % item])

    # Check that variable exists
    if task["%s_var" % item] not in nc.variables:
        print "No variable '%s' in %s file" % (task["%s_var" % item], item.upper())
        sys.exit()

    # Get data and mask invalid entries
    data[item] = np.ma.masked_invalid(nc.variables[task["%s_var" % item]][:])

    # If PFT data is given at 3D (time, lat, lon), take the first time step data
    if item in ["cvh", "cvl", "tvh", "tvl"] and len(data[item].shape) == 3: data[item] = data[item][0]

    # Apply scaling coefficient if provided
    if "%s_coef" % item in task: data[item] = data[item] * eval(task["%s_coef" % item])

    # Load latitudes/longitudes/areas if provided (obligatory for plots production)
    if "%s_lat" % item in task: lat = nc.variables[task["%s_lat" % item]][:]
    if "%s_lon" % item in task: lon = nc.variables[task["%s_lon" % item]][:]
    if "%s_area" % item in task: area = nc.variables[task["%s_area" % item]][:]
    nc.close()

    # Print some statistics on the data loaded
    print "%s data loaded : shape %s, range [%.3g, %.3g], %s data points"  % (item.upper(), data[item].shape, data[item].min(), data[item].max(), data[item].count())

    # Apply smoothing on SIF, OBS and SIM data if needed
    if item in ["sif", "obs", "sim"] and window > 0:
        data[item] = smooth(data[item])
        print "%s data smoothed : shape %s, range [%.3g, %.3g], %s data points"  % (item.upper(), data[item].shape, data[item].min(), data[item].max(), data[item].count())

# Select points where all input data streams are valid
comask = (data["sif"].mask == True) | (data["obs"].mask == True) | (data["sim"].mask == True)
SIF = np.ma.masked_where(comask, data["sif"])
OBS = np.ma.masked_where(comask, data["obs"])
SIM = np.ma.masked_where(comask, data["sim"])

# Set dimensions
npft = len(pftnames)
ntime, nlat, nlon = SIF.shape
nyears = ntime / 12
nobs = SIF.count()
nveg = np.where(data["cvh"] + data["cvl"] > 0, 1, 0).sum()
print "ntime : %s (number of time steps)" % ntime
print "nobs  : %s (number of synchronous SIF & OBS & SIM data points across all time steps)" % nobs
print "nveg  : %s (number of points with vegetation cover)" % nveg

# Create month indeces across all years for seasonal optimization
months = np.array(months * nyears)
itime = dict(winter = np.where((months == "Jan") | (months == "Feb") | (months == "Dec")),
             spring = np.where((months == "Mar") | (months == "Apr") | (months == "May")),
             summer = np.where((months == "Jun") | (months == "Jul") | (months == "Aug")),
             autumn = np.where((months == "Sep") | (months == "Oct") | (months == "Nov")),
             alldata = np.arange(ntime))

# Create vegetation cover maps for individual PFTs based on the threshold given
pftmaps = dict()
for mode in ["l", "h"]:
    mode_msk = np.where(data["cv%s" % mode] / (data["cvh"] + data["cvl"]) > float(task["veg_threshold"]), 1, 0)
    for pft in np.delete(np.unique(data["tv%s" % mode]), 0):
        pftmap = np.ma.where((mode_msk == 1) & (data["tv%s" % mode] == pft), data["cv%s" % mode], np.ma.masked)

        # if ok_zonal is set, further split PFT into bands 90N-30N, 30N-10S, 10S-90S
        if "ok_zonal" in task and task["ok_zonal"].lower() in ok_flags:
            pftmapN = np.ma.where(lat[:,None] > 30, pftmap, np.ma.masked)
            pftmapT = np.ma.where((lat[:,None] < 30) & (lat[:,None] > -10), pftmap, np.ma.masked)
            pftmapS = np.ma.where(lat[:,None] < -10, pftmap, np.ma.masked)
            if pftmapN.count() != 0: pftmaps["%02d_N" % pft] = pftmapN
            if pftmapT.count() != 0: pftmaps["%02d_T" % pft] = pftmapT
            if pftmapS.count() != 0: pftmaps["%02d_S" % pft] = pftmapS
        else:
            pftmaps["%02d" % pft] = pftmap

# Print PFT distribution statistics
npft = len(pftmaps)
print "PFT list (%s):" % npft
rest = nveg
for pft in sorted(pftmaps.keys()):
    title = pftnames[int(pft.split("_")[0])]
    if "_" in pft: title = title + ", " + bandnames[pft.split("_")[1]]
    print "   %2s : %s points, %.2g%% (%s)" % (pft, pftmaps[pft].count(), (pftmaps[pft].count() * 100. / nveg), title)
    rest = rest - pftmaps[pft].count()
print "   Non-attributed vegetation points : %s, %.2g%%" % (rest, rest * 100. / nveg)

# Convert prior/min/max/error of parameters into floats and create vectors in (a,b)-form
# If error has '%' sign, calculate it as percentage of variation range
print "Optimized parameters data :"
for item in ["a_prior", "a_min", "a_max", "a_error", "b_prior", "b_min", "b_max", "b_error"]:
    if "error" in item and "%" in task[item]:
        task[item] = (task[item[0] + "_max"] - task[item[0] + "_min"]) * float(task[item].split("%")[0]) / 100.
    else:
        task[item] = float(task[item])
    print "   %s = %s" % (item, task[item])
x_prior = np.array([task["a_prior"], task["b_prior"]])
x_min = np.array([task["a_min"], task["b_min"]])
x_max = np.array([task["a_max"], task["b_max"]])
x_error = np.array([task["a_error"], task["b_error"]])

# Function to smooth data (running window mean)
def smooth(data):
    print "smoothing data using running window mean [-%d;+%d]" % (window/2, window/2)
    smdata = np.ma.masked_all_like(data)
    for idx in range(len(data)):
        smdata[idx] = data[max(0, idx - window/2) : (idx + window/2 + 1)].mean(axis=0)
    return smdata

# Function to calculate the model (SIF ~ Hx = a * OBS + b)
def calcHx(xi):
    return xi[0] * OBS_pft + xi[1]

# Function to calculate root-mean-square difference between two data streams
def calcRMSD(x, y):
    return np.mean((x - y)**2)**0.5

# Function to calculate correlation coefficient between two data streams (NB: select only valid points)
def calcCORR(x, y):
    xvalid = x[x.mask == False].flatten()
    yvalid = y[y.mask == False].flatten()
    return np.corrcoef(xvalid, yvalid)[0,1]

# Function to calculate cost-function and its gradient
def calcCost(xi):
    Hxi = calcHx(xi)
    Jo = ((Hxi - SIF_pft)**2).sum() / R
    Jb = ((xi - x_prior)**2 / x_error**2).sum()
    G0 = (OBS_pft * (Hxi - SIF_pft)).sum() / R
    G1 = (Hxi - SIF_pft).sum() / R
    Go = 2 * np.array([G0, G1])
    Gb = 2 * (xi - x_prior) / x_error
    J = Jo + Jb
    G = Go + Gb
    print "   a = %.3g\tb = %.3g\tcost = %s" % (xi[0], xi[1], J)
    return J, G

# Function to do one-step inverse matrix calculation
def calcInv():
    # Calculate HRH = Had^t * R^-1 * Had
    HRH = np.zeros((2, 2))
    HRH[0, 0] = (OBS_pft * OBS_pft).sum() / R
    HRH[0, 1] = OBS_pft.sum() / R
    HRH[1, 0] = OBS_pft.sum() / R
    HRH[1, 1] = OBS_pft.count() / R

    # Create Binv matrix
    Binv = np.eye(2) * x_error ** -2

    # Calculate Bpost = [HRH + B^-1]^-1
    Bpost = np.linalg.inv(HRH + Binv)

    # Calculate xad = Had^t * R^-1 * (y - H(x))
    xad = np.zeros(2)
    xad[0] = (OBS_pft * (SIF_pft - Hx_prior)).sum() / R
    xad[1] = (SIF_pft - Hx_prior).sum() / R

    # Calculate posterior model
    x_post = x_prior + np.dot(Bpost, xad)
    x_post = np.clip(x_post, x_min, x_max)
    return x_post

# Function to calculate global mean weigthed by grid cell areas
def calcMean(data):
    loc_area = np.where(data.mask == False, area, 0)
    return (data * loc_area).sum(axis = (-2,-1)) / loc_area.sum(axis = (-2,-1))

# Function to plot map
def plotMap(item):
    stat = "mean = %.3g\nmin = %.3g\nmax = %.3g" % (calcMean(item["data"]), item["data"].min(), item["data"].max())
    m = Basemap(projection = "mill", resolution = "c")
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels([-80,-60,-30,0,30,60,80], labels=[1, 0, 1, 1])
    m.drawmeridians(np.arange(0., 420., 60.), labels=[0, 0, 0, 1])
    cmap = plt.cm.Reds
    cmesh = m.pcolormesh(meshlon, meshlat, item["data"], shading = "flat", cmap = cmap, latlon = True)
    cbar = m.colorbar(cmesh, pad = 0.08, location = "right")
    plt.annotate(stat, xy = (0.05, 0.05), xycoords = "axes fraction")
    plt.title(item["title"])
    plt.title("MEAN", loc = "right")
    plt.savefig(item["filename"], dpi = 150)
    plt.clf()

# Function to plot time-series
def plotTS(item):
    fig, axs = plt.subplots(1, 1)
    for data, label in zip(item["data"], item["label"]):
        axs.plot(range(len(data)), data, label = label)
    axs.set_ylabel(item["ylabel"])
    axs.set_title(item["title"])
    legend = axs.legend(fontsize=10, loc="upper right")
    legend.get_frame().set_alpha(0.7)
    plt.savefig(item["filename"], dpi = 150)
    plt.close()

# Function to plot scatter-plots
def plotScat(item):
    stat = "r = %.3g\nRMSD = %.3g" % (calcCORR(item["x"], item["y"]), calcRMSD(item["x"], item["y"]))
    fig, axs = plt.subplots(1, 1)
    axs.plot(item["x"].flatten(), item["y"].flatten(), ",", color = "red")
    axs.set_xlabel(item["xlabel"])
    axs.set_ylabel(item["ylabel"])
    axs.set_title(item["title"])
    axs.text(0.95, 0.95, stat, verticalalignment = "top", horizontalalignment = "right", transform = axs.transAxes, fontsize = 14)
    plt.savefig(item["filename"], dpi = 150)
    plt.close()

# Optimize SIF-OBS(GPP) relationship for each PFT
for pft in sorted(pftmaps.keys()):

    # Prepare PFT mask and label
    pftmask = np.ma.where(pftmaps[pft].mask, np.ma.masked, 1)
    pftlabel = "PFT%s (%s)%s" % (pft, pftnames[int(pft.split("_")[0])], "" if "_" not in pft else ", %s" % bandnames[pft.split("_")[-1]])

    # Set seasons to optimize
    if "ok_seasonal" in task and task["ok_seasonal"].lower() in ok_flags:
        seasons = ["winter", "spring", "summer", "autumn"]
    else:
        seasons = ["alldata"]

    a_post = np.zeros((ntime))
    b_post = np.zeros((ntime))
    for season in seasons:
        print "Optimization for PFT%s, %s" % (pft, season)

        # Select data according to season
        SIF_temp = np.ma.take(SIF, itime[season], axis=0)
        OBS_temp = np.ma.take(OBS, itime[season], axis=0)

        # Filter data according to PFT map
        SIF_pft = SIF_temp * pftmask
        OBS_pft = OBS_temp * pftmask

        # Calculate prior data
        Hx_prior = calcHx(x_prior)
        RMSD_prior = calcRMSD(SIF_pft, Hx_prior)
        R = RMSD_prior**2
        print "   nobs : %s" % SIF_pft.count()
        print "   RMSD_prior : %.3g" % RMSD_prior

        # Run optimization
        if "opt_method" in task and task["opt_method"] == "bfgs":
            print "   RUN L-BFGS-B OPTIMIZATION"
            x_post, fmin, info = fmin_l_bfgs_b(func = calcCost, x0 = x_prior, maxfun = 100, pgtol = 1E-5,
                                               bounds = [(task["a_min"], task["a_max"]), (task["b_min"], task["b_max"])])
        else:
            print "   RUN ONE-STEP OPTIMIZATION"
            x_post = calcInv()

        # Save optimal parameters according to season
        np.put(a_post, itime[season], x_post[0])
        np.put(b_post, itime[season], x_post[1])

        # Calculate posterior data
        Hx_post = calcHx(x_post)
        RMSD_post = calcRMSD(SIF_pft, Hx_post)
        reduction = (1 - RMSD_post/RMSD_prior) * 100
        print "   a_post : %.3g, b_post : %.3g" % (x_post[0], x_post[1])
        print "   RMSD_post : %.3g" % RMSD_post
        print "   RMSD reduction : %.2g%%" % reduction

    if "ok_plot" in task and task["ok_plot"].lower() in ok_flags:

        # Calculate data to plot
        SIF_sat = SIF * pftmask
        GPP_obs = OBS * pftmask
        GPP_sim = SIM * pftmask
        SIF_obs = a_post[:, np.newaxis, np.newaxis] * GPP_obs + b_post[:, np.newaxis, np.newaxis]
        SIF_sim = a_post[:, np.newaxis, np.newaxis] * GPP_sim + b_post[:, np.newaxis, np.newaxis]
        dGPP = GPP_obs - GPP_sim
        dSIF = SIF_sat - SIF_sim
        dGPP_inv = dSIF/a_post[:, np.newaxis, np.newaxis]

        if "ok_seasonal" in task and task["ok_seasonal"].lower() in ok_flags:
            formula = "$a_\mathsf{mth} \cdot \mathsf{GPP}_\mathsf{CTESSEL} + b_\mathsf{mth}"
        else:
            formula = "$%.3g \cdot \mathsf{GPP}_\mathsf{CTESSEL} %s %.3g$" % (a_post[0], ("+" if b_post[0]>0 else "-"), abs(b_post[0]))

        # Plot maps
        plots = [ dict(data = SIF_sat.mean(axis=0), filename = "sif_sat_%s" % pft, title = "%s\nSIF (SAT) [mW/m2/sr/nm]" % pftlabel),
                  dict(data = SIF_sim.mean(axis=0), filename = "sif_sim_%s" % pft, title = "%s\n%s" % (pftlabel, formula)) ]
        meshlon, meshlat = np.meshgrid(lon, lat)
        for item in plots: plotMap(item)

        # Plot time-series
        plots = [ dict(data = [calcMean(GPP_obs), calcMean(GPP_sim)], label = ["obs", "sim"], 
                               ylabel = "GPP", title = pftlabel, filename = "gpp_ts_%s.png" % pft),
                  dict(data = [calcMean(SIF_sat), calcMean(SIF_obs), calcMean(SIF_sim)], label = ["sat", "obs", "sim"],
                               ylabel = "SIF", title = pftlabel, filename = "sif_ts_%s.png" % pft),
                  dict(data = [calcMean(dGPP), calcMean(dGPP_inv)], label = ["$\Delta\mathsf{GPP}$", "$\Delta\mathsf{GPP}_\mathsf{inv}$"],
                               ylabel = "$\Delta$GPP", title = pftlabel, filename = "dgpp_ts_%s.png" % pft) ]
        for item in plots: plotTS(item)

        # Plot scatter-plots
        plots = [ dict(x = SIF_sat, y = SIF_sim, xlabel = "SIF SAT", ylabel = formula, title = pftlabel, filename = "sif_scat_%s.png" % pft),
                  dict(x = dGPP, y = dGPP_inv, xlabel = "$\Delta\mathsf{GPP}$", ylabel = "$\Delta\mathsf{GPP}_\mathsf{inv}$", title = pftlabel, filename = "dgpp_scat_%s.png" % pft) ]
        for item in plots: plotScat(item)

