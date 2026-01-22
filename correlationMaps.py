import sys
sys.path.append("/glade/u/home/geethma/phd_research_home/functions")
from imports import *
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from load_xarray import *
from lat_weight_mean import lat_weight_mean

# PNAS font sizes
title_fontsize = 8
label_fontsize = 7
tick_fontsize = 6
fig_dir = '/glade/derecho/scratch/geethma/figures_GFB/'

file = np.load('/glade/derecho/scratch/cisong/backup_FOR_wyom0124/variables_filtered/pe.npz')
PE_g_mean = xr.DataArray(file['data'], dims=['runs'], coords={'runs':file['runs']})  # Global dP/dLWP [s-1]
common_members = PE_g_mean.runs
common_members = common_members[common_members != 175]

var_list = ['dTSmap_gol', 'dLWPmap_gol', 'pe_maps_ol']  # 'dPE_maps_gol', 'dSWCREmap', 'dLWCREmap', 'wvp_maps_ol', 'w500_maps_ol', 
var_dict = {var: load_xarray(var) for var in var_list}

### Select only common runs for all datasets
var_keys = list(var_dict.keys())
for name in var_keys:
    print(name)
    common_members = np.intersect1d(common_members, var_dict[name]["runs"].values)

var_dict = {name: ds.sel(runs=common_members) for name, ds in var_dict.items()}

dGMT = lat_weight_mean(var_dict['dTSmap_gol'])[0]  # ∆GMT [K]
pe_maps_ol = var_dict['pe_maps_ol']*86400  # PD P-E [kgm-2day-1]
dLWP_g = var_dict['dLWPmap_gol']/dGMT  # ∆LWP/∆GMT [kgm-2K-1] 

# all pressure all cloud feedbacks
ds_kernel = {}
for mem_run in common_members:
    pk = 'all_pressure'
    tk = 'all'

    mem = str(mem_run).zfill(3)
    path = (
        f"/glade/derecho/scratch/travisa/CAM6_kernels/"
        f"member_number_{mem}/"
        f"CAM6_PPE_CRK_decomp_{pk}_{tk}_cloud.nc"
    )

    if not os.path.exists(path):
        print(f"Skipping member {mem}: file not found")
        continue

    ds_kernel[mem + pk + tk] = xr.open_dataset(path)

# now for kernel feedback take the global, subtropical, and midltitude means
def trav_fbvar(fbvar, region_dic, ds_kernel=ds_kernel, dGMT=dGMT, common_members=common_members):
    '''return each cloud feedback variable separated into regions.
       fbvar is the feedback component name in strings.'''
    # region_dic = {
    #     'global':[0,90],
    #     'tropics':[0,20],
    #      'subtropical':[20,50],
    #      'midlatitude':[50,90],
    #     # 'ten_thirty': [10,30]
    # }
    regions = list(region_dic.keys())
    
    SWkernel_regmean = {}
    for reg in (regions):
        SWkernel_list = []
        for mem_run in (common_members):
            mem = str(mem_run).zfill(3)
            SWkernel = ds_kernel[mem+'all_pressure'+'all'][fbvar].mean('time')
            # SWkernel_c = regional_concat(SWkernel,region_dic[reg])
            # SWkernel_reg = regional_mean(SWkernel_c,ocean_only=False)
            SWkernel_list.append(SWkernel)
    
        SWkernel_all = xr.concat(SWkernel_list,'runs').assign_coords({'runs':common_members})
        SWkernel_regmean[reg] = SWkernel_all.load()/dGMT
    return SWkernel_regmean

# cloud feedback components
region_dic = {
    'global' : [0, 90]
}
LWcld_tot_fb_global = trav_fbvar('LWcld_tot', region_dic)
SWcld_tot_fb_global = trav_fbvar('SWcld_tot', region_dic)
totdCRE_g = LWcld_tot_fb_global['global'] + SWcld_tot_fb_global['global']
totdCRE_g = totdCRE_g.sel(runs=common_members)

dVar_name = ['∆LWP/∆GMT * PD P-E', '$λ_{cld}$  * PD P-E']

PE_g_mean = PE_g_mean.sel(runs=common_members) 
grad = [PE_g_mean.data]
grad = np.array(grad)
grad_name = ['P/LWP']

lat_sh = dLWP_g.lat.shape[0]
lon_sh = dLWP_g.lon.shape[0]

nums = [[0, 0], [0, 1]]
title_num = ['(a) ', '(b) ', '(c) ']
# dVar_name = ['∆LWP/∆GMT', '$λ_{cld}$']
# statistical significance from Travis
dir = '/glade/derecho/scratch/travisa/Statistical_testing_for_Geethma/'
stipple = ['dPLWP_ddelLWP_statistical_test_0.3alph_1.nc', 'gpe_totcld_statistical_test_0.3alph.nc']
# correlations from my saved npz
sdir = "/glade/u/home/geethma/phd_research_home/globalFB/data/gpe_totcre_correlation_pvalue_maps.npz"
print(sdir)
data = np.load(sdir)
lon, lat, corr_00, pval_00, corr_01, pval_01 = (
    data["lon"],  #longitude
    data["lat"],  #latitude
    data["correlation_0_0"],  #correlation of P/LWP with ∆LWP
    data["pvalue_0_0"],  #p-value of P/LWP with ∆LWP
    data["correlation_0_1"],  #correlation of P/LWP with λCRE
    data["pvalue_0_1"],  #p-value of P/LWP with λCRE
)
sdir = "/glade/u/home/geethma/phd_research_home/globalFB/data/gpe_totcld_correlation_pvalue_maps.npz"
print(sdir)
data = np.load(sdir)
cld_correlation = data["correlation"]
corrlist = [corr_00, cld_correlation]

# -------------- Plot correlation of global precipitation eff with ∆LWP and lambda scaled by PD P-E in maps ----------------- 
fig, axs = plt.subplots(2, 1, figsize=(4.49, 5), subplot_kw={'projection': ccrs.Robinson()}, constrained_layout=True)  # Adjusted height to maintain vertical layout: horizontal layout:(4.49, 2.8))

for ax_n, ax in enumerate(axs.flat):
    ax.set_facecolor("darkgrey")
    g = nums[ax_n][0]
    d = nums[ax_n][1]
    
    correlation = corrlist[ax_n]
    stipdata = xr.open_dataset(dir+stipple[ax_n])
    stipcorr = stipdata['Wilks_binary']*correlation

    cbar_img = ax.pcolormesh(lon, lat, stipcorr,
                             transform=ccrs.PlateCarree(), cmap='seismic',
                             vmin=-1, vmax=1)

    ax.add_feature(cfeature.LAND, facecolor='white')
    ax.add_feature(cfeature.GSHHSFeature(scale='auto', edgecolor='black'))
    ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
    ax.set_title(title_num[ax_n] + 'Global mean ' + grad_name[g] + ' correlation with ' + dVar_name[d],
                 fontsize=title_fontsize, loc='left')
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='black', alpha=0.8, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'fontsize': label_fontsize}
    gl.ylabel_style = {'fontsize': label_fontsize}

# Single horizontal colorbar below
cbar = fig.colorbar(cbar_img, ax=axs, orientation='vertical', fraction=0.05, pad=0.07)
cbar.set_label('Correlation', fontsize=label_fontsize)
cbar.ax.tick_params(labelsize=tick_fontsize)

plt.savefig(fig_dir+'gPE_λcldpdPE_∆LWPpdPE_correlation_statistical.png',
            bbox_inches='tight', facecolor='white', dpi=600)

# -------------- Plot ∆LWP/∆GMT, lambda, and PD P-E in maps ----------------- 
lon = dLWP_g.lon
lat = dLWP_g.lat
dVAR_map_o = [dLWP_g[0].data, np.transpose(totdCRE_g[0].data), pe_maps_ol[0].data]
dVar_name = ['∆LWP/∆GMT', '$λ_{cld}$', 'PD P-E']
fig, axs = plt.subplots(3, 1, figsize=(4.49, 9), subplot_kw={'projection': ccrs.Robinson()}, constrained_layout=True)  # Adjusted height to maintain vertical layout: horizontal layout:(4.49, 2.8))

for ax_n, ax in enumerate(axs.flat):
    print(dVar_name[ax_n])
    ax.set_facecolor("darkgrey")
    # g = nums[0][ax_n]
    # d = nums[0][ax_n]
    plot_var = dVAR_map_o[ax_n]
    cbar_img = ax.pcolormesh(lon, lat, plot_var,
                             transform=ccrs.PlateCarree(), cmap='seismic')

    ax.add_feature(cfeature.LAND, facecolor='white')
    ax.add_feature(cfeature.GSHHSFeature(scale='auto', edgecolor='black'))
    ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
    ax.set_title(title_num[ax_n] + f'{dVar_name[ax_n]} for the CAM6 default run',
                 fontsize=title_fontsize, loc='left')
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='black', alpha=0.8, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'fontsize': label_fontsize}
    gl.ylabel_style = {'fontsize': label_fontsize}

    cbar = fig.colorbar(
        cbar_img,
        ax=ax,
        orientation='horizontal',
        fraction=0.05,
        pad=0.07
    )
    cbar.set_label(dVar_name[ax_n], fontsize=label_fontsize)
    cbar.ax.tick_params(labelsize=tick_fontsize)
    max_abs_value = np.nanmax(np.abs(plot_var))
    cbar_img.set_clim(-max_abs_value, max_abs_value)

plt.savefig(fig_dir+'λcld_pdPE_∆LWP_maps.png',
            bbox_inches='tight', facecolor='white', dpi=600)


