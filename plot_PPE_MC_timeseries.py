import sys
sys.path.append("/glade/u/home/geethma/phd_research_home/functions")
from imports import *
from lat_weight_mean import *

#### all saving functions in u/home/geethma/phd_research_home/TEST.ipynb

def plot_PPE_MC_timeseries(fnum):
    ### downwelling SW at TOA
    varnms = ['FSNTOA', 'FSUTOA']
    var_data = {}
    runval='PD'
    fnums = f"{fnum:03d}"
    WD='/glade/campaign/cgd/projects/ppe/cam_ppe/rerun_PPE_250/'+runval+'/'+runval+'_timeseries/'
    for varnm in varnms:
        fn=WD+'PPE_250_ensemble_'+runval+'.'+fnums+'/atm/hist/cc_PPE_250_ensemble_'+runval+'.'+fnums+'.h0.'+varnm+'.nc'
        if fnum==175:
            fn='/glade/campaign/cgd/projects/ppe/cam_ppe/PPE_250/control/control_timeseries/PPE_250_ensemble.175/atm/hist/cc_PPE_250_ensemble.'+fnums+'.h0.'+varnm+'.nc'
        #print(fn)
        f = xr.open_dataset(fn)
        Z = f[varnm]
        Z = Z.assign_coords(lon=(((Z.lon + 180) % 360) - 180))
        Z = Z.sortby('lon')
        
        var_data[varnm] = {'Z': Z.data}
        print(varnm)
    
    FSNTOA = var_data['FSNTOA']['Z']
    FSUTOA = var_data['FSUTOA']['Z']
    var = FSNTOA + FSUTOA
    unit = 'Wm-2'
    savefn = f'/glade/derecho/scratch/geethma/PPE_PD_inSolar/maps_{fnums}.npz'
    np.savez(savefn, 
             data=var, 
             lat=Z.lat,
             lon=Z.lon,
             time =Z.time,
             unit=unit)
    print(fnums+" saved.")




if __name__ == "__main__":    
       plot_PPE_MC_timeseries(int(sys.argv[1]))