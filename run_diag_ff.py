from matplotlib.pyplot import *
from numpy import *
import glob
#import xarray as xr

wd = '/glade/campaign/uwyo/wyom0178/dLWCF/'

def get_diag(i):
    varnms=['TGCLDLWP','ACTNL','FCTL','PRECC','PRECL','LHFLX','SWCF']
    data={'SST4Kf':[None],'PDf' :[None],'PDn':[None],'PIn':[None]}
    print(i)
    data['PDf']={}
    data['SST4Kf']={}
    data['PIn']={}
    data['PDn']={}
    for j in range(len(varnms)):
        #print(varnms[j])
        try:
            if (varnms[j]=='TGCLDLWP') | (varnms[j]=='PRECC') | (varnms[j]=='PRECL') | (varnms[j]=='LHFLX') | (varnms[j]=='SWCF'):
                data['PDf'][varnms[j]]=get_data4k(i,varnm=varnms[j],runval='PD')
                data['SST4Kf'][varnms[j]]=get_data4k(i,varnm=varnms[j],runval='SST4K')
            else:
                print('nond4k')
            data['PDn'][varnms[j]]=get_data(i,varnm=varnms[j],runval='PPE_rerun_ensemble_PD',type_run='PD')
            data['PIn'][varnms[j]]=get_data(i,varnm=varnms[j],runval='PPE_rerun_ensemble_PI',type_run='PI')
        except FileNotFoundError:
            print('FileNotFoundError')
        except OSError:
            print('OSError')
    FF_vals=calc_FF(data)
    
    ensn_str=str(i)
    if i<100:
        ensn_str='0'+ensn_str
        if i<10:
            ensn_str='0'+ensn_str
            
    savez(wd+ensn_str+'_coefs_ffv2',ensn=i,adj_lwppipd=FF_vals['adj_lwppipd'],dNdpipdadjnd=FF_vals['dNdpipdadjnd'],lon=FF_vals['lon'],lat=FF_vals['lat'],dLWP4K=FF_vals['dLWP4K'],map_dLWP4K=FF_vals['map_dLWP4K'],map_LWP_PD_f=FF_vals['LWPf_pd_map'],
          map_dSWCF4K=FF_vals['map_dSWCF4K'],map_SWCF_PD_f=FF_vals['SWCFf_pd_map'],
          map_dPE4K=FF_vals['map_dPE4K'],map_PE_PD_f=FF_vals['PEf_pd_map'],allow_pickle=True)
##SWCFf_pd_map; map_dSWCF4K
    return FF_vals
            
def calc_FF(data):
    LE=2264705. #J/kg= W/s/kg
    latr=[20,80]
    lat=data['PDf']['TGCLDLWP'][1]
    lon=data['PDf']['TGCLDLWP'][2]
    Ndpd=gm(data['PDn']['ACTNL'][0]/data['PDn']['FCTL'][0],data['PDn']['TGCLDLWP'][1],latr)
    Ndpi=gm(data['PIn']['ACTNL'][0]/data['PIn']['FCTL'][0],data['PDn']['TGCLDLWP'][1],latr)
    adjnd=Ndpd-Ndpi
    adj=gm(data['PDn']['TGCLDLWP'][0]-data['PIn']['TGCLDLWP'][0],data['PDn']['TGCLDLWP'][1],latr)
    dLWP4K=data['SST4Kf']['TGCLDLWP'][0]-data['PDf']['TGCLDLWP'][0]
    dSWCF4K=data['SST4Kf']['SWCF'][0]-data['PDf']['SWCF'][0]

    expn='PDf'
    PEpd=nanmean((data[expn]['PRECC'][0]+data[expn]['PRECL'][0])*1000.-data[expn]['LHFLX'][0]/LE,axis=0) ## m/s rain rate- 1 mm rain in a m2=kg m rain/1000=1kg LHFLX=w/m2 / W/s/kg= kg/m2/s
    ## kg/m2 = mm  1000mm/m -> m (1000mm/m) (kg/m2/mm)
    expn='SST4Kf'
    PE4k=nanmean((data[expn]['PRECC'][0]+data[expn]['PRECL'][0])*1000.-data[expn]['LHFLX'][0]/LE,axis=0)

    dPE4K=PE4k-PEpd
    latr=[-80,-40]    
    fbs=gm(dLWP4K,lat,latr)
    latr=[40,80]    
    fbn=gm(dLWP4K,lat,latr)
    return {'adj_lwppipd':adj,'dNdpipdadjnd':adjnd,'dLWP4K':(fbn+fbs)/2.,'map_dLWP4K':nanmean(dLWP4K,axis=0),'LWPf_pd_map':data['PDf']['TGCLDLWP'][0],'map_dPE4K':dPE4K,'PEf_pd_map':PEpd,'lat':lat,'lon':lon,
           'SWCFf_pd_map':data['PDf']['SWCF'][0],'map_dSWCF4K':nanmean(dSWCF4K,axis=0)}

def gm(data,lat,lr=[-90,90]):
    data=nanmean(data,axis=(0,2))
    lw=cos(pi*lat/180)
    ind=(lat<max(lr))&(lat>min(lr))
    return sum(data[ind]*lw[ind])/sum(lw[ind])

def get_data(fnum,varnm='FSDS',runval='PPE_rerun_ensemble_frT_PD.',type_run='PD'):
    fnums=str(fnum)
    if fnum<100:
        fnums='0'+fnums
        if fnum<10:
            fnums='0'+fnums
    WD='/glade/campaign/uwyo/wyom0124/daily/'+type_run+'/'+runval+'.'#001/atm/hist/PPE_rerun_ensemble_PD.001.cam.h0.2016-01.nc'
    fn=WD+fnums+'/atm/hist/'+runval+'.'+fnums+'.cam.h0.*.nc'
    if runval=='PPE_rerun_ensemble_frT_PD.':
        WD='/glade/scratch/cisong/archive/'
        fn=WD+runval+fnums+'/atm/hist/'+runval+fnums+'.cam.h0.*.nc'
    #WD='/glade/campaign/cgd/projects/ppe/cam_ppe/rerun_PPE_250/'+runval+'/'+runval+'_timeseries/'
    #fn=WD+'PPE_250_ensemble_'+runval+'.'+fnums+'/atm/hist/cc_PPE_250_ensemble_'+runval+'.'+fnums+'.h0.'+varnm+'.nc'
    #if fnum==175:
    #    fn='/glade/campaign/cgd/projects/ppe/cam_ppe/PPE_250/control/control_timeseries/PPE_250_ensemble.175/atm/hist/cc_PPE_250_ensemble.'+fnums+'.h0.'+varnm+'.nc'
    import netCDF4 as nc
    from datetime import datetime
    print(fn)
    f = nc.MFDataset(fn, 'r')
    tvar = 'time'
    tt = f.variables[tvar]
    latvar = 'lat'
    lonvar = 'lon'
    lat = f.variables[latvar][:]
    lon = f.variables[lonvar][:]
    Z = f.variables[varnm][:]
    lon[lon > 180] = lon[lon > 180]-360
    ind = argsort(lon)
    lon = lon[ind]
    Z = 1.*Z[:,:,ind]
    ind2 = argsort(lat)
    Z = 1.*Z[:,ind2,:]
    lat = lat[ind2]
    #import numpy.ma
    #Z=1.*Z.filled(fill_value=NaN)
    return [Z, lat, lon,tt]
    
def get_data4k(fnum,varnm='FSDS',runval='PD'):
    fnums=str(fnum)
    if fnum<100:
        fnums='0'+fnums
        if fnum<10:
            fnums='0'+fnums
    WD='/glade/campaign/cgd/projects/ppe/cam_ppe/rerun_PPE_250/'+runval+'/'+runval+'_timeseries/'
    fn=WD+'PPE_250_ensemble_'+runval+'.'+fnums+'/atm/hist/cc_PPE_250_ensemble_'+runval+'.'+fnums+'.h0.'+varnm+'.nc'
    if fnum==175:
        fn='/glade/campaign/cgd/projects/ppe/cam_ppe/PPE_250/control/control_timeseries/PPE_250_ensemble.175/atm/hist/cc_PPE_250_ensemble.'+fnums+'.h0.'+varnm+'.nc'
    import netCDF4 as nc
    from datetime import datetime
    #print(fn)
    f = nc.Dataset(fn, 'r')
    tvar = 'time'
    tt = f.variables[tvar]
    latvar = 'lat'
    lonvar = 'lon'
    lat = f.variables[latvar][:]
    lon = f.variables[lonvar][:]
    Z = f.variables[varnm][:]
    lon[lon > 180] = lon[lon > 180]-360
    ind = argsort(lon)
    lon = lon[ind]
    Z = 1.*Z[:,:,ind]
    ind2 = argsort(lat)
    Z = 1.*Z[:,ind2,:]
    lat = lat[ind2]
    import numpy.ma
    Z=1.*Z.filled(fill_value=NaN)
    return [Z, lat, lon,tt]

if __name__ == "__main__":    
       get_diag(int(sys.argv[1]))
