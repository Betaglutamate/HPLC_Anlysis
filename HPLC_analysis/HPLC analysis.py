import numpy as np
import scipy
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt


def printdata(filename):            #generates raw HPLC data from a .txt file
    f = open(filename) #opens data file

    #finds the start and end of the lc data file (courtesy of Brian)
    setlcstart = False
    setrtimestart = True
    setend = True
    for linenum,line in enumerate(f):
        if setlcstart == False:
            if "[LC Chromatogram" in line:
                setrtimestart = False
        if setrtimestart == False:
            if "R.Time (min)" in line:
                start_ = linenum
                setend = False
        if setend == False:
            if "[" in line:
                end_ = linenum
                break

    data = pd.read_csv(filename,index_col=None,delimiter="\t", names = ['R.Time (min)','Intensity'],skiprows=start_+1,nrows=end_-start_-2)
    data = data.rename(index=str, columns={"R.Time (min)": "t", "Intensity": "Intensity"})

    time = np.array(data.t.values)
    intensity = np.float64(data.Intensity)

    return (time,intensity)


def correct_baseline(data,tmin,tmax):       #removes baseline from HPLC data, so that a peak starts from and ends at 0
    baseline_pre = round(np.mean(data[1][np.where(data[0]==tmin)]))
    baseline_post = round(np.mean(data[1][np.where(data[0]==tmax)]))
    mb = (baseline_post - baseline_pre)/(tmax-tmin)
    qb = baseline_pre - tmin*(baseline_post - baseline_pre)/(tmax-tmin)
    sliced_times = data[0][np.where(np.logical_and(data[0]>=tmin,data[0]<=tmax))]
    sliced_data = data[1][np.where(np.logical_and(data[0]>=tmin,data[0]<=tmax))]
    corrected_data = sliced_data-[qb + mb*sliced_times[i] for i in range(len(sliced_times))]
    return [sliced_times,corrected_data]

def compute_area(data):           #computes the peak area from HPLC data
    return np.trapz(data[1],x=data[0])
#%%


calib_folder = "./2022/2022-08/08-12/2022-08-10-4B03-1A01-chitinase-activity-05um-Cm/HPLC-data/"
calib_files =    ["GlcNAc 0.0mM.txt",        #list containing the names of the files containing the HPLC data
                  "GlcNAc 0.1mM.txt",
                  "GlcNAc 0.2mM.txt",
                  "GlcNAc 0.3mM.txt",
                  "GlcNAc 0.4mM.txt",
                  "GlcNAc 0.5mM.txt",
                  "GlcNAc 0.6mM.txt",
                  "GlcNAc 0.7mM.txt",
                  "GlcNAc 0.8mM.txt",
                  "GlcNAc 0.9mM.txt",
                  "GlcNAc 1.0mM.txt"]

comp_std_concs = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]     #concentrations of the compound of interest in the standard samples

calib_raw_data = [printdata(calib_folder+file) for file in calib_files]

calib_ref_peaks = [correct_baseline(printdata(calib_folder+file),10.25,12.25) for file in calib_files]

calib_ref_peak_areas = [compute_area(calib_ref_peaks[i]) for i in range(len(calib_ref_peaks))]

calib_norm_factors = [calib_ref_peak_areas[i]/122902.72640416668 for i in range(len(calib_ref_peaks))]

data = [correct_baseline(calib_raw_data[i],19.4,21) for i in range(len(calib_ref_peaks))]

calib_peak_areas = [np.trapz(data[i][1],x=data[i][0])/calib_norm_factors[i] for i in range(len(calib_ref_peaks))]

GlcNAc_q,GlcNAc_m = scipy.optimize.curve_fit(lambda t,a,b: a+b*t, comp_std_concs, calib_peak_areas)[0]

#%%

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(12,7))
plt.title('HPLC data')
plt.xlabel('Time (minutes)')
plt.ylabel('Intensity (mV)')
plt.xlim(18,22)
# plt.xlim(25,40)
plt.ylim(-100,3000)
colormap = plt.get_cmap('Blues')
for i in range(len(calib_ref_peaks)):
    plt.scatter(calib_raw_data[i][0],calib_raw_data[i][1],s=20,facecolors=colormap(0.0+i/len(calib_ref_peaks)),linewidth=2.5)
plt.show()
#%%

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(14,7))
plt.title('HPLC calibration - GlcNAc')
plt.xlabel('GlcNAc concentration (mM)')
plt.ylabel('HPLC peak area')
colormap = plt.get_cmap('Blues')
for i in range(len(calib_files)):
	plt.scatter(comp_std_concs[i],calib_peak_areas[i],s=100,facecolors='none',edgecolors=colormap(0.25+i/len(calib_files)),linewidth=3)
plt.plot(np.linspace(comp_std_concs[0],comp_std_concs[-1],100),GlcNAc_q+GlcNAc_m*np.linspace(comp_std_concs[0],comp_std_concs[-1],100),color=colormap(0.5),label='Slope = '+str(round(GlcNAc_m)))
plt.legend()
plt.show()
#plt.savefig('HPLC-Glucose-calibration.jpg',dpi=300,bbox_inches="tight")
