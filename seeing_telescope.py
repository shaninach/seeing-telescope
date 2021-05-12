###########################
# this code calculates the mean seeing (FWHM) of a star (e.g. Polaris) in every image in a given folder. 
# define the path to the folder containing the fits images (path), and defint the date of measurements (start, end). 
##########################

# importing needed libraries
import numpy as np
from astropy.modeling import models, fitting
from astropy.io import fits
import os
import pandas as pd
from matplotlib import pyplot as plt 

# definitions
path = r"D:\Master\ZWO\trial_lone"  # image folder path
os.chdir(path)  # define the path as current directory
pscale = 0.604  # plate scale (arcsec/pixel)
darki = fits.open(r'D:\Master\ZWO\2021-04-12_19_27_48Z_dark\dark.fit')
dark = darki[0].data

# functions
def seeing(path,pscale): 
    """ path = path of folder where images are saved 
        pscale = plate scale of telescope (arcsec/pixel).
        output: list of seeing values per image """
    Seeing = []
    Date = []
    Hour = []
    Exposure = []
    for image in os.listdir(path): 
        hdul = fits.open(path+ r'\\' +image)  # Open image
        data = hdul[0].data  # data = pixel brightness
        data = data.astype(float)
        data -= dark
        data -= int(np.median(data))
        data /= np.std(data)
        info = hdul[0].header  # extracts image data
        time = info['DATE-OBS'] # extract image time 
        date = time [0:10] 
        hour = time[11:19]
        exposure = info['EXPOSURE']
        
        cents = np.where(data == np.max(data))  # Find brightest pixel
        hb = 30   # half box size of fit to gausian (around brightest pixel)
        yc = int(cents[0][0]) # index of centroid row 
        xc = int(cents[1][0]) # index of centroid column 
        fitbox = data[(yc-hb):(yc+hb), (xc-hb):(xc+hb)] # box limits: data[yc+-hb, xc+-hb]
        yp, xp = fitbox.shape
        y, x, = np.mgrid[:yp, :xp]  # Generate grid of same size like box to put the fit on
        f_init = models.Gaussian2D()  # Declare what function you want to fit to your data
        fit_f = fitting.LevMarLSQFitter()  # Declare what fitting function you want to use
        f = fit_f(f_init, x, y, np.array(fitbox))  # Fit the model to your data (box)
        Seeing.append(np.mean([f.x_fwhm, f.y_fwhm])*pscale) 
        Date.append(date)
        Hour.append(hour)
        Exposure.append(exposure)
        d = {'Date': Date , 'Hour': Hour, 'Seeing': Seeing, 'Exposure': Exposure}
        df = pd.DataFrame(data=d)
    return df

d = {'Date': [] , 'Hour': [], 'Seeing': [], 'Exposure': []}
df = pd.DataFrame(d) 
for folder in os.listdir(path):
    df = pd.concat([df,seeing(path+'\\'+folder, pscale)])

df['time'] = pd.to_datetime(df.Date + ' ' + df.Hour, format= '%Y-%m-%d %H:%M:%S')
df = df.set_index('time').drop('Date',1)
df = df[df['Seeing'] <= df.Seeing.mean()+ 1*np.std(df.Seeing)]
df = df[df['Seeing'] >= 0.2]
    
seeing_av = df.resample('1T').median()

print('%i images were taken. \nMean seeing %.2f". \nExposure time %.4f seconds' %(len(df.Seeing),df.Seeing.mean(),df.Exposure[0]))


"""
# check individual images - optional 
# if the error "The fit may be unsuccessful" comes up, check the fit_info, 'message': fit_f.fit_info

from matplotlib import pyplot as plt
path = r"D:\Master\ZWO\2021-04-11_21_56_07Z"  # image folder path
os.chdir(path)  # define the path as current directory
image = 'Light_ASIImg_0.01sec_Bin1_20.8C_gain179_2021-04-11_215618_frame0011.fit'
hdul = fits.open(image)  # Open image
data = hdul[0].data  # data = pixel brightness  
data = data.astype(float)

data -= np.median(data)
data = data/np.std(data)
plt.imshow(data)
plt.colorbar()

box = data[300:400,200:330]
plt.imshow(box)
plt.colorbar()

cents = np.where(data == np.max(data)) # Find brightest pixel
hb = 30   # half box size of fit to gausian (around brightest pixel)
yc = int(cents[0][2]) # index of centroid row 
xc = int(cents[1][2]) # index of centroid column 
fitbox = data[(yc-hb):(yc+hb), (xc-hb):(xc+hb)] # box limits: data[yc+-hb, xc+-hb]
plt.imshow(fitbox)

yp, xp = fitbox.shape
y, x, = np.mgrid[:yp, :xp]  # Generate grid of same size like box to put the fit on
f_init = models.Gaussian2D()  # Declare what function you want to fit to your data
fit_f = fitting.LevMarLSQFitter()  # Declare what fitting function you want to use
f = fit_f(f_init, x, y, np.array(fitbox))  # Fit the model to your data (box)

# Plot the data with the best-fit model
plt.subplot(1, 3, 1)
plt.imshow(fitbox)
plt.title("Data")
plt.subplot(1, 3, 2)
plt.imshow(f(x, y))
plt.title("Model")
plt.subplot(1, 3, 3)
plt.imshow(fitbox - f(x, y))
plt.title("Residual")
plt.show()
"""

hdul = fits.open(r'D:\Master\ZWO\trial_lone\2021-04-11_20_38_27Z_lone\Light_ASIImg_0.001sec_Bin1_20.0C_gain179_2021-04-11_203830_frame0002.fit')  # Open image
info = hdul[0].header  # extracts image data
info['DATE-OBS'] # extract image time 
