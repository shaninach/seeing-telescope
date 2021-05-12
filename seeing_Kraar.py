###########################
# this code calculates the mean seeing (FWHM) of a star (e.g. Polaris) in every image in a given folder. 
# correct path to the specific image folder, then run. 
##########################

# importing needed libraries
import numpy as np
import matplotlib.pyplot as plt  
from astropy.modeling import models, fitting
from astropy.io import fits
import os
import pandas as pd

# definitions
path = r"C:\Users\shaninac\Desktop\studies\weizmann\Rotations\third\DIMM\DIMM_Test_2020-09-01"  # image folder path
os.chdir(path)  # define the path as current directory
conC = 0.65  # conversion coefficient of the specific telescope
seeing_values = []  # creating an empty list for seeing values
Date = []
Hour = [] 

# functions
def seeing(f, conC):  # seeing function, finds FWHM, then by conversion constant, finds seeing in arcsec
    """ input = f: an object of fitting to 2D gaussian, conC: the conversion constant from pixel to arcsec.
    output =  seeing value in arcsec"""
    x_std = f.x_stddev[0]
    y_std = f.y_stddev[0]
    std_all = [x_std, y_std]
    return np.mean(std_all) * 2.355 * conC
    # multiply of std average by 2.355 and by conversion coefficient to get seeing

# main loop - output table of seeing per image 
for file in os.listdir(path):
    hdul = fits.open(file)  # Open image
    data = hdul[0].data  # data = pixel brightness
    cents = np.where(data == np.max(data))  # Find brightest pixel
    if len(cents[1])>1 or len(cents[0])>1: # if there's more then one brightest pixel, choose the first
        xc = cents[1][0]
        yc = cents[0][0]
    else:
        xc = int(cents[1])
        yc = int(cents[0])
    
    bb = 30  # half box size (around brightest pixel)
    box = data[yc - bb:yc + bb, xc - bb:xc + bb]
    yp, xp = box.shape # yp 1D (#rows) xp 2D (#columns)
    y, x, = np.mgrid[:yp, :xp]  # Generate grid of same size like box to put the fit on
    f_init = models.Gaussian2D()  # fit model 
    fit_f = fitting.LevMarLSQFitter()  # fitting function 
    f = fit_f(f_init, x, y, box)  # Fit the model to your data (box)

    seeing_values.append(seeing(f, conC))  
    
    info = hdul[0].header  # extracts image data
    time = info[7] # extract image time 
    date = time [0:10] 
    hour = time[11:16]
    Date.insert(len(Date),date) # add date value to the end
    Hour.insert(len(Hour),hour) # add hour value to the end
    
    
    d = {'Date': Date , 'Hour': Hour, 'seeing': seeing_values}
    df = pd.DataFrame(data=d)
    
    seeing_av = df[['Hour', 'seeing']].groupby('Hour').mean()


"""
# check individual images - optional 

hdul = fits.open("DIMM Test - 20200901-0906.fit")  # Open image
data = hdul[0].data  # data = pixel brightness
cents = np.where(data == np.max(data))  # Find brightest pixel
if len(cents[1])>1 or len(cents[0])>1:
    xc = cents[1][0]
    yc = cents[0][0]
else:
    xc = int(cents[1])
    yc = int(cents[0])
bb = 30  # half box size (around brightest pixel)
box = data[yc - bb:yc + bb, xc - bb:xc + bb]
yp, xp = box.shape
y, x, = np.mgrid[:yp, :xp]  # Generate grid of same size like box to put the fit on
f_init = models.Gaussian2D()  # Declare what function you want to fit to your data
fit_f = fitting.LevMarLSQFitter()  # Declare what fitting function you want to use
f = fit_f(f_init, x, y, box)  # Fit the model to your data (box)

std = [f.x_stddev[0], f.y_stddev[0]]
seeing = np.average(std) * 2.355 * conC
print('seeing is: ', seeing)

# if the error "The fit may be unsuccessful" comes up, check the fit_info, 'message': 
fit_f.fit_info


# Plot the data with the best-fit model
plt.figure(figsize=(8, 2.5))
plt.subplot(1, 3, 1)
plt.imshow(box)
plt.title("Data")
plt.subplot(1, 3, 2)
plt.imshow(f(x, y))
plt.title("Model")
plt.subplot(1, 3, 3)
plt.imshow(box - f(x, y))
plt.title("Residual")
plt.show()
"""

