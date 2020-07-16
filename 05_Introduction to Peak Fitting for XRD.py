#!/usr/bin/env python
# coding: utf-8

# # Introduction to Peak Fitting in Python

# So far you have been introduced to the python programming language, working with arrays, importing data, and doing simple curve fitting. In this tutorial we are going to apply those skills to peak fitting data to determine lattice spacing and peak intensity.
# 
# After this tutorial you should be able to:
# 
#     Determine lattice spacing and intensity using gaussian peak fitting
#     Find lattice parameters to fit Data to the (100) peak and the (200) peak 
# 
# We'll be using the same imports as our previous tutorials with one important addition. We'll be adding a theta to q function, that makes it very easy for us to import to look at our data with respect to Q

# In[4]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from scipy.optimize import curve_fit
import pandas as pd


# # Importing the Data
# For this tutorial I already have a datafile of xrd data. This datafile is stacked .csv where the first column is the value of 2-theta and subsequent columns are diffracted intensity for a variety of samples.
# To import the data we're going to define a function "csv_to_np" which will use the read_csv function built into pandas to pull the data, and then convert it to a numpy array since we're already familiar with working with those.

# First, given that It's worthwhile to make a function that will convert our data with respect to Q go ahead and use the two_to_q function built to pull the data, and takes in an array of 2-theta values and a X-ray wavelength 𝜆, and returns an array of Q values.

# In[5]:


def two_to_q(two_theta, wave):
    #two_theta is a 1D array of two_theta angles
    #wave is the X-ray energy in angstroms
    rad_theta = two_theta/2*np.pi/180
    q = 4*np.pi*np.sin(rad_theta)/wave
    return q


# In[6]:


def csv_to_np(filename):
    data = pd.read_csv(filename)
    return(np.array(data))

perov = csv_to_np('D1_MAPBIBr2_Xraydeg.csv')


# In[ ]:


# Exercise 5
q_1 = 2.04 # This will be the lower limit for Q we'll consider
q_2 =2.15
limit1 = find_nearest(q, q_1) #First our lower limit
limit2 = find_nearest(q, q_2) #And of our higher limit
q_sub = q[limit1:limit2] # We'll reduce the domain of Q
perov_sub = perov[limit1:limit2,1:]

q_linear = np.hstack((q_sub[0:10], q_sub[-11:-1])) #I'm taking the starting and ending values
perov_linear = np.hstack((perov_sub[0:10,0], perov_sub[-11:-1,0])) #We'll use these to fit a straight line
slope, intercept = np.polyfit(q_linear, perov_linear, 1) #Do linear fit
back = slope*q_sub+intercept #Create background array of the form Background = Ax+B
#print (back)

#plt.plot(q_sub,perov_sub[:,1], 'r-',label='$MAPbIBr_2$')# plot minus background 
#Let's begin by getting our data ready to analyze
perov_fit = perov_sub[:,0]-back #We'll begin by subtracting the background we calculated for this piece of data

#Now let's define a function we'll want to fit to - this is analagous to the "straight-line-model" from tutorial 03
#We'll call our function gaussian and it will calculate the expression described above
def gaussian(x, a, b, c): 
    return a*np.exp(-(x - b)**2/(2*c**2))

#We'll also give an initial guess for our fits based off of a visual interpretaion of our data
p0 = [45, 2.4, 2] #(height, center, width)

#Use scipy.optimize.curve_fit to fit our desired data
popt, pcov = curve_fit(gaussian, q_sub, perov_fit, p0)

#To confirm our fits it's always nice to plot our model versus our data.
plt.figure(figsize=(8,6)) #make plot larger
plt.plot(q_sub,perov_fit,'r-', label='$MAPbIBr_2$') #plot subfield of data
plt.plot(q_sub,gaussian(q_sub, *popt),'b--', label='Model') #plot best fit
plt.xlabel('Q [$\AA^{-1}$]',size=12) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=12)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner

