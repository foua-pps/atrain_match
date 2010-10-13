#!/usr/local/bin/python
#
# File:
#	radiance_tb_tables.py
# Comments:
#	Routines to make tables of the channel 3b,4,5 radiance
#       to brightness temperature tables. Using the spectral
#       response functions of the individual satellites.
#

import string
import filtfunc
#print "filter function module: ",filtfunc.__file__

import Numeric
import math

PI = math.pi

h_planck = 6.6262*1e-34 # SI-unit = [J*s]
k_bolzmann = 1.3806e-23 # SI-unit = [J/K]
c_speed = 2.99793e8 # SI-unit = [m/s]

# From NOAA POD users guide (http://www2.ncdc.noaa.gov/docs/podug/html/c3/sec3-3.htm):
C1 = 1.1910659*1e-5 # mW/(m2-sr-cm-4)
C2 = 1.438833 # cm-K

TRUE=1
FALSE=0

INFO_SFIELDS = ("satellite",)

# -------------------------------------------------------------------
#
# INPUT:
#   noaa:       Noaa satellite id
#   temp:       Brightness temperature in Kelvin
#
# RETURN:
#   retv:       Wavenumber given in cm-1
#
def get_central_wavenumber(noaa_number,temp):
    if noaa_number == 14:
	if temp > 190.0 and temp < 230.0:
	    retv = 2638.652,928.2603,834.4496
	elif temp < 270.0:
	    retv = 2642.807,928.8284,834.8066
	elif temp < 310.0:
	    retv = 2645.899,929.3323,835.1647	
	elif temp < 330.0:
	    retv = 2647.169,929.5878,835.374
	else:
	    retv = 2647.169,929.5878,835.374 #Repeat 310-330/KG
#	    retv = None
    elif noaa_number == 12:
	if temp > 190.0 and temp < 230.0:
	    retv = 2632.713,920.0158,836.6847
	elif temp < 270.0:
	    retv = 2636.669,920.5504,837.0251
	elif temp < 310.0:
	    retv = 2639.61,921.0291,837.3641
	elif temp < 330.0:
	    retv = 2640.817,921.2741,837.5612
	else:
	    retv = 2640.817,921.2741,837.5612  #Repeat 310-330/KG
#	    retv = None

    # NOAA KLM uses instead a fix central wave number but adjusts the temperature in Planck function

    elif noaa_number == 15:
        retv = 2695.9743,925.4075,839.8979
    elif noaa_number == 16:
        retv = 2700.1148,917.2289,838.1255
    elif noaa_number == 17:
        retv = 2669.3554,926.2947,839.8246
    elif noaa_number == 18:
        retv = 2659.7952,928.1460,833.2532
    elif noaa_number == 19:
        retv = (2670.0, 928.9, 831.9)

    # METOP-A interpreted as noaa_number 2
    
    elif noaa_number == 2:
        retv = 2687.0,927.2,837.7

    else:
        raise NotImplementedError("No support for NOAA%d yet." % noaa_number)
    
    return retv

# -------------------------------------------------------------------
#
def radiance2tb_using_central_wavenumber(cwnum,rad):

    temp = C2*cwnum/math.ln(1+C1*cwnum*cwnum*cwnum/rad)    

    return temp

# -------------------------------------------------------------------
#
def tb2radiance_using_central_wavenumber(cwnum,temp):

    rad = C1*cwnum*cwnum*cwnum/(math.exp(C2*cwnum/temp) - 1)
    
    return rad

# -------------------------------------------------------------------
# -------------------------------------------------------------------
#
def tb2radiance_using_central_wavenumber_klm(cwnum,temp,noaa_number,channel):

    # Same as previous calculations but adjusting temperature into an
    # "efficient temperature" following NOAA's KLM guidelines
    # Also, here we only calculate for one individual channel

    if noaa_number == 15:
	A_coeff_3b = 1.621256
	A_coeff_4 = 0.337810
	A_coeff_5 = 0.304558
        B_coeff_3b = 0.998015
        B_coeff_4 = 0.998719
        B_coeff_5 = 0.999024
    elif noaa_number == 16:
	A_coeff_3b = 1.592459
	A_coeff_4 = 0.332380
	A_coeff_5 = 0.674623
        B_coeff_3b = 0.998147
        B_coeff_4 = 0.998522
        B_coeff_5 = 0.998363
    elif noaa_number == 17:
	A_coeff_3b = 1.702380
	A_coeff_4 = 0.271683
	A_coeff_5 = 0.309180
        B_coeff_3b = 0.997378
        B_coeff_4 = 0.998794
        B_coeff_5 = 0.999012
    elif noaa_number == 18:
	A_coeff_3b = 1.698704
	A_coeff_4 = 0.436645
	A_coeff_5 = 0.253179
        B_coeff_3b = 0.996960
        B_coeff_4 = 0.998607
        B_coeff_5 = 0.999057
    elif noaa_number == 19:
        A_coeff_3b = 1.67396
        A_coeff_4 = 0.53959
        A_coeff_5 = 0.36064
        B_coeff_3b = 0.997364
        B_coeff_4 = 0.998534
        B_coeff_5 = 0.998913
    elif noaa_number == 2: # Assuming METOP-02
	A_coeff_3a = 2.06699
	A_coeff_4 = 0.55126
	A_coeff_5 = 0.34716
        B_coeff_3b = 0.996577
        B_coeff_4 = 0.998509
        B_coeff_5 = 0.998947
    else:
        raise NotImplementedError("No support for NOAA%d yet." % noaa_number)

    if channel == '3b':
        A_coeff = A_coeff_3b
        B_coeff = B_coeff_3b
    elif channel == '4':
        A_coeff = A_coeff_4
        B_coeff = B_coeff_4
    elif channel == '5':
        A_coeff = A_coeff_5
        B_coeff = B_coeff_5

    temp_eff = A_coeff + B_coeff*temp

    rad = C1*cwnum*cwnum*cwnum/(math.exp(C2*cwnum/temp_eff) - 1)
    
    return rad

# -------------------------------------------------------------------
#
def blackbody(wln, temperature, repr="wavelength"):
    c = c_speed
    h = h_planck
    k = k_bolzmann
    
    c1=2*h_planck*c_speed*c_speed
    c2=h_planck*c_speed/k_bolzmann
    
    if repr=="wavelength":
	lamda=wln
	return c1/(lamda*lamda*lamda*lamda*lamda*(math.exp(c2/(lamda*temperature))-1))
    elif repr=="wavenumber":
	nu=wln
	return c1*nu*nu*nu/(math.exp(c2*nu/temperature)-1)
	#return 2*h*c*c*wln*wln*wln/(math.exp(wln*h*c/(k*temperature))-1)
    else:
	print "Error: Unable to derive Blackbody radiation!"
	return None

#---------------------------------------------------------------
def noaasat(satnum):
    satellite = "NOAA%.2d"%satnum
    
    ff = filtfunc.filfunc()

    # Make a list of tuples containing the wavelength and the response:
    channels = []
    for i in range(2):
	ch = Numeric.array(ff['%s_AVHRR-ch%d'%(satellite,i+1)])
	channels.append((ch[:,0],ch[:,1]))
	
    if satnum >= 15:
	ch = Numeric.array(ff['%s_AVHRR-ch3a'%(satellite)])
	channels.append((ch[:,0],ch[:,1]))
	ch = Numeric.array(ff['%s_AVHRR-ch3b'%(satellite)])
	channels.append((ch[:,0],ch[:,1]))
    else:
	ch = Numeric.array(ff['%s_AVHRR-ch3b'%(satellite)])
	channels.append((ch[:,0],ch[:,1]))
	
    for i in [3,4]:
	ch = Numeric.array(ff['%s_AVHRR-ch%d'%(satellite,i+1)])
	channels.append((ch[:,0],ch[:,1]))

    return channels

# -------------------------------------------------------------------
def integral(x,y):
    #print "Make the integral..."
    ndim = len(x)
    #print x
    #print y
    
    value = y[0] / 2.0 # "y[-1]" = 0
    dx = x[1] - x[0] # assumed to be equal to "x[0] - x[-1]"
    sum = value*dx
    for i in range(0,ndim-1):
	#print "i,dx,value,sum: ",i,dx,value,sum
	value = (y[i+1] + y[i]) / 2.0
	dx = x[i+1] - x[i]
	sum = sum + value*dx

    return sum

# -------------------------------------------------------------------
def read_table(filename):
    info = {}    
    sfields = INFO_SFIELDS
    
    fd = open(filename, "r")
    fd.readline()
    fd.readline()
    fd.readline()

    line = fd.readline()
    sl = string.split(line)

    # Read the header lines:
    txt = line[0:3]
    while (txt != "EOH") :
	sl = string.split(string.split(line,"#")[1],":")
	if len(sl) > 1:
	    key = string.split(sl[0])[0]
	    val = string.split(sl[1])[0]
	    print key, val
	    if key in sfields:
		#info[key] = val
		info[key] = string.strip(sl[1])
	    else:
		info[key] = eval(val)
		
	line = fd.readline()
	try:
	    txt = line[0:3]
	except:
	    print "WARNING: Error in header..."
	    pass

    tb = []
    rad3b = []
    rad4 = []
    rad5 = []
    lines = fd.readlines()

    for i in range(len(lines)):
	sl = string.split(lines[i]) 
	tb.append(string.atof(sl[0]))
	rad3b.append(string.atof(sl[1]))
	rad4.append(string.atof(sl[2]))
	rad5.append(string.atof(sl[3]))
    
    fd.close()

    info["tb_start"] = string.atof(string.split(lines[0])[0])
    info["tb_incr"] = (string.atof(string.split(lines[1])[0])-
		       string.atof(string.split(lines[0])[0]))
    info["tb_end"] = string.atof(string.split(lines[len(lines)-1])[0])
    info["ndim"] = len(lines)
    
    tb = Numeric.array(tb)
    rad3b = Numeric.array(rad3b)
    rad4 = Numeric.array(rad4)
    rad5 = Numeric.array(rad5)
    rad = (rad3b,rad4,rad5)
    
    return info,tb,rad

# -------------------------------------------------------------------
# Calculate the tb to radiance conversion tables for IR channels:
def make_table(filename, satnum, repr="wavelength"):

    tb, radiance = prepare_table(satnum) #Used in original code/KG
    #tb, radiance = fold(satnum,repr)      #Now tested to see difference

    lines = []
    lines.append("# AVHRR channel 3b,4,5 radiance as a function of Tb\n")
    lines.append("# 1. column = brightness temperature (K)\n")
    if repr=="wavenumber":
	lines.append("# 2.-4. column = radiance [W/(m2*sr) * m] ch3b;ch4;ch5\n")
    else:
	lines.append("# 2.-4. column = radiance [W/(m2*sr) / m] ch3b;ch4;ch5\n")
    lines.append("# satellite: NOAA %d\n"%satnum)
    lines.append("EOH\n")
    for i in range(len(tb)):
	lines.append("%5.2f %e %e %e\n"%(tb[i],
					 radiance[i][0],
					 radiance[i][1],
					 radiance[i][2]))
	
    # Write table:
    fd = open(filename, "w")
    fd.writelines(lines)
    fd.close()
    
    return 

# -------------------------------------------------------------------
def wavelength_integration(tb,wl,response):
    y=[]
    for i in range(len(wl)):
	y.append(blackbody(wl[i], tb) * response[i])
    y = Numeric.array(y)    
    return integral(wl,y)

# -------------------------------------------------------------------
def wavenumber_integration(tb,wl,response):
    x=[]
    y=[]    
    for i in range(len(wl)):	
	x.append(1./wl[i])
	#print "wavelength,wavenumber: ",wl[i],1./wl[i]
	yval_length = blackbody(wl[i], tb) * response[i] * (wl[i]*wl[i])
	yval = blackbody(1./wl[i], tb, "wavenumber") * response[i]
	#print "YVAL: ",yval,yval_length
	y.append(yval)
    
    x = Numeric.array(x)
    y = Numeric.array(y)
    y = Numeric.choose(Numeric.argsort(x),y)
    x = Numeric.sort(x)

    response = Numeric.array(response)
    
    res = integral(x,y)
    # Normalisation:
    norm = integral(x,response)
    #print "Norm: ",norm
    
    return res/norm

# -------------------------------------------------------------------
def fold(satnum, repr="wavelength"):

    # The unit representation:
    # if wavenumber representation then multiply by 1/sqr(lambda)
    # otherwise dont do anything
    rnumber = FALSE
    if repr=="wavenumber":
	rnumber = TRUE
	print "Wavenumber space!"
    else:
	print "Wavelength space!"
	
    # Go through every possible brightnes temperature:
    #tb_start = 200.0         #Not good for Antarctica!/KG
    tb_start = 150.0         #Changed for Antarctica!/KG 16 May 2008
    #tb_start = 274.8
    tb_incr = 0.1
    #tb_end = 275.0
    tb_end = 350.0
    
    result = []
    tb = tb_start
    tblist = []
    while (tb < tb_end):
	tblist.append(tb)
	res = []
	for ch in (3,4,5):
	    # Fold the Planck curve with the filter function:
	    if satnum < 15:
		(wl,response) = noaasat(satnum)[ch-1]
	    else:
		(wl,response) = noaasat(satnum)[ch]
	    #print wl,response
	    
	    wl = wl * 1.0e-6 # microns -> meters

	    if rnumber:
		val = wavenumber_integration(tb,wl,response)
	    else:
		# Wavenumber given in m-1:
		val = wavelength_integration(tb,wl,response)
		
	    res.append(val)

	#print tb,val
	
	result.append(res)
	tb = tb + tb_incr
    
    return tblist, result

# -------------------------------------------------------------------
def prepare_table(satnumber):
	
    # Go through every possible brightnes temperature:
    #tb_start = 200.0
    tb_start = 150.0         #Changed for Antarctica!/KG 16 May 2008
    tb_incr = 0.1
    tb_end = 350.0
    
    result = []
    tb = tb_start
    tblist = []
    while (tb < tb_end):
	tblist.append(tb)
	res = []
	cwnum = get_central_wavenumber(satnumber,tb)
	if cwnum==None:
	    tb = tb + tb_incr
	    result.append([-9,-9,-9])
	    continue
			  
	for chi in range(len(cwnum)):
	    rad = tb2radiance_using_central_wavenumber(cwnum[chi],tb)
	    # To get W/(m2*sr*m-1): multiply by 1e-5
	    rad=rad*1e-5
	    res.append(rad)
	    
	result.append(res)	    
	tb = tb + tb_incr

    
    return tblist, result

# -------------------------------------------------------------------
if __name__ == "__main__":

    # Testing:
    #make_table("ir_radiance_tab.noaa12", 12, "wavenumber")
    #make_table("ir_radiance_tab.noaa14", 14, "wavenumber")
    make_table("ir_radiance_tab.noaa12_extended", 12, "wavenumber")
    make_table("ir_radiance_tab.noaa14_extended", 14, "wavenumber")
    #make_table("ir_radiance_tab.noaa15_kg", 15, "wavenumber")
    #make_table("ir_radiance_tab_X.noaa14", 14, "wavenumber")
    #make_table("ir_radiance_tab.noaa15", 15, "wavenumber")

    #info,tb,rad = read_table("ir_radiance_tab.noaa15")

