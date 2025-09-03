import numpy as np

##########################################################
# Finley & Matt (2017/2018) Torque Formulations ##########
# Version 1.0 : 12th April 2018 ##########################
##########################################################
# Function:                                              #
# FM18_surf-Takes polar strengths of the lowest harmonic #
#	    components and predicts an average Alfv\'en  #
#	    radius which produces a braking torque on 	 #
#	    the star.                                    #
##########################################################
# Input Parameters:                                      #
# Bdip- Polar Dipole Strength in Gauss                   #
# Bquad- Polar Quadrupole Strength in Gauss              #
# Boct- Polar Octupole Strength in Gauss                 #
# Mdot- Mass Loss Rate in Solar per Year                 #
# Period- Rotation Period in Days                        #
# Mass- Stellar Mass in Solar Units                      #
# Radius- Stellar Radius in Solar Units                  #
# Info- Print Topology and Rotation Information          #
##########################################################

def FM18_surf(Bdip, Bquad, Boct, Mdot, Period, Mass, Radius, Info):

	#########################################################
	#Stellar Parameters:                                    #
	#########################################################
	R_star = 6.96e10*(Radius)
	M_star = 1.99e33*(Mass)
	Omega_star = 2*np.pi/Period/24/3600
	v_esc = np.sqrt( 2*(6.67e-8)*M_star/R_star )
	f = Omega_star*R_star/v_esc*np.sqrt(2)
	if (f < 0.1 and Info==True):
		print('SLOW ROTATOR REGIME, f = '+str(f))
	elif (Info==True):
		print('FAST ROTATOR REGIME, f = '+str(f))
	Btot = abs(Bdip) + abs(Bquad) + abs(Boct)
	Mdot = Mdot*( (2e-14)*(1.99e33)/(365.25*24*3600) )
	Upsilon = Btot**2*R_star**2/Mdot/v_esc
	Rdip = Bdip/Btot
	Rquad = Bquad/Btot
	Roct = Boct/Btot


	#########################################################
	# Derived Simulation Parameters:                        #
	#########################################################
	K1 = 1.53
	m1 = 0.229
	K2 = 1.70
	m2 = 0.134
	K3 = 1.80
	m3 = 0.087
	Kvr = 0.2

	#########################################################
	# Calculate Average Alfv\'en Radius using equation (17) #
	#########################################################
	RA_dip = K1*(Rdip)**(2*m1)*(Upsilon/np.sqrt(1+f**2/Kvr**2))**(m1)
	RA_quad = K2*(Rquad+Rdip)**(2*m2)*(Upsilon/np.sqrt(1+f**2/Kvr**2))**(m2)
	RA_oct = K3*(Roct+Rquad+Rdip)**(2*m3)*(Upsilon/np.sqrt(1+f**2/Kvr**2))**(m3)
	##################### Max Function ######################
	if ( RA_dip>RA_quad and RA_dip>RA_oct and Rdip>0. ):
		RA = RA_dip
		if (Info==True):
			print('DIPOLE DOMINATED')
	elif (RA_quad>RA_dip and RA_quad>RA_oct and Rquad>0.):
		RA=(RA_quad)
		if (Info==True):
			print('MULTIPOLAR DOMINATED')
	else:
		RA=(RA_oct)
		if (Info==True):
			print('MULTIPOLAR DOMINATED')
	#########################################################

	#########################################################
	#Calculate Angular Momentum Loss Rate, equation (14)    #
	#########################################################
	Torque = Mdot*Omega_star*R_star**2*np.power(RA,2)

	return Upsilon, RA, Torque


##########################################################
# Metcalfe et al. (2025) Torque Uncertainty###############
# Version 1.0 :  3rd June 2025  ##########################
##########################################################
# Function:                                              #
# torqueUncert - Calculates the FM18 torque for Star     #
#       and evaluates the uncertainty due to variation   #
#		of the input parameters.		 #
##########################################################
# Input Parameters:                                      #
# star - Dictionary* of Star Properties including        #
#        uncertainties for each parameter                #
# Dictionary* - Name, Bdip, Bquad, Boct, Mdot, Period,   #
#       Mass, and Radius.                                #
##########################################################

def torqueUncert(star,ndec=5):

	#########################################################
	#Star Wind Torque:                                     #
	#########################################################
	print('[%s, Bd=%.2f, Bq=%.2f, Bo=%.2f]'%(star["Name"],star["Bdip"],star["Bquad"],star["Boct"]))
	Upsilon, RA, Torque = FM18_surf(Bdip=star["Bdip"], Bquad=star["Bquad"], Boct=star["Boct"], Mdot=star["Mdot"], Period=star["Period"], Mass=star["Mass"], Radius=star["Radius"], Info=False)
	#########################################################

	if ('N/A' in star["un_Bdip"]):
		star["Bdip_min"] = star['Bdip']
		star["Bdip_max"] = star['Bdip']
	else:
		star["Bdip_min"] = star['Bdip'] + star["un_Bdip"][0]
		star["Bdip_max"] = star['Bdip'] + star["un_Bdip"][1]

	if ('N/A' in star["un_Bquad"]):
		star["Bquad_min"] = star['Bquad']
		star["Bquad_max"] = star['Bquad']
	else:
		star["Bquad_min"] = star['Bquad'] + star["un_Bquad"][0]
		star["Bquad_max"] = star['Bquad'] + star["un_Bquad"][1]

	if ('N/A' in star["un_Boct"]):
		star["Boct_min"] = star['Boct']
		star["Boct_max"] = star['Boct']
	else:
		star["Boct_min"] = star['Boct'] + star["un_Boct"][0]
		star["Boct_max"] = star['Boct'] + star["un_Boct"][1]

	if ('N/A' in star["un_Mdot"]):
		star["Mdot_min"] = star['Mdot']
		star["Mdot_max"] = star['Mdot']
	else:
		star["Mdot_min"] = star['Mdot'] + star["un_Mdot"][0]
		star["Mdot_max"] = star['Mdot'] + star["un_Mdot"][1]

	if ('N/A' in star["un_Period"]):
		star["Period_min"] = star['Period']
		star["Period_max"] = star['Period']
	else:
		star["Period_min"] = star['Period'] + star["un_Period"][0]
		star["Period_max"] = star['Period'] + star["un_Period"][1]

	if ('N/A' in star["un_Mass"]):
		star["Mass_min"] = star['Mass']
		star["Mass_max"] = star['Mass']
	else:
		star["Mass_min"] = star['Mass'] + star["un_Mass"][0]
		star["Mass_max"] = star['Mass'] + star["un_Mass"][1]

	if ('N/A' in star["un_Radius"]):
		star["Radius_min"] = star['Radius']
		star["Radius_max"] = star['Radius']
	else:
		star["Radius_min"] = star['Radius'] + star["un_Radius"][0]
		star["Radius_max"] = star['Radius'] + star["un_Radius"][1]

	Upsilon, RA, Torque_min = FM18_surf(Bdip=star["Bdip_min"], Bquad=star["Bquad_min"], Boct=star["Boct_min"], Mdot=star["Mdot_min"], Period=star["Period_max"], Mass=star["Mass_max"], Radius=star["Radius_min"], Info=False)
	#########################################################
	Upsilon, RA, Torque_max = FM18_surf(Bdip=star["Bdip_max"], Bquad=star["Bquad_max"], Boct=star["Boct_max"], Mdot=star["Mdot_max"], Period=star["Period_min"], Mass=star["Mass_min"], Radius=star["Radius_max"], Info=False)
	#########################################################
	print('Torque:',round(Torque/1e30,ndec),'+',round((Torque_max-Torque)/1e30,ndec),'-',round((Torque-Torque_min)/1e30,ndec),'x10^30erg') #value with errors
	print(' ')
	#########################################################
        
ndec=5 # Number of decimal places for the outputs

# F STARS #

_iotHor = {"Name":'iotHor', "Bdip":2.13, "Bquad":3.89, "Boct":3.91, "Mdot":89.5, "Period":7.7, "Mass":1.19, "Radius":1.161,
"un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-66.5,144], "un_Period":[-0.67,0.18], "un_Mass":[-0.07,0.07], "un_Radius":[-0.017,0.017]}

_88Leo = {"Name":'88Leo', "Bdip":4.03, "Bquad":0, "Boct":0, "Mdot":6.44, "Period":14, "Mass":1.12, "Radius":1.132,
"un_Bdip":[-0.00,0.53], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-4.49,8.68], "un_Period":[-0.5,0.5], "un_Mass":[-0.07,0.07], "un_Radius":[-0.032,0.032]}

_rhoCrB = {"Name":'rhoCrB', "Bdip":1.28, "Bquad":0, "Boct":0, "Mdot":0.24, "Period":17, "Mass":1.05, "Radius":1.300,
"un_Bdip":[-0.46,0.46], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-0.20,0.32], "un_Period":[-0.5,0.5], "un_Mass":[-0.06,0.06], "un_Radius":[-0.025,0.025]}

# SOLAR ANALOGS #

_kap1Cet = {"Name":'kap1Cet', "Bdip":16.0, "Bquad":15.6, "Boct":11.1, "Mdot":118, "Period":9, "Mass":1.05, "Radius":0.914,
"un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-88.9,198], "un_Period":[-0.5,0.5], "un_Mass":[-0.06,0.06], "un_Radius":[-0.014,0.014]}

_HD76151 = {"Name":'HD76151', "Bdip":5.98, "Bquad":2.15, "Boct":0.28, "Mdot":29.2, "Period":15, "Mass":1.05, "Radius":0.964,
"un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-19.8,39.4], "un_Period":[-0.5,0.5], "un_Mass":[-0.06,0.06], "un_Radius":[-0.018,0.018]}

_18Sco = {"Name":'18Sco', "Bdip":1.34, "Bquad":2.01, "Boct":0.86, "Mdot":0.36, "Period":22.7, "Mass":1.07, "Radius":1.009,
"un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-0.23,0.39], "un_Period":[-0.5,0.5], "un_Mass":[-0.06,0.06], "un_Radius":[-0.019,0.019]}

_16CygA = {"Name":'16CygA', "Bdip":0.46, "Bquad":0, "Boct":0, "Mdot":1.30, "Period":20.5, "Mass":1.10, "Radius":1.231,
"un_Bdip":[-0.44,0.45], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-0.78,1.37], "un_Period":[-1.1,2.0], "un_Mass":[-0.07,0.07], "un_Radius":[-0.024,0.024]}

_16CygB = {"Name":'16CygB', "Bdip":0.88, "Bquad":0, "Boct":0, "Mdot":0.42, "Period":21.2, "Mass":1.05, "Radius":1.157,
"un_Bdip":[-0.73,1.04], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-0.26,0.46], "un_Period":[-1.5,1.8], "un_Mass":[-0.06,0.06], "un_Radius":[-0.019,0.019]}

_51Peg = {"Name":'51Peg', "Bdip":0.77, "Bquad":0.44, "Boct":0.65, "Mdot":0.20, "Period":21.9, "Mass":1.10, "Radius":1.174,
"un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-0.17,0.27], "un_Period":[-0.4,0.4], "un_Mass":[-0.07,0.07], "un_Radius":[-0.023,0.023]}

# COOL G STARS #

_61UMa = {"Name":'61UMa', "Bdip":11.5, "Bquad":12.0, "Boct":6.12, "Mdot":29.8, "Period":17, "Mass":0.97, "Radius":0.855,
"un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-19.0,36.3], "un_Period":[-0.5,0.5], "un_Mass":[-0.06,0.06], "un_Radius":[-0.014,0.014]}

_tauCet = {"Name":'tauCet', "Bdip":0.86, "Bquad":0, "Boct":0, "Mdot":0.10, "Period":34, "Mass":0.82, "Radius":0.836,
"un_Bdip":[-0.34,0.34], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-0.08,0.13], "un_Period":[-0.5,0.5], "un_Mass":[-0.05,0.05], "un_Radius":[-0.020,0.020]}

# K STARS #

_epsEri = {"Name":'epsEri', "Bdip":14.6, "Bquad":8.78, "Boct":5.90, "Mdot":22.4, "Period":12, "Mass":0.86, "Radius":0.694,
"un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-16.3,34.4], "un_Period":[-0.5,0.5], "un_Mass":[-0.05,0.05], "un_Radius":[-0.014,0.014]}

_sigDra = {"Name":'sigDra', "Bdip":5.68, "Bquad":4.82, "Boct":4.76, "Mdot":4.17, "Period":27, "Mass":0.86, "Radius":0.769,
"un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-3.01,5.91], "un_Period":[-0.5,0.5], "un_Mass":[-0.05,0.05], "un_Radius":[-0.013,0.013]}

_107Psc = {"Name":'107Psc', "Bdip":4.24, "Bquad":2.77, "Boct":1.37, "Mdot":0.64, "Period":35, "Mass":0.86, "Radius":0.811,
"un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-0.46,0.76], "un_Period":[-0.5,0.5], "un_Mass":[-0.05,0.05], "un_Radius":[-0.017,0.017]}

_HD103095 = {"Name":'HD103095', "Bdip":0.61, "Bquad":0, "Boct":0, "Mdot":0.05, "Period":31, "Mass":0.60, "Radius":0.641,
"un_Bdip":[-0.04,0.03], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-0.04,0.07], "un_Period":[-0.5,0.5], "un_Mass":[-0.04,0.04], "un_Radius":[-0.016,0.016]}

_HD219134 = {"Name":'HD219134', "Bdip":2.39, "Bquad":4.05, "Boct":1.19, "Mdot":0.31, "Period":42.2, "Mass":0.80, "Radius":0.724,
"un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-0.16,0.31], "un_Period":[-0.9,0.9], "un_Mass":[-0.05,0.05], "un_Radius":[-0.014,0.014]}

_HD166620 = {"Name":'HD166620', "Bdip":2.81, "Bquad":0, "Boct":0, "Mdot":0.53, "Period":43, "Mass":0.83, "Radius":0.771,
"un_Bdip":[-0.94,0.95], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-0.34,0.58], "un_Period":[-0.5,0.5], "un_Mass":[-0.05,0.05], "un_Radius":[-0.019,0.019]}

print('# TORQUE UNCERTAINTIES #')
print(' ')

torqueUncert(_107Psc)
torqueUncert(_tauCet)
torqueUncert(_iotHor)
torqueUncert(_kap1Cet)
torqueUncert(_epsEri)
torqueUncert(_HD76151)
torqueUncert(_88Leo)
torqueUncert(_61UMa)
torqueUncert(_HD103095)
torqueUncert(_rhoCrB)
torqueUncert(_18Sco)
torqueUncert(_HD166620)
torqueUncert(_sigDra)
torqueUncert(_16CygA)
torqueUncert(_16CygB)
torqueUncert(_51Peg)
torqueUncert(_HD219134)
