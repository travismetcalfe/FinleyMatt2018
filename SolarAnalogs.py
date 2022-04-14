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
# Metcalfe et al. (2022) Torque Analysis #################
# Version 1.0 :  8th April 2022 ##########################
##########################################################
# Function:                                              #
# torqueDiagnostic-Calculates the FM18 torque for StarA  #
#	    and StarB, and evaluates the relative        #
#	    influence of changing each stellar parameter #
#	    in the braking law.                          #
##########################################################
# Input Parameters:                                      #
# fromStarA- Dictionary* of StarA Properties             #
# toStarB- Dictionary* of StarB Properties               #
# RhkFieldChange- Change in <B> given difference in R'HK #
# Dictionary*- Name, Bdip, Bquad, Boct, Mdot, Period,    #
#	    Mass, and Radius.                            #
##########################################################


def torqueDiagnostic(fromStarA,toStarB,RhkFieldChange,ndec=5):

	#########################################################
	#StarA Wind Torque:                                     #
	#########################################################
	print('['+fromStarA["Name"]+', Bd='+str(fromStarA["Bdip"])+', Bq='+str(fromStarA["Bquad"])+', Bo='+str(fromStarA["Boct"])+']')
	Upsilon, RA, Torque1 = FM18_surf(Bdip=fromStarA["Bdip"], Bquad=fromStarA["Bquad"], Boct=fromStarA["Boct"], Mdot=fromStarA["Mdot"], Period=fromStarA["Period"], Mass=fromStarA["Mass"], Radius=fromStarA["Radius"], Info=False)
	print('Torque:',round(Torque1/1e30,ndec),'x10^30erg') #fiducial model StarA
	print(' ')
	#########################################################

	# COMPARE TO StarA MODEL #
	print('# COMPARE TO '+fromStarA["Name"]+' MODEL #')
	print(' ')

	#########################################################
	#StarB Wind Torque:                                     #
	#########################################################
	print('['+toStarB["Name"]+', Bd='+str(toStarB["Bdip"])+', Bq='+str(toStarB["Bquad"])+', Bo='+str(toStarB["Boct"])+']')
	Upsilon, RA, Torque2 = FM18_surf(Bdip=toStarB["Bdip"], Bquad=toStarB["Bquad"], Boct=toStarB["Boct"], Mdot=toStarB["Mdot"], Period=toStarB["Period"], Mass=toStarB["Mass"], Radius=toStarB["Radius"], Info=False)
	print('Torque:',round(Torque2/1e30,ndec),'x10^30erg') #fiducial model StarB
	print('['+fromStarA["Name"]+' -> '+toStarB["Name"]+']:',round(Torque2/Torque1-1.,ndec))
	#########################################################

	#########################################################
	#Evaluate the Influence of Changing Each Parameter:     #
	#########################################################
	print('#parameter study#')
	Upsilon, RA, Torque3 = FM18_surf(Bdip=RhkFieldChange*fromStarA["Bdip"], Bquad=RhkFieldChange*fromStarA["Bquad"], Boct=RhkFieldChange*fromStarA["Boct"], Mdot=fromStarA["Mdot"], Period=fromStarA["Period"], Mass=fromStarA["Mass"], Radius=fromStarA["Radius"], Info=False)
	print('<B>:',round(Torque3/Torque1-1.,ndec))
	Upsilon, RA, Torque = FM18_surf(Bdip=toStarB["Bdip"], Bquad=toStarB["Bquad"], Boct=toStarB["Boct"], Mdot=fromStarA["Mdot"], Period=fromStarA["Period"], Mass=fromStarA["Mass"], Radius=fromStarA["Radius"], Info=False)
	print('Morph:',round(Torque/Torque3-1.,ndec)) # StarB morphology at same <B>
	print('both:',round(Torque/Torque1-1.,ndec))
	Upsilon, RA, Torque = FM18_surf(Bdip=fromStarA["Bdip"], Bquad=fromStarA["Bquad"], Boct=fromStarA["Boct"], Mdot=toStarB["Mdot"], Period=fromStarA["Period"], Mass=fromStarA["Mass"], Radius=fromStarA["Radius"], Info=False)
	print('Mdot:',round(Torque/Torque1-1.,ndec))
	Upsilon, RA, Torque = FM18_surf(Bdip=fromStarA["Bdip"], Bquad=fromStarA["Bquad"], Boct=fromStarA["Boct"], Mdot=fromStarA["Mdot"], Period=toStarB["Period"], Mass=fromStarA["Mass"], Radius=fromStarA["Radius"], Info=False)
	print('Prot:',round(Torque/Torque1-1.,ndec))
	Upsilon, RA, Torque = FM18_surf(Bdip=fromStarA["Bdip"], Bquad=fromStarA["Bquad"], Boct=fromStarA["Boct"], Mdot=fromStarA["Mdot"], Period=fromStarA["Period"], Mass=toStarB["Mass"], Radius=fromStarA["Radius"], Info=False)
	print('M:',round(Torque/Torque1-1.,ndec))
	Upsilon, RA, Torque = FM18_surf(Bdip=fromStarA["Bdip"], Bquad=fromStarA["Bquad"], Boct=fromStarA["Boct"], Mdot=fromStarA["Mdot"], Period=fromStarA["Period"], Mass=fromStarA["Mass"], Radius=toStarB["Radius"], Info=False)
	print('R:',round(Torque/Torque1-1.,ndec))
	Upsilon, RA, Torque = FM18_surf(Bdip=RhkFieldChange*fromStarA["Bdip"], Bquad=RhkFieldChange*fromStarA["Bquad"], Boct=RhkFieldChange*fromStarA["Boct"], Mdot=toStarB["Mdot"], Period=toStarB["Period"], Mass=toStarB["Mass"], Radius=toStarB["Radius"], Info=False)
	print('all non-Morph:',round(Torque/Torque1-1.,ndec))
	Upsilon, RA, Torque = FM18_surf(Bdip=toStarB["Bdip"], Bquad=toStarB["Bquad"], Boct=toStarB["Boct"], Mdot=toStarB["Mdot"], Period=toStarB["Period"], Mass=toStarB["Mass"], Radius=toStarB["Radius"], Info=False)
	print('all plus Morph:',round(Torque/Torque1-1.,ndec))
	print(' ')
	#########################################################


def uncertaintyDiagnostic(Star,Quantities,ndec=5):

	#########################################################
	#Star Wind Torque:                                      #
	#########################################################
	print('# '+Star["Name"]+' #')
	print(' ')

	print('['+Star["Name"]+', Bd='+str(Star["Bdip"])+', Bq='+str(Star["Bquad"])+', Bo='+str(Star["Boct"])+']')
	Upsilon, RA, Torque1 = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
	print('Torque:',round(Torque1/1e30,ndec),'x10^30erg') #fiducial model Star
	print(' ')
	#########################################################
	Bdip_n = 0
	Bdip_p = 0
	Bquad_n = 0
	Bquad_p = 0
	Boct_n = 0
	Boct_p = 0
	Mdot_n = 0
	Mdot_p = 0
	Period_n = 0
	Period_p = 0
	Mass_n = 0
	Mass_p = 0
	Period_n = 0
	Period_p = 0
	#########################################################
	#Influence of Unceratinty in Each Parameter:            #
	#########################################################
	print('# uncertainties #')

	for Quantity in Quantities:
		if Quantity in list(Star.keys()):
			print('['+Star["Name"]+', '+Quantity+'='+str(Star[Quantity])+'-'+str(Star["un_"+Quantity][0])+'/+'+str(Star["un_"+Quantity][1])+']')
			if Quantity == "Bdip":
                                Upsilon, RA, Torque1 = FM18_surf(Bdip=Star["Bdip"], Bquad=0, Boct=0, Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
                                Upsilon, RA, Torque_n = FM18_surf(Bdip=abs(Star["Bdip"]-Star["un_"+Quantity][0]), Bquad=0, Boct=0, Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
                                Upsilon, RA, Torque_p = FM18_surf(Bdip=abs(Star["Bdip"]+Star["un_"+Quantity][1]), Bquad=0, Boct=0, Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
                                Bdip_n = Star["un_"+Quantity][0]
                                Bdip_p = Star["un_"+Quantity][1]
			elif Quantity == "Bquad":
                                Upsilon, RA, Torque1 = FM18_surf(Bdip=0, Bquad=Star["Bquad"], Boct=0, Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
                                Upsilon, RA, Torque_n = FM18_surf(Bdip=0, Bquad=abs(Star["Bquad"]-Star["un_"+Quantity][0]), Boct=0, Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
                                Upsilon, RA, Torque_p = FM18_surf(Bdip=0, Bquad=abs(Star["Bquad"]+Star["un_"+Quantity][1]), Boct=0, Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
                                Bquad_n = Star["un_"+Quantity][0]
                                Bquad_p = Star["un_"+Quantity][1]
			elif Quantity == "Boct":
                                Upsilon, RA, Torque1 = FM18_surf(Bdip=0, Bquad=0, Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
                                Upsilon, RA, Torque_n = FM18_surf(Bdip=0, Bquad=0, Boct=abs(Star["Boct"]-Star["un_"+Quantity][0]), Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
                                Upsilon, RA, Torque_p = FM18_surf(Bdip=0, Bquad=0, Boct=abs(Star["Boct"]+Star["un_"+Quantity][1]), Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
                                Boct_n = Star["un_"+Quantity][0]
                                Boct_p = Star["un_"+Quantity][1]
			elif Quantity == "Mdot":
				Upsilon, RA, Torque_n = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"]-Star["un_"+Quantity][0], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Upsilon, RA, Torque_p = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"]+Star["un_"+Quantity][1], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Mdot_n = Star["un_"+Quantity][0]
				Mdot_p = Star["un_"+Quantity][1]
			elif Quantity == "Period":
				Upsilon, RA, Torque_n = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"]-Star["un_"+Quantity][0], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Upsilon, RA, Torque_p = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"]+Star["un_"+Quantity][1], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Period_n = Star["un_"+Quantity][0]
				Period_p = Star["un_"+Quantity][1]
			elif Quantity == "Mass":
				Upsilon, RA, Torque_n = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"]-Star["un_"+Quantity][0], Radius=Star["Radius"], Info=False)
				Upsilon, RA, Torque_p = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"]+Star["un_"+Quantity][1], Radius=Star["Radius"], Info=False)
				Mass_n = Star["un_"+Quantity][0]
				Mass_p = Star["un_"+Quantity][1]
			elif Quantity == "Radius":
				Upsilon, RA, Torque_n = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"]-Star["un_"+Quantity][0], Info=False)
				Upsilon, RA, Torque_p = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"]+Star["un_"+Quantity][1], Info=False)
				Radius_n = Star["un_"+Quantity][0]
				Radius_p = Star["un_"+Quantity][1]
			print(Quantity+'(-):',round(Torque_n/Torque1-1.,ndec))
			print(Quantity+'(+):',round(Torque_p/Torque1-1.,ndec))
			print(' ')
		else:
			print('No Quantity {'+Quantity+'} in '+Star["Name"]+' ...')


# FIDUCIAL MODELS #

_HD76151 = {"Name":'HD76151', "Bdip":5.13, "Bquad":2.88, "Boct":1.34, "Mdot":7.5, "Period":20.5, "Mass":1.04, "Radius":1.023,
            "un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[0.6,0.6], "un_Period":[0.3,0.3], "un_Mass":[0.06,0.06], "un_Radius":[0.026,0.026]}

_18Sco = {"Name":'18 Sco', "Bdip":1.34, "Bquad":2.01, "Boct":0.864, "Mdot":0.87, "Period":22.7, "Mass":1.02, "Radius":1.01,
	  "un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[0.32,0.32], "un_Period":[0.5,0.5], "un_Mass":[0.03,0.03], "un_Radius":[0.009,0.009]}

_16CygA = {"Name":'16 Cyg A', "Bdip":0.5, "Bquad":2.5, "Boct":5.8, "Mdot":0.92, "Period":20.5, "Mass":1.072, "Radius":1.223,
	      "un_Bdip":[0.4,0.4], "un_Bquad":[2.3,2.3], "un_Boct":[5.2,5.2], "un_Mdot":[0.16,0.16], "un_Period":[1.1,2.0], "un_Mass":[0.013,0.013], "un_Radius":[0.005,0.005]}
_16CygB = {"Name":'16 Cyg B', "Bdip":0.9, "Bquad":4.5, "Boct":31.3, "Mdot":0.57, "Period":21.2, "Mass":1.038, "Radius":1.113,
	      "un_Bdip":[1.0,1.0], "un_Bquad":[5.2,5.2], "un_Boct":[39.2,39.2], "un_Mdot":[0.13,0.13], "un_Period":[1.5,1.8], "un_Mass":[0.047,0.047], "un_Radius":[0.016,0.016]}

_16CygA_HD = {"Name":'16 Cyg A', "Bdip":0, "Bquad":2.19, "Boct":2.53, "Mdot":0.92, "Period":20.5, "Mass":1.072, "Radius":1.223,
	      "un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[0.16,0.16], "un_Period":[1.1,2.0], "un_Mass":[0.013,0.013], "un_Radius":[0.005,0.005]}
_16CygB_HD = {"Name":'16 Cyg B', "Bdip":0, "Bquad":2.36, "Boct":2.73, "Mdot":0.57, "Period":21.2, "Mass":1.038, "Radius":1.113,
	      "un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[0.13,0.13], "un_Period":[1.5,1.8], "un_Mass":[0.047,0.047], "un_Radius":[0.016,0.016]}

_16CygA_18 = {"Name":'16 Cyg A', "Bdip":0, "Bquad":1.34, "Boct":1.55, "Mdot":0.92, "Period":20.5, "Mass":1.072, "Radius":1.223,
	      "un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[0.16,0.16], "un_Period":[1.1,2.0], "un_Mass":[0.013,0.013], "un_Radius":[0.005,0.005]}
_16CygB_18 = {"Name":'16 Cyg B', "Bdip":0, "Bquad":1.45, "Boct":1.67, "Mdot":0.57, "Period":21.2, "Mass":1.038, "Radius":1.113,
	      "un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[0.13,0.13], "un_Period":[1.5,1.8], "un_Mass":[0.047,0.047], "un_Radius":[0.016,0.016]}


ndec=5 # Number of decimal places for the outputs

print('# EVOLUTIONARY SEQUENCE #')
print(' ')

torqueDiagnostic(_HD76151,_18Sco,0.425) # HD 76151 x (2.464/5.797) => 0.425
torqueDiagnostic(_HD76151,_16CygA_HD,0.4207) # HD76151 x Delta<logR'hk>=-0.376 => 0.4207
torqueDiagnostic(_HD76151,_16CygB_HD,0.4539) # HD76151 x Delta<logR'hk>=-0.343 => 0.4539
print(' ')

torqueDiagnostic(_18Sco,_16CygA_18,0.6077) # 18 Sco x Delta<logR'hk>=-0.217 => 0.6077
torqueDiagnostic(_18Sco,_16CygB_18,0.6549) # 18 Sco x Delta<logR'hk>=-0.184 => 0.6549
print(' ')

torqueDiagnostic(_18Sco,_16CygA_HD,0.6077) # 18 Sco x Delta<logR'hk>=-0.217 => 0.6077
torqueDiagnostic(_18Sco,_16CygB_HD,0.6549) # 18 Sco x Delta<logR'hk>=-0.184 => 0.6549
print(' ')

print('# UNCERTAINTIES #')
print(' ')

# HD 76151 #

uncertaintyDiagnostic(_HD76151,['Mdot','Period','Mass','Radius'])

# combined uncertainty for fiducial model
Upsilon, RA, Torque1 = FM18_surf(Bdip=5.13, Bquad=2.88, Boct=1.34, Mdot=7.5, Period=20.5, Mass=1.04, Radius=1.023, Info=False)

Upsilon, RA, Torque = FM18_surf(Bdip=5.13, Bquad=2.88, Boct=1.34, Mdot=6.9, Period=20.8, Mass=1.10, Radius=0.997, Info=False)
print('all(-):',round(Torque/Torque1-1.,ndec))
Upsilon, RA, Torque = FM18_surf(Bdip=5.13, Bquad=2.88, Boct=1.34, Mdot=8.1, Period=20.2, Mass=0.98, Radius=1.049, Info=False)
print('all(+):',round(Torque/Torque1-1.,ndec))
print(' ')
print(' ')


# 18 Sco #

uncertaintyDiagnostic(_18Sco,['Mdot','Period','Mass','Radius'])

# combined uncertainty for fiducial model
Upsilon, RA, Torque1 = FM18_surf(Bdip=1.34, Bquad=2.01, Boct=0.864, Mdot=0.87, Period=22.7, Mass=1.02, Radius=1.01, Info=False)

Upsilon, RA, Torque = FM18_surf(Bdip=1.34, Bquad=2.01, Boct=0.864, Mdot=0.55, Period=23.2, Mass=1.05, Radius=1.001, Info=False)
print('all(-):',round(Torque/Torque1-1.,ndec))
Upsilon, RA, Torque = FM18_surf(Bdip=1.34, Bquad=2.01, Boct=0.864, Mdot=1.19, Period=22.2, Mass=0.99, Radius=1.019, Info=False)
print('all(+):',round(Torque/Torque1-1.,ndec))
print(' ')
print(' ')


# 16 Cyg A #

uncertaintyDiagnostic(_16CygA,['Bdip','Bquad','Boct'])
uncertaintyDiagnostic(_16CygA_18,['Mdot','Period','Mass','Radius'])

# combined uncertainty for pure quadrupole
Upsilon, RA, Torque1 = FM18_surf(Bdip=0, Bquad=2.5, Boct=0, Mdot=0.92, Period=20.5, Mass=1.072, Radius=1.223, Info=False)

Upsilon, RA, Torque = FM18_surf(Bdip=0, Bquad=0.2, Boct=0, Mdot=0.76, Period=22.5, Mass=1.085, Radius=1.218, Info=False)
print('all(-):',round(Torque/Torque1-1.,ndec))
Upsilon, RA, Torque = FM18_surf(Bdip=0, Bquad=4.8, Boct=0, Mdot=1.08, Period=19.4, Mass=1.059, Radius=1.228, Info=False)
print('all(+):',round(Torque/Torque1-1.,ndec))
print(' ')
print(' ')


# 16 Cyg B #

uncertaintyDiagnostic(_16CygB,['Bdip','Bquad','Boct'])
uncertaintyDiagnostic(_16CygB_HD,['Mdot','Period','Mass','Radius'])

# combined uncertainty for pure quadrupole
Upsilon, RA, Torque1 = FM18_surf(Bdip=0, Bquad=4.5, Boct=0, Mdot=0.57, Period=21.2, Mass=1.038, Radius=1.113, Info=False)

Upsilon, RA, Torque = FM18_surf(Bdip=0, Bquad=0.7, Boct=0, Mdot=0.44, Period=23.0, Mass=1.085, Radius=1.097, Info=False)
print('all(-):',round(Torque/Torque1-1.,ndec))
Upsilon, RA, Torque = FM18_surf(Bdip=0, Bquad=9.7, Boct=0, Mdot=0.70, Period=19.7, Mass=0.991, Radius=1.129, Info=False)
print('all(+):',round(Torque/Torque1-1.,ndec))
print(' ')
