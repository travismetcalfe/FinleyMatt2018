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
				Upsilon, RA, Torque_n = FM18_surf(Bdip=Star["Bdip"]-Star["un_"+Quantity][0], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Upsilon, RA, Torque_p = FM18_surf(Bdip=Star["Bdip"]+Star["un_"+Quantity][1], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Bdip_n = Star["un_"+Quantity][0]
				Bdip_p = Star["un_"+Quantity][1]
			elif Quantity == "Bquad":
				Upsilon, RA, Torque_n = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"]-Star["un_"+Quantity][0], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Upsilon, RA, Torque_p = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"]+Star["un_"+Quantity][1], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Bquad_n = Star["un_"+Quantity][0]
				Bquad_p = Star["un_"+Quantity][1]
			elif Quantity == "Boct":
				Upsilon, RA, Torque_n = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"]-Star["un_"+Quantity][0], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Upsilon, RA, Torque_p = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"]+Star["un_"+Quantity][1], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Boct_n = Star["un_"+Quantity][0]
				Boct_p = Star["un_"+Quantity][1]
			elif Quantity == "Mdot":
				Upsilon, RA, Torque_n = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"]-Star["un_"+Quantity][0], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Upsilon, RA, Torque_p = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"]+Star["un_"+Quantity][1], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Mdot_n = Star["un_"+Quantity][0]
				Mdot_p = Star["un_"+Quantity][1]
			elif Quantity == "Period":
				Upsilon, RA, Torque_n = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"]+Star["un_"+Quantity][1], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Upsilon, RA, Torque_p = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"]-Star["un_"+Quantity][0], Mass=Star["Mass"], Radius=Star["Radius"], Info=False)
				Period_n = Star["un_"+Quantity][0]
				Period_p = Star["un_"+Quantity][1]
			elif Quantity == "Mass":
				Upsilon, RA, Torque_n = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"]+Star["un_"+Quantity][1], Radius=Star["Radius"], Info=False)
				Upsilon, RA, Torque_p = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"]-Star["un_"+Quantity][0], Radius=Star["Radius"], Info=False)
				Mass_n = Star["un_"+Quantity][0]
				Mass_p = Star["un_"+Quantity][1]
			elif Quantity == "Radius":
				Upsilon, RA, Torque_n = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"]-Star["un_"+Quantity][0], Info=False)
				Upsilon, RA, Torque_p = FM18_surf(Bdip=Star["Bdip"], Bquad=Star["Bquad"], Boct=Star["Boct"], Mdot=Star["Mdot"], Period=Star["Period"], Mass=Star["Mass"], Radius=Star["Radius"]+Star["un_"+Quantity][1], Info=False)
				Radius_n = Star["un_"+Quantity][0]
				Radius_p = Star["un_"+Quantity][1]
			print(Quantity+'(-):',round(1.-Torque_n/Torque1,ndec))
			print(Quantity+'(+):',round(Torque_p/Torque1-1.,ndec))
			print(' ')
		else:
			print('No Quantity {'+Quantity+'} in '+Star["Name"]+' ...')


	print('['+Star["Name"]+', +/- all:'+str(Quantities)+']')
	Upsilon, RA, Torque_n = FM18_surf(Bdip=Star["Bdip"]-Bdip_n, Bquad=Star["Bquad"]-Bquad_n, Boct=Star["Boct"]-Boct_n, Mdot=Star["Mdot"]-Mdot_n, Period=Star["Period"]+Period_p, Mass=Star["Mass"]+Mass_p, Radius=Star["Radius"]-Radius_n, Info=False)
	print('all(-):',round(1.-Torque_n/Torque1,ndec))
	Upsilon, RA, Torque_p = FM18_surf(Bdip=Star["Bdip"]+Bdip_p, Bquad=Star["Bquad"]+Bquad_p, Boct=Star["Boct"]+Boct_p, Mdot=Star["Mdot"]+Mdot_p, Period=Star["Period"]-Period_n, Mass=Star["Mass"]-Mass_n, Radius=Star["Radius"]+Radius_p, Info=False)
	print('all(+):',round(Torque_p/Torque1-1.,ndec))
	print(' ')


def torqueDiagnostic(fromStarA,toStarB,FieldChange,ndec=5):

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
	Upsilon, RA, Torque3 = FM18_surf(Bdip=FieldChange*fromStarA["Bdip"], Bquad=FieldChange*fromStarA["Bquad"], Boct=FieldChange*fromStarA["Boct"], Mdot=fromStarA["Mdot"], Period=fromStarA["Period"], Mass=fromStarA["Mass"], Radius=fromStarA["Radius"], Info=False)
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
	print(' ')
	#########################################################

_16CygA_Bd = {"Name":'16 Cyg A', "Bdip":0.5, "Bquad":0, "Boct":0, "Mdot":0.92, "Period":20.5, "Mass":1.072, "Radius":1.223,
              "un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[0.16,0.16], "un_Period":[1.1,2.0], "un_Mass":[0.013,0.013], "un_Radius":[0.005,0.005]}

_alt31Aql = {"Name":'31 Aql', "Bdip":8.5, "Bquad":0.0, "Boct":0.0, "Mdot":1.05, "Period":40.3, "Mass":1.07, "Radius":1.358,
          "un_Bdip":[-1.1,1.1], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-0.58,1.08], "un_Period":[-1.5,1.5], "un_Mass":[-0.04,0.04], "un_Radius":[-0.016,0.016]}

_31Aql = {"Name":'31 Aql', "Bdip":8.5, "Bquad":0.0, "Boct":0.0, "Mdot":1.05, "Period":40.3, "Mass":1.07, "Radius":1.358,
"un_Bdip":[1.1,1.1], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[0.58,1.08], "un_Period":[1.5,1.5], "un_Mass":[0.04,0.04], "un_Radius":[0.016,0.016]}

_31Aql_Bd = {"Name":'31 Aql', "Bdip":8.5, "Bquad":0.0, "Boct":0.0, "Mdot":1.05, "Period":40.3, "Mass":1.07, "Radius":1.358,
"un_Bdip":[-1.1,1.1], "un_Bquad":['N/A','N/A'], "un_Boct":['N/A','N/A'], "un_Mdot":[-0.58,1.08], "un_Period":[-1.5,1.5], "un_Mass":[-0.04,0.04], "un_Radius":[-0.016,0.016]}

_31Aql_Bq = {"Name":'31 Aql', "Bdip":0.0, "Bquad":133, "Boct":0.0, "Mdot":1.05, "Period":40.3, "Mass":1.07, "Radius":1.358,
"un_Bdip":['N/A','N/A'], "un_Bquad":[-13,13], "un_Boct":['N/A','N/A'], "un_Mdot":[-0.58,1.08], "un_Period":[-1.5,1.5], "un_Mass":[-0.04,0.04], "un_Radius":[-0.016,0.016]}

_31Aql_Bo = {"Name":'31 Aql', "Bdip":0.0, "Bquad":0.0, "Boct":256, "Mdot":1.05, "Period":40.3, "Mass":1.07, "Radius":1.358,
"un_Bdip":['N/A','N/A'], "un_Bquad":['N/A','N/A'], "un_Boct":[-74,74], "un_Mdot":[-0.58,1.08], "un_Period":[-1.5,1.5], "un_Mass":[-0.04,0.04], "un_Radius":[-0.016,0.016]}


print('# TORQUE #')
print(' ')

torqueUncert(_31Aql_Bd)
print(' ')
torqueUncert(_31Aql_Bq)
print(' ')
torqueUncert(_31Aql_Bo)
print(' ')

print('# UNCERTAINTIES #')
print(' ')

uncertaintyDiagnostic(_31Aql,['Bdip','Mdot','Period','Mass','Radius'])

torqueDiagnostic(_16CygA_Bd,_31Aql,0.944) # Delta logR'HK=-0.025 => 0.944
print(' ')

