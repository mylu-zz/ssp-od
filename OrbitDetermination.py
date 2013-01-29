from visual import *
#Orbit Determination
#Data file used: Data.txt
#Max Lu

from math import *

#Data-Replaced by File Input/Output
#solar vectors
'''
R1 = vector(-2.462502561048408E-01,9.864109051578982E-01,-8.141041345977961E-06) # R-1 vector
R2 = vector(-3.113766307131727E-01,9.677830701688845E-01,-8.433405262724878E-06) # R0 vector
R3 = vector(-3.593272602834906E-01,9.509197039355376E-01,-9.109124036217700E-06) # R+1 vector
R = [R1,R2,R3]
RA= [(21+18./60+14.0557342423/3600)*15*pi/180,(21+28./60+43.535/3600)*15*pi/180,(21+37./60+44.1338932742/3600)*15*pi/180]
Dec =[(-10-6./60-37.9037427425/3600) * pi/180,(-15-30./60-51.454/3600)*pi/180,(-20-17./60-45.4055367936/3600)*pi/180]
t[0] = 2455383.8308
t[1] = 2455387.83443
t[2] = 2455390.83781
'''
##############################################################################################################################################################################################################
#constants and instance variable declarations
k = 0.01720209895
e = 23.439281 * pi/180
mu = 1
RA= [0,0,0]
Dec =[0,0,0]
t = [0,0,0]
R = [vector(0,0,0),vector(0,0,0),vector(0,0,0)]

##############################################################################################################################################################################################################
#function definitions
#finds f series value for a given r, rdot, and tau value
def fseries(r,rdot,i):
    f = 1 - tau[i]**2/(2 * mag(r)**3) + tau[i]**3 * dot(r,rdot)/(2 * mag(r)**5) + tau[i]**4/24 * (3/mag(r)**3 * (dot(rdot,rdot)/(mag(r)**2) - 1/(mag(r)**3)) - 15/mag(r)**7* (dot(r,rdot))**2 + 1/mag(r)**6)
    return f

#finds g series value for a given r, rdot, and tau value
def gseries(r,rdot,i):
    g = tau[i] - tau[i]**3/(6 * mag(r)**3) + tau[i]**4 * dot(r,rdot)/(4 * mag(r)**5)
    return g

#reads data file and assigns values to RA, Dec, t, and R lists
def readdata():
    f = open ('Data.txt','r')
    #observation 1
    line1 = f.readline()
    line1 = line1.split()
    t[0] = float(line1[0])
    '\n'
    line2 = f.readline()
    line2 = line2.split()
    RA[0] = (float(line2[0]) + float(line2[1])/60. + float(line2[2])/3600.) * 15 * pi/180
    Dec[0] = (float(line2[3]) - float(line2[4])/60. - float(line2[5])/3600.) * pi/180
    '\n'
    line3 = f.readline()
    line3 = line3.split()
    R[0] = vector(float(line3[0]),float(line3[1]),float(line3[2]))
    '\n'

    #observation 2
    line4 = f.readline()
    line4 = line4.split()
    t[1] = float(line4[0])
    '\n'
    line5 = f.readline()
    line5 = line5.split()
    RA[1] = (float(line5[0]) + float(line5[1])/60. + float(line5[2])/3600.) * 15 * pi/180
    Dec[1] = (float(line5[3]) - float(line5[4])/60. - float(line5[5])/3600.) * pi/180
    '\n'
    line6 = f.readline()
    line6 = line6.split()
    R[1] = vector(float(line6[0]),float(line6[1]),float(line6[2]))
    '\n'
    '\n'

    #observation 3
    line7 = f.readline()
    line7 = line7.split()
    t[2] = float(line7[0])
    '\n'
    line8 = f.readline()
    line8 = line8.split()
    RA[2] = (float(line8[0]) + float(line8[1])/60. + float(line8[2])/3600.) * 15 * pi/180
    Dec[2] = (float(line8[3]) - float(line8[4])/60. - float(line8[5])/3600.) * pi/180
    '\n'
    line9 = f.readline()
    line9 = line9.split()
    R[2] = vector(float(line9[0]),float(line9[1]),float(line9[2]))
    return

#corrects angle ambiguity based on its quadrant
def angleambiguity(sin,cos):
    if sin > 0 and cos > 0:#quadrant 1
        return acos(cos)
    if cos < 0 and sin > 0:#quadrant 2
        return acos(cos)
    if sin < 0 and cos < 0:#quadrant 3
        return 2 * pi + acos(cos)
    if cos > 0 and sin < 0:#quadrant 4
        return 2 * pi - acos(cos)

##############################################################################################################################################################################################################
#main program

readdata()

#initializes rhohat for observations[0,1,2] and converts to ecliptic coordinates
i = 0
rhohat = [vector(0,0,0),vector(0,0,0),vector(0,0,0)]
while i <3:
    rhohat[i] = vector(cos(RA[i]) * cos(Dec[i]) , sin(RA[i]) * cos(Dec[i]) , sin(Dec[i]) )
    temp = cos(e) * rhohat[i].y + sin(e) * rhohat[i].z
    rhohat[i].z = -1 * sin(e) * rhohat[i].y + cos(e) * rhohat[i].z
    rhohat[i].y = temp
    i = i + 1

#finds taus for time interval
tau1 = k * (t[0] - t[1])
tau2 = k * (t[2] - t[0])
tau3 = k * (t[2] - t[1])
tau = [tau1,tau2,tau3]
a3 = -1 * tau1/tau2
a1 = tau3/tau2

#initializes rho magnitudes based on a3 and a1 values
rho = [0,0,0]
rho[0] = (a1 * dot(cross(R[0], rhohat[1]),rhohat[2]) - dot(cross(R[1],rhohat[1]),rhohat[2]) + a3 * dot(cross(R[2],rhohat[1]),rhohat[2]))/(a1 * dot(cross(rhohat[0],rhohat[1]),rhohat[2]))
rho[1] = (a1 * dot(cross(rhohat[0], R[0]),rhohat[2]) - dot(cross(rhohat[0],R[1]),rhohat[2]) + a3 * dot(cross(rhohat[0],R[2]),rhohat[2]))/(-1 * dot(cross(rhohat[0],rhohat[1]),rhohat[2]))
rho[2] = (a1 * dot(cross(rhohat[1], R[0]),rhohat[0]) - dot(cross(rhohat[1],R[1]),rhohat[0]) + a3 * dot(cross(rhohat[1],R[2]),rhohat[0]))/(a3 * dot(cross(rhohat[1],rhohat[2]),rhohat[0]))

#finds r using rho and rhohat
r = [vector(0,0,0),vector(0,0,0),vector(0,0,0)]
r[0] = rho[0] * rhohat[0] - R[0]
r[2] = rho[2] * rhohat[2]- R[2]
r[1] = a1 * r[0] + a3 * r[2]
#finds rdots for corresponding values
rdot0 = ((r[1] - r[0])/tau[1] + (r[0] - r[2])/(-1 *tau[2]))/2
rdot1 = ((r[2] - r[1])/tau[2] + (r[1] - r[0])/(-1 *tau[0]))/2
rdot2 = ((r[0] - r[2])/tau[0] + (r[2] - r[1])/(-1 *tau[1]))/2
rdot = [rdot0,rdot1,rdot2]

#declares f and g lists
f = [0,0,0]
g = [0,0,0]


#index variables
i=0
count = 0
notconverging = True
temp = vector(0,0,0)
while(notconverging):
    print r[1]
    #previous a1 and a3 values
    a1temp = a1
    a3temp = a3

    #finds f and g series using r[1] and rdot[1] values for tau[0,1,2]
    i = 0
    while(i <3):
        f[i] = fseries(r[1],rdot[1],i)
        g[i] = gseries(r[1],rdot[1],i)
        i = i + 1
    
    #updates r[0,1,2], a1, and a3 using f and g series
    r[0] = rho[0] * rhohat[0] - R[0]
    r[2] = rho[2] * rhohat[2] - R[2]
    a1 = g[2]/(f[0] * g[2] - f[2]*g[0])
    a3 = -1 * g[0]/(f[0] * g[2] - f[2]*g[0])
    rdot[1] = ((r[2] - f[2] * r[1])/g[2] + (r[0] - f[0] * r[1])/g[0])/2.
    
    #updates rho values using new a1 and a3 values
    rho[0] = (a1 * dot(cross(R[0], rhohat[1]),rhohat[2]) - dot(cross(R[1],rhohat[1]),rhohat[2]) + a3 * dot(cross(R[2],rhohat[1]),rhohat[2]))/(a1 * dot(cross(rhohat[0],rhohat[1]),rhohat[2]))
    rho[1] = (a1 * dot(cross(rhohat[0], R[0]),rhohat[2]) - dot(cross(rhohat[0],R[1]),rhohat[2]) + a3 * dot(cross(rhohat[0],R[2]),rhohat[2]))/(-1 * dot(cross(rhohat[0],rhohat[1]),rhohat[2]))
    rho[2] = (a1 * dot(cross(rhohat[1], R[0]),rhohat[0]) - dot(cross(rhohat[1],R[1]),rhohat[0]) + a3 * dot(cross(rhohat[1],R[2]),rhohat[0]))/(a3 * dot(cross(rhohat[1],rhohat[2]),rhohat[0]))
    r[1] = rho[1] * rhohat[1] - R[1]
    
    #breaks out of loop once a1 converges
    if abs(a1 - a1temp) <= 1e-12:
        notconverging = False

    #counts number of iterations
    count = count + 1
    

print "Number of iterations: ", count
print ""
print "r0 vector: ",r[1]
print "rdot0 vector: ", rdot[1]
print "rho0 vector: ", rho[0] * rhohat[0]

#Orbital Elements
#Semi Major Axis
a = 1/(2/mag(r[1]) - dot(rdot[1],rdot[1])/mu)

#Eccentricity
e = (1-mag(cross(r[1],rdot[1]))**2/(mu * a))**(1./2)

#Orbit Inclination
h = cross(r[1],rdot[1])
i = arctan((h.x**2 + h.y**2)**(1./2)/h.z)

#Longitude of Ascending Node
sinOmega = h.x/(mag(h) * sin(i))
cosOmega = -1 * h.y/(mag(h) * sin(i))
Omega = angleambiguity(sinOmega,cosOmega)

#Argument of Perihelion
sinnu = (a * (1-e**2) * dot(r[1],rdot[1])/(e * mag(h) * mag(r[1])))
cosnu = (a * (1-e**2)/mag(r[1]) - 1)/e
nu = angleambiguity(sinnu,cosnu)
cosU = (r[1].x * cos(Omega) + r[1].y * sin(Omega))/mag(r[1])
sinU = r[1].z/(mag(r[1]) * sin(i))
U = angleambiguity(sinU,cosU)
omega = U - nu
while omega < 0:
    omega = omega + 360
omega = fmod(omega,360)

#Mean Anomaly
E = acos(1/e * (1 - mag(r[0])/a))
M = E - e * sin(E)

n = k/(a**(1.5))
te = 2455400.5
if(t[1] > te):
    Mte = n*(te - t[1]) + M
else:
    Mte = n*(te - t[1]) - M

while Mte < 0:
    Mte = Mte + 2 * pi
Mte = fmod(Mte,2*pi)
Mte = Mte * 180/pi
Omega = Omega * 180/pi
omega = omega * 180/pi
i = i * 180/pi

print ""
print "Orbital Elements:"
print "Semi Major Axis: ",a, " A.U."
print "Eccentricity: ",e
print "Orbit Inclination: ",i, " degrees"
print "Longitude of Ascending Node: ",Omega, " degrees"
print "Argument of Perihelion: ",omega, " degrees"
print "Mean Anomaly: ", Mte, " degrees"

##############################################################################################################################################################################################################


