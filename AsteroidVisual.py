# Orbit Determination
# Max Lu

from visual import *
from math import*
from numpy import*

# Function Definitions
# Corrects angle ambiguity based on its quadrant
def angleambiguity(sin,cos):
    if sin > 0 and cos > 0:#quadrant 1
        return acos(cos)
    if cos < 0 and sin > 0:#quadrant 2
        return acos(cos)
    if sin < 0 and cos < 0:#quadrant 3
        return 2 * pi + acos(cos)
    if cos > 0 and sin < 0:#quadrant 4
        return 2 * pi - acos(cos)

# Finds r vector at certain Julian Date 
def EphemerisGen(t):
    # Orbital Elements
    a = 1.88235825136
    ecc = 0.403064386884
    i = 29.1487264055*pi/180
    W = 114.531246122*pi/180
    w = 234.372957239 * pi/180
    Mte = 341.842495564 * pi/180

    # Constants
    k = 0.01720209894
    n = k/(a**(1.5))
    TE = 2455400.5
    e = 0.40909262775
    
    # Convergence of E
    E = Mte
    M = n * (t - TE) + Mte
    lastE = 0
    while abs(E - lastE) > 1e-12:
        lastE = E
        E = E - (M - (E - ecc  * sin(E)))/(ecc * cos(E) - 1)

    # Spinning the vector
    r = matrix([[a* cos(E) - a* ecc],
                [a * (1 -ecc**2)**.5 * sin(E)],
                [0]])
    M1 = matrix([[cos(W),-sin(W),0],
                [sin(W),cos(W),0],
                [0,0,1]])
    M2 = matrix([[1,0,0],
                 [0,cos(i),-sin(i)],
                 [0,sin(i),cos(i)]])
    M3 = matrix([[cos(w),-sin(w),0],
                 [sin(w),cos(w),0],
                 [0,0,1]])
    r1 = M1* M2 * M3 * r

    # Into equatorial conversion
    M4 = ([[1,0,0],
           [0,cos(e),-sin(e)],
           [0,sin(e),cos(e)]])
    r2 = M4 * r1

    return vector(r2[0,0] * 1.496e11,r2[1,0] * 1.496e11,r2[2,0] * 1.496e11)

curve(pos=[(-3 * 10**11,0,0),(3 * 10**11,0,0)]) # draws x axis
curve(pos=[(0,-3 * 10**11,0),(0,3 * 10**11,0)]) # draws y axis
curve(pos=[(0,0,-3 * 10**11),(0,0,3 * 10**11)]) # draws z axis

# constants
G = 6.67428e-11
t = 2455389.82888
dt = .005

sun = sphere(pos = vector(0,0,0), radius = 6960000000, color = color.yellow, mass = 1.9891e30)

# a = semi major axis
# b = semi minor axis
# v = velocity (initial velocity determined by formula)
# accel = acceleration due to gravity (determined by formula)

# draws earth and assigns attributes
earth = sphere(pos = vector(-3.432887344502621E-01* 1.496e11,  9.568883454208035E-01* 1.496e11, -2.209784931791160E-05* 1.496e11), radius = 637100000, eccentricity = .01671123, pe = 31557600, color = color.blue)
earth.a = 149598261000
earth.b = (earth.a**2 - (earth.eccentricity * earth.a)**2)**(1./2)
earth.v = vector(-1.592075835161679E-02* 1.496e11/86400, -5.748661396619696E-03* 1.496e11/86400, -4.842887735454858E-07* 1.496e11/86400)
earthtrail = curve()
earth.accel = -1 * G * sun.mass * earth.pos/(mag(earth.pos)**3)


ast = sphere(pos = EphemerisGen(2455389.82888),radius = 1e9,color = color.red)
asttrail = curve()

mindis = mag(ast.pos- earth.pos)
time = 0

# animates asteroid and earth
i = 0
track = false
label1 = label()
while(1==1):
    
    earth.pos = earth.pos + earth.v * dt * 86400 # multiplies dt by 86400 seconds/Julian day
    earth.v = earth.v + earth.accel * dt * 86400
    earth.accel = -1 * G * sun.mass * earth.pos/(mag(earth.pos)**3)
    earthtrail.append(pos=earth.pos)
    
    ast.pos = EphemerisGen(t)
    asttrail.append(pos = ast.pos)
    # centers camera on an object if clicked
    if(scene.mouse.clicked):
        m = scene.mouse.getclick()
        for obj in scene.objects:
            if m.pick == obj:
                scene.center = obj.pos
                track = true
                break

    # tracks object
    if(track):
        scene.center = obj.pos
    # counts number of elapsed Julian dates, converts to years, and prints out
    t = t + dt
    label1.visible = 0
    label1 = label(pos=label1.pos,text = (str((t-2455389.82888)/365)+ " years"))

  
