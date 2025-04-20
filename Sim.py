##########################################################################################################################################################################################################

import vpython as vp

BackgroundColor = vp.vec(241/255,242/255,244/255)
vp.scene.background = BackgroundColor

vp.canvas.get_selected().align = "left"

# // Utils // 

def FormatNumber(Number, Digits):
    String = str(round(Number, Digits))
    DecimalSplit = String.split(".")
    DecimalSplit.append("") # Making sure index 1 is not empty
    Decimal = DecimalSplit[1] + "0"*Digits # Adding 0s at the end, instead of nothingness

    return DecimalSplit[0] + "." + Decimal[:Digits]

# // Constants

SIM_RATE = 240 # hz

g = 9.81
l = 2.5

r = 0.5 # Distance from the seat to the center of mass of BEATRICE
I = 10 # Moment of inertia of BEATRICE (need real mesurments)

TAILLE_BEATRICE = vp.vec(0.2, 1.55, 0.35)
MASSE_BEATRICE = 46.72

# Constants you SHALL NOT TOUCH

BLOCK_OFF = -1/2 - TAILLE_BEATRICE.y/2
art_l = I/MASSE_BEATRICE # Artificial length, of an equivalent point mass system

# // Coordinates

# The life saver:
# https://web.mit.edu/jorloff/www/chaosTalk/double-pendulum/double-pendulum-en.html

def P1_upt_func():
    m1 = P1["m"]
    m2 = P2["m"]

    Num = -g*(2*m1 + m2)*vp.sin(P1["x"]) - m2*g*vp.sin(P1["x"] - 2*P2["x"]) - 2*vp.sin(P1["x"] - P2["x"])*m2*((P2["v"]**2)*art_l + (P1["v"]**2)*l*vp.cos(P1["x"] - P2["x"]))
    Denum = l*(2*m1 + m2 - m2*vp.cos(2*P1["x"] - 2*P2["x"]))
    return Num/Denum
    #P1["a"] = (-1/((l**2)*(m1 + m2))) * (m2*l*art_l*P2["a"]*vp.cos(P2["x"] - P1["x"]) - m2*l*art_l*P2["v"]*vp.sin(P2["x"] - P1["x"]) + (m1 + m2)*g*l*vp.sin(P1["x"]))

def P2_upt_func():
    m1 = P1["m"]
    m2 = P2["m"]

    Num = 2*vp.sin(P1["x"] - P2["x"])*((P1["v"]**2)*l*(m1 + m2) + g*(m1 + m2)*vp.cos(P1["x"]) + (P2["v"]**2)*art_l*m2*vp.cos(P1["x"] - P2["x"]))
    Denum = l*(2*m1 + m2 - m2*vp.cos(2*P1["x"] - 2*P2["x"]))
    return Num/Denum
	#(-1/(m2*art_l^2)) * (m2*l*art_l*P1["a"]*math.cos(P2["x"] - P1["x"]) + m2*l*art_l*P1["v"]*math.sin(P2["x"] - P1["x"]) + m2*g*art_l*math.sin(P2["x"]))

#P2 = None

obj1 = vp.box(color = vp.color.red, size = vp.vec(.5,.1,1))
obj2 = vp.box(color = vp.color.yellow, size = TAILLE_BEATRICE) # texture = "test.png"

P1 = {
    "a": 0,
    "v": 0,
    "x": 1,
    "m" : 1,
    "obj": obj1,
    "upd_func": P1_upt_func,
}

P2 = {
    "a": 0,
    "v": 0,
    "x": 0,
    "m": MASSE_BEATRICE,
    "obj": obj2,
    "upd_func": P2_upt_func,
}

# // Events & UI

Running = True

def Exit():
    global Running
    Running = False

vp.button(text = "<b>Exit</b>", pos = vp.scene.title_anchor, bind = Exit)

# // Loop //

# Label

TextLabel = vp.wtext(text="")

# Loop variables
sim_dt = 0
dt = 1/SIM_RATE
i = 0

while Running:
    i += 1

    vp.rate(SIM_RATE)

    start_tick = vp.clock()

    # Updating the block's acceleration
    a1 = P1["upd_func"]()
    a2 = P2["upd_func"]()

    P1["a"] = a1
    P2["a"] = a2

    P1["v"] += P1["a"]*dt
    P1["x"] += P1["v"]*dt

    P2["v"] += P2["a"]*dt
    P2["x"] += P2["v"]*dt

    angle_1 = P1["x"] - vp.pi/2 # -pi to have 0 at (1,0) on the unit circle
    P1["obj"].pos = l*vp.vec(vp.cos(angle_1),vp.sin(angle_1),0)
    P1["obj"].axis = vp.vec(vp.cos(angle_1 + vp.pi/2), vp.sin(angle_1 + vp.pi/2), 0) * P1["obj"].length

    angle_2 = P2["x"] - vp.pi/2 # -pi to have 0 at (1,0) on the unit circle
    P2["obj"].pos = (P1["obj"].pos + vp.vec(vp.cos(angle_2),vp.sin(angle_2),0) 
        #+ r*vp.vec(vp.cos(angle_1 + vp.pi/2), vp.sin(angle_1 + vp.pi/2), 0)
        #+ r*vp.vec(vp.cos(angle_2 + vp.pi/2), vp.sin(angle_2 + vp.pi/2), 0)
    )
    P2["obj"].axis = vp.vec(vp.cos(angle_2 + vp.pi/2), vp.sin(angle_2 + vp.pi/2), 0) * P2["obj"].length


    #MiscAnimation()

    t = i*dt

    #PlotGraph(t)

    # Debug labels

    sim_time = FormatNumber(t, 3) + " secondes"
    ms_str = FormatNumber(sim_dt*1000, 3) + " ms"
    load_string = FormatNumber((sim_dt/dt)*100, 2) + "% load"

    # Tab char in between each one
    TextLabel.text = "\u0009".join([sim_time,ms_str,load_string])

    # Debug calculations are included in the sim_dt
    sim_dt = vp.clock() - start_tick