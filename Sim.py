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

def CFrameBox(box, CFrame):
    box.pos = CFrame.Position

    box.up = CFrame.YVector * box.height
    box.axis = -CFrame.ZVector * box.length

def MatrixMul(Matrix1, Matrix2):
    Matrix = []

    for i in range(len(Matrix1)):
        Matrix.append([])
        for j in range(len(Matrix2[i])):

            Sum = 0
            for n in range(len(Matrix1[i])):
                Sum += Matrix1[i][n] * Matrix2[n][j]

            Matrix[i].append(Sum)

    return Matrix

class CFrame:
    def __init__(self, x, y, z, R00, R01, R02, R10, R11, R12, R20, R21, R22):
        self.Position = vp.vec(x, y, z)
        self.RotMatrix = [
            [R00, R01, R02],
            [R10, R11, R12],
            [R20, R21, R22],
        ]

        self.XVector = vp.vec(R00, R10, R20)
        self.YVector = vp.vec(R01, R11, R21)
        self.ZVector = vp.vec(R02, R12, R22)

    def __mul__(self, other):
        if isinstance(other, CFrame):
            pos = MatrixMul(self.RotMatrix, [[other.Position.x], [other.Position.y], [other.Position.z]])
            RotMatrix = MatrixMul(self.RotMatrix, other.RotMatrix)
            return self.fromMatrix(self.Position.x + pos[0][0], self.Position.y + pos[1][0], self.Position.z + pos[2][0], RotMatrix)
        elif isinstance(other, vp.vector):
            pos = MatrixMul(self.RotMatrix, [[other.x], [other.y], [other.z]])
            return vp.vec(pos[0][0], pos[1][0], pos[2][0])
        
        raise Exception(f"Unknown operation \"mul\" on {type(self)} and {type(other)}")

    def __str__(self):
        return f"CFrame:\n{self.Position}\n{self.XVector}\n{self.YVector}\n{self.ZVector}"

    @classmethod
    def new(self, x, y, z):
        return self(x, y, z, 
                    1, 0, 0, 
                    0, 1, 0, 
                    0, 0, 1)

    @classmethod
    def AngleX(self, rx):
        return self(0,0,0,
                    1, 0, 0,
                    0, vp.cos(rx), vp.sin(rx), 
                    0, -vp.sin(rx), vp.cos(rx))

    @classmethod
    def AngleY(self, ry):
        return self(0,0,0, 
                    vp.cos(ry), 0, -vp.sin(ry), 
                    0, 1, 0,
                    vp.sin(ry), 0,  vp.cos(ry))

    # The angle for the rotation matricies is flipped from the traditional, for parity with Roblox CFrames
    @classmethod
    def AngleZ(self, rz):
        return self(0,0,0, 
                    vp.cos(rz), vp.sin(rz), 0, 
                    -vp.sin(rz), vp.cos(rz), 0, 
                    0, 0, 1)
    
    @classmethod
    def fromVectors(self, x, y, z, XVector, YVector, ZVector):
        return self(x, y, z, 
                    XVector[0], YVector[0], ZVector[0], 
                    XVector[1], YVector[1], ZVector[1], 
                    XVector[2], YVector[2], ZVector[2])
    
    @classmethod
    def fromMatrix(self, x, y, z, RotMatrix):
        return self(x, y, z,
                    RotMatrix[0][0], RotMatrix[0][1], RotMatrix[0][2], 
                    RotMatrix[1][0], RotMatrix[1][1], RotMatrix[1][2],  
                    RotMatrix[2][0], RotMatrix[2][1], RotMatrix[2][2])

    def Invert(self):
        print("adfsdfasdfkghasdkjfhgads")

#print(CFrame.AngleY(vp.pi/16))
#print(CFrame.AngleY(vp.pi/4))
#print(CFrame.AngleY(vp.pi/4) * CFrame.AngleY(vp.pi/4))
#print(CFrame.AngleY(vp.pi/2))

Matrix1 = [
            [2, 0, 2],
            [0, 2, 0],
            [0, 0, 2],
        ]

Matrix2 = [
            [0, 0, 0],
            [3, 0, 3],
            [3, 0, 0],
        ]

#print(MatrixMul(Matrix1, Matrix2))

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

HingeCFrame = CFrame.new(0, 2, 0)

def P1_upt_func():
    m1 = P1["m"]
    m2 = P2["m"]

    Num = -g*(2*m1 + m2)*vp.sin(P1["x"]) - m2*g*vp.sin(P1["x"] - 2*P2["x"]) - 2*vp.sin(P1["x"] - P2["x"])*m2*((P2["v"]**2)*art_l + (P1["v"]**2)*l*vp.cos(P1["x"] - P2["x"]))
    Denum = l*(2*m1 + m2 - m2*vp.cos(2*P1["x"] - 2*P2["x"]))
    return Num/Denum
    #P1["a"] = (-1/((l**2)*(m1 + m2))) * (m2*l*art_l*P2["a"]*vp.cos(P2["x"] - P1["x"]) - m2*l*art_l*P2["v"]*vp.sin(P2["x"] - P1["x"]) + (m1 + m2)*g*l*vp.sin(P1["x"]))

def P1_pos_func():
	return HingeCFrame * CFrame.AngleX(P1["x"]) * CFrame.new(0,-l,0)

def P2_upt_func():
    m1 = P1["m"]
    m2 = P2["m"]

    Num = 2*vp.sin(P1["x"] - P2["x"])*((P1["v"]**2)*l*(m1 + m2) + g*(m1 + m2)*vp.cos(P1["x"]) + (P2["v"]**2)*art_l*m2*vp.cos(P1["x"] - P2["x"]))
    Denum = l*(2*m1 + m2 - m2*vp.cos(2*P1["x"] - 2*P2["x"]))
    return Num/Denum
	#(-1/(m2*art_l^2)) * (m2*l*art_l*P1["a"]*math.cos(P2["x"] - P1["x"]) + m2*l*art_l*P1["v"]*math.sin(P2["x"] - P1["x"]) + m2*g*art_l*math.sin(P2["x"]))

def P2_pos_func():
	return P1_pos_func() * CFrame.AngleX(P2["x"] - P1["x"]) * CFrame.new(0,-r,0)

#P2 = None

obj1 = vp.box(color = vp.color.red, size = vp.vec(.5,.1,1))
obj2 = vp.box(color = vp.color.yellow, size = TAILLE_BEATRICE) # texture = "test.png"

P1 = {
    "a": 0,
    "v": 0,
    "x": 1,
    "m" : 20,
    "obj": obj1,
    "upd_func": P1_upt_func,
    "pos_func": P1_pos_func,
}

P2 = {
    "a": 0,
    "v": 0,
    "x": 0,
    "m": MASSE_BEATRICE,
    "obj": obj2,
    "upd_func": P2_upt_func,
    "pos_func": P2_pos_func,
}

objs = [P1, P2]

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

    # Calculate acceleration
    Accelerations = [0] * len(objs)
    for i, v in enumerate(objs):
        Accelerations[i] = v["upd_func"]()

    # Apply the acceleration to velocity and position
    for i, v in enumerate(objs):
        v["a"] = Accelerations[i]
        v["v"] += v["a"]*dt
        v["x"] += v["v"]*dt

    # Update the object's CFrames
    for i, v in enumerate(objs):
        CFrameBox(v["obj"], v["pos_func"]())

    """
    angle_1 = P1["x"] - vp.pi/2 # -pi to have 0 at (1,0) on the unit circle
    P1["obj"].pos = l*vp.vec(vp.cos(angle_1),vp.sin(angle_1),0)
    P1["obj"].axis = vp.vec(vp.cos(angle_1 + vp.pi/2), vp.sin(angle_1 + vp.pi/2), 0) * P1["obj"].length

    angle_2 = P2["x"] - vp.pi/2 # -pi to have 0 at (1,0) on the unit circle
    P2["obj"].pos = (P1["obj"].pos + r*vp.vec(vp.cos(angle_2),r*vp.sin(angle_2),0) 
        #+ r*vp.vec(vp.cos(angle_1 + vp.pi/2), vp.sin(angle_1 + vp.pi/2), 0)
        #+ r*vp.vec(vp.cos(angle_2 + vp.pi/2), vp.sin(angle_2 + vp.pi/2), 0)
    )
    P2["obj"].axis = vp.vec(vp.cos(angle_2 + vp.pi/2), vp.sin(angle_2 + vp.pi/2), 0) * P2["obj"].length
    """


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