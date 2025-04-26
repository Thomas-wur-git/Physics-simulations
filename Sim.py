##########################################################################################################################################################################################################

import vpython as vp
import math

BackgroundColor = vp.vec(241/255,242/255,244/255)
vp.scene.background = BackgroundColor

vp.canvas.get_selected().align = "left"

# // Utils //

Colors = {
    "Permisson": vp.vec(1, 0.34902, 0.34902),
    "Electric blue": vp.vec(0.0352941, 0.537255, 0.811765),
    "Br. yellowish green": vp.vec(0.643137, 0.741176, 0.278431),
    "Baby blue": vp.vec(0.596078, 0.760784, 0.858824),
    "Flint": vp.vec(0.411765, 0.4, 0.360784),
}

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
        raise Exception(f":Invert() has not yet been implemented")

# // Constants

FRICTION = 0.01 # Air resistance

K = 80000 # Spring constant links between objects
DAMPNING = 80 # Spring dampning, to stabilize the simulation

# TODO - For some reason, the DAMPNING is slowing down the angular oscillations of the pendulum, why is that?

SIM_RATE = 240 # hz

g = vp.vec(0,-9.81,0)

SWING_SIZE = vp.vec(.5,.1,.9)

SWING_LENGTH = 2.5
SWING_MASS = 5
SWING_INERTIA = SWING_MASS * SWING_LENGTH**2

BEATRICE_SIZE = vp.vec(0.2, 1.55, 0.35)

BEATRICE_ROT_OFFSET = 0.3 # Distance from the seat to the center of mass of BEATRICE
BEATRICE_INERTIA = 40 # Moment of inertia of BEATRICE (need real mesurments)
BEATRICE_MASS = 46.72

# // Coordinates

HingeCFrame = CFrame.new(0, 2, 0)

def P1_force_func(self):

    # Spring from hinge to swing

    RotCFrame = CFrame.new(self["x"].x,self["x"].y,self["x"].z) * CFrame.AngleX(P1["r_x"])
    PoleTopCFrame = RotCFrame * CFrame.new(0,(SWING_LENGTH+add_l),0)

    DistanceVec = PoleTopCFrame.Position - HingeCFrame.Position
    DistanceVecNorm = vp.norm(DistanceVec)

    HingeForce = - vp.norm(DistanceVec)*(vp.mag(DistanceVec))*K # Attachment to the hinge

    # Dampning for spring 1

    Spring_V1 = vp.dot(self["v"], DistanceVecNorm)*DistanceVecNorm

    # Spring from swing to Hooman

    RotCFrame = CFrame.new(P2["x"].x,P2["x"].y,P2["x"].z) * CFrame.AngleX(P2["r_x"])
    Attachment = RotCFrame * CFrame.new(0,(BEATRICE_ROT_OFFSET),0)

    DistanceVec = Attachment.Position - P1["CFrame"].Position
    DistanceVecNorm = vp.norm(DistanceVec)

    SpringForce = + vp.norm(DistanceVec)*(vp.mag(DistanceVec))*K # Attachment to Hooman

    # Dampning for spring 2

    Spring_V2 = -vp.dot(self["v"] - P1["v"], DistanceVecNorm)*DistanceVecNorm

    Force = (
        g * self["m"] 
        + HingeForce
        + SpringForce
        - self["v"]*vp.mag(self["v"])*FRICTION # Using v squared for more realistic drag
        - Spring_V1*DAMPNING # Dampning in the direction of the spring
        - Spring_V2*DAMPNING # Dampning in the direction of the spring
    )

    # τ = F x r
    Tau = (
        vp.cross(CFrame.AngleX(P1["r_x"]).YVector * SWING_LENGTH, HingeForce)
    )

    return Force, -Tau.x

def P1_pos_func(self):
    add_l = BEATRICE_SIZE.x/2 + SWING_SIZE.y/2
    
    RotCFrame = CFrame.new(self["x"].x,self["x"].y,self["x"].z) * CFrame.AngleX(P1["r_x"])
    PoleCFrame = RotCFrame * CFrame.new(0,(SWING_LENGTH+add_l)/2,0)

    self["CFrame"] = RotCFrame

    CFrameBox(P1["obj"], RotCFrame)
    CFrameBox(P1["ExtraData"]["Poles"][0], PoleCFrame * CFrame.new(-SWING_SIZE.z/2 + .1, 0, 0))
    CFrameBox(P1["ExtraData"]["Poles"][1], PoleCFrame * CFrame.new(SWING_SIZE.z/2 - .1, 0, 0))

def P2_force_func(self):

    RotCFrame = CFrame.new(self["x"].x,self["x"].y,self["x"].z) * CFrame.AngleX(P2["r_x"])
    Attachment = RotCFrame * CFrame.new(0,(BEATRICE_ROT_OFFSET),0)

    DistanceVec = Attachment.Position - P1["CFrame"].Position
    DistanceVecNorm = vp.norm(DistanceVec)

    SpringForce = - vp.norm(DistanceVec)*(vp.mag(DistanceVec))*K # Attachment to the swing

    Spring_V = vp.dot(self["v"] - P1["v"], DistanceVecNorm)*DistanceVecNorm

    Force = (
        g * self["m"] 
        + SpringForce
        - self["v"]*vp.mag(self["v"])*FRICTION # Using v squared for more realistic drag
        - Spring_V*DAMPNING # Dampning in the direction of the spring
    )

    # τ = F x r
    Tau = (
        vp.cross(CFrame.AngleX(P2["r_x"]).YVector * BEATRICE_ROT_OFFSET, SpringForce)
    )

    return Force, -Tau.x


def P2_pos_func(self):
    HumanCFrame = CFrame.new(self["x"].x,self["x"].y,self["x"].z) * CFrame.AngleX(P2["r_x"])

    self["CFrame"] = HumanCFrame

    CFrameBox(P2["obj"], HumanCFrame)
    CFrameBox(P2["ExtraData"]["Face"], HumanCFrame * CFrame.new(0, -BEATRICE_SIZE.y/2 + 0.3/2 + 0.05/2, BEATRICE_SIZE.x/2 - .08/2 + 10**-3) * CFrame.AngleZ(vp.pi))

#P2 = None

add_l = BEATRICE_SIZE.x/2 + SWING_SIZE.y/2

Swing = vp.box(color = Colors["Permisson"], size = SWING_SIZE)
SwingPole1 = vp.box(color = Colors["Flint"], size = vp.vec(.1,SWING_LENGTH + add_l,.1))
SwingPole2 = vp.box(color = Colors["Flint"], size = vp.vec(.1,SWING_LENGTH + add_l,.1))

Hooman = vp.box(color = Colors["Br. yellowish green"], size = BEATRICE_SIZE)
HoomanFace = vp.box(color = vp.color.white, size = vp.vec(.08,.3,.3), texture = "BeatriceEvil.png")

INITIAL_ANGLE = 1

P1 = {
    "a": vp.vec(0,0,0),
    "v": vp.vec(0,0,0),
    "x": HingeCFrame.Position + vp.vec(0,-SWING_LENGTH*vp.cos(INITIAL_ANGLE),SWING_LENGTH*vp.sin(INITIAL_ANGLE)),
    "r_a": 0,
    "r_v": 0,
    "r_x": INITIAL_ANGLE,
    "m" : SWING_MASS,
    "I" : SWING_INERTIA,
    "obj": Swing,
    "CFrame" : CFrame.new(0,0,0),
    "force_func": P1_force_func,
    "pos_func": P1_pos_func,
    "ExtraData": {
        "Poles": [SwingPole1, SwingPole2]
    }
}

INITIAL_ANGLE = 2

P2 = {
    "a": vp.vec(0,0,0),
    "v": vp.vec(0,0,0),
    "x": P1["x"] + vp.vec(0,-BEATRICE_ROT_OFFSET*vp.cos(INITIAL_ANGLE),BEATRICE_ROT_OFFSET*vp.sin(INITIAL_ANGLE)),
    "r_a": 0,
    "r_v": 0,
    "r_x": INITIAL_ANGLE,
    "m" : BEATRICE_MASS,
    "I" : BEATRICE_INERTIA,
    "obj": Hooman,
    "CFrame" : CFrame.new(0,0,0),
    "force_func": P2_force_func,
    "pos_func": P2_pos_func,
    "ExtraData": {
        "Face": HoomanFace
    }
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
step = 0

while Running:
    step += 1

    start_tick = vp.clock()

    # Calculate acceleration
    Forces = [0] * len(objs)
    MomentForces = [0] * len(objs)
    for i, v in enumerate(objs):
        Forces[i], MomentForces[i] = v["force_func"](v)

    # Apply the acceleration to velocity and position
    for i, v in enumerate(objs):
        v["a"] = Forces[i]/v["m"]
        v["v"] += v["a"]*dt
        v["x"] += v["v"]*dt

        v["r_a"] = MomentForces[i]/v["I"]
        v["r_v"] += v["r_a"]*dt
        v["r_x"] += v["r_v"]*dt

    # Update the object's CFrames
    for i, v in enumerate(objs):
        v["pos_func"](v)

    vp.rate(SIM_RATE)

    #MiscAnimation()

    t = step*dt

    #PlotGraph(t)

    # Debug labels

    sim_time = FormatNumber(t, 3) + " secondes"
    ms_str = FormatNumber(sim_dt*1000, 3) + " ms"
    load_string = FormatNumber((sim_dt/dt)*100, 2) + "% load"

    # Tab char in between each one
    TextLabel.text = "\u0009".join([sim_time,ms_str,load_string])

    # Debug calculations are included in the sim_dt
    sim_dt = vp.clock() - start_tick