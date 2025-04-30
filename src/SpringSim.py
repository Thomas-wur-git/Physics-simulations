##########################################################################################################################################################################################################

import vpython as vp

BackgroundColor = vp.vec(241/255,242/255,244/255)
vp.scene.background = BackgroundColor

#vp.canvas.get_selected().align = "left"

vp.distant_light(direction = vp.vec(0, 1, 0), color = vp.vec(1.0, 0.9, 0.5)*0.8)
#vp.scene.ambien = vp.color.white*0.6

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

Ids = {}
def PrintAtAReasonablePace(Id, *args):
    LastPrint = Ids.get(Id, 0)
    if LastPrint + 1 > vp.clock():
        return
    
    Ids[Id] = vp.clock()
    print(args)

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

    @classmethod
    def lookAt(self, Origin, Target, UpVector):

        ZVector = vp.norm(Target - Origin)
        XVector = vp.norm(vp.cross(ZVector, UpVector))
        YVector = vp.cross(XVector, ZVector)

        return self.fromVectors(Origin.x, Origin.y, Origin.z,
                    [XVector.x, XVector.y, XVector.z], 
                    [YVector.x, YVector.y, YVector.z],  
                    [ZVector.x, ZVector.y, ZVector.z])

    def Invert(self):
        raise Exception(f":Invert() has not yet been implemented")

# // Constants

FPS = 60
SIM_RATE = FPS*12 # Should be a multiple of FPS
SLOWDOWN_MOD = 1 # Slowing down the simulation for debuging resons

g = vp.vec(0,-9.81,0)
AirResistance = .5 # Air resistance
TargetArmLength = .45

K = 80000 # Spring constant links between objects
DAMPNING = 100 # Spring dampning, to stabilize the simulation

SWING_SIZE = vp.vec(.5,.1,.9)

# Inertia for a rectangular prism:
# https://www.physicsblog.in/2025/01/moment-of-inertia-of-rectangle-prism.html

SWING_LENGTH = 2.5
SWING_MASS = 5 # Arbitrairy value
SWING_INERTIA = 5 # I around the center of mass? I think that's what it should be. Unstable if too low

BEATRICE_SIZE = vp.vec(0.2, 1.55, 0.35)

BEATRICE_ROT_OFFSET = 0.3 # Distance from the seat to the center of mass of BEATRICE
BEATRICE_MASS = 46.72
BEATRICE_INERTIA = 30  # Moment of inertia of BEATRICE 

BEATRICE_ARM_OFFSET = -0.35
SWING_ARM_ATTACHMENT = 0.5
ARM_SPRING_K = 3000
ARM_SPRING_DAMPNING = 800

# // Force functions

# Calculates the force between two points, from 2 objects
def CalculateSpringForceOnObjs(obj1, obj2, Attachment1Offset, Attachment2Offset, K, D, RestLength):

    # // Spring force

    Attachment1Pos = (obj1["CFrame"] * Attachment1Offset).Position
    Attachment2Pos = (obj2["CFrame"] * Attachment2Offset).Position

    # Vector from 2 to 1
    DistanceVector = Attachment1Pos - Attachment2Pos
    UnitVector = vp.norm(DistanceVector)

    # Force on attachment 1 (equal and opposite for attachment 2)
    SpringForce = - K * (vp.mag(DistanceVector) - RestLength)*UnitVector

    # // Dampning Forces

    Attachment1Velocity = obj1["v"] - vp.cross(obj1["r_v"], (obj1["CFrame"] * Attachment1Offset).Position - obj1["CFrame"].Position)
    Attachment2Velocity = obj2["v"] - vp.cross(obj2["r_v"], (obj2["CFrame"] * Attachment2Offset).Position - obj2["CFrame"].Position)

    # Velocity of Attachment1 relative to Attachment2
    VelocityVector = Attachment1Velocity - Attachment2Velocity
    ProjVelocityVector = vp.proj(VelocityVector, DistanceVector)

    # Force on attachment 1 (equal and opposite for attachment 2)
    DampningForce = - ProjVelocityVector*D

    # // Applying the forces

    Force = SpringForce + DampningForce

    obj1["Forces"] += Force
    obj2["Forces"] -= Force

    obj1["MomentForces"] += vp.cross(obj1["CFrame"].Position - Attachment1Pos, Force)
    obj2["MomentForces"] += vp.cross(obj2["CFrame"].Position - Attachment2Pos, -Force)

def CalculateAirResistance(obj):
    # Does not consider angular speed, cuz that would probably be a pain
    obj["Forces"] -= obj["v"]*vp.mag(obj["v"])*AirResistance

def CalculateGravity(obj):
    obj["Forces"] += obj["m"]*g

# // Coordinates

# add_l is a visual offset to the swing to make the green blob on top of the swing
# but with the springs, there is some sag, and so the green bloc is still not quite on the swing...
add_l = BEATRICE_SIZE.x/2 + SWING_SIZE.y/2

HingeCFrame = CFrame.new(0, 2, 0)

def P1_pos_func(MoveParts):
    RotCFrame = CFrame.new(P1["x"].x,P1["x"].y,P1["x"].z) * CFrame.AngleX(P1["r_x"])
    PoleCFrame = RotCFrame * CFrame.new(0,(SWING_LENGTH+add_l)/2,0)

    P1["CFrame"] = RotCFrame

    if MoveParts:
        CFrameBox(P1["obj"], RotCFrame)
        CFrameBox(P1["ExtraData"]["Poles"][0], PoleCFrame * CFrame.new(-SWING_SIZE.z/2 + .1, 0, 0))
        CFrameBox(P1["ExtraData"]["Poles"][1], PoleCFrame * CFrame.new(SWING_SIZE.z/2 - .1, 0, 0))

def P2_pos_func(MoveParts):
    HumanCFrame = CFrame.new(P2["x"].x,P2["x"].y,P2["x"].z) * CFrame.AngleX(P2["r_x"])

    P2["CFrame"] = HumanCFrame

    if MoveParts:
        CFrameBox(P2["obj"], HumanCFrame)
        CFrameBox(P2["ExtraData"]["Face"], HumanCFrame * CFrame.new(0, -BEATRICE_SIZE.y/2 + 0.3/2 + 0.05/2, BEATRICE_SIZE.x/2 - .08/2 + 10**-2) * CFrame.AngleZ(vp.pi))

        ArmHeightCFrame = HumanCFrame * CFrame.new(0,BEATRICE_ARM_OFFSET,0)
        SwingHeightCFrame = P1["CFrame"] * CFrame.new(0,SWING_ARM_ATTACHMENT + add_l,0)

        Arm1StartPos = (ArmHeightCFrame * CFrame.new(-BEATRICE_SIZE.x/2,0,0)).Position
        Arm1EndPos = (SwingHeightCFrame * CFrame.new(-SWING_SIZE.z/2 + .1, 0, +.05)).Position

        Arm2StartPos = (ArmHeightCFrame * CFrame.new(BEATRICE_SIZE.x/2,0,0)).Position
        Arm2EndPos = (SwingHeightCFrame * CFrame.new(SWING_SIZE.z/2 - .1, 0, +.05)).Position

        ArmLength = vp.mag(Arm1EndPos - Arm1StartPos)

        P2["ExtraData"]["Arms"][0].length = ArmLength
        P2["ExtraData"]["Arms"][1].length = ArmLength

        CFrameBox(P2["ExtraData"]["Arms"][0], CFrame.lookAt(Arm1StartPos, Arm1EndPos, vp.vec(0,1,0)) * CFrame.new(0,0,ArmLength/2))
        CFrameBox(P2["ExtraData"]["Arms"][1], CFrame.lookAt(Arm2StartPos, Arm2EndPos, vp.vec(0,1,0)) * CFrame.new(0,0,ArmLength/2))

Swing = vp.box(color = Colors["Permisson"], size = SWING_SIZE)
SwingPole1 = vp.box(color = Colors["Flint"], size = vp.vec(.1,SWING_LENGTH + add_l,.1))
SwingPole2 = vp.box(color = Colors["Flint"], size = vp.vec(.1,SWING_LENGTH + add_l,.1))

Hooman = vp.box(color = Colors["Br. yellowish green"], size = BEATRICE_SIZE)
HoomanFace = vp.box(color = vp.color.white, size = vp.vec(.08,.3,.3), texture = "BeatriceEvil.png")
Arm1 = vp.box(color = Colors["Br. yellowish green"], size = vp.vec(1,1,1)*0.12)
Arm2 = vp.box(color = Colors["Br. yellowish green"], size = vp.vec(1,1,1)*0.12)

INITIAL_ANGLE = 0

P1 = {
    "a": vp.vec(0,0,0),
    "v": vp.vec(0,0,0),
    "x": HingeCFrame.Position + vp.vec(0,-SWING_LENGTH*vp.cos(INITIAL_ANGLE),SWING_LENGTH*vp.sin(INITIAL_ANGLE)),
    "r_a": vp.vec(0,0,0),
    "r_v": vp.vec(0,0,0),
    "r_x": INITIAL_ANGLE,
    "m" : SWING_MASS,
    "I" : SWING_INERTIA,
    "obj": Swing,
    "Forces" : vp.vec(0,0,0),
    "MomentForces" : vp.vec(0,0,0),
    "CFrame" : CFrame.new(0,0,0),
    "ExtraData": {
        "Poles": [SwingPole1, SwingPole2]
    }
}

INITIAL_ANGLE = vp.pi*0.6

P2 = {
    "a": vp.vec(0,0,0),
    "v": vp.vec(0,0,0),
    "x": P1["x"] + vp.vec(0,-BEATRICE_ROT_OFFSET*vp.cos(INITIAL_ANGLE),BEATRICE_ROT_OFFSET*vp.sin(INITIAL_ANGLE)),
    "r_a": vp.vec(0,0,0),
    "r_v": vp.vec(0,0,0),
    "r_x": INITIAL_ANGLE,
    "m" : BEATRICE_MASS,
    "I" : BEATRICE_INERTIA,
    "obj": Hooman,
    "Forces" : vp.vec(0,0,0),
    "MomentForces" : vp.vec(0,0,0),
    "CFrame" : CFrame.new(0,0,0),
    "ExtraData": {
        "Face": HoomanFace,
        "Arms": [Arm1, Arm2]
    }
}

# Running these to update "CFrame"
P1_pos_func(True)
P2_pos_func(True)

HingeObj = {
    "a": vp.vec(0,0,0),
    "v": vp.vec(0,0,0),
    "x": HingeCFrame.Position,
    "r_a": vp.vec(0,0,0),
    "r_v": vp.vec(0,0,0),
    "r_x": 0,
    "m" : 0,
    "I" : 0,
    "obj": None,
    "Forces" : vp.vec(0,0,0),
    "MomentForces" : vp.vec(0,0,0),
    "CFrame" : HingeCFrame,
    "ExtraData": {}
}

objs = [P1, P2]

ForceFunctions = [
    # Gravity
    lambda: CalculateGravity(P1),
    lambda: CalculateGravity(P2),

    # AirResistance
    lambda: CalculateAirResistance(P1),
    lambda: CalculateAirResistance(P2),

    # Springs
    lambda: CalculateSpringForceOnObjs(P1, HingeObj, CFrame.new(0,(SWING_LENGTH+add_l),0), CFrame.new(0,0,0), K, DAMPNING, 0),
    lambda: CalculateSpringForceOnObjs(P1, P2, CFrame.new(0,0,0), CFrame.new(0,(BEATRICE_ROT_OFFSET),0), K, DAMPNING, 0),
    lambda: CalculateSpringForceOnObjs(P1, P2, CFrame.new(0,SWING_ARM_ATTACHMENT,0), CFrame.new(0,(BEATRICE_ARM_OFFSET),0), ARM_SPRING_K, ARM_SPRING_DAMPNING, TargetArmLength),
]

PositionFunctions = [
    P1_pos_func,
    P2_pos_func,
]

# // Oscillations per second counter

ThresoldSpeed = 0.1
Debounce = 0.5

TimestampStart = 0
Timestamps = [0,0,0]

# Number of frames. This is decreased until it hits 0, when it isn't stopped
IsStopped = 0
IsAtThreshold = False

def AnalyzeOscillations(t):
    global Threshold, TimestampStart, IsStopped, IsAtThreshold

    IsAtThreshold = False

    if vp.mag(P1["v"]) < ThresoldSpeed:
        IsStopped = vp.floor(SIM_RATE*Debounce)
        IsAtThreshold = True

    elif IsStopped > 0:

        IsStopped -= 1
        if IsStopped == 0:
            Timestamps.append(t)
            Timestamps.pop(0)


# // Events & UI

Running = True

def Exit():
    global Running
    Running = False

vp.button(text = "<b>Exit</b>", pos = vp.scene.title_anchor, bind = Exit)

# // Loop //

# Label

TextLabel = vp.wtext(text="", pos = vp.scene.title_anchor)

ArmLengthLabel = vp.wtext(text = "Longeur des bras: " + FormatNumber(TargetArmLength,2))

def ArmLenght(event):
    global TargetArmLength
    TargetArmLength = event.value
    ArmLengthLabel.text = "Longeur des bras: " + FormatNumber(TargetArmLength,2)

vp.slider(bind = ArmLenght, max = 1, min = 0, step = 0.01, value = TargetArmLength)

GravityLabel = vp.wtext(text = "Gravité: 1.0g")

def GravityModifier(event):
    global g
    g = vp.vec(0,-9.81,0)*event.value
    GravityLabel.text = "Gravité: " + FormatNumber(event.value,1) + "g"

vp.slider(bind = GravityModifier, max = 3, min = 0, step = 0.1, value = 1)

AirResistanceLabel = vp.wtext(text = "Résistance de l'air: " + FormatNumber(AirResistance,1))

def AirResistanceModifier(event):
    global AirResistance
    AirResistance = event.value
    AirResistanceLabel.text = "Résistance de l'air: " + FormatNumber(AirResistance,1)

vp.slider(bind = AirResistanceModifier, max = 100, min = 0, step = 0.1, value = AirResistance)

# Loop variables
sim_dt = 0
step = 0

n = SIM_RATE//60
dt = 1/(60*n)

while Running:

    vp.rate(FPS)

    start_tick = vp.clock()

    for i in range(n):

        step += 1

        # Force calculations
        for func in ForceFunctions:
            func()

        # Apply the effects of the force
        for i, v in enumerate(objs):
            v["a"] = v["Forces"]/v["m"]
            v["v"] += v["a"]*dt
            v["x"] += v["v"]*dt

            v["r_a"] = v["MomentForces"]/v["I"]
            v["r_v"] += v["r_a"]*dt
            v["r_x"] += v["r_v"].x*dt # TODO - Make this a vector as well maybe

            v["Forces"] = vp.vec(0,0,0)
            v["MomentForces"] = vp.vec(0,0,0)

        # Reposition objects
        for func in PositionFunctions:
            func(False)

        AnalyzeOscillations(step*dt)

    # Debug labels

    # Reposition objects
    for func in PositionFunctions:
        func(True)

    sim_time = FormatNumber(step*dt, 3) + " secondes"
    ms_str = FormatNumber(sim_dt*1000, 3) + " ms"
    load_string = FormatNumber((sim_dt*FPS)*100, 2) + "% load"

    sec_per_osc = ("*" if IsAtThreshold else " ") + FormatNumber(Timestamps[-1] - Timestamps[-3], 2) + " sec/osc"

    TextLabel.text = "  " + "  ".join([sim_time,ms_str,load_string,sec_per_osc])

    # Debug calculations are included in the sim_dt
    sim_dt = vp.clock() - start_tick