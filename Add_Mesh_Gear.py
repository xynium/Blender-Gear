#
# ____  ________.___._______  .___ ____ ___  _____
# \   \/  /\__  |   |\      \ |   |    |   \/     \
#  \     /  /   |   |/   |   \|   |    |   /  \ /  \
#  /     \  \____   /    |    \   |    |  /    Y    \
# /___/\  \ / ______\____|__  /___|______/\____|__  /
#       \_/ \/              \/                    \/
#
# last mod 20 dec 2019 change to 2.8
#
# Feedback xynium@laposte.net


# notes :
#  Backlash should be arond .99
#  With 20deg of pressure angle minimal tooth number is 13 check pinion
#
#  For crown train add backlash 0.99

# TODO For rack close mesh has to close the end of the rack
# TODO Has we dublicate a whole tooth there is vertices in common -> need a cleanup for duplicate vertice

import bpy
# from mathutils import *
from math import sin, cos, tan, pi, atan, sqrt
# from bpy.props import *

bl_info = {
    "name": "Gear",
    "description": "Add a gear mesh.",
    "author": "XYNIUM JP Lathuile",
    "version": (1, 0),
    "blender": (2, 80, 0),
    "location": "View3D > Add > Mesh",
    "warning": "",
    "support": "TESTING",
    "category": "Add Mesh",
}

global PitchRadius, TeethN, PressureAng, Addendum, Dedendum, Fillet, Resolution, Thickness
global HelicalAng, Width, LongRes, prop, Alpha, toDel, N, Pinion, CloseObj, jeu
global nIntoth, ax


#####################################
#   Main Add Mesh Gear + Pinion
#####################################


class ObjectGear(bpy.types.Operator):
    bl_idname = "mesh.primitive_gear"
    bl_label = "Add a Gear"
    bl_options = {'REGISTER',  'UNDO',  'PRESET'}

    global toDel

    toDel = 0
    # default blender is mm
    UseMod = bpy.props.BoolProperty(name="Use module", default=True)
    Pinions = bpy.props.BoolProperty(name="Pinion", default=False)
    Adv = bpy.props.BoolProperty(name="Advenced setting", default=False)
    CloseObjs = bpy.props.BoolProperty(name="Close Mesh", default=False)
    TeethNs = bpy.props.IntProperty(name="Teeth Number", default=40, min=9, soft_max=200)
    PitchRadiuss = bpy.props.FloatProperty(name="Reference Radius", min=0.5, soft_max=200,  default=20, unit='LENGTH',
                                           description="Pitch Radius")
    PressureAngs = bpy.props.IntProperty(name="Pressure angle",  default=20, min=10, max=40)
    HelicalAngs = bpy.props.IntProperty(name="Tooth angle", default=0, min=-45, max=45,
                                        description="Angle of the tooth in an Helical (Spiral) gear tooth")
    Addendums = bpy.props.FloatProperty(name="Addendum ", min=0.001, soft_max=30, default=1, precision=1, unit='LENGTH',
                                        description="Addendum")
    Dedendums = bpy.props.FloatProperty(name="Dedendum", min=0.001, soft_max=30, default=1, precision=1, unit='LENGTH',
                                        description="Dedendum")
    Fillets = bpy.props.FloatProperty(name="Fillet", min=0.0001, soft_max=1, default=0.05, precision=2, unit='LENGTH',
                                      description="Fillet")
    Resolutions = bpy.props.IntProperty(name="Resolution", default=3, min=1, max=7,
                                        description="Resolution of the tooth profile")
    LongRess = bpy.props.IntProperty(name="Long. Res.", default=5, min=1, max=25,
                                     description="Resolution along the tooth extrusion")
    Thicknesss = bpy.props.FloatProperty(name="Thickness", default=5,  min=0.05,  soft_max=200, unit='LENGTH',
                                         description="Thickness of the tooth extrusion")
    Widths = bpy.props.FloatProperty(name="Axis", min=0.01, soft_max=200, default=6, unit='LENGTH', precision=1,
                                     description="Concentrical inner extrusion from Inner tooth")
    props = bpy.props.EnumProperty(name="Type", items=[("nor", "Normal", "Spur gear"), ("con", "Connical",
                                   "Connical gear"), ("cro", "Crown", "Crown"), ("rac", "Rack", "Rack")], default="nor")
    Alphas = bpy.props.IntProperty(name="Cone Angle", default=90, min=10, max=180,
                                   description="Full Angle at vertex cone")
    Ns = bpy.props.FloatProperty(name="Multiply Rate", default=2, min=0.01, soft_max=30, description="Multipling rate")
    jeus = bpy.props.FloatProperty(name="Backlash", default=1, min=0.6, max=1, precision=3, description="Backlash")

    def draw(self,  context):
        layout = self.layout

        box = layout.box()
        box.prop(self,  'UseMod')
        box.prop(self,  'Pinions')
        box.prop(self,  'Adv')
        box.prop(self,  'CloseObjs')
        box.prop(self,  'props')
        box.prop(self,  'TeethNs')
        box.prop(self,  'PitchRadiuss')
        box.prop(self,  'Thicknesss')
        box.prop(self,  'Widths')

        if self.Adv is True:
            boxa = layout.box()
            boxa.prop(self, 'PressureAngs')
            boxa.prop(self, 'HelicalAngs')
            boxa.prop(self, 'jeus')
            boxa.prop(self, 'Resolutions')
            boxa.prop(self, 'LongRess')

        if self.UseMod is False:
            boxb = layout.box()
            boxb.prop(self, 'Addendums')
            boxb.prop(self, 'Dedendums')
            boxb.prop(self, 'Fillets')

        if self.props is "con":
            boxc = layout.box()
            boxc.prop(self, 'Alphas')

        if self.Pinions is True:
            boxd = layout.box()
            boxd.prop(self,  'Ns')

    def execute(self,  context):
        global PitchRadius, TeethN, PressureAng, Addendum, Dedendum, Fillet, Bevel, Resolution, Thickness
        global prop, Alpha, toDel, N, Pinion, CloseObj, jeu, LongRes, HelicalAng, Width

        # create mesh
        PitchRadius = self.PitchRadiuss
        TeethN = self.TeethNs
        PressureAng = self.PressureAngs
        Addendum = self.Addendums
        Dedendum = self.Dedendums
        Fillet = self.Fillets
        Bevel = 0
        Resolution = self.Resolutions
        Thickness = self.Thicknesss
        LongRes = self.LongRess
        HelicalAng = self.HelicalAngs
        Width = self.Widths
        Alpha = self.Alphas
        prop = self.props
        N = self.Ns
        Pinion = self.Pinions
        CloseObj = self.CloseObjs
        jeu = self.jeus

        if self.UseMod is True:
            module = 2.0 * PitchRadius / TeethN
            Addendum = module
            Dedendum = 1.25 * module
            Fillet = 0.3 * module
            Bevel = 0.25 * module
        if (prop == "con") and (Pinion is True):
            # Compute cone angle
            pTeethN = int(TeethN / N)
            N1 = TeethN/pTeethN
            Alpha = Alpha*pi/180.0
            Alpha1 = atan(sin(Alpha)/(N1+cos(Alpha)))
            Alpha = 2*(Alpha-Alpha1)*180/pi
        verts, faces = add_gear()
        if toDel == 1:
            #bpy.ops.object.mode_set(mode='OBJECT')   #add 1712
            bpy.ops.object.delete()

        obj = create_mesh_object(context,  verts,  [],  faces,  "Gear")
        toDel = 1
        # duplique les dents spin
        bpy.ops.object.mode_set(mode='EDIT')
        if (prop == "rac"):
            DiametralPitch = 2*pi*PitchRadius/TeethN
            for k in range(TeethN):
                bpy.ops.mesh.duplicate_move(TRANSFORM_OT_translate={"value": (0, DiametralPitch, 0)})
        else:
            bpy.ops.mesh.spin(steps=TeethN, dupli=True, angle=2*pi, center=bpy.context.scene.cursor.location,axis=(0.0, 0.0, 1.0))
            bpy.ops.mesh.delete(type='VERT')
        bpy.ops.object.mode_set(mode='OBJECT')
        if Pinion is True:
                Rg = PitchRadius
                pTeethN = int(TeethN / N)
                Rp = PitchRadius * pTeethN / TeethN
                if (prop == "con"):
                    Alpha = 2*Alpha1*180/pi
                PitchRadius = Rp
                TeethN = pTeethN
                HelicalAng = -HelicalAng
                verts, faces = add_gear(pinion=True)
                objp = create_mesh_object(context, verts, [], faces, "Pinion")
                # Move center
                paxis = (0.0, 0.0, 1.0)
                #changed obj to bpy.context. le 19/12
                #objectToSelect = bpy.data.objects["objectName"]
                #objectToSelect.select_set(True)
                #bpy.context.view_layer.objects.active = obj #2012 selcted object suposed to be pinion
                print('Debug test')   # to be displayed launch from terminal

                if (prop == "nor" or prop == "con"):
                    pcenter = (bpy.context.object.location[0]+Rp+Rg, bpy.context.object.location[1], bpy.context.object.location[2])
                if (prop == "cro"):
                    pcenter = (bpy.context.object.location[0]+Rg-Rp, bpy.context.object.location[1], bpy.context.object.location[2])
                if (prop == "rac"):
                    pcenter = (bpy.context.object.location[0]+Rp, bpy.context.object.location[1], bpy.context.object.location[2])
                #bpy.context.view_layer.objects.active = objp    #2012
                bpy.context.object.location = pcenter  #changed objp to bpy.context. le 19/12
                bpy.ops.object.mode_set(mode='EDIT')
                bpy.ops.mesh.spin(steps=TeethN,  dupli=True, angle=2*pi, center=pcenter,  axis=paxis)
                bpy.ops.mesh.delete(type='VERT')
                bpy.ops.object.mode_set(mode='OBJECT')
                if (prop == "con"):
                    bpy.ops.transform.rotate(value=-self.Alphas*pi/180.0, orient_axis='Y')
                    bpy.context.object.location[0] = -Rg  # bpy.context.object.location[0] + Rg
                    bpy.context.object.location[2] = Rp  # bpy.context.object.location[2] + Rp
        return {'FINISHED'}


class INFO_MT_mesh_gear_add(bpy.types.Menu):
    # Define the "Gears" menu
    bl_idname = "INFO_MT_mesh_gear_add"
    bl_label = "Gear"

    def draw(self,   context):
        layout = self.layout
        layout.operator_context = 'INVOKE_REGION_WIN'
        layout.operator("mesh.primitive_gear",    text="Mesh Gear")
        # layout.operator("mesh.primitive_animate_gear",   text="Animate Gear")


# Define "Extras" menu
def menu_func(self,   context):
        self.layout.operator(ObjectGear.bl_idname, text="Gear",   icon="PLUGIN")


def register():
        bpy.utils.register_class(ObjectGear)
        bpy.types.VIEW3D_MT_mesh_add.append(menu_func)


def unregister():
        bpy.utils.unregister_class(ObjectGear)
        bpy.types.VIEW3D_MT_mesh_add.remove(menu_func)

if __name__ == "__main__":
        register()


# Create a new mesh (object) from verts/edges/faces.
# verts/edges/faces ... List of vertices/edges/faces for the
#                       new mesh (as used in from_pydata).
# name ... Name of the new mesh (& object).

def create_mesh_object(context,   verts,  edges,  faces,  name):

    # Create new mesh
    mesh = bpy.data.meshes.new(name)

    # Make a mesh from a list of verts/edges/faces.
    mesh.from_pydata(verts,  edges,  faces)

    # Update mesh geometry after adding stuff.
    mesh.update()

    from bpy_extras import object_utils
    return object_utils.object_data_add(context,  mesh,  operator=None)


####################################################################
# Do The Gear
####################################################################

def add_gear(pinion=False):
    global PitchRadius, Alpha, prop, N, Pinion, TeethN, HelicalAng

    ####################################################################
    # Main Gear
    ####################################################################

    if (prop == "nor" or pinion is True or prop == "con"):
        (verts, faces) = TheTooth(0, 0, piniont=pinion)
        if (prop == "con"):
            AlphaG = Alpha*pi/360.0
            verts = ConifyTooth(verts,  AlphaG)

    if (prop == "cro" and pinion is False):
        (verts, faces) = TheTooth(0, 1)

    if (prop == "rac" and pinion is False):
        (verts, faces) = TheTooth(1, 0)

    return verts, faces


####################################################################
# CREATES THE BASE PROFILE
####################################################################

def TheTooth(Rack, Crown, piniont=False):

    verts = []
    faces = []

    if (Rack == 1):
        (Vert,  Norm) = RackOutline()
    else:
        (Vert) = ToothOutline(piniont)

    b = 0
    Thicn = Thickness / 2
    Psi = HelicalAng*pi/180.0
    # Special po les engrens droit le nombre de vertice est special
    if prop == "nor" or prop == "con" or piniont is True:
        ###################################################################
        # make thickness
        for i in range(LongRes):
            z = Thicn - (Thickness)*(i)/(LongRes-1)
            sq = z*tan(Psi)
            Phi = sq/PitchRadius

            for j in range(len(Vert)):
                p = Vert[j]

                x = p[0]
                y = p[1]
                verts.append((x*cos(Phi)+y*sin(Phi), -x*sin(Phi)+y*cos(Phi), z))

        ####################################################################
        # make Faces
        NV = len(Vert)
        NL = LongRes-1
        N2 = Resolution
        FVC = len(verts)

        for i in range(NL):
            for j in range(NV-1):
                faces.append((i*(NV)+j,  i*(NV)+j+1,  (i+1)*(NV)+j+1,  (i+1)*(NV)+j))
        ####################################################################
        # make top
        p0 = verts[0]
        p1 = verts[NV-1]
        th0 = atan(p0[1]/p0[0])
        th1 = atan(p1[1]/p1[0])
        R = PitchRadius
        N = 2*Resolution
        Wid = PitchRadius-Width

        for i in range(N+1):
            th = (th1-th0)*i/(N)+th0
            verts.append(((R-Wid)*cos(th), (R-Wid)*sin(th), Thickness/2.0))
        verts.append((R*cos((th1+th0)/2.0), R*sin((th1+th0)/2.0), Thickness/2.0))

        for i in range(N2-1):
            faces.append((i, i+1, FVC+i+1, FVC+i))

        faces.append((N2-1, FVC+N+1, FVC+N2, FVC+N2-1))

        for i in range(2*nIntoth+N-1):
            faces.append((N2+i-1, N2+i, FVC+N+1))

        faces.append((FVC+N2, FVC+N+1, NV-N2, FVC+N2+1))

        for i in range(N2-1):
            faces.append((NV-N2+i, NV-N2+i+1, FVC+N2+i+2, FVC+N2+i+1))

        # make bottom
        F1 = NL*NV
        FVT = FVC
        FVC = len(verts)
        tha = th1
        th1 = -th0
        th0 = -tha

        for i in range(N+1):
            th = ((th1-th0)*i/(N)+th0)
            verts.append(((R-Wid)*cos(th), (R-Wid)*sin(th), -Thickness/2.0))
        verts.append((R*cos((th1+th0)/2.0), R*sin((th1+th0)/2.0), -Thickness/2.0))

        for i in range(N2-1):
            faces.append((F1+i, F1+i+1, FVC+i+1, FVC+i))

        faces.append((F1+N2-1, FVC+N+1, FVC+N2, FVC+N2-1))

        for i in range(2*nIntoth+N-1):
            faces.append((F1+N2+i-1, F1+N2+i, FVC+N+1))

        faces.append((FVC+N2, FVC+N+1, F1+NV-N2, FVC+N2+1))

        for i in range(N2-1):
            faces.append((F1+NV-N2+i, F1+NV-N2+i+1, FVC+N2+i+2, FVC+N2+i+1))

        # Close the mesh
        if CloseObj is True:
            for i in range(N):
                faces.append((FVC+i, FVC+i+1, FVT+i+1, FVT+i))

    # case crown and rack
    else:
        if prop == "cro":
                (verts,  faces) = SolCrown(Vert)
        elif piniont is False:

            ####################################################################
            # First Vertices

            for i in range(Resolution+1):
                chi = (i/(1.0*Resolution))*pi/2.0
                z = Thicn-b*(1-cos(chi))
                sq = z*tan(Psi)

                for j in range(len(Vert)):
                    p = Vert[j]
                    m = Norm[j]
                    x = p[0]+b*(1-sin(chi))*m[0]
                    y = p[1]+b*(1-sin(chi))*m[1]
                    verts.append((x, y - sq, z))

            for i in range(LongRes):
                z = Thicn - b - (Thickness-2.0*b) * (i+1) / LongRes
                sq = z*tan(Psi)

                for j in range(len(Vert)):
                    p = Vert[j]
                    x = p[0]
                    y = p[1]
                    verts.append((x, y - sq, z))

            for i in range(Resolution):
                chi = (1.0-(i+1) / (1.0*Resolution)) * pi / 2.0
                z = - Thicn + b * (1-cos(chi))
                sq = z*tan(Psi)

                for j in range(len(Vert)):
                    p = Vert[j]
                    m = Norm[j]
                    x = p[0]+b*(1-sin(chi))*m[0]
                    y = p[1]+b*(1-sin(chi))*m[1]
                    verts.append((x, y - sq, z))

            ####################################################################
            # Then Faces
            NV = len(Vert)
            NL = 2*Resolution+LongRes
            FVC = len(verts)

            for i in range(NL):
                for j in range(NV-1):
                    faces.append((i*(NV)+j,  i*(NV)+j+1,  (i+1)*(NV)+j+1,  (i+1)*(NV)+j))

            ####################################################################
            # Add Width
            #
            # TOP

            p0 = verts[0]
            p1 = verts[NV-1]
            y0 = p0[1]
            y1 = p1[1]
            x0 = p0[0]
            N = 2*Resolution
            Wid = Width

            for i in range(N):
                y = (y1-y0)*i/(N-1)+y0
                verts.append(((x0-Wid),  y,  Thickness/2.0))
            verts.append((x0,  (y1+y0)/2.0,  Thickness/2.0))

            for i in range(N-1):
                faces.append((i,  i+1,  FVC))

            for i in range(Resolution-1):
                faces.append((i+N-1,  i+N,  FVC+i+1,  FVC+i))

            faces.append((int(N+N/2-2),  int(N+N/2-1),  FVC+N,  int(FVC+N/2-1)))
            faces.append((int(FVC+N/2-1),  FVC+N,  int(FVC+N/2)))

            for i in range(6*N):
                faces.append((int(N+N/2+i-1),  int(N+N/2+i),  FVC+N))

            faces.append((int(FVC+N/2),  FVC+N,  int(NV-N-N/2),  int(NV-N-N/2+1)))

            for i in range(Resolution-1):
                faces.append((int(NV-N-N/2+i+1), int(NV-N-N/2+i+2), int(FVC+N/2+i+1), int(FVC+N/2+i)))

            for i in range(N-1):
                faces.append((NV-i-2,  NV-i-1,  FVC+N-1))

            # BOTTOM
            F1 = NL*NV
            FVT = FVC
            FVC = len(verts)
            ya = y1
            y1 = - y0
            y0 = - ya
            Wid = Width
            for i in range(N):
                y = (y1-y0)*i/(N-1)+y0
                verts.append(((x0-Wid),   y,  -Thickness/2.0))
            verts.append((x0,  (y1+y0)/2.0,  -Thickness/2.0))

            for i in range(N-1):
                faces.append((F1+i,  F1+i+1,  FVC))

            for i in range(Resolution-1):
                faces.append((F1+i+N-1,  F1+i+N,  FVC+i+1,  FVC+i))

            faces.append((int(F1+N+N/2-2),  int(F1+N+N/2-1),  FVC+N, int(FVC+N/2-1)))
            faces.append((int(FVC+N/2-1),  FVC+N,  int(FVC+N/2)))

            for i in range(6*N):
                faces.append((int(F1+N+N/2+i-1),  int(F1+N+N/2+i), FVC+N))

            faces.append((int(FVC+N/2), FVC+N, int(F1+NV-N-N/2), int(F1+NV-N-N/2+1)))

            for i in range(Resolution-1):
                faces.append((int(F1+NV-N-N/2+i+1), int(F1+NV-N-N/2+i+2), int(FVC+N/2+i+1), int(FVC+N/2+i)))

            for i in range(N-1):
                faces.append((F1+NV-i-2, F1+NV-i-1, FVC+N-1))

            # Close the mesh
            if CloseObj is True:
                for i in range(N):
                    faces.append((FVC+i, FVC+i+1, FVT+i+1, FVT+i))
                # if(Rack==1):
                    # faces.append((F1+NV-2, FVC, FVT, NV-2))
                    # faces.append((F1+NV-N, FVC+N, FVT+N, NV-N))

    return verts, faces


####################################################################
# Deform base profile to cone
####################################################################

def ConifyTooth(verts, Alpha):

    vertsr = []

    R = PitchRadius
    Zo = Thickness/2.0
    h = R/tan(Alpha)
    for v in verts:
        x = v[0]
        y = v[1]
        z = v[2]+Zo

        r = sqrt(x*x+y*y)
        phi = atan(y/x)

        Rr = R + (r-R)*cos(Alpha)
        Rz = (r-R)*sin(Alpha)

        Mod = sqrt(Rr*Rr+(Rz-h)*(Rz-h))

        rc = Rr - z * Rr/Mod
        zc = Rz + z * (h-Rz)/Mod

        vertsr.append((rc*cos(phi), rc*sin(phi), zc))

    return vertsr


####################################################################
# CREATES THE BASE INVOLUTE PROFILE
####################################################################

def ToothOutline(pinion=False):
    global nIntoth, ax, Fillet

    ####################################################################
    # Compute ;ù@j_{[#
    #

    Ag = PressureAng * pi / 180.0
    if prop == "cro" and pinion is False:
        # page T23
        module = 2.0 * PitchRadius / TeethN
        # Cas la distance entraxe est fixé on corrige la longueur dent
        # Pteeth=int(TeethN/N)
        # x=PitchRadius*(1-Pteeth/TeethN)
        # y=(ax/module)-(TeethN-Pteeth)/2
        # alphaw=acos((TeethN-Pteeth)*cos(Ag)/(2*y+TeethN-Pteeth))
        # x=(TeethN-Pteeth)*(tan(alphaw)-alphaw-tan(Ag)+Ag)/(2*tan(Ag))
        x = 0
        Bottom = PitchRadius - (1-x) * module
        Ded = Bottom
        Add = Bottom + (2.25-x) * module
        bk = 1 / jeu
        Fillet = 0
    else:
        Bottom = PitchRadius - Dedendum - Fillet
        Ded = PitchRadius - Dedendum
        Add = PitchRadius + Addendum
        bk = jeu
    Base = PitchRadius * cos(Ag)
    DiametralPitch = TeethN/(2.0*PitchRadius)
    # ToothThickness = 1.5708/DiametralPitch
    CircularPitch = pi / DiametralPitch
    Theta0 = CircularPitch/(PitchRadius*2.0)

    # solve for involute thickness
    csphi = bk/cos(Ag)
    csphi = csphi*csphi
    csded = (Ded/PitchRadius)*(Ded/PitchRadius)/(cos(Ag)*cos(Ag))
    # csbevel = (Add/PitchRadius)*(Add/PitchRadius)/(cos(Ag)*cos(Ag))
    th = 0
    thd = 0
    aeqx = 0
    aeqDed = 0
    Precision = 400
    for i in range(Precision):
        tx = (i/Precision) * (pi/4)
        ft = ((cos(tx)+tx*sin(tx)))*((cos(tx)+tx*sin(tx))) + ((sin(tx)-tx*cos(tx))) * ((sin(tx)-tx*cos(tx)))
        aeqxl = aeqx
        aeqDedl = aeqDed
        aeqx = ft - csphi
        aeqDed = ft - csded
        if i > 1:
            if ((aeqx*aeqxl) <= 0):
                th = tx
                break
            if ((aeqDed*aeqDedl) <= 0):
                thd = tx

    txx = Theta0/2 - atan((-sin(th)+th*cos(th))/(cos(th)+th*sin(th)))
    if thd == 0:
        txd = atan(tan(txx)*Base/Ded)
    else:
        txd = txx + atan((-sin(thd)+thd*cos(thd))/(cos(thd)+thd*sin(thd)))
    # fillet angle
    thf = atan(Fillet/Bottom)
    txf = txd+thf

    '''if prop=="cro" and pinion==False :
        txx = - txx'''

    ####################################################################
    # Mesh
    #
    Nr = Resolution
    points = []

    # bottom of tooth
    points.append([Bottom*cos(Theta0), Bottom*sin(Theta0)])

    # Fillet
    xc = Ded*cos(txf)
    yc = Ded*sin(txf)
    Aw = pi/2.0 + txd - txf
    for i in range(Nr-1):
        thv = (Aw)*(i+1)/(Nr) + pi + txf
        points.append([xc + Fillet*cos(thv), yc + Fillet*sin(thv)])

    # Tooth Involute
    Theta5 = 0
    nIntoth = 0
    dteta = pi/(TeethN*50*Resolution)
    dref = (Add-Ded)/(3*Nr)
    az = xc + Fillet*cos(thv)
    bz = yc + Fillet*sin(thv)
    Rl = sqrt(az*az+bz*bz)
    thv = txx+Theta0/2
    for i in range(5000):
        sx = Base*(cos(thv)+thv*sin(thv))
        sy = -Base*(sin(thv)-thv*cos(thv))
        if sx > Ded:
            apointsx = sy * sin(txx)+sx*cos(txx)
            apointsy = sx * sin(txx)+sy*cos(txx)
            R = sqrt(apointsx*apointsx+apointsy*apointsy)
            if R > Add or apointsy <= 0:
                break
            if (R-Rl) > dref:
                Rl = R
                points.append([apointsx, apointsy])
                nIntoth = nIntoth + 1
            Theta5 = atan(apointsy/apointsx)
        thv = thv + dteta

    # Tooth Top
    for i in range(Nr-1):
        thv = Theta5 * (1-(i+1)/Nr)
        points.append([Add*cos(thv), Add*sin(thv)])

    # Mirrors
    Nr = len(points)
    for i in range(Nr):
        P = points[Nr-1-i]
        points.append([P[0], -P[1]])

    return points


####################################################################
# CREATES THE BASE RACK PROFILE
####################################################################

def RackOutline():

    ####################################################################
    # Basic Math computations: QUotes
    #
    X = {
        'Bottom': - Dedendum - Fillet,
        'Ded': - Dedendum,
        'Bevel': Addendum - Bevel,
        'Add': Addendum
    }

    ####################################################################
    # Basic Math computations: Angles
    #
    DiametralPitch = TeethN/(2*PitchRadius)
    ToothThickness = 1.5708/DiametralPitch
    CircularPitch = pi / DiametralPitch
    Pa = PressureAng*pi/180.0
    yA1 = ToothThickness/2.0
    yA2 = (-X['Ded']+Fillet*sin(Pa))*tan(Pa)
    yA3 = Fillet*cos(Pa)

    A = {
        'y0': CircularPitch/2.0,
        'y1': yA1 + yA2 + yA3,
        'y2': yA1 + yA2,
        'y3': yA1 - (X['Add']-Bevel) * tan(Pa),
        'y4': yA1 - (X['Add']-Bevel) * tan(Pa) - cos(Pa) / (1-sin(Pa)) * Bevel
    }

    ####################################################################
    # Profiling
    #
    N = Resolution
    points = []
    normals = []
    # Top half bottom of tooth
    for i in range(2*N):
        y = (A['y1'] - A['y0'])*i/(2*N-1) + A['y0']
        points.append([X['Bottom'],  y])
        normals.append([-1.0,  -0.0])

    # Bottom Fillet
    xc = X['Ded']
    yc = A['y1']
    Aw = pi/2.0 - Pa
    for i in range(N):
        th = (Aw)*(i+1)/(N) + pi
        points.append([xc + Fillet*cos(th),  yc + Fillet*sin(th)])
        normals.append([cos(th),  sin(th)])

    # Straight part
    Xded = X['Ded'] - Fillet*sin(Pa)
    for i in range(4*N):
        x = (X['Bevel']-Xded)*(i+1)/(4*N) + Xded
        points.append([x,  yA1-tan(Pa)*x])
        normals.append([-sin(Pa), -cos(Pa)])

    # Tooth Bevel
    rA = Bevel/(1-sin(Pa))
    xc = X['Add'] - rA
    yc = A['y4']
    for i in range(N):
        th = (-pi/2.0+Pa)*(i+1)/(N) + pi/2.0-Pa
        points.append([xc + rA*cos(th),  yc + rA*sin(th)])
        normals.append([-cos(th),  -sin(th)])

    # Tooth Top
    for i in range(N):
        y = -A['y4']*(i+1)/(N) + A['y4']
        points.append([X['Add'],  y])
        normals.append([-1.0,  0.0])

    # Mirrors this!
    N = len(points)
    for i in range(N-1):
        P = points[N-2-i]
        points.append([P[0], -P[1]])
        V = normals[N-2-i]
        normals.append([V[0], -V[1]])

    return points, normals


####################################################################
# CREATES THE BASE CROWN INVOLUTE
####################################################################

def CrownOutline():

    ####################################################################
    # Basic Math computations: Radii
    #
    R = {
        'Bottom': PitchRadius * cos(PressureAng*pi/180.0),
        'Base': PitchRadius * cos(PressureAng*pi/180.0) + Fillet,
        'Ded': PitchRadius + Dedendum
    }

    ####################################################################
    # Basic Math computations: Angles
    #
    DiametralPitch = TeethN/(2*PitchRadius)
    ToothThickness = 1.5708/DiametralPitch
    CircularPitch = pi / DiametralPitch

    U1 = sqrt((1-cos(PressureAng*pi/180.0))/cos(PressureAng*pi/180.0))
    U2 = sqrt(R['Ded']*R['Ded']/(R['Base']*R['Base'])-1)

    ThetaA1 = atan((sin(U1)-U1*cos(U1))/(cos(U1)+U1*sin(U1)))
    ThetaA2 = atan((sin(U2)-U2*cos(U2))/(cos(U2)+U2*sin(U2)))
    ThetaA3 = ThetaA1 + ToothThickness/(PitchRadius*2.0)

    A = {
        'Theta0': CircularPitch/(PitchRadius*2.0),
        'Theta1': (ThetaA3 + Fillet/R['Base']),
        'Theta2': ThetaA3,
        'Theta3': ThetaA3 - ThetaA2,
        'Theta4': ThetaA3 - ThetaA2 - Bevel/R['Ded']
    }

    M = A['Theta0']
    A['Theta0'] = 0
    A['Theta1'] = A['Theta1']-M
    A['Theta2'] = A['Theta2']-M
    A['Theta3'] = A['Theta3']-M
    A['Theta4'] = A['Theta4']-M

    ####################################################################
    # Profiling
    #
    N = Resolution
    apoints = []
    anormals = []

    # Top half top of tooth
    for i in range(2*N):
        th = (A['Theta1'] - A['Theta0'])*i/(2*N-1) + A['Theta0']
        apoints.append([R['Bottom']*cos(th), R['Bottom']*sin(th)])
        anormals.append([cos(th), sin(th)])

    # Bottom Bevel
    xc = R['Base']*cos(A['Theta1'])
    yc = R['Base']*sin(A['Theta1'])
    Aw = pi/2.0 + A['Theta2'] - A['Theta1']
    for i in range(N):
        th = (Aw)*(i+1)/(N) + pi + A['Theta1']
        apoints.append([xc + Fillet*cos(th), yc + Fillet*sin(th)])
        anormals.append([-cos(th), -sin(th)])

    # Tooth Involute
    for i in range(4*N):
        r = (R['Ded'] - R['Base'])*(i+1)/(4*N) + R['Base']
        u = sqrt(r*r/(R['Base']*R['Base'])-1)
        xp = R['Base']*(cos(u)+u*sin(u))
        yp = - R['Base']*(sin(u)-u*cos(u))
        apoints.append([xp*cos(A['Theta2'])-yp*sin(A['Theta2']), +xp*sin(A['Theta2'])+yp*cos(A['Theta2'])])
        anormals.append([sin(u), cos(u)])

    # Tooth Bevel
    auxth = -u
    auxth = auxth + ThetaA3 + pi/2.0
    # m = tan(auxth)
    P0 = apoints[len(apoints)-1]
    rA = Bevel/(1-cos(auxth-A['Theta4']))
    xc = P0[0] - rA*cos(auxth)
    yc = P0[1] - rA*sin(auxth)
    for i in range(N):
        th = (A['Theta4'] - auxth)*(i+1)/(N) + auxth
        apoints.append([xc + rA*cos(th), yc + rA*sin(th)])
        anormals.append([cos(th), sin(th)])

    # Tooth Top
    P0 = apoints[len(apoints)-1]
    A['Theta4'] = atan(P0[1]/P0[0])
    Ra = sqrt(P0[0]*P0[0]+P0[1]*P0[1])
    for i in range(N):
        th = (-M - A['Theta4'])*(i+1)/(N) + A['Theta4']
        apoints.append([Ra*cos(th), Ra*sin(th)])
        anormals.append([cos(th), sin(th)])

    points = []
    normals = []
    N = len(apoints)
    for i in range(N):
        points.append(apoints[N-1-i])
        normals.append(anormals[N-1-i])

    # Mirrors this!
    N = len(points)
    for i in range(N-1):
        P = points[N-2-i]
        points.append([P[0], -P[1]])
        V = normals[N-2-i]
        normals.append([V[0], -V[1]])

    return points, normals


######################################################
# Mesh the crown
######################################################
def SolCrown(Vert):

    vertcs = []
    facecs = []

    Thicn = Thickness / 2
    Psi = HelicalAng*pi / 180.0

    ###################################################################
    # make thickness
    for i in range(LongRes):
            z = Thicn - (Thickness)*(i)/(LongRes-1)
            sq = z*tan(Psi)
            Phi = sq/PitchRadius

            for j in range(len(Vert)):
                p = Vert[j]

                x = p[0]
                y = p[1]
                vertcs.append((x*cos(Phi)+y*sin(Phi), -x*sin(Phi)+y*cos(Phi), z))

    ####################################################################
    # make Faces
    NV = len(Vert)
    NL = LongRes-1
    N2 = Resolution
    FVC = len(vertcs)

    for i in range(NL):
        for j in range(NV-1):
            facecs.append((i*(NV)+j, i*(NV)+j+1, (i+1)*(NV)+j+1, (i+1)*(NV)+j))
    ####################################################################
    # make top
    p0 = vertcs[0]
    p1 = vertcs[NV-1]
    th0 = atan(p0[1]/p0[0])
    th1 = atan(p1[1]/p1[0])
    R = PitchRadius
    N = 2*Resolution
    Wid = -Width

    for i in range(N-2):
        th = (th1-th0)*i/(N-3)+th0
        vertcs.append(((R-Wid)*cos(th), (R-Wid)*sin(th), Thickness/2.0))

    for i in range(nIntoth+N2):
        facecs.append((i, i+1, FVC))

    for i in range(N-3):
        facecs.append((nIntoth+N2+i, nIntoth+N2+i+1, FVC+i+1, FVC+i))

    NN = NV-nIntoth-N2-1
    for i in range(nIntoth+N2):
        facecs.append((NN+i, NN+i+1, FVC+N-3))

    # make bottom
    F1 = NL*NV
    FVT = FVC
    FVC = len(vertcs)
    tha = th1
    th1 = -th0
    th0 = -tha

    for i in range(N-2):
        th = ((th1-th0)*i/(N-3)+th0)
        vertcs.append(((R-Wid)*cos(th), (R-Wid)*sin(th),  -Thickness/2.0))

    for i in range(nIntoth+N2):
        facecs.append((F1+i, F1+i+1,  FVC))

    for i in range(N-3):
        facecs.append((F1+nIntoth+N2+i,  F1+nIntoth+N2+i+1,  FVC+i+1,  FVC+i))

    NN = F1+NV-nIntoth-N2-1
    for i in range(nIntoth+N2):
        facecs.append((NN+i, NN+i+1,  FVC+N-3))

    # Close the mesh
    if CloseObj is True:
        for i in range(N-3):
            facecs.append((FVC+i, FVC+i+1, FVT+i+1, FVT+i))

    return vertcs, facecs
