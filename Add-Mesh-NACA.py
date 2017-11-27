#
# ____  ________.___._______  .___ ____ ___  _____
# \   \/  /\__  |   |\      \ |   |    |   \/     \
#  \     /  /   |   |/   |   \|   |    |   /  \ /  \
#  /     \  \____   /    |    \   |    |  /    Y    \
# /___/\  \ / ______\____|__  /___|______/\____|__  /
#       \_/ \/              \/                    \/
#

# notes :
# NACA cord from : naca.py
# Copyright (C) 2011 by Dirk Gorissen <dgorissen@gmail.com>



import bpy
# from mathutils import *
from math import sin, cos, tan, pi, atan, sqrt
# from bpy.props import *

bl_info = {
    "name": "NACA",
    "description": "Add a NACA foil mesh.",
    "author": "XYNIUM JP Lathuile",
    "version": (1, 0),
    "blender": (2, 79, 0),
    "location": "View3D > Add > Mesh",
    "warning": "",
    "wiki_url": "http://wiki.blender.org/index.php/Extensions:2.6/Py/Scripts/"
    "Add-Mesh-NACA.py",
    "tracker_url": "https://developer.blender.org/maniphest/task/edit/form/2/",
    "support": "TESTING",
    "category": "Add Mesh"
}

#####################################
#   Main Add NACA
#####################################

class ObjectNACA(bpy.types.Operator):
    bl_idname = "mesh.primitive_naca"
    bl_label = "Add a NACA Foil"
    bl_options = {'REGISTER',  'UNDO',  'PRESET'}

    # default blender is mm
    sNACAtype = bpy.props.StringProperty(name="NACA", default="0018", maxlen=4)
    snbrPts = bpy.props.IntProperty(name="Nbr Points", default=40, min=5, soft_max=400)
    sLength = bpy.props.IntProperty(name="Length foil [mm]", default=15, min=4)
    sChord = bpy.props.IntProperty(name="Chord [mm]", default=45, min=4)
      
    def draw(self,  context):
        layout = self.layout

        box = layout.box()
        box.prop(self, 'sNACAtype')
        box.prop(self, 'snbrPts')
        box.prop(self, 'sLength')
        box.prop(self, 'sChord')
               
    def execute(self,  context):
     
        degrad = 3.1416/180.0
        scene = bpy.context.scene

        #Create mesh
        NACAtype = self.sNACAtype
        nbrPts = self.snbrPts
        Length = self.sLength
        Larg = self.sChord
       
        v,e, f = NACAMesh(NACAtype, nbrPts, Length, Larg)
        #bpy.ops.object.mode_set(mode='OBJECT')
        '''bpy.data.screens['Default'].scene.cursor_location[0]=0
        bpy.data.screens['Default'].scene.cursor_location[1]=0
        bpy.data.screens['Default'].scene.cursor_location[2]=0
        bpy.ops.object.origin_set(type='ORIGIN_CURSOR')'''
        obj = create_mesh_object(context, v, e,  f,  "NACA")
                   
        return {'FINISHED'}


class INFO_MT_mesh_NACA_add(bpy.types.Menu):
    # Define the "NACA" menu
    bl_idname = "INFO_MT_mesh_NACA_add"
    bl_label = "NACA Foil"

    def draw(self,   context):
        layout = self.layout
        layout.operator_context = 'INVOKE_REGION_WIN'
        layout.operator("mesh.primitive_naca",    text="Mesh NACA Foil")

# Define "Extras" menu
def menu_func(self,   context):
        layout = self.layout
        layout.operator_context = 'INVOKE_REGION_WIN'
        layout.menu("INFO_MT_mesh_NACA_add",   text="NACA Foil",   icon="PLUGIN")


def register():
        bpy.utils.register_module(__name__)
        bpy.types.INFO_MT_mesh_add.append(menu_func)


def unregister():
        bpy.utils.unregister_module(__name__)
        bpy.types.INFO_MT_mesh_add.remove(menu_func)

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
# Do The Foil
# define profile
# make face
####################################################################

def NACAMesh(NACAtype, nbrPts, Length, Larg):
    v = []
    e = []
    f = []
    
    finite_TE = True
    half_cosine_spacing = True
    
    # Returns 2*n+1 points in [0 1] for the given 4 digit NACA number string
    m = float(NACAtype[0])/100.0
    p = float(NACAtype[1])/10.0
    t = float(NACAtype[2:])/100.0

    a0 = +0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = +0.2843

    if finite_TE is True:
        a4 = -0.1015 # For finite thick TE
    else:
        a4 = -0.1036 # For zero thick TE

    if half_cosine_spacing is True:
        beta = linspace(0.0,pi,nbrPts+1)
        x = [(0.5*(1.0-cos(xx))) for xx in beta]  # Half cosine based spacing
    else:
        x = linspace(0.0,1.0,nbrPts+1)

    yt = [5*t*(a0*sqrt(xx)+a1*xx+a2*pow(xx,2)+a3*pow(xx,3)+a4*pow(xx,4)) for xx in x]

    xc1 = [xx for xx in x if xx <= p]
    xc2 = [xx for xx in x if xx > p]

    if p == 0:
        xu = x
        yu = yt

        xl = x
        yl = [-xx for xx in yt]

        xc = xc1 + xc2
        zc = [0]*len(xc)
    else:
        yc1 = [m/pow(p,2)*xx*(2*p-xx) for xx in xc1]
        yc2 = [m/pow(1-p,2)*(1-2*p+xx)*(1-xx) for xx in xc2]
        zc = yc1 + yc2

        dyc1_dx = [m/pow(p,2)*(2*p-2*xx) for xx in xc1]
        dyc2_dx = [m/pow(1-p,2)*(2*p-2*xx) for xx in xc2]
        dyc_dx = dyc1_dx + dyc2_dx

        theta = [atan(xx) for xx in dyc_dx]

        xu = [xx - yy * sin(zz) for xx,yy,zz in zip(x,yt,theta)]
        yu = [xx + yy * cos(zz) for xx,yy,zz in zip(zc,yt,theta)]

        xl = [xx + yy * sin(zz) for xx,yy,zz in zip(x,yt,theta)]
        yl = [xx - yy * cos(zz) for xx,yy,zz in zip(zc,yt,theta)]

    X = xu[::-1] + xl[1:]
    Z = yu[::-1] + yl[1:]
       
    #mesh 
    isFirst = True
    for iN in range(0,2*nbrPts):  #one end of the foil
        v.append((X[iN]*Larg,Z[iN]*Larg, 0))  
    for iN in range(0,2*nbrPts):  #the other end
        v.append((X[iN]*Larg,Z[iN]*Larg, Length))
    for iN in range(0,2*nbrPts-1):  # face longeur 
        f.append((iN ,iN+1, 2*nbrPts+iN+1 ,2*nbrPts+iN)) 
    for iN in range(0,nbrPts-1):  # face  profil
        f.append((iN , 2*nbrPts-1-iN ,iN+1) ) 
        f.append((2*nbrPts -1- iN , 2*nbrPts-2-iN,iN+1 ) ) 
        f.append((2*nbrPts+iN ,2*nbrPts+iN+1, 4*nbrPts-1-iN ) ) 
        f.append((4*nbrPts -1- iN ,2*nbrPts+iN+1, 4*nbrPts-2-iN ) ) 
        
    f.append((2*nbrPts-1,0,2*nbrPts,4*nbrPts-1)) #rebouclage face bord de fuite
  
    return v,e, f 
    
def linspace(start,stop,np):
    #Emulate Matlab linspace
        
    return [start+(stop-start)*i/(np-1) for i in range(np)]    


