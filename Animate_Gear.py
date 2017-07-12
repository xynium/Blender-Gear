# ____  ________.___._______  .___ ____ ___  _____
# \   \/  /\__  |   |\      \ |   |    |   \/     \
#  \     /  /   |   |/   |   \|   |    |   /  \ /  \
#  /     \  \____   /    |    \   |    |  /    Y    \
# /___/\  \ / ______\____|__  /___|______/\____|__  /
#       \_/ \/              \/                    \/
#
# This script should be executed only one time
# before a second execution erase rigid body constraint on gear and pinion
# For conic train rotate (r y 90) cylinder pinion empty pinion and empty motor
# With crown train set in gear rigid body margin to 0.01

import bpy

bl_info = {
    "name": "Animate_Gear",
    "description": "Animate a gear train.",
    "author": "XYNIUM JP Lathuile",
    "version": (1, 0),
    "blender": (2, 65, 0),
    "location": "View3D > Add > Mesh",
    "warning": "",
    "wiki_url": "http://wiki.blender.org/index.php/Extensions:2.6/Py/Scripts/Add-Mesh-Gear.py",
    "tracker_url": "https://developer.blender.org/maniphest/task/edit/form/2/",
    "support": "TESTING",
    "category": "Add Mesh"
}


#####################################
#   Main
#   Create additonnal object to animate
#            the stuff
#####################################

class AnimateGear(bpy.types.Operator):
    bl_idname = "mesh.primitive_animate_gear"
    bl_label = "Animate a Gear"
    bl_options = {'REGISTER', 'UNDO', 'PRESET'}
    bl_description = "Animate a gear/pinion"

    def execute(self, context):
        scene = bpy.context.scene
        k = 1
        EltLst = list(bpy.data.objects)
        for Elt in EltLst:
            if Elt.name.find("Gear") == 0:
                Gear = Elt.name
                k = k + 1
            if Elt.name.find("Pinion") == 0:
                Pinion = Elt.name
                k = k + 1
        if k == 3:
            # Pinion
            bpy.ops.object.mode_set(mode='OBJECT')
            bpy.ops.object.select_all(action='DESELECT')
            bpy.data.objects[Pinion].select = True
            scene.objects.active = bpy.data.objects[Pinion]
            bpy.ops.rigidbody.object_add(type='ACTIVE')
            bpy.ops.rigidbody.constraint_add()
            bpy.data.objects[Pinion].rigid_body.collision_shape = 'MESH'
            bpy.data.objects[Pinion].rigid_body_constraint.disable_collisions = False
            bpy.data.objects[Pinion].rigid_body.collision_margin = 0

            bpy.ops.object.select_all(action='DESELECT')
            bpy.ops.mesh.primitive_cylinder_add()
            CylPin = bpy.context.selected_objects
            CylPin[0].location = bpy.data.objects[Pinion].location
            bpy.ops.rigidbody.object_add(type='PASSIVE')

            bpy.ops.object.select_all(action='DESELECT')
            bpy.ops.object.empty_add(type='PLAIN_AXES')
            EptyPin = bpy.context.selected_objects
            EptyPin[0].location = bpy.data.objects[Pinion].location
            bpy.ops.rigidbody.constraint_add()
            EptyPin[0].rigid_body_constraint.object1 = bpy.data.objects[Pinion]
            EptyPin[0].rigid_body_constraint.object2 = CylPin[0]
            EptyPin[0].rigid_body_constraint.type = "HINGE"

            # Gear
            bpy.ops.object.mode_set(mode='OBJECT')
            bpy.ops.object.mode_set(mode='OBJECT')
            bpy.ops.object.select_all(action='DESELECT')
            bpy.data.objects[Gear].select = True
            scene.objects.active = bpy.data.objects[Gear]
            bpy.ops.rigidbody.object_add(type='ACTIVE')
            bpy.ops.rigidbody.constraint_add()
            bpy.data.objects[Gear].rigid_body.collision_shape = 'MESH'
            bpy.data.objects[Gear].rigid_body_constraint.disable_collisions = False
            bpy.data.objects[Gear].rigid_body.collision_margin = 0

            bpy.ops.object.select_all(action='DESELECT')
            bpy.ops.mesh.primitive_cylinder_add()
            CylGear = bpy.context.selected_objects
            CylGear[0].location = bpy.data.objects[Gear].location
            bpy.ops.rigidbody.object_add(type='PASSIVE')

            bpy.ops.object.select_all(action='DESELECT')
            bpy.ops.object.empty_add(type='PLAIN_AXES')
            EptyGear = bpy.context.selected_objects
            EptyGear[0].location = bpy.data.objects[Gear].location
            bpy.ops.rigidbody.constraint_add()
            EptyGear[0].rigid_body_constraint.object1 = bpy.data.objects[Gear]
            EptyGear[0].rigid_body_constraint.object2 = CylGear[0]
            EptyGear[0].rigid_body_constraint.type = "HINGE"

            # Motor
            bpy.ops.object.select_all(action='DESELECT')
            bpy.ops.object.empty_add(type='PLAIN_AXES')
            EptyMotor = bpy.context.selected_objects
            EptyMotor[0].location = bpy.data.objects[Pinion].location
            bpy.ops.transform.rotate(value=90.0, axis=(0.0, 1.0, 0.0))
            bpy.ops.rigidbody.constraint_add()
            EptyMotor[0].rigid_body_constraint.type = 'MOTOR'
            EptyMotor[0].rigid_body_constraint.object1 = bpy.data.objects[Pinion]
            EptyMotor[0].rigid_body_constraint.object2 = CylPin[0]
            EptyMotor[0].rigid_body_constraint.use_motor_ang = True
            bpy.ops.screen.animation_play()
        return {'FINISHED'}


class INFO_MT_mesh_gear_add2(bpy.types.Menu):
    # Define the "Gears" menu
    bl_idname = "INFO_MT_mesh_gear_add2"
    bl_label = "Gear"

    def draw(self, context):
        layout = self.layout
        layout.operator_context = 'INVOKE_REGION_WIN'
        layout.operator("mesh.primitive_animate_gear",  text="Animate Gear")


# Define "Extras" menu
def menu_func(self, context):
    layout = self.layout
    layout.operator_context = 'INVOKE_REGION_WIN'
    # layout.menu("INFO_MT_mesh_gear_add2", text="Gear", icon="PLUGIN")
    layout.operator("Gear.mesh.primitive_animate_gear", text="Animate Gear")
    layout.operator("Gear.mesh.primitive_gear",  text="Mesh Gear")


def register():
    bpy.utils.register_module(__name__)
    bpy.types.INFO_MT_mesh_add.append(menu_func)


def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.INFO_MT_mesh_add.remove(menu_func)


if __name__ == "__main__":
    register()
