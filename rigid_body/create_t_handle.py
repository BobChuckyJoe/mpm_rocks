import bpy
import math

bpy.ops.mesh.primitive_cylinder_add(enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
bpy.ops.mesh.primitive_cylinder_add(enter_editmode=False, align='WORLD', location=(0, 0, 1), scale=(1, 1, 2), rotation=(0, math.pi/2, 0))
bpy.data.objects["Cylinder"].select_set(True)
bpy.data.objects["Cylinder.001"].select_set(True)
bpy.ops.object.join()
bpy.ops.object.origin_set(type="GEOMETRY_ORIGIN", center="MEDIAN")