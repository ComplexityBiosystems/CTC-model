import bpy  

# select object manually
def select_current():
    obj = bpy.context.object
    scene  = bpy.context.scene
    return obj, scene

def set_trajectory(scene, obj, frame_start, frame_increment, still_frames, positions):

    pos = positions[0]
    for frame_idx in range(frame_start, frame_start + still_frames):
        scene.frame_set(frame_idx)
        obj.location = pos
        obj.keyframe_insert(data_path="location", index=-1)

    frame_idx = frame_start + still_frames
    for pos in positions:
        scene.frame_set(frame_idx)
        obj.location = pos
        obj.keyframe_insert(data_path="location", index=-1)
        frame_idx += frame_increment
