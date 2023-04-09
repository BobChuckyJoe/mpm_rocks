#!/usr/bin/env python3
import math
import os
from pathlib import Path
import sys

BLENDER_PATH = "/scratch/cluster/bcj/blender-3.4.0-linux-x64/blender"
FFMPEG_PATH = "/scratch/cluster/bcj/ffmpeg"

def get_frame_names(start_frame, end_frame):
    return "|".join([f"frame_{i:06d}.png" for i in range(start_frame, end_frame + 1)])

def get_files_path(blend_file_path):
    path = Path(blend_file_path)
    if path.is_absolute():
        return path.parent.absolute()
    else:
        return os.getcwd()
# For [a,b], there will be b - a + 1 frames total, so each worker should do
# ceil((b - a + 1) / tot_workers) frames
def get_start_and_end(worker_id, tot_workers, start_frame, end_frame):
    work_size = math.ceil((end_frame - start_frame + 1) / tot_workers)
    return (start_frame + worker_id * work_size, min(start_frame + (worker_id + 1) * work_size - 1, end_frame))

def run_blender(blend_file_path, output_relative_path, start_frame, end_frame, worker_id):
    # TODO Work on this more:
    # "../blender-3.4.0-linux-x64/blender -b light_rock.blend -o //light_rock_ -F PNG -x 1 -f 1..2"
    
    # This works
    # cmd = f"{BLENDER_PATH} -b {blend_file_path} -o //{output_relative_path}/frame_###### -F JPEG -x 1 -f {start_frame}..{end_frame}"

    cmd = f"{BLENDER_PATH} -b {blend_file_path} -o //{output_relative_path}/animation_ -F AVIJPEG -x 1 -s {start_frame} -e {end_frame} -a"

    print(cmd, file=sys.stderr)
    os.system(cmd)
if __name__ == "__main__":
    # The 0th arg is the file name
    worker_id = int(sys.argv[1])
    tot_workers = int(sys.argv[2])
    start_frame = int(sys.argv[3])
    end_frame = int(sys.argv[4])
    blender_file_path = sys.argv[5]
    output_relative_path = sys.argv[6] # Will output png images 

    worker_start, worker_end = get_start_and_end(worker_id, tot_workers, start_frame, end_frame)
    run_blender(blender_file_path, output_relative_path, worker_start, worker_end, worker_id)

    # Merge the individual frames outputted ACROSS ALL RUNS into a single video using another
    # os.system(f"cd {get_files_path(blender_file_path)} && {FFMPEG_PATH} -r 60 -f image2 -s 1920x1080 -i {get_frame_names(worker_start, worker_end)} -vcodec libx264 -crf 25  -pix_fmt yuv420p out_{worker_id}.mp4")