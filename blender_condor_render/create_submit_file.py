# Creates the script and the condor submit file

# Configuration values
NUM_QUEUE = 10 # Number of jobs to queue
START_FRAME = 900 
END_FRAME = 1500
BLEND_FILE_PATH = "rock_on_dirt.blend"
OUTPUT_RELATIVE_PATH = "rock_on_dirt"

if __name__ == "__main__":
    output = f"""
Requirements = InMastodon
Rank=Memory

+Group="UNDER"
+Project="OTHER"
+ProjectDescription="Renders blender animation for thesis"
+GPUjob=true
request_gpus=1

executable=render_worker.py
arguments=$(Process) {NUM_QUEUE} {START_FRAME} {END_FRAME} {BLEND_FILE_PATH} {OUTPUT_RELATIVE_PATH}

output=last_run/render_$(Process).output
error=last_run/render_$(Process).error
log =last_run/render.log

should_transfer_file=false
stream_output=true

Queue {NUM_QUEUE}
"""

    with open("render.submit", "w") as f:
        f.write(output)
    