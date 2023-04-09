#!/usr/local/bin/bash

# For images
# NEED -start_number flag for image sequence
# export START_NUMBER=900
# /scratch/cluster/bcj/ffmpeg -r 24 -f image2 -start_number $START_NUMBER -s 1920x1080 -i frame_%06d.jpg -vcodec libx264 -crf 25  -pix_fmt yuv420p out.mp4

# For video
for f in animation_*; do echo "file '$f'" >> videos.txt; done
/scratch/cluster/bcj/ffmpeg -f concat -i videos.txt -c copy out.mp4

# Compress the video, so it takes less space
/scratch/cluster/bcj/ffmpeg -i out.mp4 -vcodec libx264 -crf 16 compressed.mp4
