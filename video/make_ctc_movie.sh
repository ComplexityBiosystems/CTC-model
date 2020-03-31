ffmpeg \
    -start_number 1 \
    -i frames/frame%04d.png \
    -c:v libx264 \
    -preset veryslow  \
    -crf 10 \
    -pix_fmt yuv420p \
    -vf scale=3840:2160 \
    out_very_large.mp4

ffmpeg \
    -start_number 1 \
    -i frames/frame%04d.png \
    -c:v libx264 \
    -preset veryslow  \
    -crf 10 \
    -pix_fmt yuv420p \
    -vf scale=1920:1080 \
    out_large.mp4

ffmpeg \
    -start_number 1 \
    -i frames/frame%04d.png \
    -c:v libx264 \
    -preset veryslow \
    -crf 0 \
    -vf scale=960:540 \
    out_medium.mkv

ffmpeg \
    -start_number 1 \
    -i frames/frame%04d.png \
    -c:v libx264 \
    -preset veryslow \
    -crf 0 \
    -vf scale=480:270 \
    out_small.mkv
