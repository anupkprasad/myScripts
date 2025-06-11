#!/bin/bash
ffmpeg -framerate 10 -pattern_type glob -i '*.png' -vf "palettegen" palette.png && \
ffmpeg -framerate 10 -pattern_type glob -i '*.png' -i palette.png -lavfi "paletteuse" output.gif