#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import cv2
from PIL import Image, ImageDraw, ImageFont
from natsort import natsorted
import numpy as np

# Function to add timestamp text to images
def add_timestamp(image, frame_index, font_size=30, position=(10, 10), font_path=None):
    """
    Adds a timestamp text (5 ns per frame) to the image.
    
    Parameters:
        image (PIL.Image): Input image.
        frame_index (int): Frame number.
        font_size (int): Font size of the text.
        position (tuple): (x, y) position of text.
        font_path (str): Path to custom font file.
    
    Returns:
        PIL.Image: Image with timestamp overlay.
    """
    draw = ImageDraw.Draw(image)
    if font_path is None:
        font = ImageFont.load_default()
    else:
        font = ImageFont.truetype(font_path, font_size)

    timestamp_text = f"Frame {frame_index} ({frame_index * 5} ns)"
    draw.text(position, timestamp_text, font=font, fill="white")

    return image

# Function to resize while keeping transparency
def process_image(image, new_size=(512, 512)):
    """
    Resizes the image and ensures transparency is maintained.

    Parameters:
        image (PIL.Image): Input image.
        new_size (tuple): New size (width, height).

    Returns:
        PIL.Image: Processed image in RGBA mode.
    """
    image = image.convert("RGBA")  # Ensure transparency
    image = image.resize(new_size, Image.LANCZOS)
    return image

# Function to create transparent GIF
def create_gif_with_annotations(image_folder, output_gif, duration=0.30, resize_dims=(512, 512), font_path=None):
    im_list = []
    
    for idx, file in enumerate(natsorted(os.listdir(image_folder))):
        if file.endswith(".png"):
            img = Image.open(os.path.join(image_folder, file))

            # Resize and keep transparency
            img = process_image(img, resize_dims)

            # Add timestamp
            img = add_timestamp(img, idx + 1, font_size=30, position=(10, 10), font_path=font_path)

            im_list.append(img)

    if im_list:
        im_list[0].save(
            output_gif,
            save_all=True,
            append_images=im_list[1:],
            duration=int(duration * 1000),  # Convert seconds to ms
            loop=0,
            disposal=2,  # Ensures no frame blending
            transparency=0  # Keep transparency
        )
        print(f"Transparent GIF saved at: {output_gif}")
    else:
        print("No valid images found in the directory!")

# Function to create a transparent video
def create_video(image_folder, output_video, fps=10, resize_dims=(512, 512), font_path=None):
    frame_array = []
    
    for idx, file in enumerate(natsorted(os.listdir(image_folder))):
        if file.endswith(".png"):
            img = Image.open(os.path.join(image_folder, file))

            # Resize and keep transparency
            img = process_image(img, resize_dims)

            # Add timestamp
            img = add_timestamp(img, idx + 1, font_size=30, position=(10, 10), font_path=font_path)

            # Convert to OpenCV format (BGRA for transparency)
            image = np.array(img)
            frame = cv2.cvtColor(image, cv2.COLOR_RGBA2BGRA)
            frame_array.append(frame)

    if frame_array:
        height, width, _ = frame_array[0].shape
        fourcc = cv2.VideoWriter_fourcc(*"vp09")  # WebM format (supports transparency)
        out = cv2.VideoWriter(output_video, fourcc, fps, (width, height), isColor=True)

        for frame in frame_array:
            out.write(frame)

        out.release()
        print(f"Transparent Video saved at: {output_video}")
    else:
        print("No valid images found for video creation!")

# Define input path and output files
input_path = "/mnt/mydrive/alpfaFold3/dynamics/Custom_MSA/dclk3/dclk3_dclk3_fmsaatemp/snaps/"
output_gif_path = os.path.join(input_path, "0_output.gif")
output_video_path = os.path.join(input_path, "0_output.webm")  # WebM for transparency

# Run GIF and Video generation
create_gif_with_annotations(input_path, output_gif_path, duration=0.30, resize_dims=(512, 512))
create_video(input_path, output_video_path, fps=10, resize_dims=(512, 512))
