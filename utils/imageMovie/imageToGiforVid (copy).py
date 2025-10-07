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

# Function to resize while keeping transparency and aspect ratio
def process_image(image, new_size=(512, 512)):
    """
    Resizes the image while maintaining aspect ratio and ensures transparency is maintained.

    Parameters:
        image (PIL.Image): Input image.
        new_size (tuple): Maximum size (width, height).

    Returns:
        PIL.Image: Processed image in RGBA mode with preserved aspect ratio.
    """
    image = image.convert("RGBA")  # Ensure transparency
    
    # Calculate aspect ratio preserving resize
    original_width, original_height = image.size
    target_width, target_height = new_size
    
    # Calculate scaling factor to fit within target size while preserving aspect ratio
    scale_factor = min(target_width / original_width, target_height / original_height)
    
    # Calculate new dimensions
    new_width = int(original_width * scale_factor)
    new_height = int(original_height * scale_factor)
    
    # Resize with preserved aspect ratio
    image = image.resize((new_width, new_height), Image.LANCZOS)
    
    # Create a new image with target dimensions and transparent background
    final_image = Image.new("RGBA", new_size, (0, 0, 0, 0))  # Transparent background
    
    # Center the resized image in the new canvas
    x_offset = (target_width - new_width) // 2
    y_offset = (target_height - new_height) // 2
    final_image.paste(image, (x_offset, y_offset), image)
    
    return final_image

# Function to create transparent GIF without annotations and with original quality
def create_gif_with_annotations(image_folder, output_gif, duration=0.30, resize_dims=None, font_path=None):
    im_list = []
    
    for idx, file in enumerate(natsorted(os.listdir(image_folder))):
        if file.endswith(".tga") or file.endswith(".png"):
            img = Image.open(os.path.join(image_folder, file))

            # Convert to RGBA to ensure transparency support
            img = img.convert("RGBA")
            
            # Only resize if resize_dims is specified, otherwise keep original size
            if resize_dims is not None:
                img = process_image(img, resize_dims)

            # Skip adding timestamp annotation
            # img = add_timestamp(img, idx + 1, font_size=30, position=(10, 10), font_path=font_path)

            im_list.append(img)

    if im_list:
        im_list[0].save(
            output_gif,
            save_all=True,
            append_images=im_list[1:],
            duration=int(duration * 1000),  # Convert seconds to ms
            loop=0,
            disposal=2,  # Ensures no frame blending
            transparency=0,  # Keep transparency
            optimize=False  # Disable optimization to preserve quality
        )
        print(f"High-quality transparent GIF saved at: {output_gif}")
    else:
        print("No valid images found in the directory!")

# Function to create a transparent video without annotations and with original quality
def create_video(image_folder, output_video, fps=10, resize_dims=None, font_path=None):
    frame_array = []
    
    for idx, file in enumerate(natsorted(os.listdir(image_folder))):
        if file.endswith(".tga") or file.endswith(".png"):
            img = Image.open(os.path.join(image_folder, file))

            # Convert to RGBA to ensure transparency support
            img = img.convert("RGBA")
            
            # Only resize if resize_dims is specified, otherwise keep original size
            if resize_dims is not None:
                img = process_image(img, resize_dims)

            # Skip adding timestamp annotation
            # img = add_timestamp(img, idx + 1, font_size=30, position=(10, 10), font_path=font_path)

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
        print(f"High-quality transparent Video saved at: {output_video}")
    else:
        print("No valid images found for video creation!")

# Define input path and output files
input_path = "/mnt/mydrive/others/dr_kannan/brady/test/"
output_gif_path = os.path.join(input_path, "0_output.gif")
output_video_path = os.path.join(input_path, "0_output.webm")  # WebM for transparency

# Run GIF and Video generation with original quality (no resizing, no annotations)
create_gif_with_annotations(input_path, output_gif_path, duration=0.20, resize_dims=None)
##create_video(input_path, output_video_path, fps=10, resize_dims=None)
