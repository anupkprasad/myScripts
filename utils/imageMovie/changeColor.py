"""
Created on Sun Aug  8 12:03:54 2021

@author: anup
"""


from PIL import Image

img = Image.open("fig_1.png")
img = img.convert("RGB")

datas = img.getdata()

new_image_data = []
for item in datas:
    # change all white (also shades of whites) pixels to yellow
    if (item[0] in list(range(50,256))) and (item[1] in list(range(50,256))) and (item[2] in list(range(0,50))):                 ###item[0] in list(range(100, 255)):
        new_image_data.append((230, 0, 255))
    else:
        new_image_data.append(item)
        
# update image data
img.putdata(new_image_data)

# save new image
img.save("test_image_altered_background.png", quality = 100)
