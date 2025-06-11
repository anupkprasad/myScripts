"""
Created on Sun Aug  8 19:29:46 2021

@author: anup
"""

from PIL import Image
  
def convertImage():
    img = Image.open("./test.png")
    img = img.convert("RGBA")
  
    datas = img.getdata()
  
    newData = []
  
    for items in datas:
        if items[0] == 255 and items[1] == 255 and items[2] == 255:
            newData.append((255, 255, 255, 0))
        else:
            newData.append(items)
  
    img.putdata(newData)
    img.save("./New.png", "PNG")
    print("Successful")
  
convertImage()