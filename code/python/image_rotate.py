from PIL import Image
import numpy as np
import os

def im_rotate(img_path):
    img = Image.open(img_path)
    # rotate image from -20~20
    angles = [-20, -15, -10, -5, 5, 10, 15, 20]
    width, height = img.size
    n = 20
    count = 1
    for i in angles:
        temp = 255*np.ones((height + 2 * n, width + 2 * n))
        temp[n:(height + n), n:(width + n)] = img
        temp = Image.fromarray(temp)
        temp_rotate = temp.rotate(i, Image.BICUBIC)
        temp_rotate = np.asarray(temp_rotate)
        im = temp_rotate[n:(height + n), n:(width + n)]
        im = Image.fromarray(im)
        if im.mode != 'RGB':
            im = im.convert('RGB')
        f, e = os.path.split(img_path)
        save_path = '{0}/{1}_{2}.png'.format(f, e.split('.')[0], count)
        print save_path
        im.save(save_path)
        count += 1
    return

img_path = '../../data/20150728212540_001_fin_002_017.jpg'
im_rotate(img_path)
