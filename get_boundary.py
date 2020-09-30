import os
import cv2
import sys
import numpy as np

img = cv2.imread(sys.argv[1])
# set pix higher than 127 to 255
ret, thresh = cv2.threshold(cv2.cvtColor(img.copy(), cv2.COLOR_BGR2GRAY), 127, 255, 0)
# find contours
cnts, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

c = np.array(cnts[0])
c = c.reshape(c.shape[0], c.shape[2])
c = c[:, :: -1]

np.savetxt('boundary_point.txt', c, fmt='%d', delimiter = ' ')