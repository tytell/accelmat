# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 13:45:42 2014

@author: etytel01
"""

import numpy as np
import cv2
from matplotlib import pyplot as plt

cap = cv2.VideoCapture('rawdata/Bg5-26_5hz_V_008.avi')

ret, img = cap.read()
cap.release()

gray = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)

sift = cv2.SIFT()
kp = sift.detect(gray,None)

img=cv2.drawKeypoints(gray,kp,flags=cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)

cv2.imwrite('sift_keypoints.jpg',img)