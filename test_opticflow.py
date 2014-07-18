# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 12:04:07 2014

@author: etytel01
"""

import cv2
import numpy as np
cap = cv2.VideoCapture('rawdata/Bg5-26_5hz_V_008.avi')

ret, frame1 = cap.read()
prvs = cv2.cvtColor(frame1,cv2.COLOR_BGR2GRAY)
prvs = np.float64(prvs)/255
hsv = np.zeros(prvs.shape + (3,), np.uint8)
hsv[...,1] = 255

while(1):
    ret, frame2 = cap.read()
    nxt = cv2.cvtColor(frame2,cv2.COLOR_BGR2GRAY)
    nxt = np.float64(nxt)/255

    flow = cv2.calcOpticalFlowFarneback(prvs,nxt, None, 0.5, 3, 15, 3, 5, 1.2, 0)

    mag, ang = cv2.cartToPolar(flow[...,0], flow[...,1])
    hsv[...,0] = ang*180/np.pi/2
    hsv[...,2] = cv2.normalize(mag,None,0,255,cv2.NORM_MINMAX)
    rgb = cv2.cvtColor(hsv,cv2.COLOR_HSV2BGR)

    cv2.imshow('frame2',rgb)
    k = cv2.waitKey(30) & 0xff
    if k == 27:
        break
    elif k == ord('s'):
        cv2.imwrite('opticalfb.png',frame2)
        cv2.imwrite('opticalhsv.png',rgb)
    prvs = nxt

cap.release()
cv2.destroyAllWindows()