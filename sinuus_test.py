#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 14:37:28 2022

@author: mayk
"""

import matplotlib.pyplot as plt
import numpy as np




# x = np.arange(0,4*np.pi,0.1)   # start,stop,step
# y = np.sin(x)



# plt.plot(x,y)
# plt.show()



# x = np.arange(0,4*np.pi,0.1)   # start,stop,step
# y = np.sin(x)
# z = np.cos(x)



# plt.plot(x,y,x,z)
# plt.show()



x = np.arange(0,4*np.pi-1,0.1)   # start,stop,step
y = np.sin(x)
z = np.cos(x)

#plt.figure(dpi=1200)


plt.plot(x,y,x,z)
plt.xlabel('x values from 0 to 4pi')  # string must be enclosed with quotes '  '
plt.ylabel('sin(x) and cos(x)')
plt.title('Plot of sin and cos from 0 to 4pi')
plt.legend(['sin(x)', 'cos(x)'])      # legend entries as seperate strings in a list
plt.show()

plt.savefig('filename.pdf' )

