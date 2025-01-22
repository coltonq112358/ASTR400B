import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class node:
    def __init__(self, x, y, px, py, m):
        self.m = m
        self.pos = np.array([x, y])
        self.mom = np.array([px, py])
        self.child = None

    
    
