import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Body:
    def __init__(self, mass, position, velocity):
        self.m = mass
        self.p = position
        self.v = velocity
    

