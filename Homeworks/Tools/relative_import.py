# Add "../Tools" to the path so that I can import ReadFile
# and ParticleProperties without copying them into the 
# directory each time.
import sys
import os

tools_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../Tools"))

sys.path.append(tools_path)


