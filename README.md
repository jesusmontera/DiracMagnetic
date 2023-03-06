# Dirac-Pauli-Schrod Magnetic Python 3.6 OpenGL QT5
Apply magnetic fields to Dirac equation spinors in python(under construction).
Based on marl0ny Dirac examples from Github.https://github.com/marl0ny/split-operator-simulations 

Works by making a numpy.dot between magnetic B vectors and the spin extracted from the Dirac psi (or a specific spin for Schroduinger equation) or applying B as potential(only in Pauli mode). 
<BR>The .ui files from qt5 designer, were converterd to py files with command line: pyuic5 xxxx.ui -o xxxx-py
Each dialog have it's own myxxxx that import and implement functionality to the  orginal dialog's.

For a good B Gerlach field load it from a npy file, that can be maked with magpylib(in that folder there is the script gerlach.py for make the B npy file. Use the B dialog, menu file open from the main aplication to load the npy file and press ok button. 




