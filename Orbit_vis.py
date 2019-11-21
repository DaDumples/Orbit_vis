from numpy import *
import os
import sys


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

PARAMETERS = ['Semimajor Axis', 'Eccentricity', 'True Anomaly', 'Inclination', 'Right Ascention', 'Argument of Perigee']

def plot_earth(ax, color = 'b', radius = 6378, alpha = .25, resolution = 100):
    # Make data
    u = linspace(0,2*pi, resolution)
    v = linspace(0,pi, resolution)
    x = radius * outer(cos(u), sin(v))
    y = radius * outer(sin(u), sin(v))
    z = radius * outer(ones(size(u)), cos(v))

    # Plot the surface
    ax.plot_surface(x, y, z, color=color,alpha = alpha)

def coes2rv(a,ecc,TA,inc,RAAN,Arg_p,mu=398600):
    
    radTA = radians(TA)
    
    rp = a*(1-ecc**2)/(1+ecc)
    h = sqrt(rp*mu*(1+ecc))
    r_perifocal = (h**2/mu)*(1/(1 + ecc*cos(radTA)))*array([cos(radTA),sin(radTA),0]).T
    v_perifocal = (mu/h)*array([-sin(radTA),ecc+cos(radTA),0]).T
    ECI2P = Cz(Arg_p)@Cx(inc)@Cz(RAAN)
    P2ECI = ECI2P.T
    r = P2ECI@r_perifocal
    v = P2ECI@v_perifocal
    
    return r, v

def Cz(ang):
    ang = radians(ang)
    c = cos(ang)
    s = sin(ang)
    
    C = array([[c,s,0],
                   [-s,c,0],
                   [0,0,1]])
    return C

def Cy(ang):
    ang = radians(ang)
    c = cos(ang)
    s = sin(ang)
    
    C = array([[c,0,-s],
                   [0,1,0],
                   [s,0,c]])
    
    return C

def Cx(ang):
    ang = radians(ang)
    c = cos(ang)
    s = sin(ang)
         
    C = array([[1,0,0],
                       [0,c,s],
                       [0,-s,c]])
    
    return C

class Window(QMainWindow):

    def __init__(self, parent = None):
        QMainWindow.__init__(self)
        self.wid = QWidget()
        self.setCentralWidget(self.wid)

        self.setFixedSize(780, 480)
        

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.ax = Axes3D(self.figure)

        self.layout = QHBoxLayout()
        self.slider_menu = QVBoxLayout()
        self.sliders = []
        self.labels = []
        self.inputs = zeros(len(PARAMETERS))

        for param in PARAMETERS:
            s = QSlider(Qt.Horizontal)
            s.setMaximum(100)
            s.setMinimum(0)
            s.setTickInterval(1)
            s.setValue(0)
            s.valueChanged.connect(self.update_inputs)

            l = QLabel(param)
            self.sliders.append(s)
            self.labels.append(l)

            self.slider_menu.addWidget(l)
            self.slider_menu.addWidget(s)

        self.update_inputs()

            
        self.layout.addWidget(self.canvas)
        self.layout.addLayout(self.slider_menu)

        self.wid.setLayout(self.layout)
        self.show()

    def update_inputs(self):
        for i in range(len(PARAMETERS)):
            if i == 0:
                x = self.sliders[i].value()
                semi_a = 6378+400 + x*10000/100
                self.inputs[i] = semi_a
                self.labels[i].setText(PARAMETERS[i] + ' : '+str(semi_a))
            elif i == 1:
                x = self.sliders[i].value()
                ecc = x/100*.3
                self.inputs[i] = ecc
                self.labels[i].setText(PARAMETERS[i] + ' : {0:.2f}'.format(ecc))
            elif i >= 2:
                x = self.sliders[i].value()
                ang = x/100*360
                self.inputs[i] = ang
                self.labels[i].setText(PARAMETERS[i] + ' : {0:.2f}'.format(ang))


        points = []
        num_pts = 100
        for i in range(360):
            TA = i/num_pts*360
            r, _ = coes2rv(self.inputs[0], self.inputs[1], TA, self.inputs[3], self.inputs[4], self.inputs[5])
            points.append(r)
        points = vstack(points)

        peri, v_peri = coes2rv(self.inputs[0], self.inputs[1], 0, self.inputs[3], self.inputs[4], self.inputs[5])
        apo, _ = coes2rv(self.inputs[0], self.inputs[1], 180, self.inputs[3], self.inputs[4], self.inputs[5])
        p_90, _ = coes2rv(self.inputs[0], self.inputs[1], 90, self.inputs[3], self.inputs[4], self.inputs[5])
        n_90 = -p_90

        h = cross(peri, v_peri)/linalg.norm(cross(peri, v_peri))*6378*2
        node_a, _ = coes2rv(self.inputs[0], self.inputs[1], self.inputs[4], 0,0,0)
        node_d, _ = coes2rv(self.inputs[0], self.inputs[1], self.inputs[4]+180, 0, 0, 0)

        pos, _ = coes2rv(self.inputs[0], self.inputs[1], self.inputs[2], self.inputs[3], self.inputs[4], self.inputs[5])

        self.ax.clear()
        self.ax.plot(points[:,0], points[:,1], points[:,2], 'k-')
        self.ax.plot([0, 0], [0,0], [-2*6378, 2*6378],'k')
        self.ax.plot([0, 0], [-2*6378, 2*6378], [0,0],'k')
        self.ax.plot([-2*6378, 2*6378], [0, 0], [0,0], 'k')

        self.ax.plot([0, h[0]], [0, h[1]], [0, h[2]], 'orange')

        self.ax.plot([apo[0], peri[0]], [apo[1], peri[1]], [apo[2], peri[2]], 'red')
        self.ax.plot([p_90[0], n_90[0]], [p_90[1], n_90[1]], [p_90[2], n_90[2]], 'red')

        self.ax.plot([node_a[0], node_d[0]], [node_a[1], node_d[1]], [node_a[2], node_d[2]], 'g--')

        self.ax.scatter(pos[0], pos[1], pos[2], 'ro')


        plot_earth(self.ax, resolution = 10)

        scale = linalg.norm(apo)
        self.ax.set_xlim(-scale, scale)
        self.ax.set_ylim(-scale, scale)
        self.ax.set_zlim(-scale, scale)


        self.ax.set_axis_off()

        self.canvas.draw()




if __name__ == '__main__':
    app = QApplication(sys.argv)

    main = Window()
    main.show()
    sys.exit(app.exec_())