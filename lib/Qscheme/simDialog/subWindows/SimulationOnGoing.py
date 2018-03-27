from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

class SimulationOnGoingWindow(QtWidgets.QWidget):
    def __init__(self,parent=None):
        super(SimulationOnGoingWindow,self).__init__(parent, Qt.Window)
        self.setWindowTitle("simulation report")
        self.init_GUI()
        
    def init_GUI(self):        
        self.vlayout = QtWidgets.QVBoxLayout()
        self.setLayout(self.vlayout)        
        
        self.progress_bar = QtWidgets.QProgressBar(self)
        self.progress_bar.setValue(0)
        self.vlayout.addWidget(self.progress_bar)
        
        self.time_left_label = QtWidgets.QLabel(self)
        self.time_left_label.setText("time left: ")
        self.vlayout.addWidget(self.time_left_label)
    
    def update_progress(self,simulator,percent=0):
        print("update progress: ",percent)
        self.progress_bar.setValue(percent*100)
        self.time_left_label.setText("time left: {}".format(simulator.simulation_start_time))
        self.repaint()
        print(self.parent())