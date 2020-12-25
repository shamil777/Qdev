from PyQt5.QtWidgets import QWidget,QListWidget
from PyQt5.QtCore import Qt

from .. SLBase import SLBaseWidget

class SimulationDatasetsWidget(QWidget,SLBaseWidget):
    def __init__(self, parent=None, flags=Qt.WindowFlags()):
        super(SimulationDatasetsWidget,self).__init__(parent,flags)
        SLBaseWidget.__init__(self)        
        self.fill_SL_names()

        self.ref_to_parent = parent        
        self.datasets_choose_list = QListWidget()
        self.datasetsFrame = self.ref_to_parent.simulator.simulation_datasets
    
    def init_GUI(self):
        names = self.datasetsFrame
    
    def fill_SL_names(self):
        pass
