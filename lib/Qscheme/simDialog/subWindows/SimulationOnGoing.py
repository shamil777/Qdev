from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

class SimulationOnGoingWindow(QtWidgets.QWidget):
    def __init__(self,parent=None):
        super(SimulationOnGoingWindow,self).__init__(parent, Qt.Window)
        self.ref_to_parent= self.parent()
        self.setWindowTitle("simulation report")
        self.init_GUI()
        
    def init_GUI(self):
        self.vlayout = QtWidgets.QVBoxLayout()
        self.setLayout(self.vlayout)        
        
        self.progress_bar = QtWidgets.QProgressBar(self)
        self.progress_bar.setValue(0)
        self.vlayout.addWidget(self.progress_bar)

        self.start_time_label = QtWidgets.QLabel(self)
        self.start_time_label.setText("start time: ")
        self.vlayout.addWidget(self.start_time_label)
        
        self.time_left_label = QtWidgets.QLabel(self)
        self.time_left_label.setText("approx time left: ")
        self.vlayout.addWidget(self.time_left_label)
        
        self.end_time_label = QtWidgets.QLabel(self)
        self.end_time_label.setText("approx end time: ")
        self.vlayout.addWidget(self.end_time_label)
        
        self.cancel_btn = QtWidgets.QPushButton(self)
        self.cancel_btn.setText("Cancel")
        self.cancel_btn.clicked.connect(self.cancel_btn_clicked_handler)
        self.vlayout.addWidget(self.cancel_btn)
        
        self.vlayout.setAlignment(Qt.AlignTop)
        
        self.setFixedWidth(400)
        
    def cancel_btn_clicked_handler(self):
        self.ref_to_parent.simulator.cancel_calculation()
        
    datetime_format_str = "%Y-%m-%d %H:%M:%S"
    
    def _deltatime_format_text(self,time_d):
        total_seconds = int(time_d.total_seconds())
        days, rem_d = divmod( total_seconds, 24*60*60 )
        hours, rem_h = divmod(rem_d, 60*60)
        minutes, seconds = divmod(rem_h,60)
        
        return "{}d  {}h {}m {}s".format( days,hours,minutes,seconds )
    
    def update_progress(self,progress_timer=None):
        if( progress_timer is not None ):
            self.progress_bar.setValue(progress_timer.percentage*100)
            self.start_time_label.setText("start time: {}".format(progress_timer.start_time.strftime(self.datetime_format_str)))
            self.time_left_label.setText("approx time left: " + self._deltatime_format_text(progress_timer.approx_process_duration))
            self.end_time_label.setText("approx end time: {}".format(progress_timer.approx_process_end.strftime(self.datetime_format_str)))
        
        QtWidgets.QApplication.instance().processEvents()
        