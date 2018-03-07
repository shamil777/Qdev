from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

from collections import OrderedDict

from .SLBase import SLBase
from ..variables import Var

class SymNum_LineEdit(QtWidgets.QLineEdit):
    def __init__(self,param_sym=None,col_idx=0,parent=None):
        super(SymNum_LineEdit,self).__init__(parent)
        self.param_sym = param_sym
        self.col_idx = col_idx
        
class Sym_ComboBox(QtWidgets.QComboBox):
    def __init__(self,param_sym=None,parent=None):
        super(Sym_ComboBox,self).__init__(parent)
        self.param_sym = param_sym


class ParametersSetupWidget(QtWidgets.QWidget,SLBase):
    param_kinds_strs = ["fixed","sweep","equation"]
    sweep_param_fill_strs = ["start","stop","N points"]
    
    def __init__(self, parent=None, flags=Qt.WindowFlags()):
        super(ParametersSetupWidget,self).__init__(parent,flags)
        self.param_sym_labels = OrderedDict()
        
        self.param_intervals_lineEdits = OrderedDict()        
        self.scheme_param_settings = OrderedDict()
        self.aux_param_settings = OrderedDict()
        
        self.param_kinds_comboBoxes = OrderedDict()
        self.param_kinds = OrderedDict()
        
        self.fill_SL_names()
        
        self.scheme_grid_layout = QtWidgets.QGridLayout()
        self.scheme_grid_layout.setSpacing(10)
        self.aux_grid_layout = QtWidgets.QGridLayout()
        self.aux_grid_layout.setSpacing(10)
        
        
        self.ref_to_parent = self.parent()
        
        self.init_GUI()
        
        
        
    def init_GUI(self):
        v_layout = QtWidgets.QVBoxLayout()
        self.setLayout(v_layout)
        
        v_layout.addLayout(self.scheme_grid_layout)
        v_layout.addLayout(self.aux_grid_layout)
        self.refresh_grid()
        
        add_param_label= QtWidgets.QLabel(self)
        add_param_label.setText("add parameter")
        v_layout.addWidget(add_param_label)
        
        add_param_row_hLayout = QtWidgets.QHBoxLayout()
        
        add_param_lineEdit_label = QtWidgets.QLabel()
        add_param_lineEdit_label.setText("parameter name:")
        add_param_row_hLayout.addWidget(add_param_lineEdit_label)
        
        self.add_param_sym_lineEdit = QtWidgets.QLineEdit()
        add_param_row_hLayout.addWidget(self.add_param_sym_lineEdit)
        
        add_var_button = QtWidgets.QPushButton("add",self)
        add_var_button.clicked.connect(self.add_var_button_clicked_handler)
        add_param_row_hLayout.addWidget(add_var_button)
        
        v_layout.addLayout(add_param_row_hLayout)        
        self.show()
    
    def fill_SL_names(self):
        self.SL_attributes_names = ["scheme_param_settings",
                                    "aux_param_settings",
                                    "param_kinds"]

    
    def add_var_button_clicked_handler(self):
        text = self.add_param_sym_lineEdit.text()
        if( text == "" ):
            return
        new_var = Var(text)
        
        new_param_sym_label = QtWidgets.QLabel()
        new_param_sym_label.setText(str(new_var.sym))
        new_param_start_lineEdit = SymNum_LineEdit(new_var.sym,self)
        new_param_stop_lineEdit = SymNum_LineEdit(new_var.sym,self)
        new_param_nPts_lineEdit = SymNum_LineEdit(new_var.sym,self)
        self.param_intervals_lineEdits[new_var.sym] = [new_param_start_lineEdit, 
                                                     new_param_stop_lineEdit,
                                                     new_param_nPts_lineEdit]
        self.scheme_param_settings[new_var.sym] = None
        new_row = self.scheme_grid_layout.rowCount()
        self.scheme_grid_layout.addWidget(new_param_sym_label,new_row,0)
        self.scheme_grid_layout.addWidget(new_param_start_lineEdit,new_row,1)
        self.scheme_grid_layout.addWidget(new_param_stop_lineEdit,new_row,2)
        self.scheme_grid_layout.addWidget(new_param_nPts_lineEdit,new_row,3)
    
    
    
    def set_row_by_sym(self,row_i,param_sym):
        grid_layout = None
        param_settings = None
        scheme_params = list(self.ref_to_parent.simulator.scheme.params.values())
        scheme_syms = [param.sym for param in scheme_params]
        
        if( param_sym in scheme_syms ):
            grid_layout = self.scheme_grid_layout
            param_settings = self.scheme_param_settings
        else:
            grid_layout = self.aux_grid_layout
            param_settings = self.aux_param_settings
        
        param_sym_label = QtWidgets.QLabel(self)
        param_sym_label.setText(str(param_sym))
        self.param_sym_labels[param_sym] = param_sym_label
        grid_layout.addWidget(param_sym_label,row_i+1,0)
        
        param_kind_ComboBox = Sym_ComboBox(param_sym,self)
        param_kind_ComboBox.addItems(self.param_kinds_strs)
        self.param_kinds_comboBoxes[param_sym] = param_kind_ComboBox
        if( param_sym in self.param_kinds and self.param_kinds[param_sym] is not None ):
            param_kind_ComboBox.setCurrentIndex(
                    self.param_kinds_strs.index(self.param_kinds[param_sym])
                    )
        else:
            self.param_kinds[param_sym] = param_kind_ComboBox.currentText()
        param_kind_ComboBox.currentIndexChanged.connect(self.param_kind_ComboBox_CIC_handler)
            
        grid_layout.addWidget(param_kind_ComboBox,row_i+1,1)
        
        param_kind = self.param_kinds[param_sym]
        if( param_kind == "fixed" ):
            fixed_param_lineEdit = SymNum_LineEdit(param_sym,self)
            fixed_param_lineEdit.editingFinished.connect(self.fixed_param_lineEdit_editingFinished_handler)
            self.param_intervals_lineEdits[param_sym] = fixed_param_lineEdit
            param_settings[param_sym] = None
            grid_layout.addWidget(fixed_param_lineEdit,row_i+1,2)
            
        elif( param_kind == "sweep" ):
            self.param_intervals_lineEdits[param_sym] = [SymNum_LineEdit(param_sym,0,self),
                                                         SymNum_LineEdit(param_sym,1,self), 
                                                         SymNum_LineEdit(param_sym,2,self)]
            for j,lineEdit in enumerate(self.param_intervals_lineEdits[param_sym]):
                lineEdit.editingFinished.connect(self.intervals_lineEdit_editingFinished_handler)
                param_settings[param_sym] = [None]*3
                grid_layout.addWidget(lineEdit,row_i+1,j+2)
            
        elif( param_kind == "equation" ):
            equation_param_lineEdit = SymNum_LineEdit(param_sym,self)
            equation_param_lineEdit.editingFinished.connect(self.equation_LineEdit_editingFinished_handler)
            self.param_intervals_lineEdits[param_sym] = equation_param_lineEdit
            param_settings[param_sym] = ""
            
            grid_layout.addWidget(equation_param_lineEdit,row_i+1,2)
    
    def delete_sym(self,param_sym):
        if( param_sym in self.scheme_param_settings ):
            del self.scheme_param_settings[param_sym]
        elif( param_sym in self.aux_param_settings ):
            del self.aux_param_settings[param_sym]
        else:
            return
        # finally
        param_kind = self.param_kinds_comboBoxes[param_sym]
        if( param_kind in ["fixed","equation"] ):
            self.param_intervals_lineEdits[param_sym].setParent(None)
        elif( param_kind == "sweep" ):
            for lineEdit in self.param_intervals_lineEdits[param_sym]:
                lineEdit.setParent(None)
                
        self.param_sym_labels[param_sym].setParent(None)
        self.param_kinds_comboBoxes[param_sym].setParent(None)
    
    def refresh_scheme_grid(self):
        params = self.ref_to_parent.simulator.scheme.params
        
        for i,param in enumerate(params.values()):
            param_sym = param.sym
            self.delete_sym(param_sym)            
            self.set_row_by_sym(i,param.sym)
           
            
    def refresh_aux_params_grid(self):
        aux_params = list(self.aux_param_settings.keys())
        for i,aux_param_sym in enumerate(aux_params):
            self.delete_sym(aux_param_sym)
            self.set_row_by_sym(i,aux_param_sym)
            
    def refresh_grid(self):
        self.refresh_scheme_grid()
        self.refresh_aux_params_grid()
        
    def subscripts_changed(self):       
        self.refresh_scheme_grid()
            
    
    def param_kind_ComboBox_CIC_handler(self):
        source = self.sender()
        print("combo Box handler ",str(source.param_sym))
        self.delete_sym(source.param_sym)

        self.param_kinds[source.param_sym] = source.currentText()
        row_i = list(self.param_kinds.keys()).index(source.param_sym)
        self.set_row_by_sym(row_i,source.param_sym)
        
        
        
    
    def fixed_param_lineEdit_editingFinished_handler(self):
        source = self.sender()
        try:
            val = float(source.text())
        except ValueError:
            return
        if( source.sym in self.scheme_param_settings ):
            self.scheme_param_settings[source.param_sym] = val
        elif( source.sym in self.aux_param_settings ):
            self.aux_param_settings[source.param_sym] = val
    
    def intervals_lineEdit_editingFinished_handler(self):
        signal_source = self.sender()
        param_sym = signal_source.param_sym
        col_idx = signal_source.col_idx
        
        text_to_parse = signal_source.text()
        try:
            if( col_idx == 0 or col_idx == 1 ):
                    val = float(text_to_parse)
            elif( col_idx == 2 ):
                    val = int(text_to_parse)
        except ValueError:
            self.scheme_param_settings[param_sym][col_idx] = None
            return
        
        if( param_sym in self.scheme_param_settings ):
            self.scheme_param_settings[param_sym][col_idx] = val
        elif( param_sym in self.aux_param_settings ):
            self.aux_param_settings[param_sym][col_idx] = val
            
    def equation_LineEdit_editingFinished_handler(self):
        source = self.sender()
        if( source.param_sym in self.scheme_param_settings ):
            self.scheme_param_settings[source.param_sym] = source.text()
        elif( source.param_sym in self.aux_param_settings ):
            self.aux_param_settings[source.param_sym] = source.text()
                
    def transfer_internal_to_widget(self):
        params_settings_tmp = self.scheme_param_settings.copy()
        self.refresh_grid() # refreshing based on simulator class
        self.scheme_param_settings = params_settings_tmp
        for param_sym,param_vals in self.scheme_param_settings.items():
            for j in range(0,3):
                text = str(param_vals[j])
                if( text == "None" ):
                    text = ""
                self.param_intervals_lineEdits[param_sym][j].setText(text)
        