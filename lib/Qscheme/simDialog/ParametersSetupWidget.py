from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

from collections import OrderedDict

from .SLBase import SLBaseWidget
from ..variables import Var
from .._KEYHASHABLE import VAR_KIND

class SymNum_LineEdit(QtWidgets.QLineEdit):
    def __init__(self,var_sym=None,col_idx=0,parent=None):
        super(SymNum_LineEdit,self).__init__(parent)
        self.var_sym = var_sym
        self.col_idx = col_idx
        
class Sym_ComboBox(QtWidgets.QComboBox):
    def __init__(self,var_sym=None,row_i=None,parent=None):
        super(Sym_ComboBox,self).__init__(parent)
        self.var_sym = var_sym
        self.row_i = row_i

class VarsGridWidget(QtWidgets.QWidget,SLBaseWidget):
    var_kinds_strs = ["fixed","sweep","equation"]
    sweep_var_fill_strs = ["start","stop","N points"]
    
    def __init__(self,var_list=[],parent=None,flags=Qt.WindowFlags()):
        super(VarsGridWidget,self).__init__(parent,flags)
        self._var_list = var_list
        
        self.var_rows_i = OrderedDict()
        self.var_sym_labels = OrderedDict()
        
        # combo Box related structures
        self.var_kinds_comboBoxes = OrderedDict()
        self.var_kinds = OrderedDict()        
        
        # editable text fields related structures
        self.var_intervals_lineEdits = OrderedDict()        
        self.var_settings = OrderedDict()
        
        # Save/Load attributes fields filling for SLBase class
        self.fill_SL_names()
        
        self.grid_layout = QtWidgets.QGridLayout()
        self.grid_layout.setSpacing(10)
        self.setLayout(self.grid_layout)
        
        self.ref_to_parent = self.parent()
        self.init_GUI()
        
    def init_GUI(self):
        for i,var in enumerate(self._var_list):
            self.set_row_variable_with_SL(i,var.sym)
            
    def delete_row_save_SL(self,row_i):
        # if grid widgets are not existing, no need to delete anything
        if( self.grid_layout.isEmpty() ):
            return
        
        # if this row is empty (this is verified via the first item in row,
        # because row is either fully filled either the whole row is nonEmpty )
        # then no need to delete anything
        grid_label_object = self.grid_layout.itemAtPosition(row_i,0)
        if( grid_label_object is None ):
            return
        
        row_label = self.grid_layout.itemAtPosition(row_i,0).widget()
        comboBox = self.grid_layout.itemAtPosition(row_i,1).widget()
        var_sym = Var(row_label.text()).sym

        for j in range(self.grid_layout.columnCount()):
            lineEditObj = self.grid_layout.itemAtPosition(row_i,j)
            if( lineEditObj is not None ):
                lineEditObj.widget().setParent(None)
        del self.var_intervals_lineEdits[var_sym]
        
        
        # delete comboBox
        comboBox.setParent(None)
        del self.var_kinds_comboBoxes[var_sym]
        # delete label
        row_label.setParent(None)
        del self.var_sym_labels[var_sym]
        
        return var_sym
            
    def delete_row_with_SL(self,row_i,var_sym):
        self.delete_row_save_SL(row_i)
        
        ## erasing entries from Save/Load dictionaries
        for attr_name_SL in self.SL_attributes_names:
            del getattr(self,attr_name_SL)[var_sym]
            
    def clear_grid_with_SL(self):
        for row_i in reversed(range(self.grid_layout.rowCount())):
            self.delete_row_save_SL(row_i)
        
        for attr_name_SL in self.SL_attributes_names:
            setattr(self,attr_name_SL,OrderedDict())
        
    def reinit(self,var_list):
        # clear all rows
        
        self.clear_grid_with_SL()
            
        for i,var in enumerate(var_list):
            self.set_row_variable_with_SL(i,var.sym)
        
        
    def set_row_variable_with_SL(self,row_i,var_sym):
        '''
        @description: 
            function that adds new widgets row into grid.
                If variable with 'var_sym' symbol is already in Save/Load structures
            then editable values in row are filled with respect to data stored in 
            this structures.
                Otherwise the row has default format and Save/Load structures are filled
            with default values.
                If this row is already occupied with another variable,
            this variable will be permamently deleted from all structures like it
            had never existed.
            
        @arguments:
            row_i - number of row you wish to insert your new variable
            var_sym - symbol of the variable you wish to insert
        @return: None
        '''        

        row_label_object = self.grid_layout.itemAtPosition(row_i,0)
        
        # if this row is occupied already
        if( row_label_object is not None ):
            row_label = row_label_object.widget()
            row_var_sym = Var(row_label.text()).sym
            # row occupied with the same variable
            if( var_sym == row_var_sym ):
                self.delete_row_save_SL(row_i)
            else: # if this row is occupied with another variable
                self.delete_row_with_SL(row_i,row_var_sym)
                
        
        self.var_rows_i[var_sym] = row_i
        
        var_sym_label = QtWidgets.QLabel(self)
        var_sym_label.setText(str(var_sym))
        self.var_sym_labels[var_sym] = var_sym_label
        self.grid_layout.addWidget(var_sym_label,row_i,0)
        
        var_kind_ComboBox = Sym_ComboBox(var_sym,row_i,self)
        var_kind_ComboBox.addItems([VAR_KIND.FIXED,
                                   VAR_KIND.SWEEP,
                                   VAR_KIND.EQUATION])
        self.var_kinds_comboBoxes[var_sym] = var_kind_ComboBox
        if( var_sym in self.var_kinds and self.var_kinds[var_sym] is not None ):
            var_kind_ComboBox.setCurrentText(self.var_kinds[var_sym])
        else:
            self.var_kinds[var_sym] = var_kind_ComboBox.currentText()
            
        var_kind_ComboBox.currentIndexChanged.connect(self.var_kind_ComboBox_CIC_handler)
            
        self.grid_layout.addWidget(var_kind_ComboBox,row_i,1)
        
        var_kind = self.var_kinds[var_sym]
        
        if( var_kind == VAR_KIND.FIXED ):
            fixed_var_lineEdit = SymNum_LineEdit(var_sym,self)
            fixed_var_lineEdit.editingFinished.connect(self.fixed_var_lineEdit_editingFinished_handler)
            self.var_intervals_lineEdits[var_sym] = [fixed_var_lineEdit]
            if( var_sym in self.var_settings and self.var_settings[var_sym] is not None):
                fixed_var_lineEdit.setText(str(self.var_settings[var_sym]))
            else:
                self.var_settings[var_sym] = None
            self.grid_layout.addWidget(fixed_var_lineEdit,row_i,2)
            
        elif( var_kind == VAR_KIND.SWEEP ):
            self.var_intervals_lineEdits[var_sym] = [SymNum_LineEdit(var_sym,0,self),
                                                         SymNum_LineEdit(var_sym,1,self), 
                                                         SymNum_LineEdit(var_sym,2,self)]
            for j,lineEdit in enumerate(self.var_intervals_lineEdits[var_sym]):
                lineEdit.editingFinished.connect(self.intervals_lineEdit_editingFinished_handler)
                if( var_sym in self.var_settings and self.var_settings[var_sym][j] is not None ):
                    lineEdit.setText(str(self.var_settings[var_sym][j]))
                else:
                    self.var_settings[var_sym] = [None]*3
                self.grid_layout.addWidget(lineEdit,row_i,j+2)
            
        elif( var_kind == VAR_KIND.EQUATION ):
            equation_var_lineEdit = SymNum_LineEdit(var_sym,self)
            equation_var_lineEdit.editingFinished.connect(self.equation_LineEdit_editingFinished_handler)
            self.var_intervals_lineEdits[var_sym] = [equation_var_lineEdit]
            if( var_sym in self.var_settings and self.var_settings[var_sym] is not None ):
                equation_var_lineEdit.setText(self.var_settings[var_sym])
            else:
                self.var_settings[var_sym] = None
            
            self.grid_layout.addWidget(equation_var_lineEdit,row_i,2)
    
    # Save/Load attributes fields filling for SLBase class        
    def fill_SL_names(self):
        self.SL_attributes_names = ["var_kinds",
                                    "var_settings",
                                    "var_rows_i"]
    
    def transfer_internal_to_widget(self):
        for row_i in reversed(range(self.grid_layout.rowCount())):
            self.delete_row_save_SL(row_i)
        
        for i,var_sym in enumerate(self.var_kinds):
            self.set_row_variable_with_SL(i,var_sym)
            
            
    ## SLOTS SECTION START ##
    def var_kind_ComboBox_CIC_handler(self):
        CB_source = self.sender()
        current_kind = CB_source.currentText()
        var_sym = CB_source.var_sym
        row_i = CB_source.row_i
        # erasing row like it never existed
        self.delete_row_with_SL(row_i,var_sym)
        
        # setting SL variable to be used in set_row_variable_with_SL
        self.var_kinds[var_sym] = current_kind
        self.set_row_variable_with_SL(row_i,var_sym)        
    
    def fixed_var_lineEdit_editingFinished_handler(self):
        lineEdit = self.sender()
        try:
            val = float(lineEdit.text())
        except ValueError:
            self.var_settings[lineEdit.var_sym] = None
            return
        self.var_settings[lineEdit.var_sym] = val
        
    
    def intervals_lineEdit_editingFinished_handler(self):
        lineEdit = self.sender()
        var_sym = lineEdit.var_sym
        col_idx = lineEdit.col_idx
        text_to_parse = lineEdit.text()
        try:
            if( col_idx == 0 or col_idx == 1 ):
                    val = float(text_to_parse)
            elif( col_idx == 2 ):
                    val = int(text_to_parse)
        except ValueError:
            self.var_settings[var_sym][col_idx] = None
            return
        
        self.var_settings[var_sym][col_idx] = val
            
    def equation_LineEdit_editingFinished_handler(self):
        lineEdit = self.sender()
        self.var_settings[lineEdit.var_sym] = lineEdit.text()
    ## SLOTS SECTION END ##
    

class ParametersSetupWidget(QtWidgets.QWidget,SLBaseWidget):
    def __init__(self, parent=None, flags=Qt.WindowFlags()):
        super(ParametersSetupWidget,self).__init__(parent,flags)
        self.ref_to_parent = self.parent()
        
        self.fill_SL_names()        
    
    def fill_SL_names(self):
        self.SL_attributes_names = []
        self.SL_children_names = ["scheme_vars_grid",
                                  "aux_vars_grid"] 
        
    def init_GUI(self):
        v_layout = QtWidgets.QVBoxLayout()
        self.setLayout(v_layout)
        
        vars_list = self.ref_to_parent.simulator.scheme.params
        self.scheme_vars_grid = VarsGridWidget(vars_list)
        self.aux_vars_grid = VarsGridWidget()
        
        ## adding 2 grid for scheme and auxillary variables
        v_layout.addWidget(self.scheme_vars_grid)
        v_layout.addWidget(self.aux_vars_grid)
        
        # adding label for parameter adding/deleting section
        add_param_label= QtWidgets.QLabel(self)
        add_param_label.setText("add parameter")
        v_layout.addWidget(add_param_label)
        
        # adding/deleting interface occupies the single 
        # row with layout 'add_param_row_hLayout'
        add_param_row_hLayout = QtWidgets.QHBoxLayout()
        v_layout.addLayout(add_param_row_hLayout)  
        
        # user navigation label
        add_param_lineEdit_label = QtWidgets.QLabel()
        add_param_lineEdit_label.setText("parameter name:")
        add_param_row_hLayout.addWidget(add_param_lineEdit_label)
        
        # name of the variable that is needed to be added/deleted to/from
        # auxillary variables grid
        self.add_param_sym_lineEdit = QtWidgets.QLineEdit()
        add_param_row_hLayout.addWidget(self.add_param_sym_lineEdit)
        
        # add button
        add_var_button = QtWidgets.QPushButton("add",self)
        add_var_button.clicked.connect(self.add_var_button_clicked_handler)
        add_param_row_hLayout.addWidget(add_var_button)
        
        # delete button
        delete_var_button = QtWidgets.QPushButton("delete",self)
        delete_var_button.clicked.connect(self.delete_var_button_clicked_handler)
        add_param_row_hLayout.addWidget(delete_var_button) 
        self.show()

    
    def add_var_button_clicked_handler(self):
        text = self.add_param_sym_lineEdit.text()
        if( text == "" ):
            return
        new_var = Var(text)
        
        row_i = self.aux_vars_grid.grid_layout.rowCount()
        self.aux_vars_grid.set_row_variable_with_SL(row_i,new_var.sym)
    
    def delete_var_button_clicked_handler(self):
        text = self.add_param_sym_lineEdit.text()
        if( text == "" ):
            return
        
        var_to_del = Var(text)
        
        row_i = self.aux_vars_grid.var_rows_i[var_to_del.sym]
        self.aux_vars_grid.delete_row_with_SL(row_i,var_to_del.sym)    
    
    def scheme_subscripts_changed_handler(self):
        var_list = self.ref_to_parent.simulator.scheme.params.values()
        self.scheme_vars_grid.reinit(list(var_list))
    
    def transfer_internal_to_widget(self):
        pass
        