# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 21:23:10 2017

@author: user
"""
import sys
from PyQt5 import QtGui,QtCore,QtWidgets

class MWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MWindow,self).__init__()
        self.setGeometry(100,100,500,300)
        self.setWindowTitle("Qcheme")
        
        extractAction = QtWidgets.QAction("&GET TO THE CHOPPA",self)
        extractAction.setShortcut("Ctrl+Q")
        extractAction.setStatusTip("Leave application")
        extractAction.triggered.connect(self.close_application) 
        
        openEditor = QtWidgets.QAction("&editor",self)
        openEditor.setShortcut("Ctrl+E")
        openEditor.setStatusTip("Open editor")
        openEditor.triggered.connect(self.editor)
        
        openFile = QtWidgets.QAction("&open file",self)
        openFile.setShortcut("Ctrl+O")
        openFile.setStatusTip("file open dialog pops up")
        openFile.triggered.connect(self.file_open)
        
        self.statusBar()
        
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu("&File")
        fileMenu.addAction(extractAction)
        
        editorMenu = mainMenu.addMenu("Editor")
        editorMenu.addAction(openEditor)
        editorMenu.addAction(openFile)
        
        self.home()
        
    def home(self):
        btn = QtWidgets.QPushButton("Quit",self)
        btn.clicked.connect(self.close_application)
        
        btn.resize(btn.minimumSizeHint())
        btn.move(100,100)
        
        extractAction = QtWidgets.QAction(QtGui.QIcon("icon.png"),"action in home",self)
        extractAction.triggered.connect(self.close_application)
        self.toolbar = self.addToolBar("Exctraction")
        self.toolbar.addAction(extractAction)

        extractAction = QtWidgets.QAction("Font1",self)
        extractAction.triggered.connect(self.font_choice)        
        self.toolbar = self.addToolBar("Font2")
        self.toolbar.addAction(extractAction)
        
        checkBox = QtWidgets.QCheckBox("Enlarge window",self)
        checkBox.move(200,200)
        checkBox.stateChanged.connect(self.enlarge_window)
        checkBox.toggle()

        self.progress = QtWidgets.QProgressBar(self)
        self.progress.setGeometry(200,400,250,20 )
        
        self.btn = QtWidgets.QPushButton("Download",self)
        self.btn.move(550,400)
        self.btn.clicked.connect(self.download_btn)
        
        print(self.style().objectName())
        self.styleChoice = QtWidgets.QLabel("Windows", self)
        
        comboBox = QtWidgets.QComboBox(self)
        comboBox.addItem("motif")
        comboBox.addItem("Windows")
        comboBox.addItem("cde")
        comboBox.addItem("Plastique")
        comboBox.addItem("Cleanlooks")
        comboBox.addItem("windowsvista")
        
        comboBox.move(500,500)
        self.styleChoice.move(220,220)
        comboBox.activated[str].connect(self.combo_style)
        
        color = QtGui.QColor(0,0,0)
        fontColor = QtWidgets.QAction("font bg color",self)
        fontColor.triggered.connect(self.color_picker)
        self.toolbar.addAction(fontColor)

        cal = QtWidgets.QCalendarWidget(self)
        cal.resize(200,200)
        cal.move(600,600)
        
        self.show()

    def file_open(self):
        name = QtWidgets.QFileDialog.getOpenFileName(self,"Open File Name")[0]
        print(name)
        file = open(name,"r")
        
        self.editor()
        
        with file:
            text = file.read()
            self.textEditor.setText(text)

    def editor(self):
        self.textEditor = QtWidgets.QTextEdit()
        self.setCentralWidget(self.textEditor)

    def color_picker(self):
        color = QtWidgets.QColorDialog.getColor()
        self.styleChoice.setStyleSheet("QWidget { background-color: %s}" % color.name())

    def font_choice(self):
        font, valid = QtWidgets.QFontDialog.getFont()
        if valid:
            self.styleChoice.setFont(font)

    def combo_style(self, text):
        self.styleChoice.setText(text)
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create(text))
        
    def download_btn(self):
        self.completed = 0
        
        while( self.completed < 100 ):
            self.completed += 0.00001
            self.progress.setValue(self.completed)
    
        
    def enlarge_window(self, state):
        if( state == QtCore.Qt.Checked ):
            self.setGeometry(100,100,1000,600)
        else:
            self.setGeometry(100,100,500,300)
    
    def close_application(self):
        choice = QtWidgets.QMessageBox.question(self, 
                                                "Exctract!", 
                                                "Are you shure?",
                                                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if( choice == QtWidgets.QMessageBox.Yes ):
            print("quitting")
            sys.exit()
        else:
            print("lol")
        
        


def run():
    app = QtCore.QCoreApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)
        
    
    window = MWindow()
    app.exec_()
    
run()