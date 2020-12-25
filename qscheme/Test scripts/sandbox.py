import appdirs
import os

if ( __name__ == "__main__" ):
    
    programData_dir = appdirs.user_data_dir()
    programFolder_name = "Qdev"
    programData_fileName = "Qdev_data.txt"
    
    print(appdirs.user_data_dir())
    if( not os.path.isdir(programData_dir + "\\" + programFolder_name) ):
        os.mkdir(programData_dir + "\\" + programFolder_name)
    
    