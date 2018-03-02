class SLBase():
    def __init__(self):
        self.SL_children = []
    
    def load_from_dict(self,loaded_attributes_dictionary):
        safe_dict = self.__dict__.copy()
        del safe_dict["SL_children"]
        
        for attr_key in safe_dict:
            try:
                setattr(self,attr_key,loaded_attributes_dictionary[attr_key])
            except KeyError:
                print("key error while loading in class ",self.__class__)
                print("error key: ",attr_key)
            
        for child_class in self.SL_children:
            print(child_class)
            child_class.load_from_dict(loaded_attributes_dictionary)
            
        
    def transfer_internal_to_widget(self):
        raise NotImplementedError
    
    def load_from_dict_and_visualize(self,loaded_attributes_dictionary):
        safe_dict = self.__dict__.copy()
        del safe_dict["SL_children"]
        
        for attr_key in safe_dict:
            try:
                setattr(self,attr_key,loaded_attributes_dictionary[attr_key])
            except KeyError:
                print("key error while loading in class ",self.__class__)
                print("error key: ",attr_key) 
            
        for child_class in self.SL_children:
            child_class.load_from_dict_and_visualize(loaded_attributes_dictionary)
            
        self.transfer_internal_to_widget()
        