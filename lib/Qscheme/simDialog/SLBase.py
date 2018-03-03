from collections import OrderedDict

class SLBase():
    def __init__(self):
        self.SL_children_dict = OrderedDict()
        self.SL_attributes_dict = OrderedDict()
    
    def _fill_SL_dicts(self):
        raise NotImplementedError
    
    def transfer_internal_to_widget(self):
        raise NotImplementedError

    def transfer_internal_to_widget_tree(self):
        self.transfer_internal_to_widget()
        for child in self.SL_children_dict.values():
            child.transfer_internal_to_widget_tree()
    
    def load_from_dict_tree(self,loaded_attributes_dictionary):
   
        for load_attr_key,load_attr_val in loaded_attributes_dictionary.items():
            if( load_attr_key in self.SL_attributes_dict ):
                self.SL_attributes_dict[load_attr_key] = load_attr_val
                setattr(self,load_attr_key,load_attr_val)
            elif( load_attr_key in self.SL_children_dict ):
                self.SL_children_dict[load_attr_key].load_from_dict_tree(load_attr_val)
        

    def return_save_dict(self):
        child_safe_dict = OrderedDict()

        self._fill_SL_dicts()
        for child_key,child_class in self.SL_children_dict.items():
            child_safe_dict[child_key] = child_class.return_save_dict()

        return dict(self.SL_attributes_dict, **child_safe_dict)
            