from collections import OrderedDict

class SLBase():
    def __init__(self):
        # this values will contain OrderedDicts after the use of _fill_SL_dicts()
        self.SL_children_names = []
        self.SL_attributes_names = []
        
    def fill_SL_names(self):
        raise NotImplementedError        

    def transfer_internal_to_widget_tree(self):
        self.transfer_internal_to_widget()
        for child_name in self.SL_children_names:
            getattr(self,child_name).transfer_internal_to_widget_tree()
    
    def load_from_dict_tree(self,loaded_attributes_dictionary):
        for name in self.SL_attributes_names:
            setattr(self,name,loaded_attributes_dictionary[name])
        for child_name in self.SL_children_names:
            getattr(self,child_name).load_from_dict_tree(loaded_attributes_dictionary[child_name])
        

    def return_save_dict(self):
        attrs_dict = OrderedDict()
        for attr_name in self.SL_attributes_names:
            attrs_dict[attr_name] = getattr(self,attr_name)
        child_dict = OrderedDict()
        for child_name in self.SL_children_names:
            child_dict[child_name] = getattr(self,child_name).return_save_dict()
        return dict(attrs_dict, **child_dict)


class SLBaseWidget(SLBase):
    def transfer_internal_to_widget(self):
        raise NotImplementedError