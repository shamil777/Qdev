import re
import csv

import networkx as nx

from Qscheme.schematic.elements import SIS_JJ, Cap, Battery

def _get_refDes_groupNames_PADS( file_rows ):
    ''' 
    @description:
        parses file and return edge types, represented by electrical elements
    @parameters:
        file_rows - list of lists containing csv file (e.g. list( csv.reader(file)) 
    @return: 
        "reference designators" and "group names" of objects used
        in file that represents schematic
    '''
    refDes_list = []
    group_name_list = []
    
    element_i = -1
    for i,row in enumerate( file_rows ):
        if( row[0] == "*PART*" ):
            element_i = i+1
            break
    
    row = file_rows[element_i]
    while( len(row) > 0 ):
        refDes = row[0]
        group_name = row[5].split(',')[1]
        
        refDes_list.append(refDes)
        group_name_list.append(group_name)
        
        
        element_i += 1
        row = file_rows[element_i]
        
    return refDes_list, group_name_list
	
def build_graph_from_netlist_PADS( file_rows ):
    ''' 
    @description:
        parses file and MultiGraph object that represents the schematic
        Nodes are numbers of nets in the electrical schematic, 
        enumeration is kept the same as in the file
        Edges are electrical elements (capacitors, SIS josephson junctions etc.)
    @parameters:
        file_rows - list of lists containing csv file, e.g. list( csv.reader(file) ) 
    @return: MultiGraph object from networkx module, containing electrical schematic
    '''
    
    # getting reference designators and group names of the 
    # elements presented in file
    refDes_list, group_name_list = _get_refDes_groupNames_PADS( file_rows )        
    
    # contains pairs of nodes that corresponds to the elements
    # stored in refDes_list
    nets_list = [[] for x in refDes_list] 
    graph = nx.MultiGraph() # structure that will be returned by this function
    
    # filling the nets_list structure
    for i,row in enumerate( file_rows ):
        if( len(row) < 2 ):
            continue
            
        net_id = row[1].split('_') 
        
        # found row with net information (net is the node of the graph)
        if( net_id[0] == "Net" ): #checking if this is a net
            graph.add_node( int(net_id[1]) ) # adding net with corresponding number
            element_row = file_rows[i+1]
            for element in element_row[:-1]: # last element is empty due to DipTrace export...
                # parsing CSV row
                refDes,pinDesc = element.split('.')
                i = refDes_list.index( refDes )
                
                # in case that corrent element has positive and
                # negative terminals, the negative will always be
                # stored firstly in nets_list[i]
                if( pinDesc == "NEG" ):
                    nets_list[i].insert( 0,int(net_id[1]) )
                else:
                    nets_list[i].append( int(net_id[1]) )
    #nets_list is filled
    
    # iterating through nets_list and ref_des
    # and fulfilling MultiGraph structure that shall
    # be returned after
    for i,edge in enumerate(nets_list):
        # getting refDes of the current element
        refDes = refDes_list[i]
        element_arg = None
        
        # getting all the letters from refDes
        # they correspond to the element type 
        # element type part of refDes is predefined and unchanged
        # during the software developement process
        refName = re.compile("[A-Za-z]*").match(refDes)[0]
        
        # getting the user-defined group name
        group_name = group_name_list[i]
        
        # parsing the refDes element type and adding this to graph
        if( refName == "J" ): # Josephson junction
            element_arg = SIS_JJ(refDes,group_name,edge)
        elif( refName == "C" ): #capacitor
            element_arg = Cap(refDes,group_name, edge)
        elif( refName == "B" ): # battery
            element_arg = Battery(refDes,group_name,edge)
          
        graph.add_edge( *edge, element=element_arg )
        
    return graph

def build_graph_from_file_PADS( file_path ):
    with open( file_path, "r" ) as file:
        rows = list(csv.reader(file, delimiter=' '))
    return build_graph_from_netlist_PADS( rows )