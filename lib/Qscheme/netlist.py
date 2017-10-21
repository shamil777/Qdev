import re
import networkx as nx

class Element():
    def __init__(self, element_type, refDes, group_name=None ):
        if( group_name == None ):
            group_name = RefDes
        self.element_type = element_type
        self.group_name = group_name
        self.refDes = refDes
        
    def __eq__( self, Element2 ):
        if( self.RefDes == Element2.RefDes ):
            return True
        else:
            return False
        
    def __nq__( self, Element2 ):
        return (not self.__eq__( Element2 ))
        
class SIS_JJ( Element ):
    def __init__(self, refDes, group_name=None, capacitance=4.5e-15, E_J = 38e9 ):
        super( SIS_JJ,self ).__init__( "SIS_JJ", refDes, group_name )
        
        self.capacitance = capacitance
        self.E_J = E_J
        
class Cap( Element ):
    def __init__(self, refDes, group_name=None, capacitance=51e-15 ):
        super( Cap,self ).__init__( "SIS_JJ", refDes, group_name )
        
        self.capacitance = capacitance

class Battery( Element ):
    def __init__(self, refDes, group_name=None ):
        super( Battery,self ).__init__( "Battery", refDes, group_name )
		

def _get_elements_PADS( file_rows ):
    ''' 
    @description:
        parses file and return edge types, represented by electrical elements
    @parameters:
        file_rows - list of lists containing csv file (e.g. list( csv.reader(file)) 
    @return: list of Element() objects used in this netlist file
    '''
    elements = []
    
    element_i = -1
    for i,row in enumerate( file_rows ):
        if( row[0] == "*PART*" ):
            element_i = i+1
            break
    
    row = file_rows[element_i]
    while( len(row) > 0 ):
        refDes = row[0]
        refName = re.compile("[A-Za-z]*").match(refDes)[0]
        group_name = row[5].split(',')[1]
        #print( group_name )
        if( refName == "J" ): # Josephson junction
            elements.append( SIS_JJ(refDes,group_name) )
        elif( refName == "C" ): #capacitor
            elements.append( Cap( refDes,group_name) )
        elif( refName == "B" ): # battery
            elements.append( Battery(refDes,group_name) )
        element_i += 1
        row = file_rows[element_i]
        
    return elements
	
def build_graph_from_netlist_PAD( file_rows ):
    ''' 
    @description:
        parses file and MultiGraph object that represents the schematic
        Nodes are nets in electrical schematic, enumeration is kept the same as in the file
        Edges are electrical elements (capacitors, SIS josephson junctions etc.)
    @parameters:
        file_rows - list of lists containing csv file (e.g. list( csv.reader(file)) 
    @return: MultiGraph object from networkx module, containing electrical schematic
    '''
    elements = _get_elements_PADS( file_rows )
    graph = nx.MultiGraph()
    refDes_list = [x.refDes for x in elements]
    print(refDes_list)
    nets_list = [[] for x in elements]
    
    for i,row in enumerate( file_rows ):
        if( len(row) < 2 ):
            continue
            
        net_id = row[1].split('_') 
        
        # found row with net information (net == node of the graph)
        if( net_id[0] == "Net" ):
            graph.add_node( net_id[1] ) # adding net with corresponding number
            
            element_row = file_rows[i+1]
            for element in element_row[:-1]: # last element is empty due to DipTrace export...
                refDes = element.split('.')[0]
                i = refDes_list.index( refDes )
                nets_list[i].append( net_id[1] )
    print( nets_list )
    for i,edge in enumerate(nets_list):
        graph.add_edge( *edge, element=elements[i] )
        
    return graph