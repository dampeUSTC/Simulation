#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#@
#@ Converts single the single solid, containting all the silicon detectors [output_file_name],
#@ into the true array of silicon detectors [output_file_name]  
#@ Usage:
#@   python conver_silicon.py
#@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#@ Options:
input_file_name   = "./converter_coordinates_tmp.xml"
#input_file_name   = "./stk_foils_tmp.gdml"
#input_file_name   = "./DAM_EQM_LshapeSIMPL.gdml"
#output_file_name  = "./stk_si_advanced.gdml"
#output_file_name  = "./stk_foils.gdml"
output_file_name  = "./converter_coordinates.xml"

#@
#@ Auxiliary gdml parsing functions   ---NOT USED YET!!!!
#@
def parce_vertices_for_gdml_file(filelines):
    vertices = {}
    for line in filelines:
        if "<position" not in line: continue
        name = line.split('name=')[1].split('"')[1]
        
        #print "line: ", line
        
        x = line.split('x=')[1].split('"')[1]
        #x = x.split(".")[0] + "." + x.split(".")[1][0:5]
        x = float(x)        
         
        y = line.split('y=')[1].split('"')[1]
        #y = y.split(".")[0] + "." + y.split(".")[1][0:5]
        y = float(y)
        
        z = line.split('z=')[1].split('"')[1]
        #z = z.split(".")[0] + "." + z.split(".")[1][0:5]
        z = float(z)
        
        if name in vertices: continue
        vertices[name] = [x,y,z] 
        
    return vertices



f = open(input_file_name,'r')
lines = f.readlines()
z = map(lambda x: x[2], parce_vertices_for_gdml_file(lines).values())
z = sorted(z)
all_z = []
for value in z:
    if value in all_z: continue
    all_z.append(value)


def create_output():
    f = open(output_file_name,'w')
    for line in lines:
        if 'z="' in line:
            z = float(line.split('z="')[1].split('"')[0])
            index = all_z.index(z)
            #if index < 0: raise Exception("z position not found in the vertex list")
            
            """
            #@ For foils
            index = 6 - index / 4
            z = z + 0.05935668945309658 - index * 0.25
            """
            
            """
            #@ For ladders
            index = 6 - (index +3) / 6
            z = z + 0.05935668945309658 - index * 0.25
            """
            
            """
            #@ For hcomb
            index = 6 - index  / 2
            z = z + 0.05935668945309658 - index * 0.25
            """
            
            """
            #@ For silicon
            index = 6 - (index +2) / 4
            z = z + 0.05935668945309658 - index * 0.25
            """
            
            #@ For tumgsten
            index = 6 - (index +2*3) / 2
            z = z + 0.05935668945309658 - index * 0.25
            
            line = line.split('z="')[0]+'z="'+str(z)+'"'+line.split('"')[-1]        
        f.write(line)
         
create_output()



        