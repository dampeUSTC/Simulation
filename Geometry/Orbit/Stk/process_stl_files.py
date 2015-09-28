import os

input  = "input.stl"
output = "output.gdml"

pathinput  = "/home/andrii/tmp/input"
pathoutput = "/home/andrii/tmp/output"
#inputfiles = ["AMS_LDR_1"]
inputfiles = map(lambda x: x.split(".stl")[0], os.listdir(pathinput))
#inputfiles = [C1_endcap_1]
cadmesh_command = '/home/andrii/cadmesh/cadmesh/examples/basic/build/cadmesh_example'



structure = '''      
  <structure>
    <volume name="%s">
    <materialref ref="%s"/>
    <solidref ref="%s"/>
  </volume>
'''
    
"""
materials_tag = '''
<!DOCTYPE gdml[
<!ENTITY materials SYSTEM "../Parameters/materials.xml">
]>
'''
"""
materials_tag = '''
<!DOCTYPE gdml[
<!ENTITY materials SYSTEM "../materials.xml">
]>
'''



#materialline ='      <materialref ref="%s"/>'

materials = {'AMS_LDR'    : 'Silicon',
             'B_field'    : 'Vacuum',
             'C1_endcap_' : 'Plexiglass',
             'C2_endcap_' : 'Plexiglass',
             'C1_env'     : 'Plexiglass',
             'C2_env'     : 'Plexiglass',
             'C1.'        : 'CO2',
             'C2.'        : 'CO2',
             'MWPC'       : 'ArCO2',
             'Pb_glass'   : 'LeadGlass',
             'S'          : 'PLasticScintillatorMaterialBeamTest',
             'Iron_Dump'  : 'Iron',
             'Concrete'   : 'Concrete'
             }  
             
def guess_material(volumename):
    for key in materials:
        if volumename.find(key) !=0: continue
        return materials[key]
    raise Exception("No material found for the volume: "+volumename)
             


def process_stl_files(guessMaterial=False):
    
    writeline = True
    for file in inputfiles:
        
        #volName = None      
       
        print "\n\n\n Processing file: %s \n\n\n"%file
    
        os.system("cp %s/%s.stl %s/%s"%(pathinput,file,pathinput,input))
        #os.system("./build/cadmesh_example")
        os.system(cadmesh_command)
        #os.system("mv %s/%s %s/%s.gdml"%(pathoutput,output,pathoutput,file))
        
        fin  = open(pathoutput+"/"+output)
        fout = open(pathoutput+"/"+file+".gdml","w")
        lines = fin.readlines()
        material = "Vacuum"
        if(guessMaterial): material = guess_material(file+".stl")
            
        for line in lines:        
            line = line.replace(pathinput+"/"+input,file)
            #line = line.replace('unit="mm"','unit="cm"')
            #line = line.replace("http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd","../GDMLSchema/gdml.xsd")
            line = line.replace("http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd","../../GDMLSchema/gdml.xsd")
    
            if '<tessellated' in line: solidname = line.split('name="')[-1].split('">')[0]
    
            if '<world' in line: line = '    <world ref="'+file+'"/>\n'
    
            if '<structure>' in line: writeline = False
            if '<structure>' in line: fout.write(structure%(file,material,solidname))
            if '</structure>' in line: writeline = True 
              
            if '<materials>' in line: writeline = False
            if '</materials>' in line: writeline = True; line = "  &materials;"#line =""      
            
            if '<box ' in line: line="" 
            
            if '<gdml' in line: fout.write(materials_tag) 
            
            #if '<volume name=' in line:
            #    if volName: raise Exception("Volume name already exists for this gdml file!")
            #    volName = line.split('<volume name="')[1].split('>')[0]    
                
            #if(guessMaterial): 
            #    #if '<materialref ref="' in line: line = materialline%guess_material(volName)
            #    if '<materialref ref="' in line: line = materialline%guess_material(file+".stl")
    
            if writeline: fout.write(line)    
        
        os.system("rm %s/%s"%(pathinput,input))  
        os.system("rm %s/%s"%(pathoutput,output))  


"""
bt_detector_gdml_tag = '''
    <physvol>
      <file name="./Geometry/BeamTestOctober2014/Components/%s.gdml"/> 
    </physvol>
'''
"""

bt_detector_gdml_tag = '''
    <physvol>
      <file name="./Geometry/BeamTestNov2014_SPS/Components/%s.gdml"/> 
    </physvol>
'''

def print_gdml_bt_detectors():
    for file in inputfiles:
        print bt_detector_gdml_tag%file
    
    

