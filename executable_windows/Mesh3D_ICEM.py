# import gmsh
import numpy as np
import os



# Class capable of generating a 3D mesh suitable for axisymmetric BFM analysis in SU2.
class ICEM3D:
    wedge =         1    # 3D wedge angle in degrees. Should be lower than 180.
    n_sec =         1    # Number of sections in tangential direction in the wedge.
    n_point =       20   # Number of points in axial direction for each blade row.
    inlet_fac =     1.0  # Inlet cell size factor. High numbers mean a more coarse inlet grid.
    outlet_fac =    1.0  # Outlet cell size factor.
    fileName =      "3DBFM.su2"  # Mesh file name.
    pointID =       None  # List of point identification numbers.
    xCoords =       None  # List of x-coordinates of the mesh points.
    yCoords =       None  # List of y-coordinates of the mesh points.
    zCoords =       None  # List of z-coordinates of the mesh points.
    lC =            None  # List of mesh cell sizes.
    rot_axis =      [0, 0, 1]
    lineID =        None  # List of line identification numbers.
    lineStart =     None  # List of line start point identification numbers.
    lineEnd =       None  # List of line end point identification numbers.

    Rev =           None  # Revolution class.


    def __init__(self, Meangen, IN):
        # Storing Meangen class.
        self.M = Meangen

        # Reading mesh parameters from input file.

        #self.wedge = IN["WEDGE"][0]     # Importing mesh wedge angle
        self.wedge = 30
        self.n_sec = int(IN["SECTIONS_PER_DEGREE"][0]) * self.wedge     # Calculating number of tangential nodes
        self.n_point = int(IN["AXIAL_POINTS"][0])       # Importing axial node count
        self.inlet_fac = IN["INLET_FACTOR"][0]          # Importing inlet cell size factor
        self.outlet_fac = IN["OUTLET_FACTOR"][0]        # Importing outlet cell size factor
        self.rot_axis = IN["Rotation_axis"]
        self.BL_thick = IN["BOUNDARY_LAYER_THICKNESS"][0]

        self.Coords = type('', (), {})()
        # # Initiating Gmesh.
        # self.model = gmsh.model
        # self.factory = self.model.geo
        # gmsh.initialize()
        # gmsh.option.setNumber("General.Terminal", 1)
        # self.model.add("3DBFM")

        # Write the initial paragraph for the replay file
        f = open("ICEM_input.txt", "w")

        #set project file directory
        #f.write("ic_chdir C:/Users/GMrx1/Desktop/ANSYS_WORKS/lecture1\n")

        # set the settings
        f.write("ic_set_global geo_cad 0 toptol_userset\n")
        f.write("ic_set_global geo_cad 0.0 toler\n")
        f.write("ic_set_meshing_params global 0\n")
        f.write("ic_set_global geo_cad 1 toptol_userset\n")
        f.write("ic_set_meshing_params global 0 gttol 1E-9 gtrel 1\n")
        f.write("ic_regenerate_tris\n")
        f.write("ic_geo_set_units {}\n")
        f.write("ic_set_global geo_cad 1 toptol_userset\n")
        f.write("ic_set_meshing_params global 0 gttol 9.9999997e-10 gtrel 1\n")
        f.write("ic_geo_set_units {}\n")
        
        f.close()

        # Creating all the mesh points.
        self.makePoints()

        #  #Connecting the mesh points to create lines.
        self.makeLines()

        # Making lateral periodic faces.
        self.makeSymFaces()

        # Revolving hub and shroud around rotation axis to create walls.
        self.makeWalls()

        # Revolving radial lines around rotation axis to create inlet and outlet.
        self.makeInOut()

        # Setting names to boundaries.
        self.nameBoundaries()

        # Fix points that were overwritten by ICEM (GEOM problem is not fixable)
        self.fixPoints()

        # Make the blocks
        self.blocking()

        # Creating 3D mesh.
        self.mesh()
        

        # # Applying mesh refinement along the lines in the mesh
        # self.refineLines()

        # # Generate the 3D mesh geometry
        # gmsh.model.mesh.generate(3)

        # # Saving mesh as su2 file
        # gmsh.write(self.fileName)

        # Transforming mesh to periodic mesh.
        # print("Building periodic mesh...")
        # self.makePerio()
        # print("Done!")
        # HOME = os.environ["M2BFM"]
        ICEMDIR = os.environ["ICEMDIR"]
        os.system("copy ICEM_input.txt " + ICEMDIR +"\\ICEM_input.rpl")

        

    #     # In case of selection of mesh visualization option, the gmesh GUI will display the 3D BFM mesh.
    #     if IN["PLOT_MESH"] == 'YES':
    #         self.plotmesh()
    #     gmsh.finalize()

    # def plotmesh(self):
    #     # This function plots the mesh in the GMesh GUI. It highlights the mesh boundaries and gives them a distinct
    #     # color.

    #     # Allow for addition of text comments in the GMesh GUI
    #     v = gmsh.view.add("comments")

    #     # Calculating the average coordinates of the mesh boundaries. This puts the boundary tag in the middle of the
    #     # boundary patch.
    #     inlet_coords = np.sum(self.Coords.Coords_inlet, axis=0)/np.shape(self.Coords.Coords_inlet)[0]
    #     outlet_coords = np.sum(self.Coords.Coords_outlet, axis=0) / np.shape(self.Coords.Coords_outlet)[0]
    #     hub_coords = np.sum(self.Coords.Coords_hub, axis=0) / np.shape(self.Coords.Coords_hub)[0]
    #     shroud_coords = np.sum(self.Coords.Coords_shroud, axis=0) / np.shape(self.Coords.Coords_shroud)[0]
    #     perio_1_coords = np.sum(self.Coords.Coords_periodic_1, axis=0) / np.shape(self.Coords.Coords_periodic_1)[0]
    #     perio_2_coords = np.sum(self.Coords.Coords_periodic_2, axis=0) / np.shape(self.Coords.Coords_periodic_2)[0]

    #     # Adding the boundary patch tags to each of the boundary patches.
    #     gmsh.view.addListDataString(v, [i for i in inlet_coords], ["inlet"], ["Align", "Center", "Font", "Helvetica"])
    #     gmsh.view.addListDataString(v, [i for i in outlet_coords], ["outlet"], ["Align", "Center", "Font", "Helvetica"])
    #     gmsh.view.addListDataString(v, [i for i in hub_coords], ["hub"], ["Align", "Center", "Font", "Helvetica"])
    #     gmsh.view.addListDataString(v, [i for i in shroud_coords], ["shroud"], ["Align", "Center", "Font", "Helvetica"])
    #     gmsh.view.addListDataString(v, [i for i in perio_2_coords], ["periodic_2"],
    #                                 ["Align", "Center", "Font", "Helvetica"])
    #     gmsh.view.addListDataString(v, [i for i in perio_1_coords], ["periodic_1"], ["Align", "Center", "Font", "Helvetica"])

    #     # Running the GMesh GUI
    #     gmsh.fltk.run()

    # def makePerio(self):
    #     # This function modifies the mesh to be periodic through SU2_PERIO.
    #     # Creating SU2_PERIO configuration file.
    #     file = open("createPerio.cfg", "w+")

    #     # Writing periodic boundary command and specifying wedge angle.
    #     file.write("MARKER_PERIODIC= (periodic_1, periodic_2, 0.0, 0.0, 0.0, 0.0, 0.0, " + str(
    #         float(self.wedge)) + ", 0.0, 0.0, 0.0)\n")
    #     file.write("MESH_FILENAME= "+self.fileName+"\n")
    #     file.write("MESH_FORMAT= SU2\n")
    #     file.write("MESH_OUT_FILENAME= "+self.fileName[:-4]+"_perio.su2\n")
    #     file.close()

    #     # Executing SU2_PERIO to create periodic mesh and storing output in output file.
    #     os.system("SU2_PERIO createPerio.cfg > SU2_PERIO.out")


    def mesh(self):
         f = open("ICEM_input.txt", "a")

         # Seed the edges in axial direction
         f.write("ic_hex_set_mesh 37 70 n "+str(self.n_point)+" h1rel "+str(0.05)+" h2rel "+str(self.BL_thick/self.length_hub[0])+" r1 1.2 r2 1.2 lmax 0 exp2 copy_to_parallel unlocked\n")
         f.write("ic_hex_set_mesh 70 86 n "+str(self.n_point)+" h1rel "+str(self.BL_thick/self.length_hub[1])+" h2rel "+str(self.BL_thick/self.length_hub[1])+" r1 1.2 r2 1.2 lmax 0 biexponential copy_to_parallel unlocked\n")
         f.write("ic_hex_set_mesh 86 102 n "+str(self.n_point)+" h1rel "+str(self.BL_thick/self.length_hub[2])+" h2rel "+str(self.BL_thick/self.length_hub[2])+" r1 1.2 r2 1.2 lmax 0 uniform copy_to_parallel unlocked\n")
         f.write("ic_hex_set_mesh 102 118 n "+str(self.n_point)+" h1rel "+str(self.BL_thick/self.length_hub[3])+" h2rel "+str(self.BL_thick/self.length_hub[3])+" r1 1.2 r2 1.2 lmax 0 biexponential copy_to_parallel unlocked\n")
         f.write("ic_hex_set_mesh 118 38 n "+str(self.n_point)+" h1rel "+str(self.BL_thick/self.length_hub[4])+" h2rel "+str(0.05)+" r1 1.2 r2 1.2 lmax 0 exp1 copy_to_parallel unlocked\n")
        
         # Seed the edges in radial direction
         n_radial = round(self.n_point*self.length_rad[0]/self.length_hub[0])
         f.write("ic_hex_set_mesh 37 41 n "+str(n_radial)+" h1rel "+str(0.01)+" h2rel "+str(0.01)+" r1 1.2 r2 1.2 lmax 0 uniform copy_to_parallel unlocked\n")


         f.write("ic_hex_create_mesh GEOM RIGHT_SYM LEFT_SYM HUB_WALL SHROUD_WALL INLET OUTLET SOLID proj 2 dim_to_mesh 3\n")
         f.close()






    def blocking(self):

        f = open("ICEM_input.txt", "a")

         # # Create the initial big block
        f.write("ic_geo_new_family SOLID\n")
        f.write("ic_boco_set_part_color SOLID\n")
        loop = "{"
        for i in self.points_hub:
            loop = loop + " point pnt." + str(i)
        for i in self.points_shroud:
            loop = loop + " point pnt." + str(i)
        for i in self.points_shroud2:
            loop = loop + " point pnt." + str(i)
        for i in self.points_hub2:
            loop = loop + " point pnt." + str(i)
        f.write("ic_hex_initialize_blocking " + loop + "} SOLID 0 101\n")

        # f.write("ic_hex_initialize_blocking {point pnt." +str(self.points_hub[0])+ " point pnt." +str(self.points_hub[-1])+
        #         " point pnt." +str(self.points_shroud[0])+ " point pnt." +str(self.points_shroud[-1])+
        #         " point pnt." +str(self.points_shroud2[0])+ " point pnt." +str(self.points_shroud2[-1])+
        #         " point pnt." +str(self.points_hub2[0])+ " point pnt." +str(self.points_hub2[-1])+"} SOLID 0 101\n")
        f.write("ic_hex_unblank_blocks \n")
        f.write("ic_hex_multi_grid_level 0\n")
        f.write("ic_hex_projection_limit 0\n")
        f.write("ic_hex_default_bunching_law default 2.0\n")
        f.write("ic_hex_floating_grid off\n")
        f.write("ic_hex_transfinite_degree 1\n")
        f.write("ic_hex_unstruct_face_type one_tri\n")
        f.write("ic_hex_set_unstruct_face_method uniform_quad\n")
        f.write("ic_hex_set_n_tetra_smoothing_steps 20\n")
        f.write("ic_hex_error_messages off_minor\n")

        # Split block
        f.write("ic_hex_split_grid 41 42 pnt."+str(self.points_shroud[1])+" m GEOM RIGHT_SYM LEFT_SYM HUB_WALL SHROUD_WALL INLET OUTLET SOLID\n")
        f.write("ic_hex_split_grid 74 42 pnt."+str(self.points_shroud[2])+" m GEOM RIGHT_SYM LEFT_SYM HUB_WALL SHROUD_WALL INLET OUTLET SOLID\n")
        f.write("ic_hex_split_grid 90 42 pnt."+str(self.points_shroud[3])+" m GEOM RIGHT_SYM LEFT_SYM HUB_WALL SHROUD_WALL INLET OUTLET SOLID\n")
        f.write("ic_hex_split_grid 106 42 pnt."+str(self.points_shroud[4])+" m GEOM RIGHT_SYM LEFT_SYM HUB_WALL SHROUD_WALL INLET OUTLET SOLID\n")
  
        # Assign nodes to vertices
        nodes_hub = [37,70,86,102,118,38]
        nodes_shroud = [41,74,90,106,122,42]
        nodes_hub2 = [21,69,85,101,117,22]
        nodes_shroud2 = [25,73,89,105,121,26] 
        for i in range(len(self.points_hub)):
            f.write("ic_hex_move_node "+str(nodes_hub[i])+" pnt." +str(self.points_hub[i])+"\n")
            f.write("ic_hex_move_node "+str(nodes_shroud[i])+" pnt." +str(self.points_shroud[i])+"\n")
            f.write("ic_hex_move_node "+str(nodes_hub2[i])+" pnt." +str(self.points_hub2[i])+"\n")
            f.write("ic_hex_move_node "+str(nodes_shroud2[i])+" pnt." +str(self.points_shroud2[i])+"\n")

        # Assign edges to curves
        for i in range(len(self.points_hub)):
            f.write("ic_hex_set_edge_projection "+str(nodes_hub[i])+" "+str(nodes_hub2[i])+" 0 1 crv."+str(self.lines_wall_hub[i])+"\n")
            f.write("ic_hex_set_edge_projection "+str(nodes_shroud[i])+" "+str(nodes_shroud2[i])+" 0 1 crv."+str(self.lines_wall_shroud[i])+"\n")
        
       
        # O-Grid split
        blocks = [13,27,28,29,30]
        for i in blocks:
            f.write("ic_hex_mark_blocks superblock "+str(i)+"\n")
        loop =""
        for i in range(len(nodes_hub)):
            loop =loop+"{"+str(nodes_hub2[i])+" "+str(nodes_hub[i])+" "+str(nodes_shroud2[i])+" "+str(nodes_shroud[i])+"} "
        for i in range(len(nodes_shroud)-1):
            loop =loop+"{"+str(nodes_shroud2[i])+" "+str(nodes_shroud2[i+1])+" "+str(nodes_shroud[i])+" "+str(nodes_shroud[i+1])+"} "
        f.write("ic_hex_mark_blocks face_neighbors corners "+loop+"\n")
        f.write("ic_hex_ogrid 1 m GEOM RIGHT_SYM LEFT_SYM HUB_WALL SHROUD_WALL INLET OUTLET SOLID -version 50\n")
        f.write("ic_hex_mark_blocks unmark\n")
        
        # Assign splitted nodes to vertices
        nodes_hub3 = [153,154,155,156,157,158]
        nodes_shroud3 = [161,162,163,164,165,166]
        nodes_hub23 = [129,130,131,132,133,134]
        nodes_shroud23 = [137,138,139,140,141,142]
        for i in range(len(self.points_hub)):
            f.write("ic_hex_move_node "+str(nodes_hub3[i])+" pnt." +str(self.points_lines_wall_hub3[i])+"\n")
            f.write("ic_hex_move_node "+str(nodes_shroud3[i])+" pnt." +str(self.points_lines_wall_shroud3[i])+"\n")
            f.write("ic_hex_move_node "+str(nodes_hub23[i])+" pnt." +str(self.points_lines_wall_hub23[i])+"\n")
            f.write("ic_hex_move_node "+str(nodes_shroud23[i])+" pnt." +str(self.points_lines_wall_shroud23[i])+"\n")

         # Assign edges to curves
        for i in range(len(self.points_hub)):
            f.write("ic_hex_set_edge_projection "+str(nodes_hub[i])+" "+str(nodes_hub3[i])+" 0 1 crv."+str(self.lines_wall_hub[i])+"\n")
            f.write("ic_hex_set_edge_projection "+str(nodes_hub3[i])+" "+str(nodes_hub23[i])+" 0 1 crv."+str(self.lines_wall_hub[i])+"\n")
            f.write("ic_hex_set_edge_projection "+str(nodes_hub23[i])+" "+str(nodes_hub2[i])+" 0 1 crv."+str(self.lines_wall_hub[i])+"\n")
            
            f.write("ic_hex_set_edge_projection "+str(nodes_shroud[i])+" "+str(nodes_shroud3[i])+" 0 1 crv."+str(self.lines_wall_shroud[i])+"\n")
            f.write("ic_hex_set_edge_projection "+str(nodes_shroud3[i])+" "+str(nodes_shroud23[i])+" 0 1 crv."+str(self.lines_wall_shroud[i])+"\n")
            f.write("ic_hex_set_edge_projection "+str(nodes_shroud23[i])+" "+str(nodes_shroud2[i])+" 0 1 crv."+str(self.lines_wall_shroud[i])+"\n")
        
        # Remove bottom blocks that are not needed
        low_blocks=[36,39,42,45,48]
        for i in low_blocks:
            f.write("ic_hex_mark_blocks superblock "+str(i)+"\n")
            f.write("ic_hex_change_element_id VORFN\n")
        
        f.close()






    def fixPoints(self):
        f = open("ICEM_input.txt", "a")
        
        # Delete the reference points of the rotation axis
        for i in self.points_rot_axis:
             f.write("ic_geo_incident point pnt." + str(i) + " 1\n")
             f.write("ic_delete_geometry point names pnt." + str(i) + " 0 1\n")
             f.write("ic_set_dormant_pickable point 0 {}\n")

    #   # Rewrite the 4 external vertex of the wedge (ICEM problem: it overwite stuff when creating revolution surfaces)
        
        f.write("ic_point {} GEOM pnt."+str(self.points_hub[0])+" "+str(self.coords_hub [0,0])+","+str(self.coords_hub [1,0])+","+str(self.coords_hub [2,0])+"\n")
        f.write("ic_point {} GEOM pnt."+str(self.points_hub[-1])+" "+str(self.coords_hub [0,-1])+","+str(self.coords_hub [1,-1])+","+str(self.coords_hub [2,-1])+"\n")
        f.write("ic_point {} GEOM pnt."+str(self.points_shroud[0])+" "+str(self.coords_shroud [0,0])+","+str(self.coords_shroud [1,0])+","+str(self.coords_shroud [2,0])+"\n")
        f.write("ic_point {} GEOM pnt."+str(self.points_shroud[-1])+" "+str(self.coords_shroud [0,-1])+","+str(self.coords_shroud [1,-1])+","+str(self.coords_shroud [2,-1])+"\n")
        


        # Rotate the reference points for future blocking set
        lines_wall_hub = self.lines_wall_hub
        lines_wall_shroud = self.lines_wall_shroud
        points_count = self.points_count
        points_count = self.points_count
        points_lines_wall_hub3 = []
        points_lines_wall_hub23 = []
        points_lines_wall_shroud3 = []
        points_lines_wall_shroud23 = []

        f.write("ic_set_global geo_cad 0.0006 toler\n")
       
        for i in range(len(lines_wall_hub)):
            points_count += 1
            f.write("ic_geo_duplicate_set_fam_and_data point pnt." + str(self.points_hub[i]) + " pnt." + str(points_count) + " {} _0\n")
            points_lines_wall_hub3.append(points_count)
            points_count += 1
            f.write("ic_geo_duplicate_set_fam_and_data point pnt." + str(self.points_hub[i]) + " pnt." + str(points_count) + " {} _0\n")
            points_lines_wall_hub23.append(points_count)
        
        for i in range(len(lines_wall_shroud)):
            points_count += 1
            f.write("ic_geo_duplicate_set_fam_and_data point pnt." + str(self.points_shroud[i]) + " pnt." + str(points_count) + " {} _0\n")
            points_lines_wall_shroud3.append(points_count)
            points_count += 1
            f.write("ic_geo_duplicate_set_fam_and_data point pnt." + str(self.points_shroud[i]) + " pnt." + str(points_count) + " {} _0\n")
            points_lines_wall_shroud23.append(points_count)
            
        loop3 = "{"     
        for i in range(len(lines_wall_hub)):
            loop3 = loop3 + " pnt."+ str( points_lines_wall_hub3[i])
        for i in range(len(lines_wall_shroud)):
            loop3 = loop3 + " pnt."+str(points_lines_wall_shroud3[i])
        loop3 = loop3 + " }"

        loop23 = "{"     
        for i in points_lines_wall_hub23:
            loop23 = loop23 + " pnt."+str(i)
        for i in points_lines_wall_shroud23:
            loop23 = loop23 + " pnt."+str(i) 
        loop23 = loop23 + " }"


        f.write("ic_move_geometry point names " + loop3 + " rotate " + str(self.wedge/3) + " rotate_axis {0 0 1} cent {0 0 0}\n")
        f.write("ic_move_geometry point names " + loop23 + " rotate " + str(self.wedge*2/3) + " rotate_axis {0 0 1} cent {0 0 0}\n")
        f.write("ic_geo_reset_data_structures\n")
        f.write("ic_geo_configure_one_attribute surface shade wire\n")

        # # Storing the reference points identifiers into the class
        self.points_lines_wall_hub3 = points_lines_wall_hub3
        self.points_lines_wall_hub23 = points_lines_wall_hub23
        self.points_lines_wall_shroud3 = points_lines_wall_shroud3
        self.points_lines_wall_shroud23 = points_lines_wall_shroud23
    #     f.write("ic_set_global geo_cad 0.0006 toler\n")
    #     f.write("ic_geo_duplicate_set_fam_and_data point pnt." + str(self.points_hub[0]) + " pnt." + str(self.points_hub2[0]) + " {} _0\n")
    #     f.write("ic_geo_duplicate_set_fam_and_data point pnt." + str(self.points_hub[-1]) + " pnt." + str(self.points_hub2[-1]) + " {} _0\n")
    #     f.write("ic_geo_duplicate_set_fam_and_data point pnt." + str(self.points_shroud[0]) + " pnt." + str(self.points_shroud2[0]) + " {} _0\n")
    #     f.write("ic_geo_duplicate_set_fam_and_data point pnt." + str(self.points_shroud[-1]) + " pnt." + str(self.points_shroud2[-1]) + " {} _0\n")
        
    #     loop = "{pnt."+str(self.points_hub2[0])+" pnt." + str(self.points_hub2[-1]) +" pnt." + str(self.points_shroud2[0]) +" pnt." + str(self.points_shroud2[-1]) +"}"
    #     f.write("ic_move_geometry point names " + loop + " rotate " + str(self.wedge) + " rotate_axis {0 0 1} cent {0 0 0}\n")
    #     f.write("ic_geo_reset_data_structures\n")
    #     f.write("ic_geo_configure_one_attribute surface shade wire\n")

        f.close()





   


    def nameBoundaries(self):
        # This function gives names to all the boundaries of the 3D mesh so boundary conditions can be assigned in the
        # SU2 configuration file.
        f = open("ICEM_input.txt", "a")

        # Create left and right symmetric surfaces
        f.write("ic_geo_set_part surface srf." +str(self.faces_simmetry[0])+ " RIGHT_SYM 0\n")
        f.write("ic_delete_empty_parts\n")
        f.write("ic_geo_set_part surface srf." +str(self.faces_simmetry[-1])+ " LEFT_SYM 0\n")
        f.write("ic_delete_empty_parts\n")

        # Create walls surfaces
        f.write("ic_geo_set_part surface srf." +str(self.walls[0])+ " HUB_WALL 0\n")
        f.write("ic_delete_empty_parts\n")
        f.write("ic_geo_set_part surface srf." +str(self.walls[-1])+ " SHROUD_WALL 0\n")
        f.write("ic_delete_empty_parts\n")

        # Create inlet and outlet
        f.write("ic_geo_set_part surface srf." +str(self.inlet[0])+ " INLET 0\n")
        f.write("ic_delete_empty_parts\n")
        f.write("ic_set_family_color_for_name INLET #00ff00\n")
        f.write("ic_geo_set_part surface srf." +str(self.outlet[0])+ " OUTLET 0\n")
        f.write("ic_delete_empty_parts\n")


        f.close()

    




    

    def makeInOut(self):
        # This function revolves the radial curves around the rotation axis to create the inlet and outlet.
        # The center point of the revolution is located at the origin and the rotation axis is set to be the z-axis.
        inlet = []
        outlet = []
        i_surf = self.surfaces_count
        f = open("ICEM_input.txt", "a")

        line_inlet= str(self.lines_rad[-1])
        f.write("ic_set_global geo_cad 0.0006 toler\n")
        f.write("ic_geo_cre_srf_rev GEOM srf." + str(i_surf+1) + " crv." + line_inlet + " pnt." + str(self.points_rot_axis[0]) + " {0 0 1} 0 " + str(-self.wedge) + " c 1\n")
        f.write("ic_set_global geo_cad 0.001 toler\n")
        f.write("ic_set_dormant_pickable point 0 {}\n")
        f.write("ic_set_dormant_pickable curve 0 {}\n")
        inlet.append(i_surf+1)

        line_outlet= str(self.lines_rad[0])
        f.write("ic_set_global geo_cad 0.0006 toler\n")
        f.write("ic_geo_cre_srf_rev GEOM srf." + str(i_surf+2) + " crv." + line_outlet + " pnt." + str(self.points_rot_axis[0]) + " {0 0 1} 0 " + str(-self.wedge) + " c 1\n")
        f.write("ic_set_global geo_cad 0.001 toler\n")
        f.write("ic_set_dormant_pickable point 0 {}\n")
        f.write("ic_set_dormant_pickable curve 0 {}\n")
        outlet.append(i_surf+2)

        f.close()

        # Storing the inlet and outlet surface identifiers into the class.
        self.inlet = inlet
        self.outlet = outlet
        self.surfaces_count += (len(inlet) + len(outlet))






    def makeWalls(self):
        # This function revolves the hub and shroud curves around the rotation axis to create the solid walls.
        # The center point of the revolution is located at the origin and the rotation axis is set to be the z-axis.
        
        # Create the reference points for the rotation axis
        points_rot_axis = []
        points_count =self.points_count

        f = open("ICEM_input.txt", "a")
        f.write("ic_set_global geo_cad 0.0002 toler\n")
        for i in range(len(self.points_hub)):
            points_count += 1
            f.write("ic_point {} GEOM pnt." + str(points_count) + " 0,0," + str(self.coords_hub[2,i]) + "\n")
            points_rot_axis.append(points_count)
       
        # Revolve the hub and shroud curves 
        walls = []
        i_surf = self.surfaces_count

        line_hub2= str(self.lines_hub2[0])
        f.write("ic_set_global geo_cad 0.0006 toler\n")
        f.write("ic_geo_cre_srf_rev GEOM srf." + str(i_surf+1) + " crv." + line_hub2 + " pnt." + str(points_count-1) + " {0 0 1} 0 " + str(self.wedge) + " c 1\n")
        f.write("ic_set_global geo_cad 0.001 toler\n")
        f.write("ic_set_dormant_pickable point 0 {}\n")
        f.write("ic_set_dormant_pickable curve 0 {}\n")
        walls.append(i_surf+1)
        
        line_shroud2= str(self.lines_shroud2[0])
        f.write("ic_set_global geo_cad 0.0006 toler\n")
        f.write("ic_geo_cre_srf_rev GEOM srf." + str(i_surf+2) + " crv." + line_shroud2 + " pnt." + str(points_count-1) + " {0 0 1} 0 " + str(self.wedge) + " c 1\n")
        f.write("ic_set_global geo_cad 0.001 toler\n")
        f.write("ic_set_dormant_pickable point 0 {}\n")
        f.write("ic_set_dormant_pickable curve 0 {}\n")
        walls.append(i_surf+2)

        # Create wall lines as a reference for blocks edges association. Also need 2 point per line
        lines_count = self.lines_count
        lines_wall_hub = []
        for i in range(len(points_rot_axis)):
            lines_count += 1 
            lines_wall_hub.append(lines_count)
            f.write("ic_curve arc_ctr_rad GEOM crv."+str(lines_count)+" {pnt."+str(points_rot_axis[i])+" pnt."+str(self.points_hub[i])+" pnt."+str(self.points_hub2[i])+" 0.0 {} {} 0}\n")
        lines_wall_shroud = []
        for i in range(len(points_rot_axis)):
            lines_count += 1 
            lines_wall_shroud.append(lines_count)
            f.write("ic_curve arc_ctr_rad GEOM crv."+str(lines_count)+" {pnt."+str(points_rot_axis[i])+" pnt."+str(self.points_shroud[i])+" pnt."+str(self.points_shroud2[i])+" 0.0 {} {} 0}\n")
        
        


        f.close()


        # Storing rotational axis into the class
        self.points_rot_axis = points_rot_axis
        self.points_count = points_count
        
        # Storing the hub and shroud wall surface identifiers into the class.
        self.lines_wall_hub = lines_wall_hub
        self.lines_wall_shroud = lines_wall_shroud
        self.lines_count = lines_count

        # Storing the hub and shroud wall surface identifiers into the class.
        self.walls = walls
        self.surfaces_count += len(walls)
        


    





    def makeSymFaces(self):
        # This function takes the lines defining the hub, shroud, inlet and outlet and builds a plane surface bound by
        # the meridional shape.
        # Then the other face with symmetric B.C. (wedge sides) is created by rotating the previous one around z axis.
        # Last step is to rotate the reference points and radial lines from the original surface onto the new one.

        faces_simmetry =[]
        i_face = 1 #Initialize surface count

        # Defining a list containing the line identifiers to make a GMesh curve loop.
        loop = "{crv."+str(self.lines_rad[0])
        for i in self.lines_hub:
            loop = loop + " crv." + str(i)
        loop = loop + " crv."+str(self.lines_rad[-1])
        for i in self.lines_shroud[::-1]:
            loop = loop + " crv." + str(i) 
        loop = loop + "}"

        # Create surface from n lines with 0.01 tolerance
        f = open("ICEM_input.txt", "a")
        f.write("ic_set_global geo_cad 0 toptol_userset\n")
        f.write("ic_set_global geo_cad 0.0002 toler\n")
        f.write("ic_surface bsinterp GEOM srf." + str(i_face) +" " + loop + "\n")
        f.write("ic_set_global geo_cad 0.0002 toler\n")
        f.write("ic_set_dormant_pickable point 0 {}\n")
        f.write("ic_set_dormant_pickable curve 0 {}\n")
        
        faces_simmetry.append(i_face)
            
        f.write("ic_set_global geo_cad 0.0002 toler\n")
        f.write("ic_geo_duplicate_set_fam_and_data surface srf." +str(i_face) + " srf." + str(i_face+1) + " {} _0\n")
        f.write("ic_move_geometry surface names srf." + str(i_face+1) + " rotate " + str(self.wedge) + " rotate_axis {0 0 1} cent {0 0 0}\n")
        f.write("ic_geo_reset_data_structures\n")
        f.write("ic_geo_configure_one_attribute surface shade wire\n")

        i_face += 1
        faces_simmetry.append(i_face)
        

        # Rotate the curves for future hub/shroud revolution
        lines_count = self.lines_count + 2            #Trick to avoid problems in ICEM: cannot recombine a list of lines into one with a name present in the list
        lines_hub2 = []
        lines_shroud2 = []

        loop_hub = "{"
        loop_shroud = ""
        f.write("ic_set_global geo_cad 0.0006 toler\n")
        # Copy hub lines to hub2
        for i in (self.lines_hub):
            lines_count += 1
            f.write("ic_geo_duplicate_set_fam_and_data curve crv." + str(i) + " crv." + str(lines_count) + " {} _0\n")
            loop_hub = loop_hub + "crv." + str(lines_count) + " "
            lines_hub2.append(lines_count)

        # Copy shroud lines to shroud2
        for i in (self.lines_shroud):
            lines_count += 1
            f.write("ic_geo_duplicate_set_fam_and_data curve crv." + str(i) + " crv." + str(lines_count) + " {} _0\n")
            loop_shroud = loop_shroud + "crv." + str(lines_count) + " "
            lines_shroud2.append(lines_count)

        loop = loop_hub + loop_shroud + "}"
        f.write("ic_move_geometry curve names " + loop + " rotate " + str(self.wedge) + " rotate_axis {0 0 1} cent {0 0 0}\n")
        f.write("ic_geo_reset_data_structures\n")
        f.write("ic_geo_configure_one_attribute surface shade wire\n")

        # Recombine lines
        f.write("ic_set_global geo_cad 0.0006 toler\n")
        f.write("ic_curve concat GEOM crv." + str(self.lines_count+1) + " " + loop_hub + "}\n")  # Recombine hub
        f.write("ic_curve concat GEOM crv." + str(self.lines_count+2) + " {" + loop_shroud + "}\n")  # Recombine shroud
        lines_hub2 = [self.lines_count+1]
        lines_shroud2 = [self.lines_count+2]


        # Rotate the reference points for future blocking set
        points_count = self.points_count
        points_hub2 = []
        points_shroud2 = []

        loop = "{"
        f.write("ic_set_global geo_cad 0.0006 toler\n")
        # Copy hub points to hub2
        for i in (self.points_hub):
            # i_point = i + point_count
            points_count += 1
            f.write("ic_geo_duplicate_set_fam_and_data point pnt." + str(i) + " pnt." + str(points_count) + " {} _0\n")
            loop = loop + " pnt." + str(points_count)
            points_hub2.append(points_count)
        # Copy shroud points to shroud2
        for i in (self.points_shroud):
            points_count += 1
            f.write("ic_geo_duplicate_set_fam_and_data point pnt." + str(i) + " pnt." + str(points_count) + " {} _0\n")
            loop = loop + " pnt." + str(points_count)
            points_shroud2.append(points_count)
            
        loop = loop + " }"
        f.write("ic_move_geometry point names " + loop + " rotate " + str(self.wedge) + " rotate_axis {0 0 1} cent {0 0 0}\n")
        f.write("ic_geo_reset_data_structures\n")
        f.write("ic_geo_configure_one_attribute surface shade wire\n")

        f.close()


        # Storing the hub2 and shroud2 point identifiers into the class.
        self.points_hub2 = points_hub2
        self.points_shroud2 = points_shroud2
        self.points_count = points_count

        # Storing the hub2 and shroud2 line identifiers into the class.
        self.lines_hub2 = lines_hub2
        self.lines_shroud2 = lines_shroud2
        self.lines_count = self.lines_count+2

        # Storing the surface identifiers into the class.
        self.faces_simmetry = faces_simmetry
        self.surfaces_count = len(faces_simmetry)










    def makeLines(self):
        # This function connects the points defined in makePoints to form lines.

        # Empty lists for the line identifiers.
        lines_hub = []
        length_hub = []
        lines_shroud = []
        length_shroud = []
        lines_rad = []
        length_rad = []
        i_line = 1
        # i_point = 1
        f = open("ICEM_input.txt", "a")

        # Looping through the hub and shroud points to build lines defining the hub and shroud shape.
        for i in range(len(self.points_hub)-1):
            # self.factory.addLine(self.points_hub[i], self.points_hub[i+1], i_line)
            f.write("ic_delete_geometry curve names crv."+str(i_line)+" 0\n")
            f.write("ic_curve point GEOM crv."+str(i_line)+" {pnt."+str(self.points_hub[i])+" pnt."+str(self.points_hub[i+1])+"}\n")
            lines_hub.append(i_line)
            length_hub.append(np.sqrt(np.sum(np.array(self.coords_hub[:, i+1] - self.coords_hub[:, i])**2)))
            i_line += 1
            # i_point += 1
        for i in range(len(self.points_shroud)-1):
            # self.factory.addLine(self.points_shroud[i], self.points_shroud[i+1], i_line)
            f.write("ic_delete_geometry curve names crv."+str(i_line)+" 0\n")
            f.write("ic_curve point GEOM crv."+str(i_line)+" {pnt."+str(self.points_shroud[i])+" pnt."+str(self.points_shroud[i+1])+"}\n")
            lines_shroud.append(i_line)
            length_shroud.append(np.sqrt(np.sum(np.array(self.coords_shroud[:, i + 1] - self.coords_shroud[:, i]) ** 2)))
            i_line += 1

        # Looping through the points in axial direction to build lines in radial direction. These will be used for the
        # definition of the inlet and outlet and refinement along the blade leading and trailing edges.
        for i in range(len(self.points_hub)):
            # self.factory.addLine(self.points_hub[i], self.points_shroud[i], i_line)
            f.write("ic_delete_geometry curve names crv."+str(i_line)+" 0\n")
            f.write("ic_curve point GEOM crv."+str(i_line)+" {pnt."+str(self.points_hub[i])+" pnt."+str(self.points_shroud[i])+"}\n")
            lines_rad.append(i_line)
            length_rad.append(np.sqrt(np.sum(np.array(self.coords_shroud[:, i] - self.coords_hub[:, i]) ** 2)))
            i_line += 1

        f.close()

        
        # Storing the line identifiers and line lengths into the class.
        self.lines_hub = lines_hub
        self.lines_shroud = lines_shroud
        self.lines_rad = lines_rad
        self.lines_count = i_line-1

        self.length_hub = length_hub
        self.length_shroud = length_shroud
        self.length_rad = length_rad












    def makePoints(self):
        # This function defines the coordinates for the points defining the annulus shape and registers them in GMesh

        # Extracting the leading and trailing edge coordinates from Meangen.
        X_LE = self.M.X_LE
        R_LE = self.M.Z_LE
        X_TE = self.M.X_TE
        R_TE = self.M.Z_TE

        # Calculating number of rows.
        n_rows = len(X_LE[0, :])

        # The inlet and outlet patches are placed two axial chords from the first and last blade row respectively.
        # The axial coordinates of the inlet and outlet are calculated here.
        x_min = min(X_LE[:, 0] - 2*(X_TE[:, 0] - X_LE[:, 0]))
        x_max = max(X_TE[:, -1] + 2*(X_TE[:, -1] - X_LE[:, -1]))

        # Each point has an identifier. This number is progressively updated with each point and stored in the
        # respective lists for the hub and shroud patches.
        i_point = 1     # Point identifier number starts at 1

        # Empty lists for the hub and shroud coordinates
        Z_hub = []
        X_hub = []
        Y_hub = []
        Z_shroud = []
        X_shroud = []
        Y_shroud = []

        # Empty lists for the hub and shroud point identifiers
        points_hub = []
        points_shroud = []

        # The mesh wedge is oriented symmetrically around the Z-Y plane. It therefore rotates with half the wedge angle
        # in positive and negative direction around the Z-axis to set the periodic boundary patches.
        theta = 0.5*self.wedge*np.pi/180    # Converting the wedge angle to radians

        # Defining the inlet point at the hub section
        Z_hub.append(x_min)
        X_hub.append(R_LE[0, 0]*np.sin(theta))
        Y_hub.append(R_LE[0, 0] * np.cos(theta))
        points_hub.append(i_point)
        i_point += 1

        # Defining the points defining the hub shape between the inlet and outlet
        for i in range(n_rows):
            Z_hub.append(X_LE[0, i])
            X_hub.append(R_LE[0, i] * np.sin(theta))
            Y_hub.append(R_LE[0, i] * np.cos(theta))
            points_hub.append(i_point)
            i_point += 1

            Z_hub.append(X_TE[0, i])
            X_hub.append(R_TE[0, i] * np.sin(theta))
            Y_hub.append(R_TE[0, i] * np.cos(theta))
            points_hub.append(i_point)
            i_point += 1

        # Defining the outlet point at the hub section
        Z_hub.append(x_max)
        X_hub.append(R_TE[0, -1] * np.sin(theta))
        Y_hub.append(R_TE[0, -1] * np.cos(theta))
        points_hub.append(i_point)
        i_point += 1

        # Defining the inlet point at the shroud section
        Z_shroud.append(x_min)
        X_shroud.append(R_LE[-1, 0] * np.sin(theta))
        Y_shroud.append(R_LE[-1, 0] * np.cos(theta))
        points_shroud.append(i_point)
        i_point += 1

        # Defining the points defining the shroud shape between the inlet and outlet
        for i in range(n_rows):
            Z_shroud.append(X_LE[-1, i])
            X_shroud.append(R_LE[-1, i] * np.sin(theta))
            Y_shroud.append(R_LE[-1, i] * np.cos(theta))
            points_shroud.append(i_point)
            i_point += 1

            Z_shroud.append(X_TE[-1, i])
            X_shroud.append(R_TE[-1, i] * np.sin(theta))
            Y_shroud.append(R_TE[-1, i] * np.cos(theta))
            points_shroud.append(i_point)
            i_point += 1

        # Defining the outlet point on the shroud section
        Z_shroud.append(x_max)
        X_shroud.append(R_TE[-1, -1] * np.sin(theta))
        Y_shroud.append(R_TE[-1, -1] * np.cos(theta))
        points_shroud.append(i_point)
        i_point += 1

        # Storing the boundary patch coordinates in the class
        self.Coords.Coords_inlet = np.array([[X_hub[0], Y_hub[0], Z_hub[0]],
                                           [X_shroud[0], Y_shroud[0], Z_shroud[0]],
                                           [-X_hub[0], Y_hub[0], Z_hub[0]],
                                           [-X_shroud[0], Y_shroud[0], Z_shroud[0]]])
        self.Coords.Coords_outlet = np.array([[X_hub[-1], Y_hub[-1], Z_hub[-1]],
                                           [X_shroud[-1], Y_shroud[-1], Z_shroud[-1]],
                                           [-X_hub[-1], Y_hub[-1], Z_hub[-1]],
                                           [-X_shroud[-1], Y_shroud[-1], Z_shroud[-1]]])
        self.Coords.Coords_hub = np.transpose(np.array([X_hub + [-x for x in X_hub], Y_hub + Y_hub, Z_hub + Z_hub]))
        self.Coords.Coords_shroud = np.transpose(np.array([X_shroud + [-x for x in X_shroud], Y_shroud + Y_shroud, Z_shroud + Z_shroud]))
        self.Coords.Coords_periodic_1 = np.transpose(np.array([X_hub + X_shroud, Y_hub + Y_shroud, Z_hub + Z_shroud]))
        self.Coords.Coords_periodic_2 = np.transpose(np.array([[-x for x in X_hub] + [-x for x in X_shroud],
                                                               [x for x in Y_hub] + [x for x in Y_shroud],
                                                               [x for x in Z_hub] + [x for x in Z_shroud]]))

        # Storing the hub and shroud point identifier lists
        self.points_hub = points_hub
        self.coords_hub = np.mat([X_hub, Y_hub, Z_hub])
        self.points_shroud = points_shroud
        self.coords_shroud = np.mat([X_shroud, Y_shroud, Z_shroud])
        self.points_count = len(self.points_hub) + len(self.points_shroud)


        # Writing the hub and shroud points in ICEM .rep file
        f = open("ICEM_input.txt", "a")

        f.write("ic_geo_new_family GEOM\n")
        f.write("ic_boco_set_part_color GEOM\n")
        f.write("ic_empty_tetin\n")
        i_point=1
        for i in range(len(X_hub)):
            f.write("ic_point {} GEOM pnt."+str(i_point)+" "+str(X_hub[i])+","+str(Y_hub[i])+","+str(Z_hub[i])+"\n")
            i_point=i_point+1
        for j in range(len(X_hub)):
            f.write("ic_point {} GEOM pnt."+str(i_point+j)+" "+str(X_shroud[j])+","+str(Y_shroud[j])+","+str(Z_shroud[j])+"\n")
        f.close()
 
