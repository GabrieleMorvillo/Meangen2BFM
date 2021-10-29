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

        self.wedge = IN["WEDGE"][0]     # Importing mesh wedge angle
        self.n_sec = int(IN["SECTIONS_PER_DEGREE"][0]) * self.wedge     # Calculating number of tangential nodes
        self.n_point = int(IN["AXIAL_POINTS"][0])       # Importing axial node count
        self.inlet_fac = IN["INLET_FACTOR"][0]          # Importing inlet cell size factor
        self.outlet_fac = IN["OUTLET_FACTOR"][0]        # Importing outlet cell size factor
        self.rot_axis = IN["Rotation_axis"]

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
        f.write("ic_chdir C:/Users/GMrx1/Desktop/ANSYS_WORKS/lecture1\n")

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


        # Creating all the mesh points.
        self.makePoints()

        # Connecting the mesh points to create lines.
        self.makeLines()

        # # Making periodic plane.
        # self.makePlane()

        # # Revolving plane around rotation axis to create 3D wedge.
        # self.revolve()

        # # Setting names to boundaries.
        # self.nameBoundaries()

        # # Creating 3D mesh.
        # gmsh.model.geo.synchronize()

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

        f.close()

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

    # def nameBoundaries(self):
    #     # This function gives names to all the boundaries of the 3D mesh so boundary conditions can be assigned in the
    #     # SU2 configuration file.

    #     # Specifying periodic boundaries.
    #     # First periodic boundary is the plane created by makePlane().
    #     periodic_1 = self.model.addPhysicalGroup(2, [1])
    #     self.model.setPhysicalName(2, periodic_1, "periodic_1")
    #     self.model.setColor([(2, 1)], 255, 255, 0)
    #     # Second periodic boundary is the first plane in the revolution volume list.
    #     periodic_2 = self.model.addPhysicalGroup(2, [self.Rev[0][1]])
    #     self.model.setPhysicalName(2, periodic_2, "periodic_2")
    #     self.model.setColor([self.Rev[-1]], 255, 102, 0)

    #     # Setting the revolution volume as a 3D physical group.
    #     self.model.addPhysicalGroup(3, [1], 1)
    #     self.model.setPhysicalName(3, 1, "FlowField")
    #     # self.model.setColor((3, 1), 255, 102, 0)
    #     # Going over the other surfaces of the revolution volume and naming them accordingly.
    #     i = 2

    #     # Naming the inlet boundary.
    #     inlet = self.model.addPhysicalGroup(2, [self.Rev[i][1]])
    #     self.model.setPhysicalName(2, inlet, "inlet")
    #     self.model.setColor([self.Rev[i]], 255, 0, 0)
    #     i += 1

    #     # Appending all planes on the shroud side of the volume to a single shroud boundary.
    #     shroud_list = []
    #     for j in range(i, i + 3 * self.M.n_stage + 2):
    #         shroud_list.append(self.Rev[j][1])
    #         i += 1
    #     shroud = self.model.addPhysicalGroup(2, shroud_list)
    #     self.model.setPhysicalName(2, shroud, "shroud")
    #     self.model.setColor([(2, j) for j in shroud_list], 0, 255, 0)

    #     # Naming the outlet boundary
    #     outlet = self.model.addPhysicalGroup(2, [self.Rev[i][1]])
    #     self.model.setPhysicalName(2, outlet, "outlet")
    #     self.model.setColor([self.Rev[i]], 0, 0, 255)
    #     i += 1

    #     # Appending all planes on the shroud side of the volume to a single shroud boundary.
    #     hub_list = []
    #     for j in range(i, i + 3 * self.M.n_stage + 2):
    #         hub_list.append(self.Rev[j][1])
    #         i += 1
    #     hub = self.model.addPhysicalGroup(2, hub_list)
    #     self.model.setPhysicalName(2, hub, "hub")
    #     self.model.setColor([(2, j) for j in hub_list], 255, 0, 255)

    # def refineLines(self):
    #     # This function sets the number of nodes along each of the curves bounding the geometry. This allows for
    #     # refinements within the bladed regions along the leading and trailing edges of the blades, while the mesh is
    #     # kept relatively coarse near the inlet and outlet patches

    #     # The application of the mesh refinement process follows the natural progression of the machine, starting with
    #     # the inlet farfield.
    #     # Setting node count along hub and shroud boundaries.
    #     self.model.mesh.setTransfiniteCurve(self.lines_hub[0], self.n_point)
    #     self.model.mesh.setTransfiniteCurve(self.lines_shroud[0], self.n_point)

    #     # Looping over the hub and shroud lines to refine them according to the axial node count specified in the input
    #     # file. In the row and stage gaps, the cell count is adjusted according to the size of the gap with respect to
    #     # the axial chord.
    #     for i in range(1, len(self.lines_hub)-1):
    #         # In case of a row or stage gap, the cell count in axial direction is altered.
    #         if i % 2 == 0:
    #             self.model.mesh.setTransfiniteCurve(self.lines_hub[i], max([int(self.n_point*self.length_hub[i] /
    #                                                                             self.length_hub[i-1]), 3]))
    #             self.model.mesh.setTransfiniteCurve(self.lines_shroud[i], max([int(self.n_point*self.length_shroud[i] /
    #                                                                                self.length_shroud[i-1]), 3]))
    #         else:
    #             self.model.mesh.setTransfiniteCurve(self.lines_hub[i], self.n_point)
    #             self.model.mesh.setTransfiniteCurve(self.lines_shroud[i], self.n_point)
    #     self.model.mesh.setTransfiniteCurve(self.lines_hub[-1], self.n_point)
    #     self.model.mesh.setTransfiniteCurve(self.lines_shroud[-1], self.n_point)

    #     # The lines in radial direction are embedded into the periodic plane surface, which allows for their respective
    #     # refinements to apply to the plane surface. This allows for the cell size to be nearly constant within the
    #     # bladed region.
    #     for i in range(len(self.lines_rad)):
    #         # Embedding the radial line into the surface.
    #         self.model.mesh.embed(1, [self.lines_rad[i]], 2, 1)
    #         # In case of the inlet or outlet line, the cell size is increased according to the inlet and outlet
    #         # refinement factors as specified by the user.
    #         if i == 0:
    #             self.model.mesh.setTransfiniteCurve(self.lines_rad[i], int(self.n_point/self.inlet_fac))
    #         elif i == len(self.lines_rad)-1:
    #             self.model.mesh.setTransfiniteCurve(self.lines_rad[i], int(self.n_point / self.outlet_fac))
    #         else:
    #             # The number of cells in along the radial lines is adjusted according to the respective length with
    #             # respect to the local, average axial chord.
    #             if i % 2 == 0:
    #                 L_av = 0.5 * (self.length_hub[i-1] + self.length_shroud[i-1])
    #             else:
    #                 L_av = 0.5 * (self.length_hub[i] + self.length_shroud[i])
    #             self.model.mesh.setTransfiniteCurve(self.lines_rad[i], max([int(self.n_point*self.length_rad[i]/L_av), 3]))

    # def revolve(self):
    #     # This function revolves the 2D periodic plane around the rotation axis to create a wedge.
    #     # The center point of the revolution is located at the origin and the rotation axis is set to be the x-axis.
    #     self.Rev = self.factory.revolve([(2, 1)], 0, 0, 0, 0, 0, 1, self.wedge*np.pi/180, [self.n_sec])

    # def makePlane(self):
    #     # This function takes the lines defining the hub, shroud, inlet and outlet and builds a plane surface bound by
    #     # the meridional shape.

    #     # Defining a list containing the line identifiers to make a GMesh curve loop.
    #     loop = [-self.lines_rad[0]]
    #     for i in self.lines_hub:
    #         loop.append(i)
    #     loop.append(self.lines_rad[-1])
    #     for i in self.lines_shroud[::-1]:
    #         loop.append(-i)

    #     # Building the curve loop in GMesh
    #     self.factory.addCurveLoop(loop, 1)

    #     # Building the plane surface in GMesh
    #     self.factory.addPlaneSurface([1], 1)

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

        # Looping through the hub and shroud points to build lines defining the hub and shroud shape.
        for i in range(len(self.points_hub)-1):
            self.factory.addLine(self.points_hub[i], self.points_hub[i+1], i_line)
            lines_hub.append(i_line)
            length_hub.append(np.sqrt(np.sum(np.array(self.coords_hub[:, i+1] - self.coords_hub[:, i])**2)))
            i_line += 1
        for i in range(len(self.points_shroud)-1):
            self.factory.addLine(self.points_shroud[i], self.points_shroud[i+1], i_line)
            lines_shroud.append(i_line)
            length_shroud.append(np.sqrt(np.sum(np.array(self.coords_shroud[:, i + 1] - self.coords_shroud[:, i]) ** 2)))
            i_line += 1

        # Looping through the points in axial direction to build lines in radial direction. These will be used for the
        # definition of the inlet and outlet and refinement along the blade leading and trailing edges.
        for i in range(len(self.points_hub)):
            self.factory.addLine(self.points_hub[i], self.points_shroud[i], i_line)
            lines_rad.append(i_line)
            length_rad.append(np.sqrt(np.sum(np.array(self.coords_shroud[:, i] - self.coords_hub[:, i]) ** 2)))
            i_line += 1

        # Storing the line identifiers and line lengths into the class.
        self.lines_hub = lines_hub
        self.lines_shroud = lines_shroud
        self.lines_rad = lines_rad

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

        # # Building the hub and shroud points in GMesh
        # for i in range(len(X_hub)):
        #     self.factory.addPoint(X_hub[i], Y_hub[i], Z_hub[i], 0.01, points_hub[i])
        #     self.factory.addPoint(X_shroud[i], Y_shroud[i], Z_shroud[i], 0.01, points_shroud[i])

        # Write the replay procedure for the points

        # Writing the hub and shroud points in ICEM .rep file
        f.write("ic_geo_new_family GEOM\n")
        f.write("ic_boco_set_part_color GEOM\n")
        f.write("ic_empty_tetin\n")
        for i in range(len(X_hub)):
            f.write("ic_point {} GEOM pnt."+str(i)+" "+str(X_hub[i])+","+str(Y_hub[i])+","+str(Z_hub[i])+"\n")
        for i in range(len(X_hub)):
            f.write("ic_point {} GEOM pnt."+str(i)+" "+str(X_shroud[i])+","+str(Y_shroud[i])+","+str(Z_shroud[i])+"\n")

 
