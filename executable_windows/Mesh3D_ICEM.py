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
        DIR = os.getcwd() 
        if not os.path.isdir("MESHOutput"):
            os.mkdir(DIR + "\\MESHOutput")
           
        DIRMESH = DIR + "\\MESHOutput"  

        # Write the initial paragraph for the replay file
        f = open("ICEM_input.txt", "w")

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
        
        # Save mesh to ANSYS CFX input file
        self.savemesh()

        # # Applying mesh refinement along the lines in the mesh
        # self.refineLines()

        # # Saving mesh as su2 file
        # gmsh.write(self.fileName)

        os.system("copy ICEM_input.txt " + DIRMESH +"\\ICEM_input.rpl")

        

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
    
    def savemesh(self):
        DIR = os.getcwd()
        DIRMESH = DIR + "\\MESHOutput"
        DIRMESH = DIRMESH.replace("\\","/")
        # # os.chdir("\\MESHOutput")

        f = open("ICEM_input.txt", "a")
        f.write("ic_hex_write_file ./hex.uns GEOM RIGHT_SYM LEFT_SYM HUB_WALL SHROUD_WALL INLET OUTLET SOLID proj 2 dim_to_mesh 3 no_boco\n")
        f.write("ic_uns_load ./hex.uns 3 0 {} 1\n")
        f.write("ic_uns_update_family_type visible {INLET SHROUD_WALL GEOM OUTLET ORFN HUB_WALL RIGHT_SYM SOLID LEFT_SYM} {!NODE !LINE_2 QUAD_4 !HEXA_8} update 0\n")
        f.write("ic_boco_solver\n")
        f.write("ic_boco_clear_icons\n")
        f.write("ic_uns_update_family_type visible {INLET SHROUD_WALL GEOM OUTLET ORFN HUB_WALL RIGHT_SYM SOLID LEFT_SYM} {!NODE LINE_2 QUAD_4 !HEXA_8} update 0\n")
        f.write("ic_boco_solver CGNS\n")
        f.write("ic_solver_mesh_info CGNS\n")
        f.write("ic_boco_solver\n")
        f.write("ic_boco_solver CGNS\n")
        f.write("ic_solution_set_solver CGNS 1\n")
        f.write("ic_boco_save {" + DIRMESH +"\\ICEM_MESH.fbc}\n")
        f.write("ic_boco_save_atr {" + DIRMESH +"\\ICEM_MESH.atr}\n")
        f.write("ic_delete_empty_parts\n")
        f.write("ic_save_tetin project1.tin 0 0 {} {} 0 0 1\n")
        f.write("ic_uns_check_duplicate_numbers\n")
        f.write("ic_save_unstruct project1.uns 1 {} {} {}\n")
        f.write("ic_uns_set_modified 1\n")
        f.write("ic_hex_save_blocking project1.blk\n")
        f.write("ic_boco_solver\n")
        f.write("ic_boco_solver CGNS\n")
        f.write("ic_solution_set_solver CGNS 1\n")
        f.write("ic_boco_save project1.fbc\n")
        f.write("ic_boco_save_atr project1.atr\n")
        f.write("ic_save_project_file "+DIRMESH+"/project1.prj {array\ set\ file_name\ \{ {    catia_dir .} {    parts_dir .} {    domain_loaded 0} {    cart_file_loaded 0} {    cart_file {}} {    domain_saved project1.uns} {    archive {}} {    med_replay {}} {    topology_dir .} {    ugparts_dir .} {    icons {{$env(ICEM_ACN)/lib/ai_env/icons} {$env(ICEM_ACN)/lib/va/EZCAD/icons} {$env(ICEM_ACN)/lib/icons} {$env(ICEM_ACN)/lib/va/CABIN/icons}}} {    tetin project1.tin} {    family_boco project1.fbc} {    iges_dir .} {    solver_params_loaded 0} {    attributes_loaded 0} {    project_lock {}} {    attributes project1.atr} {    domain project1.uns} {    domains_dir .} {    settings_loaded 0} {    settings project1.prj} {    blocking project1.blk} {    hexa_replay {}} {    transfer_dir .} {    mesh_dir .} {    family_topo {}} {    gemsparts_dir .} {    family_boco_loaded 0} {    tetin_loaded 0} {    project_dir .} {    topo_mulcad_out {}} {    solver_params {}} \} array\ set\ options\ \{ {    expert 1} {    remote_path {}} {    tree_disp_quad 2} {    tree_disp_pyra 0} {    evaluate_diagnostic 0} {    histo_show_default 1} {    select_toggle_corners 0} {    remove_all 0} {    keep_existing_file_names 0} {    record_journal 0} {    edit_wait 0} {    face_mode all} {    select_mode all} {    med_save_emergency_tetin 1} {    user_name GMrx1} {    diag_which all} {    uns_warn_if_display 500000} {    bubble_delay 1000} {    external_num 1} {    tree_disp_tri 2} {    apply_all 0} {    default_solver {ANSYS Fluent}} {    temporary_directory {}} {    flood_select_angle 0} {    home_after_load 1} {    project_active 0} {    histo_color_by_quality_default 1} {    undo_logging 1} {    tree_disp_hexa 0} {    histo_solid_default 1} {    host_name LAPTOP-1G7TMHD3} {    xhidden_full 1} {    replay_internal_editor 1} {    editor notepad} {    mouse_color orange} {    clear_undo 1} {    remote_acn {}} {    remote_sh csh} {    tree_disp_penta 0} {    n_processors 1} {    remote_host {}} {    save_to_new 0} {    quality_info Quality} {    tree_disp_node 0} {    med_save_emergency_mesh 1} {    redtext_color red} {    tree_disp_line 0} {    select_edge_mode 0} {    use_dlremote 0} {    max_mesh_map_size 1024} {    show_tris 1} {    remote_user {}} {    enable_idle 0} {    auto_save_views 1} {    max_cad_map_size 512} {    display_origin 0} {    uns_warn_user_if_display 1000000} {    detail_info 0} {    win_java_help 0} {    show_factor 1} {    boundary_mode all} {    clean_up_tmp_files 1} {    auto_fix_uncovered_faces 1} {    med_save_emergency_blocking 1} {    max_binary_tetin 0} {    tree_disp_tetra 0} \} array\ set\ disp_options\ \{ {    uns_dualmesh 0} {    uns_warn_if_display 500000} {    uns_normals_colored 0} {    uns_icons 0} {    uns_locked_elements 0} {    uns_shrink_npos 0} {    uns_node_type None} {    uns_icons_normals_vol 0} {    uns_bcfield 0} {    backup Wire} {    uns_nodes 0} {    uns_only_edges 0} {    uns_surf_bounds 0} {    uns_wide_lines 0} {    uns_vol_bounds 0} {    uns_displ_orient Triad} {    uns_orientation 0} {    uns_directions 0} {    uns_thickness 0} {    uns_shell_diagnostic 0} {    uns_normals 0} {    uns_couplings 0} {    uns_periodicity 0} {    uns_single_surfaces 0} {    uns_midside_nodes 1} {    uns_shrink 100} {    uns_multiple_surfaces 0} {    uns_no_inner 0} {    uns_enums 0} {    uns_disp Wire} {    uns_bcfield_name {}} {    uns_color_by_quality 0} {    uns_changes 0} {    uns_cut_delay_count 1000} \} {set icon_size1 24} {set icon_size2 35} {set thickness_defined 0} {set solver_type 1} {set solver_setup -1} array\ set\ prism_values\ \{ {    n_triangle_smoothing_steps 5} {    min_smoothing_steps 6} {    first_layer_smoothing_steps 1} {    new_volume {}} {    height {}} {    prism_height_limit {}} {    interpolate_heights 0} {    n_tetra_smoothing_steps 10} {    do_checks {}} {    delete_standalone 1} {    ortho_weight 0.50} {    max_aspect_ratio {}} {    ratio_max {}} {    incremental_write 0} {    total_height {}} {    use_prism_v10 0} {    intermediate_write 1} {    delete_base_triangles {}} {    ratio_multiplier {}} {    verbosity_level 1} {    refine_prism_boundary 1} {    max_size_ratio {}} {    triangle_quality {}} {    max_prism_angle 180} {    tetra_smooth_limit 0.3} {    max_jump_factor 5} {    use_existing_quad_layers 0} {    layers 3} {    fillet 0.10} {    into_orphan 0} {    init_dir_from_prev {}} {    blayer_2d 0} {    do_not_allow_sticking {}} {    top_family {}} {    law exponential} {    min_smoothing_val 0.1} {    auto_reduction 0} {    stop_columns 1} {    stair_step 1} {    smoothing_steps 12} {    side_family {}} {    min_prism_quality 0.01} {    ratio 1.2} \} {set aie_current_flavor {}} array\ set\ vid_options\ \{ {    wb_import_mat_points 0} {    wb_NS_to_subset 0} {    wb_import_surface_bodies 1} {    wb_import_cad_att_pre {SDFEA;DDM}} {    wb_import_mix_res_line 0} {    wb_import_tritol 0.001} {    auxiliary 1} {    wb_import_cad_att_trans 1} {    wb_import_mix_res -1} {    wb_import_mix_res_surface 0} {    show_name 0} {    wb_import_solid_bodies 1} {    wb_import_delete_solids 0} {    do_intersect_self_part 1} {    wb_import_mix_res_solid 0} {    wb_import_save_pmdb {}} {    inherit 0} {    default_part GEOM} {    new_srf_topo 0} {    wb_import_associativity_model_name {}} {    DelPerFlag 0} {    show_item_name 0} {    wb_import_line_bodies 0} {    wb_import_save_partfile 0} {    composite_tolerance 1.0} {    wb_NS_to_entity_parts 0} {    wb_import_en_sym_proc 1} {    wb_import_sel_proc 1} {    wb_import_work_points 0} {    wb_import_reference_key 0} {    wb_import_mix_res_point 0} {    wb_import_pluginname {}} {    wb_NS_only 0} {    wb_import_geom 0} {    wb_import_create_solids 0} {    wb_import_refresh_pmdb 0} {    wb_import_lcs 0} {    wb_import_sel_pre {}} {    wb_import_scale_geo Default} {    wb_import_load_pmdb {}} {    replace 0} {    wb_import_cad_associativity 0} {    same_pnt_tol 1e-4} {    tdv_axes 1} {    wb_import_mesh 0} {    vid_mode 0} {    DelBlkPerFlag 0} \} {set savedTreeVisibility {geomNode 1 geom_subsetNode 2 geomPointNode 0 geomCurveNode 2 geomSurfNode 2 meshNode 1 mesh_subsetNode 2 meshPointNode 0 meshLineNode 2 meshShellNode 2 meshQuadNode 2 meshVolumeNode 0 meshHexaNode 0 blockingNode 1 block_subsetNode 2 block_vertNode 0 block_edgeNode 2 block_faceNode 0 block_blockNode 0 block_meshNode 0 topoNode 2 topo-root 2 partNode 2 part-GEOM 2 part-HUB_WALL 2 part-INLET 2 part-LEFT_SYM 2 part-OUTLET 2 part-RIGHT_SYM 2 part-SHROUD_WALL 2 part-SOLID 2 part-VORFN 0}} {set last_view {rot {0.0225677294222 -0.946480402031 -0.0425933790289 -0.319141583359} scale {857.662704449 857.662704449 857.662704449} center {0 0 0} pos {173.289871722 -283.900428141 0}}} array\ set\ cut_info\ \{ {    active 0} {    whole 1} \} array\ set\ hex_option\ \{ {    default_bunching_ratio 2.0} {    floating_grid 0} {    project_to_topo 0} {    n_tetra_smoothing_steps 20} {    sketching_mode 0} {    trfDeg 1} {    wr_hexa7 0} {    smooth_ogrid 0} {    find_worst 1-3} {    hexa_verbose_mode 0} {    old_eparams 0} {    uns_face_mesh_method uniform_quad} {    multigrid_level 0} {    uns_face_mesh one_tri} {    check_blck 0} {    proj_limit 0} {    check_inv 0} {    project_bspline 0} {    hexa_update_mode 1} {    default_bunching_law BiGeometric} {    worse_criterion Quality} \} array\ set\ saved_views\ \{ {    views {}} \}} {ICEM CFD}\n")
        f.write("ic_write_file domain_list {"+DIRMESH+"/project1.uns\n}\n")
        f.write("ic_exec {C:/Program Files/ANSYS Inc/v195/icemcfd/win64_amd/icemcfd/output-interfaces/cgns} -b project1.fbc -dom_list domain_list -unstr -scale 1.0 ./project1.cgns\n")
        f.write("exit\n")







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
 
