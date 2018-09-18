sock = tcpip("localhost",30000,'NetworkRole', 'server');
sock.InputBufferSize = 100000*8;

DATA_FILENAME = "S_DATA.csv";
SONNET_PROJ_DIRNAME = "Sonnet_projects";

while 1
    % waiting for connection
    disp(get(sock,"Status"));
    fopen(sock);
    disp(get(sock,"Status"));
    
    % cleanup before previous run
    status = rmdir(pwd + "\\" + SONNET_PROJ_DIRNAME,'s');
    status = rmdir(pwd + "\\sondata",'s');
    status = mkdir( pwd + "\\" + SONNET_PROJ_DIRNAME);
    
    proj = SonnetProject();
    proj.saveAs( 'Sonnet_projects/matlab_proj.son');

    proj.changeLengthUnit("UM");
    proj.changeDielectricLayerThickness(1,500);
    proj.changeDielectricLayerThickness(2,3000);

    metal_type_name = "Al-supercond";
    proj.defineNewResistorMetalType(metal_type_name,0);
    csv_name = pwd + "\" + SONNET_PROJ_DIRNAME + "\" + DATA_FILENAME;
    while 1
        data = fread(sock, 1,"uint16");
        if data == CMD.CLOSE
            respond_ok(sock);
            fclose(sock);
            break
        elseif data == CMD.SAY_HELLO
            respond_ok(sock);
            disp("HELLO");
        elseif data == CMD.POLYGON
            respond_ok(sock);
            polygon = receive_polygon(sock);
            
            % ATOMIC EXPRESSION START
            proj.addMetalPolygonEasy(0,polygon.points_x,polygon.points_y,1);
            poly_idx = length(proj.GeometryBlock.ArrayOfPolygons);
            if polygon.ports == FLAG.TRUE
                for i = transpose(polygon.port_edges_num_list)
                    proj.addPortToPolygon(poly_idx,i);
                end
            end
            % ATOMIC EXPRESSION END
        elseif data == CMD.BOX_PROPS
            respond_ok(sock);
            boxSettings = receive_boxProps(sock);
            proj.changeBoxSizeXY(boxSettings.dim_X_um, boxSettings.dim_Y_um);
            proj.changeNumberOfCells(boxSettings.cells_X_num,boxSettings.cells_Y_num);
        elseif data == CMD.CLEAR_POLYGONS
            respond_ok(sock);
            for i = 1:length(proj.GeometryBlock.ArrayOfPolygons)
                proj.deletePolygonUsingIndex(1);
            end
        elseif data == CMD.SET_ABS
            respond_ok(sock);
            abs_params = receive_abs_parameters(sock);
            proj.addAbsFrequencySweep(abs_params.start_freq,abs_params.stop_freq);
        elseif data == CMD.SIMULATE
            respond_ok(sock);
            % output csv file
            % de-embeded data
            % including comments
            % with high precision
            % S-data
            % real-image data pairs
            proj.addFileOutput("CSV","D","Y",DATA_FILENAME,"IC","Y","S","RI","R",50);
            proj.simulate('-c');
            % sending current working directory
            
            fwrite(sock,csv_name + newline);
        elseif data == CMD.VISUALIZE
            respond_ok(sock);
            response_data = csvread(csv_name,8);
            freq = response_data(:,1);
            s21_re = response_data(:,4);
            s21_im = response_data(:,5);
            plot(freq,s21_re);
            drawnow;
        end
    end
end

function respond_ok(sock)
    fwrite(sock,RESPONSE.OK,"uint16");
end

function result=receive_flag(sock)
    result = fread(sock,1,"uint16");
    respond_ok(sock);
end

function result=receive_uint32_x1(sock)
    result = fread(sock,1,"uint32");
    respond_ok(sock);
end

function result=receive_uint32_xnum(sock)
    num = receive_uint32_x1(sock);
    idxs_array = fread( sock, num, "uint32" );
    respond_ok(sock);
    
    result = idxs_array;
end

function result=receive_float64_x1(sock)
    result=fread(sock,1,"float64");
    respond_ok(sock);
end 

function result=receive_float64_xnum(sock)
    num = fread(sock,1,"uint32");
    respond_ok(sock);
    array_float = fread( sock, num, "float64" );
    respond_ok(sock);
    
    result = array_float;
end

function result_poly=receive_polygon(sock)
    result_poly = Polygon();
    
    result_poly.ports = receive_flag(sock);
    if result_poly.ports == FLAG.TRUE 
        result_poly.port_edges_num_list = receive_uint32_xnum(sock);
    else
        result_poly.port_edges_num_list = -1;
    end
    
    result_poly.points_x = receive_float64_xnum(sock);
    result_poly.points_y = receive_float64_xnum(sock);
end

function boxSettings=receive_boxProps(sock)
    boxSettings = BoxProps();
    boxSettings.dim_X_um = receive_float64_x1(sock);
    boxSettings.dim_Y_um = receive_float64_x1(sock);
    boxSettings.cells_X_num = receive_uint32_x1(sock);
    boxSettings.cells_Y_num = receive_uint32_x1(sock);
end

function absParams=receive_abs_parameters(sock)
    absParams = ABSparams();
    absParams.start_freq = receive_float64_x1(sock);
    absParams.stop_freq = receive_float64_x1(sock);
end