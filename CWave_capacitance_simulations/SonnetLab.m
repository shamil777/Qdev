% Superconducting Aluminia on Si 500 um substrate setup
Project = SonnetProject();
Project.saveAs('CWave.son');

% set units to um
Project.changeLengthUnit("UM");

% setup dielectric layers
substrate_layer = Project.getLayer(2);
substrate_layer.NameOfDielectricLayer = "Si";
substrate_layer.Thickness = 500;
substrate_layer.RelativeDielectricConstant = 13.6;

air_layer = Project.getLayer(1);
air_layer.NameOfDielectricLayer = "Air";
air_layer.Thickness = 3000;
air_layer.RelativeDielectricConstant = 1.0;

% adding metal layer
%Project.addMetalTechLayer('Al','DXF',0,1,'',0);
Project.defineNewResistorMetalType("Al",0);

% try to import
Project.addMetalTechLayer("Al","Metal",1,0,"Cwave_first_0.gds",0);

a = 20;
Project.addMetalPolygonEasy(0,[0,0,a,a],[0,a,a,0],0);

% open project in sonnet
Project.openInSonnet();
