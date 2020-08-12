% FUNCTION DRAW_SO3_OBJECT(ROT,OBJECT): 
% 		 Plots the SO(3) object as a 3d Object.
% ------------------------------------------------------------
% This file is part of the project "C1 bezier paths on surfaces"
% and is intended in plotting the surface of Bezier for SO3.
%
% INPUT: ROT	:	The rotation to apply to the object
% 		 OBJECT	:	The 3D object to plot
% 		 PLACE  : 	The position of the object
%
% OUTPUT: The SO(3) object as a 3Dobject.
% ------------------------------------------------------------
% Initial author: Pierre-Yves Gousenbourger
% Versions
% 	28/07/2015: first version.
% ------------------------------------------------------------
function plot_so3_object(rot,object,type,place)

	hold on;
	
	object = object - center_of_mass(object);
	Mp = place + qrot(object, dc2quat(rot));
	switch type
	case 'normal'
		color 	= [0.7,0.7,0.7]; % gris
		color   = [.2,.6,1]; % bleu
		plot(Mp,'FaceColor',color);
	case 'interp'
		color 	= [255 102 102]/255;
		plot(Mp,'FaceColor',color);
	case 'control'; 
		color 	= [102 204   0]/255;
		plot(Mp,'FaceColor',color);
	end
	
	axis off;
end
